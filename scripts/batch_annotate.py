"""
batch_annotate.py — Claude Batch API annotation for unannotated reactions.

Pipeline:
  1. Load model.xml + data/unannotated_reactions.csv
  2. For each reaction: collect metabolite info, pre-screen candidates from
     ec_to_mnxr and fingerprint_index (up to 10 per reaction)
  3. Build Batch API requests with structured prompts
  4. Submit via anthropic Message Batches API (or --dry-run to inspect JSON)
  5. Poll until complete, download results
  6. Save data/batch_annotation_results.csv with audit columns

Usage:
  python scripts/batch_annotate.py \\
      --model model.xml \\
      --csv data/unannotated_reactions.csv \\
      --mnx-dir data/metanetx \\
      --api-key sk-ant-...

  python scripts/batch_annotate.py --dry-run   # write requests JSON only
"""

import argparse
import csv
import hashlib
import json
import logging
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

MODEL_ID = "claude-sonnet-4-20250514"
MAX_TOKENS = 50
CANDIDATES_PER_REACTION = 10
_MNXM_IN_EQ = re.compile(r"(MNXM\d+)@")
_UBIQUITOUS_MNXM = frozenset({"MNXM1", "WATER", "MNXM3", "MNXM9", "MNXM5"})


# ── MetaNetX loaders (minimal, self-contained) ────────────────────────────────

def _read_tsv(path: Path, names: list[str]):
    import pandas as pd
    return pd.read_csv(
        path, sep="\t", comment="#", header=None,
        names=names, dtype=str, low_memory=False,
    ).fillna("")


def _load_reac_xref(path: Path) -> dict:
    logger.info("Loading reac_xref.tsv …")
    df = _read_tsv(path, ["source", "mnx_id", "description"])
    df = df[df["mnx_id"].str.startswith("MNXR")]

    by_mnxr: dict[str, list] = defaultdict(list)
    ec_to_mnxr: dict[str, list[str]] = defaultdict(list)
    desc_index: dict[str, str] = {}
    mnxr_to_desc: dict[str, str] = {}

    for source, mnx_id, desc in df.itertuples(index=False):
        if ":" in source:
            prefix, sid = source.split(":", 1)
            by_mnxr[mnx_id].append((prefix, sid))
            if prefix == "ec-code":
                ec_to_mnxr[sid].append(mnx_id)
        desc_head = desc.split("||")[0].strip() if "||" in desc else desc.strip()
        if desc_head.startswith("EC:"):
            ec_num = desc_head[3:]
            if ec_num:
                ec_to_mnxr[ec_num].append(mnx_id)
        else:
            short = desc_head.lower()
            if short and short not in desc_index:
                desc_index[short] = mnx_id
                if mnx_id not in mnxr_to_desc:
                    mnxr_to_desc[mnx_id] = short

    logger.info(f"  {len(by_mnxr):,} MNXR IDs, {len(ec_to_mnxr):,} EC numbers")
    return {
        "by_mnxr": dict(by_mnxr),
        "ec_to_mnxr": dict(ec_to_mnxr),
        "desc_index": desc_index,
        "mnxr_to_desc": mnxr_to_desc,
    }


def _load_reac_prop(path: Path) -> dict:
    logger.info("Loading reac_prop.tsv …")
    df = _read_tsv(path, ["mnx_id", "equation", "reference", "classifs",
                           "is_balanced", "is_transport"])
    df = df[df["mnx_id"].str.startswith("MNXR")]

    fingerprint_index: dict[frozenset, list[str]] = defaultdict(list)
    mnxr_to_equation: dict[str, str] = {}

    for mnx_id, equation, *_ in df.itertuples(index=False):
        mnxm_ids = frozenset(_MNXM_IN_EQ.findall(equation))
        if len(mnxm_ids) >= 2:
            fingerprint_index[mnxm_ids].append(mnx_id)
        mnxr_to_equation[mnx_id] = equation

    logger.info(f"  {len(fingerprint_index):,} fingerprints indexed")
    return {
        "fingerprint_index": dict(fingerprint_index),
        "mnxr_to_equation": mnxr_to_equation,
    }


# ── Candidate pre-screening ───────────────────────────────────────────────────

def _get_ec_candidates(ec_list: list[str], ec_to_mnxr: dict) -> list[str]:
    seen: list[str] = []
    dedup: set[str] = set()
    for ec in ec_list:
        for mnxr in ec_to_mnxr.get(ec, []):
            if mnxr not in dedup:
                dedup.add(mnxr)
                seen.append(mnxr)
    return seen


def _get_fingerprint_candidates(
    mnxm_ids: list[str],
    fingerprint_index: dict,
) -> list[str]:
    query = frozenset(m for m in mnxm_ids if m and m != "—") - _UBIQUITOUS_MNXM
    if len(query) < 2:
        return []
    exact = fingerprint_index.get(query, [])
    if exact:
        return list(exact)
    # Jaccard similarity fallback: find fingerprints that share ≥1 MNXM
    hits: dict[str, float] = {}
    for fp, mnxr_list in fingerprint_index.items():
        fp_filt = fp - _UBIQUITOUS_MNXM or fp
        if not (query & fp_filt):
            continue
        score = len(query & fp_filt) / len(query | fp_filt)
        if score >= 0.3:
            for mnxr in mnxr_list:
                if score > hits.get(mnxr, -1):
                    hits[mnxr] = score
    return sorted(hits, key=lambda m: hits[m], reverse=True)


# ── Prompt construction ───────────────────────────────────────────────────────

def _build_prompt(row: dict, candidates: list[str], by_mnxr: dict,
                  mnxr_to_desc: dict, mnxr_to_equation: dict) -> str:
    met_names = row.get("metabolite_names", "").split("|")
    met_mnxm  = row.get("metabolite_mnxm_ids", "").split(",")
    met_pairs = [
        f"{name.strip()} ({mnxm.strip()})"
        for name, mnxm in zip(met_names, met_mnxm)
    ]
    met_list = "; ".join(met_pairs) if met_pairs else "—"

    cand_lines = []
    for i, mnxr in enumerate(candidates[:CANDIDATES_PER_REACTION], 1):
        desc = mnxr_to_desc.get(mnxr, "")
        eq   = mnxr_to_equation.get(mnxr, "")
        cand_lines.append(f"{i}. {mnxr}: {desc or '(no description)'} | equation: {eq or '(unavailable)'}")
    cand_block = "\n".join(cand_lines) if cand_lines else "(no candidates identified)"

    return (
        "You are a metabolic pathway expert. Given this reaction from a "
        "Yarrowia lipolytica genome-scale model, select the best matching "
        "MetaNetX reaction ID (MNXR) from the candidates, or respond \"NONE\" "
        "if no candidate is appropriate.\n\n"
        f"Reaction: {row['reaction_name']}\n"
        f"Metabolites: {met_list}\n"
        f"EC: {row.get('ec_numbers', '') or '—'}\n"
        f"Compartments: {row.get('compartments', '')}\n\n"
        f"Candidates:\n{cand_block}\n\n"
        "Respond with ONLY the MNXR ID or \"NONE\". No explanation."
    )


def _prompt_hash(prompt: str) -> str:
    return hashlib.sha256(prompt.encode()).hexdigest()[:8]


# ── Model helpers ─────────────────────────────────────────────────────────────

def _get_annotation_value(ann: dict, key: str) -> str:
    if not ann:
        return ""
    raw = ann.get(key, "")
    if isinstance(raw, list):
        return raw[0] if raw else ""
    return str(raw) if raw else ""


def _extract_ec_numbers(rxn) -> list[str]:
    ecs: list[str] = []
    seen: set[str] = set()
    for gene in rxn.genes:
        g_ann = gene.annotation if isinstance(gene.annotation, dict) else {}
        raw = g_ann.get("ec-code", [])
        if isinstance(raw, str):
            raw = [raw]
        for ec in raw:
            if ec and ec not in seen:
                seen.add(ec)
                ecs.append(ec)
    return ecs


# ── Main pipeline ─────────────────────────────────────────────────────────────

def build_requests(
    model_path: Path,
    csv_path: Path,
    mnx_dir: Path,
) -> tuple[list[dict], list[dict]]:
    """
    Returns (batch_requests, metadata_rows) where:
      batch_requests  — list of dicts ready for Batch API
      metadata_rows   — parallel list with reaction context for result parsing
    """
    try:
        from cobra.io import read_sbml_model
    except ImportError:
        logger.error("COBRApy not installed: pip install cobra")
        sys.exit(1)

    logger.info(f"Loading model: {model_path}")
    model = read_sbml_model(str(model_path))

    # Load unannotated CSV (reaction_id set)
    unannotated_ids: set[str] = set()
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            unannotated_ids.add(row["reaction_id"])
    logger.info(f"  {len(unannotated_ids)} unannotated reactions in CSV")

    # Load MetaNetX indexes
    xref = _load_reac_xref(mnx_dir / "reac_xref.tsv")
    prop = _load_reac_prop(mnx_dir / "reac_prop.tsv")

    by_mnxr           = xref["by_mnxr"]
    ec_to_mnxr        = xref["ec_to_mnxr"]
    mnxr_to_desc      = xref["mnxr_to_desc"]
    fingerprint_index = prop["fingerprint_index"]
    mnxr_to_equation  = prop["mnxr_to_equation"]

    # Build reaction lookup
    rxn_by_id = {rxn.id: rxn for rxn in model.reactions}

    batch_requests: list[dict] = []
    metadata_rows:  list[dict] = []
    no_candidates = 0

    for rxn_id in sorted(unannotated_ids):
        rxn = rxn_by_id.get(rxn_id)
        if rxn is None:
            logger.warning(f"  Reaction {rxn_id} not found in model — skipping")
            continue

        mets = list(rxn.metabolites)
        met_names = [m.name or m.id for m in mets]
        met_mnxm  = [
            _get_annotation_value(
                m.annotation if isinstance(m.annotation, dict) else {}, "metanetx.chemical"
            ) or "—"
            for m in mets
        ]
        comps = sorted({m.compartment for m in mets})
        ec_list = _extract_ec_numbers(rxn)

        # Pre-screen candidates
        ec_cands   = _get_ec_candidates(ec_list, ec_to_mnxr)
        fp_cands   = _get_fingerprint_candidates(met_mnxm, fingerprint_index)
        # Merge: EC candidates first (higher precision), then fingerprint
        merged: list[str] = []
        seen_mnxr: set[str] = set()
        for mnxr in ec_cands + fp_cands:
            if mnxr not in seen_mnxr:
                seen_mnxr.add(mnxr)
                merged.append(mnxr)

        candidates = merged[:CANDIDATES_PER_REACTION]
        candidate_tag = "no_candidates" if not candidates else ",".join(candidates)

        row = {
            "reaction_id":         rxn_id,
            "reaction_name":       rxn.name or rxn_id,
            "metabolite_names":    "|".join(met_names),
            "metabolite_mnxm_ids": ",".join(met_mnxm),
            "compartments":        ",".join(comps),
            "gene_reaction_rule":  rxn.gene_reaction_rule or "",
            "ec_numbers":          ",".join(ec_list),
            "candidates":          candidate_tag,
        }

        if not candidates:
            no_candidates += 1

        prompt = _build_prompt(row, candidates, by_mnxr, mnxr_to_desc, mnxr_to_equation)
        phash  = _prompt_hash(prompt)

        batch_requests.append({
            "custom_id": rxn_id,
            "params": {
                "model":      MODEL_ID,
                "max_tokens": MAX_TOKENS,
                "messages": [{"role": "user", "content": prompt}],
            },
        })
        metadata_rows.append({
            **row,
            "prompt_hash": phash,
            "prompt":      prompt,
        })

    logger.info(
        f"Built {len(batch_requests)} batch requests  "
        f"({no_candidates} with no candidates)"
    )
    return batch_requests, metadata_rows


def submit_batch(requests: list[dict], api_key: str) -> str:
    import anthropic
    client = anthropic.Anthropic(api_key=api_key)
    logger.info(f"Submitting {len(requests)} requests to Batch API …")
    batch = client.messages.batches.create(requests=requests)
    logger.info(f"Batch created: {batch.id}  status={batch.processing_status}")
    return batch.id


def poll_batch(batch_id: str, api_key: str, poll_interval: int = 30) -> None:
    import anthropic
    client = anthropic.Anthropic(api_key=api_key)
    while True:
        batch = client.messages.batches.retrieve(batch_id)
        status = batch.processing_status
        counts = batch.request_counts
        logger.info(
            f"[{batch_id}] status={status}  "
            f"succeeded={counts.succeeded}  errored={counts.errored}  "
            f"processing={counts.processing}  canceled={counts.canceled}"
        )
        if status == "ended":
            break
        time.sleep(poll_interval)


def download_results(
    batch_id: str,
    api_key: str,
    metadata_rows: list[dict],
    out_path: Path,
) -> None:
    import anthropic
    client = anthropic.Anthropic(api_key=api_key)

    # Index metadata by reaction_id
    meta_by_id = {r["reaction_id"]: r for r in metadata_rows}

    timestamp = datetime.now(timezone.utc).isoformat()
    output_rows: list[dict] = []

    for result in client.messages.batches.results(batch_id):
        rxn_id   = result.custom_id
        meta     = meta_by_id.get(rxn_id, {})

        if result.result.type == "succeeded":
            raw_response = result.result.message.content[0].text.strip()
        elif result.result.type == "errored":
            raw_response = f"ERROR:{result.result.error.type}"
        else:
            raw_response = f"CANCELED"

        # Parse MNXR ID: accept MNXRnnnnn or NONE
        matched_mnxr = ""
        if re.match(r"^MNXR\d+$", raw_response):
            matched_mnxr = raw_response
        elif raw_response.upper() == "NONE":
            matched_mnxr = "NONE"
        else:
            matched_mnxr = "PARSE_ERROR"

        output_rows.append({
            "reaction_id":    rxn_id,
            "reaction_name":  meta.get("reaction_name", ""),
            "matched_mnxr":   matched_mnxr,
            "source":         "claude-batch",
            "model_used":     MODEL_ID,
            "timestamp":      timestamp,
            "prompt_hash":    meta.get("prompt_hash", ""),
            "raw_response":   raw_response,
            "candidates":     meta.get("candidates", ""),
            "ec_numbers":     meta.get("ec_numbers", ""),
        })

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "reaction_id", "reaction_name", "matched_mnxr",
        "source", "model_used", "timestamp", "prompt_hash",
        "raw_response", "candidates", "ec_numbers",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    matched = sum(1 for r in output_rows if r["matched_mnxr"].startswith("MNXR"))
    none_   = sum(1 for r in output_rows if r["matched_mnxr"] == "NONE")
    errors  = len(output_rows) - matched - none_
    logger.info(
        f"Results: {len(output_rows)} total | "
        f"matched={matched}  NONE={none_}  errors/parse_errors={errors}"
    )
    logger.info(f"Saved to: {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch-annotate unannotated reactions via Claude Batch API"
    )
    parser.add_argument("--model",   default="model.xml",
                        help="Path to SBML model (default: model.xml)")
    parser.add_argument("--csv",     default="data/unannotated_reactions.csv",
                        help="Unannotated reactions CSV (default: data/unannotated_reactions.csv)")
    parser.add_argument("--mnx-dir", default="data/metanetx",
                        help="MetaNetX data directory (default: data/metanetx)")
    parser.add_argument("--api-key", default=None,
                        help="Anthropic API key (or set ANTHROPIC_API_KEY env var)")
    parser.add_argument("--out",     default="data/batch_annotation_results.csv",
                        help="Output CSV path (default: data/batch_annotation_results.csv)")
    parser.add_argument("--requests-json", default="data/batch_requests.json",
                        help="Path to save request JSON (default: data/batch_requests.json)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Build requests and save JSON without submitting")
    parser.add_argument("--batch-id", default=None,
                        help="Resume: poll + download a previously submitted batch")
    parser.add_argument("--poll-interval", type=int, default=30,
                        help="Seconds between status polls (default: 30)")
    args = parser.parse_args()

    model_path = Path(args.model)
    csv_path   = Path(args.csv)
    mnx_dir    = Path(args.mnx_dir)
    out_path   = Path(args.out)
    req_path   = Path(args.requests_json)

    for p, label in [(model_path, "--model"), (csv_path, "--csv"), (mnx_dir, "--mnx-dir")]:
        if not p.exists():
            logger.error(f"{label} not found: {p}")
            sys.exit(1)

    # Resolve API key (skip for dry-run)
    api_key = args.api_key
    if not api_key and not args.dry_run and not args.batch_id:
        import os
        api_key = os.environ.get("ANTHROPIC_API_KEY")
        if not api_key:
            logger.error("Provide --api-key or set ANTHROPIC_API_KEY")
            sys.exit(1)

    # ── Resume mode: only poll + download ────────────────────────────────────
    if args.batch_id:
        if not api_key:
            import os
            api_key = os.environ.get("ANTHROPIC_API_KEY")
        # Reload metadata from saved requests JSON if available
        metadata_rows: list[dict] = []
        if req_path.exists():
            with open(req_path, encoding="utf-8") as fh:
                saved = json.load(fh)
                metadata_rows = saved.get("metadata", [])
        else:
            logger.warning(f"requests JSON not found at {req_path} — prompt_hash will be empty")
        poll_batch(args.batch_id, api_key, args.poll_interval)
        download_results(args.batch_id, api_key, metadata_rows, out_path)
        return

    # ── Normal mode ───────────────────────────────────────────────────────────
    batch_requests, metadata_rows = build_requests(model_path, csv_path, mnx_dir)

    # Always save requests JSON (for dry-run inspection and resume capability)
    req_path.parent.mkdir(parents=True, exist_ok=True)
    with open(req_path, "w", encoding="utf-8") as fh:
        json.dump(
            {"requests": batch_requests, "metadata": metadata_rows},
            fh, indent=2, ensure_ascii=False,
        )
    logger.info(f"Requests JSON written to: {req_path}")

    if args.dry_run:
        logger.info("--dry-run: skipping API submission. Inspect requests at:")
        logger.info(f"  {req_path}")
        # Print first 2 prompts for quality review
        for meta in metadata_rows[:2]:
            print("\n" + "=" * 70)
            print(f"custom_id: {meta['reaction_id']}  candidates: {meta['candidates']}")
            print("-" * 70)
            print(meta["prompt"])
        return

    batch_id = submit_batch(batch_requests, api_key)
    logger.info(f"Batch ID: {batch_id}")
    logger.info(f"To resume later: python {__file__} --batch-id {batch_id} --api-key <key>")

    poll_batch(batch_id, api_key, args.poll_interval)
    download_results(batch_id, api_key, metadata_rows, out_path)


if __name__ == "__main__":
    main()
