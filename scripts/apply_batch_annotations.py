"""
apply_batch_annotations.py — Apply Claude Batch API results to model.xml.

For each row in data/batch_annotation_results.csv where matched_mnxr starts
with "MNXR":
  - Write metanetx.reaction + cross-references from by_mnxr into the reaction
  - Add annotation_source: ["claude-batch"] to mark LLM-derived annotations

Outputs:
  model.xml (updated in-place by default, or --out for a separate path)
  data/annotation_audit.csv  — one row per reaction, all sources

Usage:
  python scripts/apply_batch_annotations.py \\
      --model model.xml \\
      --csv data/batch_annotation_results.csv \\
      --mnx-dir data/metanetx \\
      --out model.xml

  python scripts/apply_batch_annotations.py --dry-run   # report only, no writes
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ── MetaNetX loader ───────────────────────────────────────────────────────────

def _load_by_mnxr(reac_xref_path: Path) -> dict[str, list[tuple[str, str]]]:
    import pandas as pd
    df = pd.read_csv(
        reac_xref_path, sep="\t", comment="#", header=None,
        names=["source", "mnx_id", "description"], dtype=str, low_memory=False,
    ).fillna("")
    df = df[df["mnx_id"].str.startswith("MNXR")]

    by_mnxr: dict[str, list] = defaultdict(list)
    for source, mnx_id, _ in df.itertuples(index=False):
        if ":" in source:
            prefix, sid = source.split(":", 1)
            by_mnxr[mnx_id].append((prefix, sid))
    return dict(by_mnxr)


# ── Annotation helpers ────────────────────────────────────────────────────────

def _apply_annotation(rxn, mnxr_id: str, by_mnxr: dict) -> None:
    """Merge MNXR cross-refs + annotation_source into rxn.annotation."""
    new_ann: dict[str, list] = defaultdict(list)
    for db_prefix, db_id in by_mnxr.get(mnxr_id, []):
        new_ann[db_prefix].append(db_id)
    new_ann["metanetx.reaction"] = [mnxr_id]
    new_ann["annotation_source"] = ["claude-batch"]

    merged = dict(rxn.annotation)
    for key, val in new_ann.items():
        if key not in merged:
            merged[key] = val
    rxn.annotation = merged


def _get_annotation_value(ann: dict, key: str) -> str:
    if not ann:
        return ""
    raw = ann.get(key, "")
    if isinstance(raw, list):
        return raw[0] if raw else ""
    return str(raw) if raw else ""


def _annotation_source_label(rxn) -> str:
    """Classify annotation source for audit CSV."""
    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    if not {k for k in ann if k != "sbo"}:
        return "unannotated"
    src = ann.get("annotation_source", [])
    if isinstance(src, str):
        src = [src]
    if "claude-batch" in src:
        return "claude-batch"
    if "manual" in src:
        return "manual"
    if "gap-fill" in src:
        return "gap-fill"
    if ann.get("metanetx.reaction"):
        return "metanetx-auto"
    return "unannotated"


# ── Audit CSV ─────────────────────────────────────────────────────────────────

def _write_audit(model, batch_results: dict[str, dict], audit_path: Path) -> None:
    """
    Generate annotation_audit.csv with one row per reaction.

    Columns: reaction_id, annotation_source, mnxr_id, ec_code, timestamp
    """
    timestamp = datetime.now(timezone.utc).isoformat()

    fieldnames = ["reaction_id", "reaction_name", "annotation_source",
                  "mnxr_id", "ec_code", "timestamp"]
    rows: list[dict] = []

    for rxn in sorted(model.reactions, key=lambda r: r.id):
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        source = _annotation_source_label(rxn)

        # For claude-batch rows, prefer the batch result timestamp
        if source == "claude-batch" and rxn.id in batch_results:
            ts = batch_results[rxn.id].get("timestamp", timestamp)
        else:
            ts = timestamp

        mnxr_id = _get_annotation_value(ann, "metanetx.reaction")
        ec_code  = _get_annotation_value(ann, "ec-code")

        rows.append({
            "reaction_id":       rxn.id,
            "reaction_name":     rxn.name or "",
            "annotation_source": source,
            "mnxr_id":           mnxr_id,
            "ec_code":           ec_code,
            "timestamp":         ts,
        })

    audit_path.parent.mkdir(parents=True, exist_ok=True)
    with open(audit_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    by_source: dict[str, int] = {}
    for r in rows:
        by_source[r["annotation_source"]] = by_source.get(r["annotation_source"], 0) + 1
    logger.info(f"Audit written to: {audit_path}")
    for src, n in sorted(by_source.items()):
        logger.info(f"  {src:20s}: {n:5d} reactions")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply Claude batch annotation results to the GEM"
    )
    parser.add_argument("--model",   default="model.xml",
                        help="Input SBML model (default: model.xml)")
    parser.add_argument("--csv",     default="data/batch_annotation_results.csv",
                        help="Batch results CSV (default: data/batch_annotation_results.csv)")
    parser.add_argument("--mnx-dir", default="data/metanetx",
                        help="MetaNetX data directory (default: data/metanetx)")
    parser.add_argument("--out",     default=None,
                        help="Output model path (default: overwrite --model)")
    parser.add_argument("--audit",   default="data/annotation_audit.csv",
                        help="Audit CSV output (default: data/annotation_audit.csv)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would be applied without writing files")
    args = parser.parse_args()

    model_path  = Path(args.model)
    csv_path    = Path(args.csv)
    mnx_dir     = Path(args.mnx_dir)
    out_path    = Path(args.out) if args.out else model_path
    audit_path  = Path(args.audit)
    xref_path   = mnx_dir / "reac_xref.tsv"

    for p, label in [(model_path, "--model"), (csv_path, "--csv"), (xref_path, "reac_xref.tsv")]:
        if not p.exists():
            logger.error(f"Not found ({label}): {p}")
            sys.exit(1)

    try:
        from cobra.io import read_sbml_model, write_sbml_model
    except ImportError:
        logger.error("COBRApy not installed: pip install cobra")
        sys.exit(1)

    # Load batch results
    batch_results: dict[str, dict] = {}
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["matched_mnxr"].startswith("MNXR"):
                batch_results[row["reaction_id"]] = row

    logger.info(f"Batch results: {len(batch_results)} reactions with valid MNXR matches")

    # Load model
    logger.info(f"Loading model: {model_path}")
    model = read_sbml_model(str(model_path))

    # Load MetaNetX cross-refs
    logger.info("Loading reac_xref.tsv …")
    by_mnxr = _load_by_mnxr(xref_path)

    # Build reaction lookup
    rxn_by_id = {rxn.id: rxn for rxn in model.reactions}

    applied = 0
    skipped_missing = 0
    skipped_already = 0

    for rxn_id, result_row in sorted(batch_results.items()):
        rxn = rxn_by_id.get(rxn_id)
        if rxn is None:
            logger.warning(f"  {rxn_id} not found in model — skipping")
            skipped_missing += 1
            continue

        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        meaningful = {k: v for k, v in ann.items() if k not in ("sbo", "annotation_source")}
        if meaningful:
            logger.debug(f"  {rxn_id} already annotated — skipping")
            skipped_already += 1
            continue

        mnxr_id = result_row["matched_mnxr"]
        if args.dry_run:
            logger.info(f"  [dry-run] Would annotate {rxn_id} ({rxn.name!r}) → {mnxr_id}")
        else:
            _apply_annotation(rxn, mnxr_id, by_mnxr)
        applied += 1

    logger.info(
        f"Annotations: applied={applied}  "
        f"already_annotated={skipped_already}  missing_in_model={skipped_missing}"
    )

    if args.dry_run:
        logger.info("--dry-run: model NOT written. Audit will still be generated.")
        _write_audit(model, batch_results, audit_path)
        return

    # Write updated model
    logger.info(f"Writing model to: {out_path}")
    write_sbml_model(model, str(out_path))

    # Write audit CSV
    _write_audit(model, batch_results, audit_path)


if __name__ == "__main__":
    main()
