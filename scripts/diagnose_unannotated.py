"""
diagnose_unannotated.py — export all unannotated reactions to CSV for manual review.

A reaction is considered "unannotated" when its annotation dict is empty or
contains only an "sbo" key (i.e. no database cross-references have been assigned).

Output columns:
  reaction_id          model reaction ID
  reaction_name        human-readable name
  type                 exchange | transport | has_GPR | no_GPR
  n_metabolites        number of participating metabolites
  compartments         comma-separated list of unique compartments
  gene_reaction_rule   GPR string (blank if none)
  metabolite_mnxm_ids  metanetx.chemical IDs of all metabolites (comma-separated)
  metabolite_bigg_ids  bigg.metabolite IDs of all metabolites (comma-separated)

Usage:
  python scripts/diagnose_unannotated.py [--model model.xml] [--out data/unannotated_reactions.csv]
  python scripts/diagnose_unannotated.py --mnx-dir data/metanetx --diag-out data/transport_diagnosis.txt
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def _get_annotation_value(ann: dict, key: str) -> str:
    """Return the first annotation value for key, or ''."""
    if not ann:
        return ""
    raw = ann.get(key, "")
    if isinstance(raw, list):
        return raw[0] if raw else ""
    return str(raw) if raw else ""


def _is_unannotated(rxn) -> bool:
    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    meaningful = {k: v for k, v in ann.items() if k != "sbo"}
    return not meaningful


def _classify(rxn, exchanges: set, demands: set) -> str:
    if rxn in exchanges or len(rxn.metabolites) == 1:
        return "exchange"
    comps = {m.compartment for m in rxn.metabolites}
    if len(comps) >= 2:
        return "transport"
    if rxn.gene_reaction_rule and rxn.gene_reaction_rule.strip():
        return "has_GPR"
    return "no_GPR"


def diagnose_transport(
    rows: list[dict],
    mnx_dir: Path,
    diag_out: Path,
) -> None:
    """
    For the first 10 unannotated transport reactions whose name contains "transport",
    probe desc_index with three matching strategies and write results to diag_out
    (and stdout).
    """
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
        from gem_annotate.io import load_reac_xref
    except ImportError as e:
        print(f"WARNING: could not import gem_annotate.io ({e}) — skipping transport diagnosis")
        return

    reac_xref_path = mnx_dir / "reac_xref.tsv"
    if not reac_xref_path.exists():
        print(f"WARNING: {reac_xref_path} not found — skipping transport diagnosis")
        return

    reac_xref  = load_reac_xref(reac_xref_path)
    desc_index = reac_xref["desc_index"]

    # Pre-build normalised desc_index (strip all non-alphanumeric chars)
    norm_desc_index: dict[str, str] = {}
    for k, v in desc_index.items():
        nk = re.sub(r"[^a-z0-9]", "", k)
        if nk not in norm_desc_index:
            norm_desc_index[nk] = v

    # Select up to 10 transport rows whose name contains "transport"
    targets = [
        r for r in rows
        if r["type"] == "transport" and "transport" in r["reaction_name"].lower()
    ][:10]

    if not targets:
        print("\n[transport diagnosis] No matching reactions found.")
        return

    lines: list[str] = []

    def emit(text: str = "") -> None:
        lines.append(text)
        print(text)

    emit("=" * 70)
    emit("TRANSPORT REACTION DIAGNOSIS")
    emit(f"desc_index size: {len(desc_index):,}  |  targets: {len(targets)}")
    emit("=" * 70)

    for row in targets:
        rxn_id   = row["reaction_id"]
        rxn_name = row["reaction_name"]
        emit()
        emit(f"Reaction : {rxn_id}")
        emit(f"Name     : {rxn_name}")

        name_lower = rxn_name.lower().strip()

        # ── Probe 1: exact name match ────────────────────────────────────────
        exact = desc_index.get(name_lower)
        emit(f"  [B1] exact match on name_lower={name_lower!r}")
        emit(f"       → {exact or '(no match)'}")

        # ── Probe 2: normalised name match ───────────────────────────────────
        norm_name = re.sub(r"[^a-z0-9]", "", name_lower)
        norm_hit  = norm_desc_index.get(norm_name)
        emit(f"  [B2] normalised name={norm_name!r}")
        emit(f"       → {norm_hit or '(no match)'}")

        # ── Extract met_name (text before "transport") ───────────────────────
        met_name = name_lower.split("transport")[0].strip().rstrip("-").strip()
        emit(f"  met_name extracted: {met_name!r}")

        # ── Probe 3: desc_index keys that contain met_name ───────────────────
        if met_name:
            containing = [
                (k, v) for k, v in desc_index.items() if met_name in k
            ][:5]
            emit(f"  [B3a] desc_index keys containing {met_name!r} (top 5):")
            if containing:
                for k, v in containing:
                    emit(f"         {k!r:50s} → {v}")
            else:
                emit("         (none)")

            # ── Probe 4: met_name + transport/permease/diffusion ─────────────
            transport_related = [
                (k, v) for k, v in desc_index.items()
                if met_name in k and any(
                    kw in k for kw in ("transport", "permease", "diffusion")
                )
            ][:5]
            emit(f"  [B3b] keys containing {met_name!r} + transport/permease/diffusion (top 5):")
            if transport_related:
                for k, v in transport_related:
                    emit(f"         {k!r:50s} → {v}")
            else:
                emit("         (none)")
        else:
            emit("  [B3a/b] met_name is empty — skipping substring searches")

        emit("-" * 70)

    diag_out.parent.mkdir(parents=True, exist_ok=True)
    with open(diag_out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nTransport diagnosis written to: {diag_out}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Export unannotated reactions to CSV")
    parser.add_argument(
        "--model",
        default="model.xml",
        help="Path to the SBML model file (default: model.xml)",
    )
    parser.add_argument(
        "--out",
        default="data/unannotated_reactions.csv",
        help="Output CSV path (default: data/unannotated_reactions.csv)",
    )
    parser.add_argument(
        "--mnx-dir",
        default="data/metanetx",
        help="MetaNetX data directory containing reac_xref.tsv (default: data/metanetx)",
    )
    parser.add_argument(
        "--diag-out",
        default="data/transport_diagnosis.txt",
        help="Transport diagnosis output path (default: data/transport_diagnosis.txt)",
    )
    args = parser.parse_args()

    model_path = Path(args.model)
    if not model_path.exists():
        print(f"ERROR: model file not found: {model_path}", file=sys.stderr)
        sys.exit(1)

    try:
        from cobra.io import read_sbml_model
    except ImportError:
        print("ERROR: COBRApy is not installed (pip install cobra)", file=sys.stderr)
        sys.exit(1)

    print(f"Loading model: {model_path}")
    model = read_sbml_model(str(model_path))
    print(f"  {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

    exchanges = set(model.exchanges)
    demands   = set(model.demands)

    rows = []
    for rxn in model.reactions:
        if not _is_unannotated(rxn):
            continue

        mets = list(rxn.metabolites)
        comps = sorted({m.compartment for m in mets})

        mnxm_ids = []
        bigg_ids = []
        for m in mets:
            ann = m.annotation if isinstance(m.annotation, dict) else {}
            mnxm = _get_annotation_value(ann, "metanetx.chemical")
            bigg = _get_annotation_value(ann, "bigg.metabolite")
            mnxm_ids.append(mnxm if mnxm else "—")
            bigg_ids.append(bigg if bigg else "—")

        rows.append({
            "reaction_id":         rxn.id,
            "reaction_name":       rxn.name or "",
            "type":                _classify(rxn, exchanges, demands),
            "n_metabolites":       len(mets),
            "compartments":        ",".join(comps),
            "gene_reaction_rule":  rxn.gene_reaction_rule or "",
            "metabolite_mnxm_ids": ",".join(mnxm_ids),
            "metabolite_bigg_ids": ",".join(bigg_ids),
        })

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "reaction_id", "reaction_name", "type", "n_metabolites",
        "compartments", "gene_reaction_rule",
        "metabolite_mnxm_ids", "metabolite_bigg_ids",
    ]

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Summary by type
    by_type: dict[str, int] = {}
    for row in rows:
        by_type[row["type"]] = by_type.get(row["type"], 0) + 1

    print(f"\nUnannotated reactions: {len(rows)} total")
    for t, n in sorted(by_type.items()):
        print(f"  {t:12s}: {n}")
    print(f"\nWritten to: {out_path}")

    # Transport diagnosis (requires MetaNetX reac_xref.tsv)
    diagnose_transport(rows, Path(args.mnx_dir), Path(args.diag_out))


if __name__ == "__main__":
    main()
