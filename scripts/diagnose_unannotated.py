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
"""

import argparse
import csv
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


if __name__ == "__main__":
    main()
