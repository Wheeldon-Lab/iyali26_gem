"""Apply reviewed batch annotations to a new deterministic SBML output.

The previous in-place implementation remains in the external legacy archive.
This safe version never replaces canonical ``model.xml`` or an existing result.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pandas as pd
from cobra.io import read_sbml_model

from scripts.gem_annotate.config import REPO_ROOT
from scripts.gem_annotate.sbml import write_deterministic_sbml_model


def _load_by_mnxr(path: Path) -> dict[str, list[tuple[str, str]]]:
    table = pd.read_csv(path, sep="\t", comment="#", header=None, names=["source", "mnxr", "description"], dtype=str).fillna("")
    mapping: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for source, mnxr in table.loc[table["mnxr"].str.startswith("MNXR"), ["source", "mnxr"]].itertuples(index=False):
        if ":" in source:
            mapping[mnxr].append(tuple(source.split(":", 1)))
    return dict(mapping)


def apply_batch_annotations(model, rows: list[dict[str, str]], mapping: dict[str, list[tuple[str, str]]]) -> int:
    reactions = {reaction.id: reaction for reaction in model.reactions}
    applied = 0
    for row in rows:
        mnxr = row.get("matched_mnxr", "")
        reaction = reactions.get(row.get("reaction_id", ""))
        if reaction is None or not mnxr.startswith("MNXR"):
            continue
        annotation = dict(reaction.annotation) if isinstance(reaction.annotation, dict) else {}
        if {key: value for key, value in annotation.items() if key not in {"sbo", "annotation_source"}}:
            continue
        for prefix, identifier in mapping.get(mnxr, []):
            annotation.setdefault(prefix, [identifier])
        annotation["metanetx.reaction"] = [mnxr]
        annotation["annotation_source"] = ["claude-batch"]
        reaction.annotation = annotation
        applied += 1
    return applied


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Apply reviewed batch annotations safely")
    parser.add_argument("--model", type=Path, default=REPO_ROOT / "model.xml")
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--mnx-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    model_path, output_path = args.model.resolve(), args.out.resolve()
    xref_path = args.mnx_dir.resolve() / "reac_xref.tsv"
    for path in (model_path, args.csv.resolve(), xref_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")
    if not args.dry_run:
        if output_path in {model_path, (REPO_ROOT / "model.xml").resolve()}:
            raise ValueError("Batch annotations require a distinct non-canonical --out path")
        if output_path.exists():
            raise FileExistsError(f"Refusing to overwrite output: {output_path}")
    with args.csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    model = read_sbml_model(str(model_path))
    applied = apply_batch_annotations(model, rows, _load_by_mnxr(xref_path))
    if not args.dry_run:
        write_deterministic_sbml_model(model, output_path)
    print(f"applied={applied} dry_run={args.dry_run}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
