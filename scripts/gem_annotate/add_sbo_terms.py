"""Create a separate, deterministically written SBO-annotated model.

The historical in-place script is retained in the external legacy archive.
This entry point requires a distinct, new output model to preserve the baseline.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

from cobra.io import read_sbml_model

from .config import REPO_ROOT
from .sbml import write_deterministic_sbml_model

logger = logging.getLogger(__name__)
_BIOMASS_RE = re.compile(r"BIOMASS|biomass|newBiom|R1372")
_MAINTENANCE_RE = re.compile(r"MAINTENANCE|ATPM")


def _set_sbo(obj, term: str) -> bool:
    annotation = obj.annotation if isinstance(obj.annotation, dict) else {}
    if "sbo" in annotation:
        return False
    obj.annotation = {**annotation, "sbo": term}
    return True


def add_sbo_terms(model) -> dict[str, int]:
    counts = {"metabolites": 0, "genes": 0, "reactions": 0}
    for metabolite in model.metabolites:
        counts["metabolites"] += _set_sbo(metabolite, "SBO:0000247")
    for gene in model.genes:
        counts["genes"] += _set_sbo(gene, "SBO:0000243")
    exchanges, demands, sinks = set(model.exchanges), set(model.demands), set(model.sinks)
    for reaction in model.reactions:
        if reaction in exchanges:
            term = "SBO:0000627"
        elif reaction in demands:
            term = "SBO:0000628"
        elif reaction in sinks:
            term = "SBO:0000632"
        elif _BIOMASS_RE.search(reaction.id):
            term = "SBO:0000629"
        elif _MAINTENANCE_RE.search(reaction.id):
            term = "SBO:0000630"
        elif len({metabolite.compartment for metabolite in reaction.metabolites}) >= 2:
            term = "SBO:0000185"
        else:
            term = "SBO:0000176"
        counts["reactions"] += _set_sbo(reaction, term)
    return counts


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add missing SBO terms to a new SBML output")
    parser.add_argument("--model", type=Path, default=REPO_ROOT / "model.xml")
    parser.add_argument("--out", type=Path, required=True)
    return parser


def main() -> int:
    args = _parser().parse_args()
    model_path, output_path = args.model.resolve(), args.out.resolve()
    if not model_path.is_file():
        raise FileNotFoundError(f"Model not found: {model_path}")
    if output_path in {model_path, (REPO_ROOT / "model.xml").resolve()}:
        raise ValueError("SBO annotation requires a distinct non-canonical --out path")
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {output_path}")
    model = read_sbml_model(str(model_path))
    counts = add_sbo_terms(model)
    write_deterministic_sbml_model(model, output_path)
    logger.info("SBO additions: %s", counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
