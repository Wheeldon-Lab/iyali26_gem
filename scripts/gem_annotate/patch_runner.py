"""Safe, deterministic entry point for one curated model patch.

This module is deliberately the only supported way to run an individual
curated patch outside the full ``gem_annotate`` rebuild.  It never writes the
canonical ``model.xml`` and refuses an existing output, so an exploratory
patch cannot silently replace either a baseline or an earlier experiment.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Sequence

from cobra.io import read_sbml_model

from .config import REPO_ROOT
from .essentiality_evidence import sha256_file
from .patches import (
    add_isozyme_gprs,
    annotate_isozyme_genes,
    clean_ec_overload,
    extend_acyl_pool_c161,
)
from .sbml import write_deterministic_sbml_model


PatchFunction = Callable[[object], int]


def _patches(*, allow_network: bool) -> dict[str, PatchFunction]:
    return {
        "c161-pool-extension": extend_acyl_pool_c161,
        "ec-overload-cleanup": clean_ec_overload,
        "isozyme-gprs": add_isozyme_gprs,
        "isozyme-gene-annotations": lambda model: annotate_isozyme_genes(
            model, network=allow_network
        ),
    }


def run_patch(
    patch_name: str,
    *,
    input_model: Path,
    output_model: Path,
    allow_network: bool = False,
) -> dict[str, object]:
    """Apply exactly one curated patch to a new, deterministic SBML output."""

    input_model = Path(input_model).resolve()
    output_model = Path(output_model).resolve()
    canonical_model = (REPO_ROOT / "model.xml").resolve()
    if not input_model.is_file():
        raise FileNotFoundError(f"Input model does not exist: {input_model}")
    if output_model == canonical_model:
        raise ValueError("Patch runner refuses to overwrite canonical model.xml")
    if output_model == input_model:
        raise ValueError("Patch runner requires a distinct output model path")
    if output_model.exists():
        raise FileExistsError(
            f"Patch runner refuses to overwrite an existing output: {output_model}"
        )
    patch = _patches(allow_network=allow_network).get(patch_name)
    if patch is None:
        raise ValueError(f"Unknown curated patch: {patch_name}")

    model = read_sbml_model(str(input_model))
    before = {
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
    }
    changes = patch(model)
    model.optimize()  # a patch output must at least remain solvable by its configured solver
    write_deterministic_sbml_model(model, output_model)
    after = {
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
    }
    audit_path = output_model.with_suffix(output_model.suffix + ".patch-audit.json")
    if audit_path.exists():
        raise FileExistsError(f"Patch audit already exists: {audit_path}")
    audit = {
        "schema_version": 1,
        "patch": patch_name,
        "input_model": {"path": str(input_model), "sha256": sha256_file(input_model)},
        "output_model": {"path": str(output_model), "sha256": sha256_file(output_model)},
        "changes": changes,
        "counts_before": before,
        "counts_after": after,
        "deterministic_writer": "scripts.gem_annotate.sbml.write_deterministic_sbml_model",
    }
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return audit


def _parser(patch_name: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Apply one curated iYali26 patch safely")
    if patch_name is None:
        parser.add_argument("--patch", choices=sorted(_patches(allow_network=False)))
    parser.add_argument("--input-model", type=Path, required=True)
    parser.add_argument("--output-model", type=Path, required=True)
    parser.add_argument(
        "--allow-network",
        action="store_true",
        help="allow the optional UniProt lookup for isozyme gene annotations",
    )
    return parser


def main_for_legacy(patch_name: str, argv: Sequence[str] | None = None) -> int:
    """Compatibility CLI used by the former one-off scripts."""

    args = _parser(patch_name).parse_args(argv)
    audit = run_patch(
        patch_name,
        input_model=args.input_model,
        output_model=args.output_model,
        allow_network=args.allow_network,
    )
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    audit = run_patch(
        args.patch,
        input_model=args.input_model,
        output_model=args.output_model,
        allow_network=args.allow_network,
    )
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
