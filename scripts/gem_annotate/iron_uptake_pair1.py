"""Isolated R2193 Pair-1 GPR experiments.

This module intentionally does not participate in the canonical model build.
It creates a separate SBML artifact, refuses to overwrite ``model.xml``, and
records the evidence boundary in both the reaction notes and a sidecar audit.

Two modes are available:

``ftra_only``
    Replace the inherited, mismatched GPR with ``YALI1D08564g``.  This is the
    only Pair-1 gene whose exact ORF has same-species experimental support.

``pair1_and``
    Use ``YALI1D08564g and YALI1D08684g`` as a sensitivity experiment.  The
    FetC identity and the AND relationship are homology-based and must not be
    described as directly validated in Yarrowia.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from cobra.io import read_sbml_model

from .config import REPO_ROOT
from .essentiality_evidence import sha256_file
from .sbml import write_deterministic_sbml_model


REACTION_ID = "R2193"
INHERITED_GPR = "YALI1F07747g"
FTRA_GENE = "YALI1D08564g"
FETC_GENE = "YALI1D08684g"
EVIDENCE_URL = "https://doi.org/10.1016/j.synbio.2026.02.004"

MODE_GPRS = {
    "ftra_only": FTRA_GENE,
    "pair1_and": f"{FTRA_GENE} and {FETC_GENE}",
}

MODE_EVIDENCE = {
    "ftra_only": {
        "status": "same_species_experimental_support",
        "scope": (
            "The study's YlFtr1 primers identify NCBI Gene 2910500 / "
            "YALI1D08564g. VHb and YlFtr1 were co-expressed, so the study is "
            "not an Ftr1-only causal test."
        ),
        "limitation": (
            "No Ftr1-only control, knockout/complementation, uptake kinetics, "
            "or membrane-localization experiment."
        ),
    },
    "pair1_and": {
        "status": "sensitivity_hypothesis",
        "scope": (
            "YALI1D08564g has same-species experimental support; "
            "YALI1D08684g is a computationally annotated multicopper oxidase."
        ),
        "limitation": (
            "The cited study did not manipulate FetC or test the proposed "
            "FtrA-FetC interaction; the AND relationship is inferred."
        ),
    },
}

GENE_ANNOTATIONS = {
    FTRA_GENE: {
        "sbo": "SBO:0000243",
        "ncbigene": "2910500",
        "kegg.genes": "yli:2910500",
        "uniprot": "Q6CA15",
        "refseq": "XP_502497.1",
    },
    FETC_GENE: {
        "sbo": "SBO:0000243",
        "ncbigene": "2910503",
        "kegg.genes": "yli:2910503",
        "uniprot": "Q6CA12",
        "refseq": "XP_502500.2",
    },
}


def _reaction_signature(reaction) -> dict[str, object]:
    return {
        "stoichiometry": {
            metabolite.id: float(coefficient)
            for metabolite, coefficient in sorted(
                reaction.metabolites.items(), key=lambda item: item[0].id
            )
        },
        "lower_bound": float(reaction.lower_bound),
        "upper_bound": float(reaction.upper_bound),
        "gpr": str(reaction.gene_reaction_rule),
    }


def apply_pair1_overlay(model, *, mode: str = "ftra_only") -> dict[str, object]:
    """Apply one evidence-labelled R2193 GPR overlay in memory.

    The function is idempotent and refuses an unexpected pre-existing GPR.
    Stoichiometry and bounds are guarded before and after the assignment.
    """

    if mode not in MODE_GPRS:
        raise ValueError(
            f"Unsupported Pair-1 mode {mode!r}; choose one of {sorted(MODE_GPRS)}"
        )
    try:
        reaction = model.reactions.get_by_id(REACTION_ID)
    except KeyError as exc:
        raise ValueError(f"Pair-1 overlay requires reaction {REACTION_ID}") from exc

    before = _reaction_signature(reaction)
    target_gpr = MODE_GPRS[mode]
    allowed_existing = {INHERITED_GPR, *MODE_GPRS.values()}
    if before["gpr"] not in allowed_existing:
        raise ValueError(
            f"{REACTION_ID} has unexpected GPR {before['gpr']!r}; "
            "refusing to overwrite unrelated curation"
        )

    changed = before["gpr"] != target_gpr
    if changed:
        reaction.gene_reaction_rule = target_gpr

    target_genes = [FTRA_GENE] if mode == "ftra_only" else [FTRA_GENE, FETC_GENE]
    for gene_id in target_genes:
        gene = model.genes.get_by_id(gene_id)
        annotation = dict(gene.annotation) if isinstance(gene.annotation, dict) else {}
        annotation.update(GENE_ANNOTATIONS[gene_id])
        gene.annotation = annotation

    evidence = MODE_EVIDENCE[mode]
    notes = dict(reaction.notes) if isinstance(reaction.notes, dict) else {}
    notes.update(
        {
            "experimental_overlay": "r2193_pair1_v1",
            "experimental_gpr_mode": mode,
            "experimental_evidence_status": evidence["status"],
            "experimental_evidence_scope": evidence["scope"],
            "experimental_evidence_url": EVIDENCE_URL,
            "experimental_limitation": evidence["limitation"],
            "canonical_status": "not_accepted_needs_more_evidence",
        }
    )
    reaction.notes = notes

    after = _reaction_signature(reaction)
    for field in ("stoichiometry", "lower_bound", "upper_bound"):
        if before[field] != after[field]:
            raise RuntimeError(f"Pair-1 overlay changed protected R2193 field {field}")

    old_gene_orphan = False
    try:
        old_gene_orphan = not bool(model.genes.get_by_id(INHERITED_GPR).reactions)
    except KeyError:
        pass

    return {
        "reaction_id": REACTION_ID,
        "mode": mode,
        "changed": changed,
        "before": before,
        "after": after,
        "target_genes": target_genes,
        "gene_annotations": {
            gene_id: dict(model.genes.get_by_id(gene_id).annotation)
            for gene_id in target_genes
        },
        "inherited_gene_retained_as_orphan": old_gene_orphan,
        "evidence_url": EVIDENCE_URL,
        "evidence_status": evidence["status"],
        "evidence_scope": evidence["scope"],
        "limitation": evidence["limitation"],
        "canonical_model_changed": False,
    }


def run_pair1_overlay(
    *,
    input_model: Path,
    output_model: Path,
    mode: str = "ftra_only",
) -> dict[str, object]:
    """Write a deterministic, non-canonical Pair-1 experimental model."""

    input_model = Path(input_model).resolve()
    output_model = Path(output_model).resolve()
    canonical_model = (REPO_ROOT / "model.xml").resolve()
    if not input_model.is_file():
        raise FileNotFoundError(f"Input model does not exist: {input_model}")
    if output_model == canonical_model:
        raise ValueError("Pair-1 overlay refuses to overwrite canonical model.xml")
    if output_model == input_model:
        raise ValueError("Pair-1 overlay requires a distinct output model")
    if output_model.exists():
        raise FileExistsError(
            f"Pair-1 overlay refuses to overwrite existing output: {output_model}"
        )
    audit_path = output_model.with_suffix(output_model.suffix + ".pair1-audit.json")
    if audit_path.exists():
        raise FileExistsError(f"Pair-1 overlay audit already exists: {audit_path}")

    canonical_sha_before = (
        sha256_file(canonical_model) if canonical_model.is_file() else None
    )
    model = read_sbml_model(str(input_model))
    counts_before = {
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
    }
    overlay = apply_pair1_overlay(model, mode=mode)
    solution = model.optimize()
    if solution.status != "optimal":
        raise RuntimeError(
            f"Pair-1 experimental model is not solvable: {solution.status}"
        )

    output_model.parent.mkdir(parents=True, exist_ok=True)
    write_deterministic_sbml_model(model, output_model)

    roundtrip = read_sbml_model(str(output_model))
    roundtrip_reaction = roundtrip.reactions.get_by_id(REACTION_ID)
    if roundtrip_reaction.gene_reaction_rule != MODE_GPRS[mode]:
        raise RuntimeError("Pair-1 GPR did not survive deterministic SBML roundtrip")
    if (
        _reaction_signature(roundtrip_reaction)["stoichiometry"]
        != overlay["before"]["stoichiometry"]
    ):
        raise RuntimeError("R2193 stoichiometry changed during SBML roundtrip")

    canonical_sha_after = (
        sha256_file(canonical_model) if canonical_model.is_file() else None
    )
    if canonical_sha_after != canonical_sha_before:
        raise RuntimeError("canonical model.xml changed during Pair-1 experiment")

    audit = {
        "schema_version": 1,
        "workflow": "r2193_pair1_experimental_overlay",
        "mode": mode,
        "input_model": {
            "path": str(input_model),
            "sha256": sha256_file(input_model),
        },
        "output_model": {
            "path": str(output_model),
            "sha256": sha256_file(output_model),
        },
        "canonical_model": {
            "path": str(canonical_model),
            "sha256_before": canonical_sha_before,
            "sha256_after": canonical_sha_after,
            "unchanged": canonical_sha_before == canonical_sha_after,
        },
        "counts_before": counts_before,
        "counts_after": {
            "reactions": len(roundtrip.reactions),
            "metabolites": len(roundtrip.metabolites),
            "genes": len(roundtrip.genes),
        },
        "default_objective_value": float(solution.objective_value or 0.0),
        "overlay": overlay,
    }
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return audit


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create an isolated R2193 Pair-1 GPR experimental model"
    )
    parser.add_argument("--input-model", type=Path, required=True)
    parser.add_argument("--output-model", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=sorted(MODE_GPRS),
        default="ftra_only",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    audit = run_pair1_overlay(
        input_model=args.input_model,
        output_model=args.output_model,
        mode=args.mode,
    )
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
