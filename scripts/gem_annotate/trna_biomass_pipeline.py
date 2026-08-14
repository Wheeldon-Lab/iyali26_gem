"""Reproducible build-and-screen pipeline for the experimental B-group model.

The B-group profile replaces the 20 free cytosolic amino-acid biomass drains
with 20 independent, carrier-conserving AA-tRNA -> tRNA + protein-residue
reactions.  It starts from the released canonical ``model.xml`` so the
experiment is an isolated overlay, rather than an unreviewed full rebuild from
the raw source model.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys


WORKFLOW = "trna_biomass_group_b"
MODEL_NAME = "model.trna_biomass_group_b.xml"
SCREEN_DIRECTORY = "essentiality_screen"
MANIFEST_NAME = "group_b_pipeline_manifest.json"
SPLIT_REACTION_PREFIX = "TRNA_BIOMASS_"
SPLIT_RESIDUE_PREFIX = "trna_biomass_residue_"
EXPECTED_PO1F_SCREEN = {
    0.01: (57, 265),
    0.05: (63, 259),
    0.10: (67, 255),
    0.15: (79, 243),
}


def _utc_compact() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _reaction_signature(reaction) -> tuple[tuple[tuple[str, float], ...], tuple[float, float], str]:
    stoichiometry = tuple(
        sorted(
            (metabolite.id, float(coefficient))
            for metabolite, coefficient in reaction.metabolites.items()
        )
    )
    return (
        stoichiometry,
        (float(reaction.lower_bound), float(reaction.upper_bound)),
        str(reaction.gene_reaction_rule),
    )


def validate_group_b_overlay(canonical_model_path: Path, group_b_model_path: Path) -> dict[str, object]:
    """Verify that the experimental model differs only by the B-group overlay."""

    from cobra.io import read_sbml_model

    canonical = read_sbml_model(str(canonical_model_path))
    group_b = read_sbml_model(str(group_b_model_path))

    canonical_reactions = {reaction.id: reaction for reaction in canonical.reactions}
    group_b_reactions = {reaction.id: reaction for reaction in group_b.reactions}
    canonical_metabolites = {metabolite.id for metabolite in canonical.metabolites}
    group_b_metabolites = {metabolite.id for metabolite in group_b.metabolites}

    added_reactions = sorted(set(group_b_reactions) - set(canonical_reactions))
    added_metabolites = sorted(group_b_metabolites - canonical_metabolites)
    removed_reactions = sorted(set(canonical_reactions) - set(group_b_reactions))
    removed_metabolites = sorted(canonical_metabolites - group_b_metabolites)
    canonical_is_b_group = (
        canonical_reactions["biomass_C"].notes.get("experimental_trna_biomass_mode")
        == "split_v1"
    )

    if canonical_is_b_group:
        if added_reactions or added_metabolites or removed_reactions or removed_metabolites:
            raise ValueError(
                "Canonical B-group copy changed model objects: "
                f"added_reactions={added_reactions}, added_metabolites={added_metabolites}, "
                f"removed_reactions={removed_reactions}, removed_metabolites={removed_metabolites}"
            )
    else:
        if removed_reactions or removed_metabolites:
            raise ValueError(
                "B-group overlay removed canonical objects: "
                f"reactions={removed_reactions}, metabolites={removed_metabolites}"
            )
        if len(added_reactions) != 20 or any(
            not reaction_id.startswith(SPLIT_REACTION_PREFIX)
            for reaction_id in added_reactions
        ):
            raise ValueError(
                "B-group overlay must add exactly 20 split reactions; found "
                f"{added_reactions}"
            )
        if len(added_metabolites) != 20 or any(
            not metabolite_id.startswith(SPLIT_RESIDUE_PREFIX)
            for metabolite_id in added_metabolites
        ):
            raise ValueError(
                "B-group overlay must add exactly 20 private protein residues; found "
                f"{added_metabolites}"
            )

    changed_preexisting: list[str] = []
    for reaction_id, canonical_reaction in canonical_reactions.items():
        if _reaction_signature(canonical_reaction) != _reaction_signature(
            group_b_reactions[reaction_id]
        ):
            changed_preexisting.append(reaction_id)
    expected_changed = [] if canonical_is_b_group else ["biomass_C"]
    if changed_preexisting != expected_changed:
        if not canonical_is_b_group:
            raise ValueError(
                "Only biomass_C may change among pre-existing reactions; found "
                f"{changed_preexisting}"
            )
        raise ValueError(
            "Unexpected changed pre-existing reactions; expected "
            f"{expected_changed}, found {changed_preexisting}"
        )

    biomass = group_b_reactions["biomass_C"]
    if biomass.notes.get("experimental_trna_biomass_mode") != "split_v1":
        raise ValueError("B-group biomass mode metadata is missing or stale")
    objective_reactions = {
        reaction.id
        for reaction in group_b.reactions
        if abs(float(reaction.objective_coefficient)) > 0
    }
    if objective_reactions != {"biomass_C"}:
        raise ValueError(
            f"B-group objective must remain biomass_C; found {sorted(objective_reactions)}"
        )

    for reaction_id in added_reactions:
        reaction = group_b_reactions[reaction_id]
        coefficients = sorted(float(value) for value in reaction.metabolites.values())
        if coefficients != [-1.0, 1.0, 1.0]:
            raise ValueError(
                f"{reaction_id} does not conserve one tRNA carrier: "
                f"{reaction.reaction}"
            )
        if reaction.bounds != (0.0, 1000.0):
            raise ValueError(
                f"{reaction_id} has unexpected bounds {reaction.bounds}"
            )

    return {
        "canonical_reactions": len(canonical.reactions),
        "canonical_metabolites": len(canonical.metabolites),
        "group_b_reactions": len(group_b.reactions),
        "group_b_metabolites": len(group_b.metabolites),
        "added_reactions": added_reactions,
        "added_metabolites": added_metabolites,
        "changed_preexisting_reactions": changed_preexisting,
        "canonical_b_group": canonical_is_b_group,
        "objective_reaction": "biomass_C",
        "carrier_conserving_split_reactions": 20,
        "experimental_mode": "split_v1",
    }


def build_group_b_overlay(canonical_model_path: Path, output_model_path: Path) -> list[dict[str, object]]:
    """Copy the released model and apply only the carrier-conserving B overlay."""

    from cobra.io import read_sbml_model

    from .patches import split_trna_charging_from_biomass
    from .sbml import write_deterministic_sbml_model

    model = read_sbml_model(str(canonical_model_path))
    split_audit = split_trna_charging_from_biomass(model)
    if len(split_audit) not in {0, 20}:
        raise RuntimeError(
            "B-group overlay must create 20 reactions or validate a complete "
            "canonical B-group state; found "
            f"{len(split_audit)}"
        )
    write_deterministic_sbml_model(model, output_model_path)
    return split_audit


def restore_free_amino_acid_biomass(model) -> list[dict[str, object]]:
    """Exactly invert split_v1 in an in-memory model copy for A/B attribution."""
    biomass = model.reactions.get_by_id("biomass_C")
    if biomass.notes.get("experimental_trna_biomass_mode") != "split_v1":
        raise ValueError("A counterfactual requires a complete split_v1 B-group biomass")
    split_reactions = sorted(
        [reaction for reaction in model.reactions if reaction.id.startswith(SPLIT_REACTION_PREFIX)],
        key=lambda reaction: reaction.id,
    )
    if len(split_reactions) != 20:
        raise ValueError(
            f"A counterfactual requires exactly 20 B-group split reactions; found {len(split_reactions)}"
        )
    residues = []
    audit: list[dict[str, object]] = []
    for reaction in split_reactions:
        notes = reaction.notes
        required = {
            "amino_acid_id",
            "charged_trna_id",
            "uncharged_trna_id",
            "protein_residue_id",
            "biomass_coefficient",
        }
        missing = sorted(required - set(notes))
        if missing:
            raise ValueError(f"{reaction.id} has incomplete B-group notes: {missing}")
        amount = float(notes["biomass_coefficient"])
        amino_acid = model.metabolites.get_by_id(str(notes["amino_acid_id"]))
        residue = model.metabolites.get_by_id(str(notes["protein_residue_id"]))
        charged = model.metabolites.get_by_id(str(notes["charged_trna_id"]))
        uncharged = model.metabolites.get_by_id(str(notes["uncharged_trna_id"]))
        if dict(reaction.metabolites) != {charged: -1.0, uncharged: 1.0, residue: 1.0}:
            raise ValueError(f"{reaction.id} is not a carrier-conserving B-group split")
        if abs(float(biomass.metabolites.get(residue, 0.0)) + amount) > 1e-12:
            raise ValueError(f"biomass_C does not consume expected residue {residue.id}")
        biomass.add_metabolites({amino_acid: -amount, residue: amount})
        residues.append(residue)
        audit.append(
            {
                "split_reaction_id": reaction.id,
                "amino_acid_id": amino_acid.id,
                "protein_residue_id": residue.id,
                "biomass_coefficient": amount,
            }
        )
    model.remove_reactions(split_reactions, remove_orphans=False)
    model.remove_metabolites(residues, destructive=False)
    for key in (
        "canonical_trna_biomass_mode",
        "experimental_trna_biomass_mode",
        "experimental_trna_biomass_template",
        "experimental_trna_biomass_design",
    ):
        biomass.notes.pop(key, None)
    return audit


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build the fully split B-group AA-tRNA biomass overlay from the "
            "released canonical model, then run the positive-only four-cutoff screen."
        )
    )
    parser.add_argument(
        "--research-root",
        type=Path,
        help="External iYali26 research workspace",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "New result directory. Defaults to a provenance-keyed directory "
            "under artifacts/results/essentiality/trna_biomass_group_b/pipeline_runs."
        ),
    )
    parser.add_argument("--solver", default="gurobi")
    parser.add_argument("--force-rerun", action="store_true")
    parser.add_argument("--reproduction-reason")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.force_rerun and not args.reproduction_reason:
        parser.error("--force-rerun requires --reproduction-reason")
    if args.research_root is not None:
        os.environ["IYALI26_RESEARCH_ROOT"] = str(args.research_root.resolve())

    # Import after configuring IYALI26_RESEARCH_ROOT.  Validation resolves its
    # external data paths at import time.
    from .config import load_project_paths
    from .essentiality_evidence import sha256_file
    from .run_registry import (
        DuplicateRunError,
        build_run_key,
        guard_duplicate_run,
        register_run,
    )
    from .validate_essential_genes import (
        DEFAULT_CUTOFFS,
        PRIMARY_CUTOFF,
        print_summary,
        validate_essential_genes,
    )
    from .strain_overlay import load_strain_profile

    try:
        paths = load_project_paths(args.research_root, required=True)
        paths.require(
            paths.output_model,
            paths.essentiality,
            paths.media,
            paths.strain_profiles,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        parser.error(str(exc))

    experimental_path = paths.essentiality / "consensus_essential_genes.csv"
    media_path = paths.media / "sd_leu.csv"
    strain_profile_path = paths.strain_profiles / "po1f_sd_leu.json"
    for path, label in (
        (experimental_path, "experimental positive list"),
        (media_path, "SD-Leu medium"),
        (strain_profile_path, "PO1f runtime strain profile"),
    ):
        if not path.is_file():
            parser.error(f"{label} not found: {path}")
    strain_profile = load_strain_profile(strain_profile_path)

    code_paths = (
        Path(__file__),
        Path(__file__).with_name("patches.py"),
        Path(__file__).with_name("sbml.py"),
        Path(__file__).with_name("strain_overlay.py"),
        Path(__file__).with_name("validate_essential_genes.py"),
    )
    code_sources = {
        str(path.resolve()): sha256_file(path)
        for path in code_paths
    }
    inputs = {
        "canonical_model": {
            "path": str(paths.output_model.resolve()),
            "sha256": sha256_file(paths.output_model),
        },
        "experimental": {
            "path": str(experimental_path.resolve()),
            "sha256": sha256_file(experimental_path),
        },
        "medium": {
            "path": str(media_path.resolve()),
            "sha256": sha256_file(media_path),
        },
        "strain_profile": {
            "path": str(strain_profile_path.resolve()),
            "sha256": sha256_file(strain_profile_path),
        },
    }
    workspace_manifest = paths.research_root / "manifest.json"
    if workspace_manifest.is_file():
        inputs["research_workspace_manifest"] = {
            "path": str(workspace_manifest.resolve()),
            "sha256": sha256_file(workspace_manifest),
        }
    configuration = {
        "trna_biomass_mode": "split",
        "overlay_base": "released_canonical_model",
        "positive_only": True,
        "primary_cutoff": PRIMARY_CUTOFF,
        "growth_cutoffs": list(DEFAULT_CUTOFFS),
        "solver": args.solver,
        "strain_profile_id": str(strain_profile["profile_id"]),
        "canonical_model_mutation_allowed": False,
    }
    run_key = build_run_key(
        WORKFLOW,
        inputs=inputs,
        code_sources=code_sources,
        configuration=configuration,
    )

    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else (
            paths.results
            / "essentiality"
            / "trna_biomass_group_b"
            / "pipeline_runs"
            / f"{run_key[:12]}-{_utc_compact()}"
        ).resolve()
    )
    if output_dir.exists():
        parser.error(f"output directory already exists: {output_dir}")
    try:
        previous = guard_duplicate_run(
            paths.research_root,
            workflow=WORKFLOW,
            run_key=run_key,
            output_dir=output_dir,
            force_rerun=args.force_rerun,
            reproduction_reason=args.reproduction_reason,
        )
    except (DuplicateRunError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))

    canonical_sha_before = sha256_file(paths.output_model)
    output_dir.mkdir(parents=True)
    model_path = output_dir / MODEL_NAME
    screen_dir = output_dir / SCREEN_DIRECTORY
    manifest_path = output_dir / MANIFEST_NAME

    try:
        build_group_b_overlay(paths.output_model, model_path)
        canonical_sha_after_build = sha256_file(paths.output_model)
        if canonical_sha_after_build != canonical_sha_before:
            raise RuntimeError(
                "canonical model.xml changed during the B-group experimental build"
            )

        overlay_audit = validate_group_b_overlay(paths.output_model, model_path)
        summary = validate_essential_genes(
            experimental_path=experimental_path,
            model_path=model_path,
            media_path=media_path,
            output_dir=screen_dir,
            primary_cutoff=PRIMARY_CUTOFF,
            growth_cutoffs=DEFAULT_CUTOFFS,
            positive_only=True,
            solver=args.solver,
            strain_profile_path=strain_profile_path,
            run_key=run_key,
            code_sources=code_sources,
        )
        canonical_sha_after_screen = sha256_file(paths.output_model)
        if canonical_sha_after_screen != canonical_sha_before:
            raise RuntimeError(
                "canonical model.xml changed during the B-group screen"
            )
        if not summary["medium"].get("iron_uptake_open"):
            raise RuntimeError("B-group screen did not apply the iron-replete SD-Leu medium")
        if not summary["strain_overlay"].get("enabled"):
            raise RuntimeError("B-group screen did not apply the PO1f runtime overlay")
        if summary["medium"].get("uracil_uptake_bound") != 1000.0:
            raise RuntimeError(
                "PO1f screen requires permissive uracil availability "
                "(R1354 uptake bound 1000)"
            )
        observed_curve = {
            float(row["cutoff_fraction_of_wt"]): (
                int(row["TP"]),
                int(row["FN"]),
            )
            for row in summary["cutoff_curve"]
        }
        if observed_curve != EXPECTED_PO1F_SCREEN:
            raise RuntimeError(
                "PO1f+B-group four-cutoff regression changed: expected "
                f"{EXPECTED_PO1F_SCREEN}, found {observed_curve}"
            )

        manifest = {
            "schema_version": "1.1",
            "workflow": WORKFLOW,
            "status": "complete",
            "run_key": run_key,
            "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "inputs": inputs,
            "code_sources": code_sources,
            "configuration": configuration,
            "invariants": {
                "canonical_model_sha256_before": canonical_sha_before,
                "canonical_model_sha256_after": canonical_sha_after_screen,
                "canonical_model_unchanged": True,
                **overlay_audit,
            },
            "outputs": {
                "group_b_model": {
                    "path": str(model_path.resolve()),
                    "sha256": sha256_file(model_path),
                },
                "essentiality_summary": str(
                    (screen_dir / "essentiality_summary.json").resolve()
                ),
                "essentiality_per_gene": str(
                    (screen_dir / "essentiality_per_gene.tsv").resolve()
                ),
                "essentiality_run_manifest": str(
                    (screen_dir / "run_manifest.json").resolve()
                ),
            },
            "screen": {
                "wt_growth": summary["wt_growth"],
                "primary": summary["primary"],
                "cutoff_curve": summary["cutoff_curve"],
                "medium": summary["medium"],
                "strain_overlay": summary["strain_overlay"],
                "simulation_context": summary["simulation_context"],
            },
        }
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        register_run(
            paths.research_root,
            workflow=WORKFLOW,
            run_key=run_key,
            output_dir=output_dir,
            inputs=inputs,
            code_sources=code_sources,
            configuration=configuration,
            status="complete",
            manifest_path=manifest_path,
            previous=previous,
            reproduction_reason=args.reproduction_reason,
        )
    except Exception as exc:
        failure = {
            "schema_version": "1.1",
            "workflow": WORKFLOW,
            "status": "failed",
            "run_key": run_key,
            "error": f"{type(exc).__name__}: {exc}",
            "canonical_model_sha256_before": canonical_sha_before,
            "canonical_model_sha256_after": (
                sha256_file(paths.output_model)
                if paths.output_model.is_file()
                else None
            ),
        }
        manifest_path.write_text(
            json.dumps(failure, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        raise

    print_summary(summary)
    print(f"\nB-group pipeline manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
