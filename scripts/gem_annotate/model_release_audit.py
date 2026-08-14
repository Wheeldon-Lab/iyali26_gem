"""Deterministic SBML release diff and metabolic-bypass safety gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any

from cobra.io import read_sbml_model


SAFE_REACTION_STOICHIOMETRY = {"R2041"}
SAFE_REACTION_BOUNDS = {"R2041"}
SAFE_REACTION_ANNOTATIONS = {"R815", "R816"}
SAFE_METABOLITE_CHEMISTRY = {"m884[C_cy]"}
PROTECTED_REACTIONS = {"R634", "R_PGAM1_PhosHydro", "R1372"}
REVIEW_REACTIONS = PROTECTED_REACTIONS | {
    "R2041",
    "R643",
    "R815",
    "R816",
}
REVIEW_METABOLITES = {"m884[C_cy]"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalized(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _normalized(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple, set)):
        normalized = [_normalized(item) for item in value]
        return sorted(normalized, key=lambda item: json.dumps(item, sort_keys=True))
    return value


def _reaction_record(reaction) -> dict[str, Any]:
    return {
        "stoichiometry": {
            metabolite.id: float(coefficient)
            for metabolite, coefficient in sorted(
                reaction.metabolites.items(), key=lambda item: item[0].id
            )
        },
        "bounds": [float(reaction.lower_bound), float(reaction.upper_bound)],
        "gpr": str(reaction.gene_reaction_rule),
        "annotation": _normalized(dict(reaction.annotation or {})),
        "notes": _normalized(dict(reaction.notes or {})),
    }


def _metabolite_record(metabolite) -> dict[str, Any]:
    return {
        "formula": metabolite.formula,
        "charge": metabolite.charge,
        "name": metabolite.name,
        "compartment": metabolite.compartment,
        "annotation": _normalized(dict(metabolite.annotation or {})),
        "notes": _normalized(dict(metabolite.notes or {})),
    }


def _changed_ids(
    baseline: dict[str, dict[str, Any]],
    candidate: dict[str, dict[str, Any]],
    field: str,
) -> list[str]:
    return sorted(
        identifier
        for identifier in baseline.keys() & candidate.keys()
        if baseline[identifier][field] != candidate[identifier][field]
    )


def compare_models(
    baseline_path: str | Path,
    candidate_path: str | Path,
    candidate_repeat_path: str | Path | None = None,
) -> dict:
    """Compare two SBML models and evaluate the safe-stage release whitelist."""

    baseline_path = Path(baseline_path)
    candidate_path = Path(candidate_path)
    candidate_repeat_path = (
        Path(candidate_repeat_path) if candidate_repeat_path is not None else None
    )
    baseline_model = read_sbml_model(str(baseline_path))
    candidate_model = read_sbml_model(str(candidate_path))
    baseline_reactions = {
        reaction.id: _reaction_record(reaction)
        for reaction in baseline_model.reactions
    }
    candidate_reactions = {
        reaction.id: _reaction_record(reaction)
        for reaction in candidate_model.reactions
    }
    baseline_metabolites = {
        metabolite.id: _metabolite_record(metabolite)
        for metabolite in baseline_model.metabolites
    }
    candidate_metabolites = {
        metabolite.id: _metabolite_record(metabolite)
        for metabolite in candidate_model.metabolites
    }

    reaction_diff = {
        "added": sorted(candidate_reactions.keys() - baseline_reactions.keys()),
        "removed": sorted(baseline_reactions.keys() - candidate_reactions.keys()),
        "stoichiometry_changed": _changed_ids(
            baseline_reactions, candidate_reactions, "stoichiometry"
        ),
        "bounds_changed": _changed_ids(
            baseline_reactions, candidate_reactions, "bounds"
        ),
        "gpr_changed": _changed_ids(baseline_reactions, candidate_reactions, "gpr"),
        "annotation_changed": _changed_ids(
            baseline_reactions, candidate_reactions, "annotation"
        ),
        "notes_changed": _changed_ids(
            baseline_reactions, candidate_reactions, "notes"
        ),
    }
    metabolite_diff = {
        "added": sorted(candidate_metabolites.keys() - baseline_metabolites.keys()),
        "removed": sorted(baseline_metabolites.keys() - candidate_metabolites.keys()),
        "formula_changed": _changed_ids(
            baseline_metabolites, candidate_metabolites, "formula"
        ),
        "charge_changed": _changed_ids(
            baseline_metabolites, candidate_metabolites, "charge"
        ),
        "name_changed": _changed_ids(
            baseline_metabolites, candidate_metabolites, "name"
        ),
        "compartment_changed": _changed_ids(
            baseline_metabolites, candidate_metabolites, "compartment"
        ),
        "annotation_changed": _changed_ids(
            baseline_metabolites, candidate_metabolites, "annotation"
        ),
        "notes_changed": _changed_ids(
            baseline_metabolites, candidate_metabolites, "notes"
        ),
    }

    violations: list[str] = []
    for field in ("added", "removed", "gpr_changed"):
        if reaction_diff[field]:
            violations.append(f"reaction.{field}")
    for field in ("added", "removed", "name_changed", "compartment_changed"):
        if metabolite_diff[field]:
            violations.append(f"metabolite.{field}")
    whitelist_checks = (
        ("reaction.stoichiometry_changed", reaction_diff["stoichiometry_changed"], SAFE_REACTION_STOICHIOMETRY),
        ("reaction.bounds_changed", reaction_diff["bounds_changed"], SAFE_REACTION_BOUNDS),
        ("reaction.annotation_changed", reaction_diff["annotation_changed"], SAFE_REACTION_ANNOTATIONS),
        ("reaction.notes_changed", reaction_diff["notes_changed"], SAFE_REACTION_STOICHIOMETRY),
        ("metabolite.formula_changed", metabolite_diff["formula_changed"], SAFE_METABOLITE_CHEMISTRY),
        ("metabolite.charge_changed", metabolite_diff["charge_changed"], SAFE_METABOLITE_CHEMISTRY),
        ("metabolite.annotation_changed", metabolite_diff["annotation_changed"], set()),
        ("metabolite.notes_changed", metabolite_diff["notes_changed"], set()),
    )
    for label, changed, allowed in whitelist_checks:
        unexpected = sorted(set(changed) - allowed)
        if unexpected:
            violations.append(f"{label}:unexpected={len(unexpected)}")

    protected = {}
    for reaction_id in sorted(PROTECTED_REACTIONS):
        before = baseline_reactions.get(reaction_id)
        after = candidate_reactions.get(reaction_id)
        protected[reaction_id] = {
            "present_in_both": before is not None and after is not None,
            "stoichiometry_unchanged": (
                before is not None
                and after is not None
                and before["stoichiometry"] == after["stoichiometry"]
            ),
            "bounds_unchanged": (
                before is not None
                and after is not None
                and before["bounds"] == after["bounds"]
            ),
            "gpr_unchanged": (
                before is not None
                and after is not None
                and before["gpr"] == after["gpr"]
            ),
        }
        if not all(protected[reaction_id].values()):
            violations.append(f"protected_reaction:{reaction_id}")

    candidate_sha256 = _sha256(candidate_path)
    repeat_sha256 = (
        _sha256(candidate_repeat_path)
        if candidate_repeat_path is not None
        else None
    )
    reproducible = repeat_sha256 is None or repeat_sha256 == candidate_sha256
    if not reproducible:
        violations.append("candidate_reproducibility")

    return {
        "schema_version": 1,
        "baseline": {
            "path": str(baseline_path.resolve()),
            "sha256": _sha256(baseline_path),
        },
        "candidate": {
            "path": str(candidate_path.resolve()),
            "sha256": candidate_sha256,
        },
        "candidate_reproducibility": {
            "repeat_path": (
                str(candidate_repeat_path.resolve())
                if candidate_repeat_path is not None
                else None
            ),
            "repeat_sha256": repeat_sha256,
            "sha256_identical": reproducible,
        },
        "reaction_diff": reaction_diff,
        "metabolite_diff": metabolite_diff,
        "protected_reactions": protected,
        "review_records": {
            "reactions": {
                reaction_id: {
                    "baseline": baseline_reactions.get(reaction_id),
                    "candidate": candidate_reactions.get(reaction_id),
                }
                for reaction_id in sorted(REVIEW_REACTIONS)
            },
            "metabolites": {
                metabolite_id: {
                    "baseline": baseline_metabolites.get(metabolite_id),
                    "candidate": candidate_metabolites.get(metabolite_id),
                }
                for metabolite_id in sorted(REVIEW_METABOLITES)
            },
        },
        "safe_stage_whitelist": {
            "reaction_stoichiometry": sorted(SAFE_REACTION_STOICHIOMETRY),
            "reaction_bounds": sorted(SAFE_REACTION_BOUNDS),
            "reaction_annotations": sorted(SAFE_REACTION_ANNOTATIONS),
            "metabolite_chemistry": sorted(SAFE_METABOLITE_CHEMISTRY),
        },
        "release_gate": {
            "passed": not violations,
            "violations": sorted(set(violations)),
        },
    }


def write_audit_atomic(report: dict, output_path: str | Path) -> Path:
    """Write an audit without leaving a partial report on failure."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, output_path)
    except Exception:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--candidate-repeat")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = compare_models(
        args.baseline, args.candidate, args.candidate_repeat
    )
    write_audit_atomic(report, args.output)
    print(json.dumps(report["release_gate"], sort_keys=True))
    return 0 if report["release_gate"]["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
