"""Optional, globally active isozyme-capacity hypotheses.

The canonical evidence-gated essentiality patch table is intentionally not
used here.  This module builds a separately named experimental model profile
in which a shared OR reaction is partitioned into a primary reaction and a
capacity-limited backup reaction.  Every affected gene and reaction is marked
for later replacement with condition-matched protein abundance and kcat data.
"""

from __future__ import annotations

import ast
import copy
import csv
import math
from pathlib import Path
from typing import Any

from cobra import Reaction
from cobra.core.gene import GPR

from .essentiality_evidence import sha256_file, target_fingerprint


REQUIRED_COLUMNS = {
    "capacity_id",
    "status",
    "source_reaction_id",
    "expected_gpr",
    "primary_gpr",
    "backup_gpr",
    "backup_reaction_id",
    "provisional_upper_bound",
    "units",
    "parameter_basis",
    "case_id",
    "validated_model_sha256",
    "target_fingerprint",
    "requires_protein_abundance",
    "requires_kcat",
    "replacement_formula",
    "rationale",
}
ACTIVE_STATUS = "active_exploratory"
MARKER_STATUS = "provisional_capacity_requires_proteomics_and_kcat"


def _parse_bool(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Unrecognised Boolean value: {value!r}")


def _reaction_context(reaction) -> dict[str, Any]:
    return {
        "reaction_id": reaction.id,
        "stoichiometry": {
            metabolite.id: float(coefficient)
            for metabolite, coefficient in reaction.metabolites.items()
        },
        "lower_bound": float(reaction.lower_bound),
        "upper_bound": float(reaction.upper_bound),
        "gpr": reaction.gene_reaction_rule,
    }


def _genes(rule: str) -> set[str]:
    if not rule.strip():
        return set()
    return set(GPR.from_string(rule).genes)


def _is_or_only_rule(rule: str) -> bool:
    """Return whether a GPR contains only genes joined by OR operators."""
    body = GPR.from_string(rule).body

    def allowed(node: ast.AST | None) -> bool:
        if isinstance(node, ast.Name):
            return True
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.Or):
            return bool(node.values) and all(allowed(value) for value in node.values)
        return False

    return allowed(body)


def _append_marker(notes: dict[str, Any], key: str, value: str) -> None:
    existing = {
        item.strip()
        for item in str(notes.get(key, "")).split(";")
        if item.strip()
    }
    existing.add(value)
    notes[key] = ";".join(sorted(existing))


def load_provisional_capacity_table(path: Path) -> list[dict[str, Any]]:
    """Read and validate the static fields of a provisional capacity profile."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Provisional capacity table is missing columns {sorted(missing)}"
            )
        rows = [dict(row) for row in reader]

    if not rows:
        raise ValueError(f"Provisional capacity table is empty: {path}")

    seen_capacity_ids: set[str] = set()
    seen_source_reactions: set[str] = set()
    seen_backup_reactions: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        for column in REQUIRED_COLUMNS:
            row[column] = str(row[column]).strip()
        if row["status"].lower() != ACTIVE_STATUS:
            raise ValueError(
                f"Provisional capacity row {line_number} must have "
                f"status={ACTIVE_STATUS}"
            )
        for column in (
            "capacity_id",
            "source_reaction_id",
            "expected_gpr",
            "primary_gpr",
            "backup_gpr",
            "backup_reaction_id",
            "units",
            "parameter_basis",
            "case_id",
            "validated_model_sha256",
            "target_fingerprint",
            "replacement_formula",
            "rationale",
        ):
            if not row[column]:
                raise ValueError(
                    f"Provisional capacity row {line_number} has an empty {column}"
                )
        try:
            upper_bound = float(row["provisional_upper_bound"])
        except ValueError as exc:
            raise ValueError(
                f"Invalid provisional upper bound on row {line_number}"
            ) from exc
        if not math.isfinite(upper_bound) or upper_bound <= 0:
            raise ValueError(
                f"Provisional upper bound on row {line_number} must be finite and positive"
            )
        row["provisional_upper_bound"] = upper_bound
        row["requires_protein_abundance"] = _parse_bool(
            row["requires_protein_abundance"]
        )
        row["requires_kcat"] = _parse_bool(row["requires_kcat"])
        if not row["requires_protein_abundance"] or not row["requires_kcat"]:
            raise ValueError(
                "Every active provisional capacity must explicitly require both "
                "protein abundance and kcat replacement"
            )

        capacity_id = row["capacity_id"]
        source_id = row["source_reaction_id"]
        backup_id = row["backup_reaction_id"]
        if capacity_id in seen_capacity_ids:
            raise ValueError(f"Duplicate capacity_id: {capacity_id}")
        if source_id in seen_source_reactions:
            raise ValueError(f"Duplicate source reaction: {source_id}")
        if backup_id in seen_backup_reactions:
            raise ValueError(f"Duplicate backup reaction: {backup_id}")
        if source_id == backup_id:
            raise ValueError(f"Capacity {capacity_id} reuses its source reaction ID")
        seen_capacity_ids.add(capacity_id)
        seen_source_reactions.add(source_id)
        seen_backup_reactions.add(backup_id)
    return rows


def _validate_partition(row: dict[str, Any], source_reaction) -> None:
    capacity_id = row["capacity_id"]
    expected_gpr = row["expected_gpr"]
    primary_gpr = row["primary_gpr"]
    backup_gpr = row["backup_gpr"]
    if not all(
        _is_or_only_rule(rule)
        for rule in (expected_gpr, primary_gpr, backup_gpr)
    ):
        raise ValueError(
            f"Capacity {capacity_id} supports only simple OR isozyme partitions"
        )
    expected_genes = _genes(expected_gpr)
    primary_genes = _genes(primary_gpr)
    backup_genes = _genes(backup_gpr)
    live_genes = {gene.id for gene in source_reaction.genes}
    if not primary_genes or not backup_genes:
        raise ValueError(f"Capacity {capacity_id} has an empty GPR partition")
    if primary_genes & backup_genes:
        raise ValueError(f"Capacity {capacity_id} assigns a gene to both partitions")
    if primary_genes | backup_genes != expected_genes:
        raise ValueError(
            f"Capacity {capacity_id} primary/backup genes do not reconstruct expected_gpr"
        )
    if live_genes != expected_genes:
        raise ValueError(
            f"Capacity {capacity_id} live GPR genes differ from the declared partition: "
            f"live={sorted(live_genes)}, expected={sorted(expected_genes)}"
        )
    if not math.isclose(
        float(source_reaction.lower_bound),
        0.0,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise ValueError(
            f"Capacity {capacity_id} requires a zero lower bound; positive minimum "
            "flux and reversible capacity need an explicit allocation model"
        )
    if row["provisional_upper_bound"] > source_reaction.upper_bound:
        raise ValueError(
            f"Capacity {capacity_id} would expand rather than tighten the source bound"
        )
    if (
        source_reaction.upper_bound - row["provisional_upper_bound"]
        < max(0.0, float(source_reaction.lower_bound))
    ):
        raise ValueError(
            f"Capacity {capacity_id} leaves an invalid primary-reaction bound"
        )


def _mark_gene(
    model,
    gene_id: str,
    row: dict[str, Any],
    role: str,
) -> None:
    gene = model.genes.get_by_id(gene_id)
    notes = dict(gene.notes or {})
    notes["provisional_capacity_status"] = MARKER_STATUS
    notes["provisional_capacity_replacement"] = row["replacement_formula"]
    notes["provisional_capacity_requires_protein_abundance"] = "true"
    notes["provisional_capacity_requires_kcat"] = "true"
    _append_marker(notes, "provisional_capacity_ids", row["capacity_id"])
    _append_marker(
        notes,
        "provisional_capacity_roles",
        f"{row['capacity_id']}:{role}",
    )
    _append_marker(
        notes,
        "provisional_capacity_reference_model_sha256",
        row["validated_model_sha256"],
    )
    _append_marker(
        notes,
        "provisional_capacity_target_fingerprints",
        row["target_fingerprint"],
    )
    gene.notes = notes


def _mark_reaction(
    reaction,
    row: dict[str, Any],
    role: str,
) -> None:
    notes = dict(reaction.notes or {})
    notes.update(
        {
            "provisional_capacity_status": MARKER_STATUS,
            "provisional_capacity_id": row["capacity_id"],
            "provisional_capacity_role": role,
            "provisional_capacity_case_id": row["case_id"],
            "provisional_capacity_parameter_basis": row["parameter_basis"],
            "provisional_capacity_replacement_formula": row["replacement_formula"],
            "provisional_capacity_requires_protein_abundance": "true",
            "provisional_capacity_requires_kcat": "true",
            "provisional_capacity_reference_model_sha256": row[
                "validated_model_sha256"
            ],
            "provisional_capacity_target_fingerprint": row["target_fingerprint"],
            "provisional_capacity_warning": "Not a measured biological parameter",
        }
    )
    reaction.notes = notes


def _already_applied(model, row: dict[str, Any]) -> bool:
    source_id = row["source_reaction_id"]
    backup_id = row["backup_reaction_id"]
    if source_id not in model.reactions or backup_id not in model.reactions:
        return False
    source = model.reactions.get_by_id(source_id)
    backup = model.reactions.get_by_id(backup_id)
    if not math.isclose(float(source.lower_bound), 0.0, abs_tol=1e-15):
        return False
    if not math.isclose(float(backup.lower_bound), 0.0, abs_tol=1e-15):
        return False
    if source.gene_reaction_rule != row["primary_gpr"]:
        return False
    if backup.gene_reaction_rule != row["backup_gpr"]:
        return False
    if not math.isclose(
        float(backup.upper_bound),
        float(row["provisional_upper_bound"]),
        rel_tol=1e-12,
        abs_tol=1e-15,
    ):
        return False
    original_upper_text = str(
        source.notes.get("provisional_capacity_original_total_upper_bound", "")
    )
    try:
        original_upper = float(original_upper_text)
    except ValueError:
        return False
    if not math.isclose(
        float(source.upper_bound) + float(backup.upper_bound),
        original_upper,
        rel_tol=1e-12,
        abs_tol=1e-12,
    ):
        return False
    if dict(source.metabolites) != dict(backup.metabolites):
        return False
    reconstructed_original_fingerprint = target_fingerprint(
        [
            {
                "reaction_id": source_id,
                "stoichiometry": {
                    metabolite.id: float(coefficient)
                    for metabolite, coefficient in source.metabolites.items()
                },
                "lower_bound": 0.0,
                "upper_bound": original_upper,
                "gpr": row["expected_gpr"],
            }
        ]
    )
    if reconstructed_original_fingerprint != row["target_fingerprint"]:
        return False
    if source.notes.get("provisional_capacity_id") != row["capacity_id"]:
        return False
    if backup.notes.get("provisional_capacity_id") != row["capacity_id"]:
        return False
    if (
        source.notes.get("provisional_capacity_target_fingerprint")
        != row["target_fingerprint"]
        or backup.notes.get("provisional_capacity_target_fingerprint")
        != row["target_fingerprint"]
    ):
        return False
    source_groups = {group.id for group in model.groups if source in group.members}
    backup_groups = {group.id for group in model.groups if backup in group.members}
    stored_source_groups = {
        group_id
        for group_id in str(
            source.notes.get("provisional_capacity_original_group_ids", "")
        ).split(";")
        if group_id
    }
    stored_backup_groups = {
        group_id
        for group_id in str(
            backup.notes.get("provisional_capacity_original_group_ids", "")
        ).split(";")
        if group_id
    }
    if not (
        source_groups
        == backup_groups
        == stored_source_groups
        == stored_backup_groups
    ):
        return False
    live_split_fingerprint = target_fingerprint(
        [_reaction_context(source), _reaction_context(backup)]
    )
    return (
        source.notes.get("provisional_capacity_split_fingerprint")
        == live_split_fingerprint
        and backup.notes.get("provisional_capacity_split_fingerprint")
        == live_split_fingerprint
    )


def apply_provisional_isozyme_capacities(
    model,
    table_path: Path,
    *,
    reference_model_sha256: str,
) -> list[dict[str, Any]]:
    """Apply a marked, reversible capacity profile to an experimental model."""
    rows = load_provisional_capacity_table(table_path)
    profile_sha256 = sha256_file(table_path)
    audit: list[dict[str, Any]] = []

    # Validate every row before mutating anything. A stale second row must not
    # leave the first capacity split partially applied in memory.
    preflight: list[tuple[dict[str, Any], bool]] = []
    for row in rows:
        if row["validated_model_sha256"] != reference_model_sha256:
            raise ValueError(
                f"Capacity {row['capacity_id']} was calibrated against model SHA "
                f"{row['validated_model_sha256']}, not {reference_model_sha256}"
            )
        source_id = row["source_reaction_id"]
        backup_id = row["backup_reaction_id"]
        if source_id not in model.reactions:
            raise ValueError(
                f"Capacity {row['capacity_id']} source reaction not found: {source_id}"
            )
        already_applied = _already_applied(model, row)
        if already_applied:
            preflight.append((row, True))
            continue
        if backup_id in model.reactions:
            raise ValueError(
                f"Capacity {row['capacity_id']} backup ID already exists with "
                "different content"
            )
        source = model.reactions.get_by_id(source_id)
        if source.gene_reaction_rule != row["expected_gpr"]:
            raise ValueError(
                f"Capacity {row['capacity_id']} expected GPR mismatch: "
                f"{row['expected_gpr']!r} != {source.gene_reaction_rule!r}"
            )
        _validate_partition(row, source)
        live_fingerprint = target_fingerprint([_reaction_context(source)])
        if live_fingerprint != row["target_fingerprint"]:
            raise ValueError(
                f"Capacity {row['capacity_id']} target fingerprint is stale: "
                f"{row['target_fingerprint']} -> {live_fingerprint}"
            )
        preflight.append((row, False))

    for row, already_applied in preflight:
        source_id = row["source_reaction_id"]
        backup_id = row["backup_reaction_id"]
        if already_applied:
            source = model.reactions.get_by_id(source_id)
            backup = model.reactions.get_by_id(backup_id)
            outcome = "already_applied"
        else:
            source = model.reactions.get_by_id(source_id)
            original_groups = [
                group for group in model.groups if source in group.members
            ]
            original_upper_bound = float(source.upper_bound)
            source.upper_bound = (
                original_upper_bound - float(row["provisional_upper_bound"])
            )
            backup = Reaction(
                backup_id,
                name=f"{source.name} [provisional backup isozyme capacity]",
                subsystem=source.subsystem,
                lower_bound=max(0.0, float(source.lower_bound)),
                upper_bound=float(row["provisional_upper_bound"]),
            )
            backup.add_metabolites(dict(source.metabolites))
            backup.annotation = copy.deepcopy(source.annotation or {})
            backup.notes = copy.deepcopy(source.notes or {})
            model.add_reactions([backup])
            backup.gene_reaction_rule = row["backup_gpr"]
            source.gene_reaction_rule = row["primary_gpr"]
            for group in original_groups:
                group.add_members([backup])
            original_group_ids = ";".join(
                sorted(group.id for group in original_groups)
            )
            source.notes["provisional_capacity_original_total_upper_bound"] = str(
                original_upper_bound
            )
            backup.notes["provisional_capacity_original_total_upper_bound"] = str(
                original_upper_bound
            )
            source.notes["provisional_capacity_original_group_ids"] = (
                original_group_ids
            )
            backup.notes["provisional_capacity_original_group_ids"] = (
                original_group_ids
            )
            outcome = "applied"

        _mark_reaction(source, row, "primary_isozyme_residual_parent_capacity")
        _mark_reaction(backup, row, "capacity_limited_backup_pool")
        backup.notes["provisional_capacity_upper_bound"] = str(
            row["provisional_upper_bound"]
        )
        backup.notes["provisional_capacity_units"] = row["units"]
        split_fingerprint = target_fingerprint(
            [_reaction_context(source), _reaction_context(backup)]
        )
        source.notes["provisional_capacity_split_fingerprint"] = split_fingerprint
        backup.notes["provisional_capacity_split_fingerprint"] = split_fingerprint
        for gene_id in sorted(_genes(row["primary_gpr"])):
            _mark_gene(model, gene_id, row, "primary")
        for gene_id in sorted(_genes(row["backup_gpr"])):
            _mark_gene(model, gene_id, row, "backup")

        audit.append(
            {
                "capacity_id": row["capacity_id"],
                "case_id": row["case_id"],
                "source_reaction_id": source_id,
                "backup_reaction_id": backup_id,
                "primary_gpr": source.gene_reaction_rule,
                "backup_gpr": backup.gene_reaction_rule,
                "primary_upper_bound": float(source.upper_bound),
                "backup_upper_bound": float(backup.upper_bound),
                "total_upper_bound": (
                    float(source.upper_bound) + float(backup.upper_bound)
                ),
                "units": row["units"],
                "parameter_basis": row["parameter_basis"],
                "profile_sha256": profile_sha256,
                "reference_model_sha256": row["validated_model_sha256"],
                "target_fingerprint": row["target_fingerprint"],
                "replacement_formula": row["replacement_formula"],
                "outcome": outcome,
            }
        )

    model_notes = dict(model.notes or {})
    model_notes["provisional_capacity_profile"] = str(table_path)
    model_notes["provisional_capacity_profile_sha256"] = profile_sha256
    model_notes["provisional_capacity_status"] = MARKER_STATUS
    model_notes["provisional_capacity_replacement"] = (
        "sum_i(kcat_i_per_s*E_i_mmol_per_gDW*3600)"
    )
    for row in rows:
        _append_marker(
            model_notes,
            "provisional_capacity_reference_model_sha256",
            row["validated_model_sha256"],
        )
    model_notes["provisional_capacity_warning"] = (
        "Experimental profile only; replace all marked bounds with measured "
        "protein abundance and kcat before canonical promotion"
    )
    model.notes = model_notes
    return audit
