"""Auditable reference chemistry for essentiality-linked reaction cases.

The curation table records Rhea-backed target equations separately from model
application.  ``component_review`` proposals can be evaluated in memory, but
this module never activates them and never writes SBML.  A proposal is ready
only when its target reaction becomes mass/charge balanced *and* the complete
microspecies substitution breaks no previously balanced internal reaction.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from cobra.io import read_sbml_model

from .config import ESSENTIALITY_DIR, REPO_ROOT
from .microspecies import (
    DEFAULT_MICROSPECIES_TABLE,
    MicrospeciesRow,
    _fully_balanced_internal,
    _is_allowed_current_pair,
    _metabolite_pair,
    _reaction_balance_record,
    _resolve_pinned_targets,
    balance_protons_and_water,
    load_curated_microspecies,
    metabolite_base_name,
)

DEFAULT_REACTION_CHEMISTRY_TABLE = (
    ESSENTIALITY_DIR / "reaction_chemistry_curation.csv"
)
DEFAULT_DOSSIER_DIR = ESSENTIALITY_DIR / "evidence"

_ALLOWED_STATUSES = {"component_review", "active"}
_CASE_ID_RE = re.compile(r"^EGC-[0-9a-z-]+$")
_RHEA_ID_RE = re.compile(r"^RHEA:\d+$")
_REQUIRED_COLUMNS = {
    "schema_version",
    "status",
    "case_id",
    "reaction_id",
    "reference_reaction_id",
    "coefficient_updates_json",
    "required_microspecies_families",
    "evidence_url",
    "rationale",
}


@dataclass(frozen=True)
class ReactionChemistryProposal:
    schema_version: int
    status: str
    case_id: str
    reaction_id: str
    reference_reaction_id: str
    coefficient_updates: Mapping[str, float]
    required_microspecies_families: tuple[str, ...]
    evidence_url: str
    rationale: str


def _parse_coefficient_updates(raw: str) -> dict[str, float]:
    try:
        payload = json.loads(raw or "{}")
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid coefficient_updates_json: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("coefficient_updates_json must be a JSON object")
    updates: dict[str, float] = {}
    for metabolite_id, coefficient in payload.items():
        if not isinstance(metabolite_id, str) or not metabolite_id.strip():
            raise ValueError("coefficient update metabolite IDs must be nonempty strings")
        if isinstance(coefficient, bool) or not isinstance(coefficient, (int, float)):
            raise ValueError(
                f"coefficient for {metabolite_id!r} must be a finite number"
            )
        value = float(coefficient)
        if not math.isfinite(value):
            raise ValueError(
                f"coefficient for {metabolite_id!r} must be a finite number"
            )
        updates[metabolite_id.strip()] = value
    return updates


def load_reaction_chemistry_proposals(
    table_path: str | Path = DEFAULT_REACTION_CHEMISTRY_TABLE,
    *,
    microspecies_table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> list[ReactionChemistryProposal]:
    """Load and validate the durable reference-reaction proposal table."""

    path = Path(table_path)
    if not path.exists():
        raise FileNotFoundError(f"reaction chemistry table not found: {path}")

    microspecies_rows = load_curated_microspecies(microspecies_table_path)
    known_families = {row.family_id for row in microspecies_rows}
    errors: list[str] = []
    proposals: list[ReactionChemistryProposal] = []
    seen_targets: set[tuple[str, str]] = set()

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        headers = set(reader.fieldnames or [])
        missing = sorted(_REQUIRED_COLUMNS - headers)
        if missing:
            raise ValueError(
                f"reaction chemistry table lacks required columns: {missing}"
            )

        for line_number, raw in enumerate(reader, start=2):
            prefix = f"row {line_number}"
            try:
                schema_version = int(raw["schema_version"].strip())
                status = raw["status"].strip()
                case_id = raw["case_id"].strip()
                reaction_id = raw["reaction_id"].strip()
                reference_reaction_id = raw["reference_reaction_id"].strip()
                coefficient_updates = _parse_coefficient_updates(
                    raw["coefficient_updates_json"].strip()
                )
                families = tuple(
                    value.strip()
                    for value in raw["required_microspecies_families"].split(";")
                    if value.strip()
                )
                evidence_url = raw["evidence_url"].strip()
                rationale = raw["rationale"].strip()
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"{prefix}: cannot parse row ({exc})")
                continue

            if schema_version != 1:
                errors.append(f"{prefix}: unsupported schema_version {schema_version}")
            if status not in _ALLOWED_STATUSES:
                errors.append(f"{prefix}: unsupported status {status!r}")
            if not _CASE_ID_RE.fullmatch(case_id):
                errors.append(f"{prefix}: invalid case_id {case_id!r}")
            if not reaction_id:
                errors.append(f"{prefix}: reaction_id is required")
            if not _RHEA_ID_RE.fullmatch(reference_reaction_id):
                errors.append(
                    f"{prefix}: invalid reference_reaction_id "
                    f"{reference_reaction_id!r}"
                )
            if not families:
                errors.append(f"{prefix}: at least one microspecies family is required")
            if len(families) != len(set(families)):
                errors.append(f"{prefix}: duplicate microspecies family")
            if tuple(sorted(families)) != families:
                errors.append(
                    f"{prefix}: microspecies families must be sorted deterministically"
                )
            unknown = sorted(set(families) - known_families)
            if unknown:
                errors.append(f"{prefix}: unknown microspecies families: {unknown}")
            if not evidence_url.startswith(("https://", "http://")):
                errors.append(f"{prefix}: evidence_url must be HTTP(S)")
            if not rationale:
                errors.append(f"{prefix}: rationale is required")

            target = (case_id, reaction_id)
            if target in seen_targets:
                errors.append(f"{prefix}: duplicate case/reaction target {target!r}")
            seen_targets.add(target)

            proposals.append(
                ReactionChemistryProposal(
                    schema_version=schema_version,
                    status=status,
                    case_id=case_id,
                    reaction_id=reaction_id,
                    reference_reaction_id=reference_reaction_id,
                    coefficient_updates=coefficient_updates,
                    required_microspecies_families=families,
                    evidence_url=evidence_url,
                    rationale=rationale,
                )
            )

    if errors:
        raise ValueError(
            "invalid reaction chemistry curation table:\n- " + "\n- ".join(errors)
        )
    return proposals


def _row_index(
    microspecies_table_path: str | Path,
) -> dict[str, MicrospeciesRow]:
    rows = load_curated_microspecies(microspecies_table_path)
    return {row.family_id: row for row in rows}


def audit_reference_reaction_chemistry(
    model,
    case_id: str,
    *,
    table_path: str | Path = DEFAULT_REACTION_CHEMISTRY_TABLE,
    microspecies_table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
    tolerance: float = 1e-9,
) -> dict:
    """Evaluate one case's Rhea chemistry transactionally in memory.

    The target formulas/charges and final reaction coefficients are applied to
    the supplied model only long enough to calculate the audit.  They are
    restored even if validation fails.  ``ready_for_activation`` additionally
    requires zero regressions among previously balanced touched reactions and
    zero generic H+/H2O bookkeeping suggestions.  Such suggestions are useful
    diagnostics, but must first become explicit reference-backed coefficient
    proposals before they can authorize activation.
    """

    proposals = [
        proposal
        for proposal in load_reaction_chemistry_proposals(
            table_path, microspecies_table_path=microspecies_table_path
        )
        if proposal.case_id == case_id
    ]
    if not proposals:
        raise ValueError(f"no reaction chemistry proposal for {case_id!r}")
    proposals.sort(key=lambda proposal: proposal.reaction_id)

    rows_by_family = _row_index(microspecies_table_path)
    family_ids = sorted(
        {
            family_id
            for proposal in proposals
            for family_id in proposal.required_microspecies_families
        }
    )
    family_rows = [rows_by_family[family_id] for family_id in family_ids]

    targets_by_family: dict[str, list] = {}
    metabolite_snapshots: dict[str, tuple[str | None, int | None]] = {}
    for row in family_rows:
        targets = _resolve_pinned_targets(model, row)
        targets_by_family[row.family_id] = targets
        for metabolite in targets:
            current = _metabolite_pair(metabolite)
            if not _is_allowed_current_pair(row, metabolite):
                raise ValueError(
                    f"{row.family_id}: {metabolite.id} has unexpected current pair "
                    f"{current!r}"
                )
            metabolite_snapshots[metabolite.id] = (
                metabolite.formula,
                metabolite.charge,
            )

    reaction_snapshots: dict[tuple[str, str], float] = {}
    proposal_reactions = []
    for proposal in proposals:
        try:
            reaction = model.reactions.get_by_id(proposal.reaction_id)
        except KeyError as exc:
            raise ValueError(
                f"{proposal.case_id}: missing reaction {proposal.reaction_id!r}"
            ) from exc
        proposal_reactions.append(reaction)
        for metabolite_id in proposal.coefficient_updates:
            try:
                metabolite = model.metabolites.get_by_id(metabolite_id)
            except KeyError as exc:
                raise ValueError(
                    f"{proposal.reaction_id}: missing metabolite {metabolite_id!r}"
                ) from exc
            reaction_snapshots[(reaction.id, metabolite.id)] = float(
                reaction.metabolites.get(metabolite, 0.0)
            )

    touched_reactions = {
        reaction
        for targets in targets_by_family.values()
        for metabolite in targets
        for reaction in metabolite.reactions
    } | set(proposal_reactions)
    balanced_before = _fully_balanced_internal(touched_reactions)
    target_before = {
        reaction.id: _reaction_balance_record(reaction, tolerance)
        for reaction in proposal_reactions
    }
    changed_metabolites: list[str] = []
    coefficient_changes: list[dict] = []
    bookkeeping_changes: list[dict] = []
    bookkeeping_rejected: list[dict] = []
    preliminary_regressions: list[str] = []

    try:
        for row in family_rows:
            for metabolite in targets_by_family[row.family_id]:
                if _metabolite_pair(metabolite) != row.target_pair:
                    metabolite.formula = row.target_formula
                    metabolite.charge = row.target_charge
                    changed_metabolites.append(metabolite.id)

        for proposal, reaction in zip(proposals, proposal_reactions, strict=True):
            for metabolite_id, target_coefficient in sorted(
                proposal.coefficient_updates.items()
            ):
                metabolite = model.metabolites.get_by_id(metabolite_id)
                before_coefficient = float(reaction.metabolites.get(metabolite, 0.0))
                delta = target_coefficient - before_coefficient
                if not math.isclose(delta, 0.0, abs_tol=tolerance):
                    reaction.add_metabolites({metabolite: delta})
                coefficient_changes.append(
                    {
                        "reaction_id": reaction.id,
                        "metabolite_id": metabolite.id,
                        "before_coefficient": before_coefficient,
                        "after_coefficient": float(
                            reaction.metabolites.get(metabolite, 0.0)
                        ),
                    }
                )

        target_after = {
            reaction.id: _reaction_balance_record(reaction, tolerance)
            for reaction in proposal_reactions
        }
        preliminary_balanced_after = _fully_balanced_internal(touched_reactions)
        preliminary_regressions = sorted(
            balanced_before - preliminary_balanced_after
        )

        # Once every participant uses the same pH-7.3 microspecies, an
        # otherwise balanced single-compartment reaction may need a different
        # explicit H+/H2O coefficient.  The shared algebraic gate admits that
        # bookkeeping only when hydrogen, oxygen and charge can all be solved
        # simultaneously; heavy-element residuals, transport reactions and
        # boundary reactions remain rejected.  These edits are audit-only and
        # are rolled back below.
        if preliminary_regressions:
            bookkeeping_report = balance_protons_and_water(
                model,
                reaction_ids=preliminary_regressions,
                tolerance=tolerance,
                table_path=microspecies_table_path,
            )
            bookkeeping_changes = bookkeeping_report["changes"]
            bookkeeping_rejected = bookkeeping_report["rejected"]

        target_after = {
            reaction.id: _reaction_balance_record(reaction, tolerance)
            for reaction in proposal_reactions
        }
        balanced_after = _fully_balanced_internal(touched_reactions)
    finally:
        for change in reversed(bookkeeping_changes):
            reaction = model.reactions.get_by_id(change["reaction_id"])
            reaction.add_metabolites(
                {
                    model.metabolites.get_by_id(metabolite_id): -float(coefficient)
                    for metabolite_id, coefficient in change["additions"].items()
                }
            )
        for (reaction_id, metabolite_id), coefficient in reaction_snapshots.items():
            reaction = model.reactions.get_by_id(reaction_id)
            metabolite = model.metabolites.get_by_id(metabolite_id)
            current = float(reaction.metabolites.get(metabolite, 0.0))
            delta = coefficient - current
            if not math.isclose(delta, 0.0, abs_tol=tolerance):
                reaction.add_metabolites({metabolite: delta})
        for metabolite_id, (formula, charge) in metabolite_snapshots.items():
            metabolite = model.metabolites.get_by_id(metabolite_id)
            metabolite.formula = formula
            metabolite.charge = charge

    regressed = sorted(balanced_before - balanced_after)
    fixed = sorted(balanced_after - balanced_before)
    unbalanced_targets = sorted(
        reaction_id
        for reaction_id, record in target_after.items()
        if record["status"] != "balanced"
    )

    fingerprint_payload = {
        "schema": "reference_reaction_chemistry_audit_v1",
        "case_id": case_id,
        "proposals": [
            {
                "reaction_id": proposal.reaction_id,
                "reference_reaction_id": proposal.reference_reaction_id,
                "coefficient_updates": dict(proposal.coefficient_updates),
                "families": list(proposal.required_microspecies_families),
            }
            for proposal in proposals
        ],
        "current_metabolite_pairs": {
            metabolite_id: [formula, charge]
            for metabolite_id, (formula, charge) in sorted(
                metabolite_snapshots.items()
            )
        },
        "target_before": target_before,
    }
    fingerprint = "sha256:" + hashlib.sha256(
        json.dumps(
            fingerprint_payload, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()

    reported_bookkeeping_changes = [
        {
            **change,
            "status": "provisional_requires_reference_review",
            "ready_for_activation": False,
        }
        for change in bookkeeping_changes
    ]

    return {
        "schema_version": 1,
        "case_id": case_id,
        "audit_fingerprint": fingerprint,
        "status": "reference_verified",
        "reference_reaction_ids": {
            proposal.reaction_id: proposal.reference_reaction_id
            for proposal in proposals
        },
        "required_microspecies_families": family_ids,
        "changed_metabolite_ids": sorted(changed_metabolites),
        "coefficient_changes": coefficient_changes,
        "microspecies_bookkeeping_changes": reported_bookkeeping_changes,
        "microspecies_bookkeeping_rejected": bookkeeping_rejected,
        "preliminary_regressed_reaction_ids": preliminary_regressions,
        "target_reaction_balances": {
            reaction.id: {
                "before": target_before[reaction.id],
                "after": target_after[reaction.id],
            }
            for reaction in proposal_reactions
        },
        "fixed_reaction_ids": fixed,
        "regressed_reaction_ids": regressed,
        "unbalanced_target_reaction_ids": unbalanced_targets,
        "reference_equations_balanced": not unbalanced_targets,
        "ready_for_activation": (
            not unbalanced_targets and not regressed and not bookkeeping_changes
        ),
    }


def _annotation_strings(value) -> list[str]:
    if value in (None, ""):
        return []
    if isinstance(value, (list, tuple, set, frozenset)):
        return sorted({str(item) for item in value if str(item).strip()})
    return [str(value)]


def audit_global_reference_component(
    model,
    *,
    table_path: str | Path = DEFAULT_REACTION_CHEMISTRY_TABLE,
    microspecies_table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
    tolerance: float = 1e-9,
) -> dict:
    """Audit every deferred microspecies and Rhea coefficient as one batch.

    This is a read-only planning transaction.  It ranks the remaining reaction
    frontier by shared outside metabolite family so the next literature pass
    can address the largest chemically coupled cluster first.
    """

    proposals = sorted(
        (
            proposal
            for proposal in load_reaction_chemistry_proposals(
                table_path, microspecies_table_path=microspecies_table_path
            )
            if proposal.status == "component_review"
        ),
        key=lambda proposal: (proposal.reaction_id, proposal.case_id),
    )
    all_microspecies_rows = load_curated_microspecies(microspecies_table_path)
    rows = sorted(
        (
            row
            for row in all_microspecies_rows
            if row.status == "component_review"
        ),
        key=lambda row: row.family_id,
    )
    if not proposals or not rows:
        raise ValueError("global component audit requires deferred proposals and families")

    targets_by_family: dict[str, list] = {}
    metabolite_snapshots: dict[str, tuple[str | None, int | None]] = {}
    for row in rows:
        targets = _resolve_pinned_targets(model, row)
        targets_by_family[row.family_id] = targets
        for metabolite in targets:
            current = _metabolite_pair(metabolite)
            if not _is_allowed_current_pair(row, metabolite):
                raise ValueError(
                    f"{row.family_id}: {metabolite.id} has unexpected current pair "
                    f"{current!r}"
                )
            metabolite_snapshots[metabolite.id] = (
                metabolite.formula,
                metabolite.charge,
            )

    proposal_reactions = []
    reaction_snapshots: dict[tuple[str, str], float] = {}
    requested_coefficients: dict[tuple[str, str], float] = {}
    for proposal in proposals:
        try:
            reaction = model.reactions.get_by_id(proposal.reaction_id)
        except KeyError as exc:
            raise ValueError(
                f"{proposal.case_id}: missing reaction {proposal.reaction_id!r}"
            ) from exc
        proposal_reactions.append(reaction)
        for metabolite_id, coefficient in proposal.coefficient_updates.items():
            key = (reaction.id, metabolite_id)
            previous = requested_coefficients.get(key)
            if previous is not None and not math.isclose(
                previous, coefficient, abs_tol=tolerance
            ):
                raise ValueError(f"conflicting coefficient proposals for {key!r}")
            requested_coefficients[key] = coefficient
            try:
                metabolite = model.metabolites.get_by_id(metabolite_id)
            except KeyError as exc:
                raise ValueError(
                    f"{proposal.reaction_id}: missing metabolite {metabolite_id!r}"
                ) from exc
            reaction_snapshots[key] = float(
                reaction.metabolites.get(metabolite, 0.0)
            )

    selected_metabolite_ids = set(metabolite_snapshots)
    curated_metabolite_membership: dict[str, tuple[str, str]] = {}
    for row in all_microspecies_rows:
        for metabolite_id in row.expected_metabolite_ids:
            curated_metabolite_membership[metabolite_id] = (
                row.family_id,
                row.status,
            )
    touched_reactions = {
        reaction
        for targets in targets_by_family.values()
        for metabolite in targets
        for reaction in metabolite.reactions
    } | set(proposal_reactions)
    before = {
        reaction.id: _reaction_balance_record(reaction, tolerance)
        for reaction in touched_reactions
    }
    balanced_before = _fully_balanced_internal(touched_reactions)
    changed_metabolites: list[str] = []
    coefficient_changes: list[dict] = []
    bookkeeping_changes: list[dict] = []
    bookkeeping_rejected: list[dict] = []
    preliminary_regressions: list[str] = []

    try:
        for row in rows:
            for metabolite in targets_by_family[row.family_id]:
                if _metabolite_pair(metabolite) != row.target_pair:
                    metabolite.formula = row.target_formula
                    metabolite.charge = row.target_charge
                    changed_metabolites.append(metabolite.id)

        for (reaction_id, metabolite_id), target_coefficient in sorted(
            requested_coefficients.items()
        ):
            reaction = model.reactions.get_by_id(reaction_id)
            metabolite = model.metabolites.get_by_id(metabolite_id)
            before_coefficient = float(reaction.metabolites.get(metabolite, 0.0))
            delta = target_coefficient - before_coefficient
            if not math.isclose(delta, 0.0, abs_tol=tolerance):
                reaction.add_metabolites({metabolite: delta})
            coefficient_changes.append(
                {
                    "reaction_id": reaction_id,
                    "metabolite_id": metabolite_id,
                    "before_coefficient": before_coefficient,
                    "after_coefficient": float(
                        reaction.metabolites.get(metabolite, 0.0)
                    ),
                }
            )

        preliminary_after = {
            reaction.id: _reaction_balance_record(reaction, tolerance)
            for reaction in touched_reactions
        }
        preliminary_balanced_after = _fully_balanced_internal(touched_reactions)
        preliminary_regressions = sorted(
            balanced_before - preliminary_balanced_after
        )
        if preliminary_regressions:
            bookkeeping = balance_protons_and_water(
                model,
                reaction_ids=preliminary_regressions,
                tolerance=tolerance,
                table_path=microspecies_table_path,
            )
            bookkeeping_changes = bookkeeping["changes"]
            bookkeeping_rejected = bookkeeping["rejected"]
        after = {
            reaction.id: _reaction_balance_record(reaction, tolerance)
            for reaction in touched_reactions
        }
        balanced_after = _fully_balanced_internal(touched_reactions)
    finally:
        for change in reversed(bookkeeping_changes):
            reaction = model.reactions.get_by_id(change["reaction_id"])
            reaction.add_metabolites(
                {
                    model.metabolites.get_by_id(metabolite_id): -float(coefficient)
                    for metabolite_id, coefficient in change["additions"].items()
                }
            )
        for (reaction_id, metabolite_id), coefficient in reaction_snapshots.items():
            reaction = model.reactions.get_by_id(reaction_id)
            metabolite = model.metabolites.get_by_id(metabolite_id)
            current = float(reaction.metabolites.get(metabolite, 0.0))
            delta = coefficient - current
            if not math.isclose(delta, 0.0, abs_tol=tolerance):
                reaction.add_metabolites({metabolite: delta})
        for metabolite_id, (formula, charge) in metabolite_snapshots.items():
            metabolite = model.metabolites.get_by_id(metabolite_id)
            metabolite.formula = formula
            metabolite.charge = charge

    regressed = sorted(balanced_before - balanced_after)
    fixed = sorted(balanced_after - balanced_before)
    unresolved = sorted(
        reaction_id
        for reaction_id, record in after.items()
        if record["status"] in {"imbalanced", "uncheckable", "error"}
    )
    proposal_reaction_ids = sorted({proposal.reaction_id for proposal in proposals})
    unbalanced_proposals = sorted(
        reaction_id
        for reaction_id in proposal_reaction_ids
        if after[reaction_id]["status"] != "balanced"
    )

    # Generic proton/water balancing is only a provisional diagnostic.  Keep
    # every preliminary regression in the research frontier so an incomplete
    # component cannot disappear merely because a local H+ coefficient makes
    # it algebraically balanced under still-uncurated partner microspecies.
    regression_set = set(preliminary_regressions)
    unresolved_set = set(unresolved)
    family_records: dict[str, dict] = {}
    family_reactions: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"regressed": set(), "unresolved": set()}
    )
    for reaction_id in sorted(regression_set | unresolved_set):
        reaction = model.reactions.get_by_id(reaction_id)
        for metabolite in reaction.metabolites:
            if metabolite.id in selected_metabolite_ids:
                continue
            family_id = metabolite_base_name(metabolite) or metabolite.id.casefold()
            record = family_records.setdefault(
                family_id,
                {
                    "candidate_family": family_id,
                    "metabolite_ids": set(),
                    "current_pairs": set(),
                    "compartments": set(),
                    "rhea_ids": set(),
                    "kegg_reaction_context": set(),
                    "existing_family_ids": set(),
                    "existing_family_statuses": set(),
                },
            )
            record["metabolite_ids"].add(metabolite.id)
            record["current_pairs"].add(
                (metabolite.formula, metabolite.charge)
            )
            record["compartments"].add(metabolite.compartment)
            existing_membership = curated_metabolite_membership.get(metabolite.id)
            if existing_membership is not None:
                record["existing_family_ids"].add(existing_membership[0])
                record["existing_family_statuses"].add(existing_membership[1])
            annotation = metabolite.annotation or {}
            record["rhea_ids"].update(_annotation_strings(annotation.get("rhea")))
            reaction_annotation = reaction.annotation or {}
            record["kegg_reaction_context"].update(
                _annotation_strings(reaction_annotation.get("kegg.reaction"))
            )
            if reaction_id in regression_set:
                family_reactions[family_id]["regressed"].add(reaction_id)
            if reaction_id in unresolved_set:
                family_reactions[family_id]["unresolved"].add(reaction_id)

    frontier_ranking = []
    for family_id, record in family_records.items():
        pairs = sorted(
            (
                {"formula": formula, "charge": charge}
                for formula, charge in record["current_pairs"]
            ),
            key=lambda item: (str(item["formula"]), str(item["charge"])),
        )
        polymer_or_uncheckable = any(
            not pair["formula"]
            or "(" in str(pair["formula"])
            or "*" in str(pair["formula"])
            for pair in pairs
        )
        regressed_reactions = sorted(family_reactions[family_id]["regressed"])
        unresolved_reactions = sorted(family_reactions[family_id]["unresolved"])
        frontier_ranking.append(
            {
                "candidate_family": family_id,
                "metabolite_ids": sorted(record["metabolite_ids"]),
                "current_pairs": pairs,
                "compartments": sorted(record["compartments"]),
                "regressed_reaction_count": len(regressed_reactions),
                "regressed_reaction_ids": regressed_reactions,
                "unresolved_reaction_count": len(unresolved_reactions),
                "unresolved_reaction_ids": unresolved_reactions,
                "polymer_or_uncheckable": polymer_or_uncheckable,
                "existing_family_ids": sorted(record["existing_family_ids"]),
                "existing_family_statuses": sorted(
                    record["existing_family_statuses"]
                ),
                "rhea_ids": sorted(record["rhea_ids"]),
                "kegg_reaction_context": sorted(record["kegg_reaction_context"]),
            }
        )
    frontier_ranking.sort(
        key=lambda item: (
            -item["regressed_reaction_count"],
            -item["unresolved_reaction_count"],
            item["candidate_family"],
        )
    )
    actionable_frontier_ranking = [
        item
        for item in frontier_ranking
        if item["regressed_reaction_count"] > 0
        and not item["existing_family_ids"]
        and not item["polymer_or_uncheckable"]
    ]

    fingerprint_payload = {
        "schema": "global_reference_component_audit_v1",
        "families": [row.family_id for row in rows],
        "proposals": [
            {
                "case_id": proposal.case_id,
                "reaction_id": proposal.reaction_id,
                "coefficient_updates": dict(proposal.coefficient_updates),
            }
            for proposal in proposals
        ],
        "current_metabolite_pairs": {
            metabolite_id: [formula, charge]
            for metabolite_id, (formula, charge) in sorted(
                metabolite_snapshots.items()
            )
        },
        "before": {reaction_id: before[reaction_id] for reaction_id in sorted(before)},
    }
    audit_fingerprint = "sha256:" + hashlib.sha256(
        json.dumps(
            fingerprint_payload, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()

    reported_bookkeeping_changes = [
        {
            **change,
            "status": "provisional_requires_reference_review",
            "ready_for_activation": False,
        }
        for change in bookkeeping_changes
    ]

    return {
        "schema_version": 1,
        "audit_fingerprint": audit_fingerprint,
        "status": "component_review",
        "family_ids": [row.family_id for row in rows],
        "proposal_reaction_ids": proposal_reaction_ids,
        "changed_metabolite_ids": sorted(changed_metabolites),
        "coefficient_changes": coefficient_changes,
        "touched_reaction_ids": sorted(before),
        "fixed_reaction_ids": fixed,
        "preliminary_regressed_reaction_ids": preliminary_regressions,
        "frontier_regressed_reaction_ids": sorted(regression_set),
        "regressed_reaction_ids": regressed,
        "unresolved_reaction_ids": unresolved,
        "unbalanced_proposal_reaction_ids": unbalanced_proposals,
        "microspecies_bookkeeping_changes": reported_bookkeeping_changes,
        "microspecies_bookkeeping_rejected": bookkeeping_rejected,
        "frontier_metabolite_ranking": frontier_ranking,
        "actionable_frontier_ranking": actionable_frontier_ranking,
        "reaction_balances": {
            reaction_id: {
                "before": before[reaction_id],
                "preliminary_after": preliminary_after[reaction_id],
                "after": after[reaction_id],
            }
            for reaction_id in sorted(before)
        },
        "reference_equations_balanced": not unbalanced_proposals,
        "ready_for_activation": (
            not unbalanced_proposals and not regressed and not bookkeeping_changes
        ),
    }


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_reference_chemistry_audits(
    model_path: str | Path,
    output_dir: str | Path,
    *,
    case_ids: list[str] | None = None,
    table_path: str | Path = DEFAULT_REACTION_CHEMISTRY_TABLE,
    microspecies_table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> list[Path]:
    """Write provenance-rich chemistry audits atomically for durable review."""

    model_path = Path(model_path)
    output_dir = Path(output_dir)
    proposals = load_reaction_chemistry_proposals(
        table_path, microspecies_table_path=microspecies_table_path
    )
    available_case_ids = sorted({proposal.case_id for proposal in proposals})
    selected_case_ids = sorted(set(case_ids or available_case_ids))
    unknown = sorted(set(selected_case_ids) - set(available_case_ids))
    if unknown:
        raise ValueError(f"no chemistry proposals for case IDs: {unknown}")

    model = read_sbml_model(str(model_path))
    provenance = {
        "model_path": str(model_path.resolve()),
        "model_sha256": _sha256_file(model_path),
        "reaction_chemistry_table_sha256": _sha256_file(table_path),
        "microspecies_table_sha256": _sha256_file(microspecies_table_path),
        "generated_at": dt.datetime.now(dt.UTC).isoformat(),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for case_id in selected_case_ids:
        report = audit_reference_reaction_chemistry(
            model,
            case_id,
            table_path=table_path,
            microspecies_table_path=microspecies_table_path,
        )
        payload = {**provenance, **report}
        output_path = output_dir / f"{case_id}.json"
        temporary_path = output_path.with_suffix(".json.tmp")
        temporary_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary_path.replace(output_path)
        written.append(output_path)
    return written


def write_global_reference_component_audit(
    model_path: str | Path,
    output_path: str | Path,
    *,
    table_path: str | Path = DEFAULT_REACTION_CHEMISTRY_TABLE,
    microspecies_table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> Path:
    """Write one provenance-rich global component audit atomically."""

    model_path = Path(model_path)
    output_path = Path(output_path)
    model = read_sbml_model(str(model_path))
    report = audit_global_reference_component(
        model,
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )
    payload = {
        "model_path": str(model_path.resolve()),
        "model_sha256": _sha256_file(model_path),
        "reaction_chemistry_table_sha256": _sha256_file(table_path),
        "microspecies_table_sha256": _sha256_file(microspecies_table_path),
        "generated_at": dt.datetime.now(dt.UTC).isoformat(),
        **report,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary_path.replace(output_path)
    return output_path


def attach_reference_chemistry_audit(
    dossier_path: str | Path,
    audit_path: str | Path,
    *,
    model_path: str | Path = REPO_ROOT / "model.xml",
    table_path: str | Path = DEFAULT_REACTION_CHEMISTRY_TABLE,
    microspecies_table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> dict:
    """Attach a current-input audit summary to one durable evidence dossier.

    The full audit remains in ``evidence/chemistry``.  This guarded importer
    records only its fingerprint, gate results and counts in the case dossier.
    It refuses stale model/table inputs and never marks unapplied reference
    chemistry as ``balanced``.
    """

    dossier_path = Path(dossier_path)
    audit_path = Path(audit_path)
    if not dossier_path.is_file():
        raise FileNotFoundError(f"evidence dossier not found: {dossier_path}")
    if not audit_path.is_file():
        raise FileNotFoundError(f"reference chemistry audit not found: {audit_path}")

    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    case_id = str(dossier.get("case_id", "")).strip()
    if not case_id or case_id != str(audit.get("case_id", "")).strip():
        raise ValueError("chemistry audit case_id does not match evidence dossier")

    expected_hashes = {
        "model_sha256": _sha256_file(model_path),
        "reaction_chemistry_table_sha256": _sha256_file(table_path),
        "microspecies_table_sha256": _sha256_file(microspecies_table_path),
    }
    for field, expected in expected_hashes.items():
        if str(audit.get(field, "")).strip() != expected:
            raise ValueError(f"stale chemistry audit: {field} does not match")
    if str(dossier.get("model_sha256", "")).strip() != expected_hashes[
        "model_sha256"
    ]:
        raise ValueError("evidence dossier model SHA does not match current model")
    if audit.get("status") != "reference_verified":
        raise ValueError("chemistry audit status is not reference_verified")

    chemistry_review = dossier.setdefault("chemistry_review", {})
    if not isinstance(chemistry_review, dict):
        raise ValueError("evidence dossier chemistry_review must be an object")
    try:
        stored_path = audit_path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()
    except ValueError:
        stored_path = str(audit_path.resolve())

    rejected = audit.get("microspecies_bookkeeping_rejected", [])
    rejection_counts: dict[str, int] = {}
    for row in rejected if isinstance(rejected, list) else []:
        reason = str(row.get("reason", "unknown"))
        rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
    ready = audit.get("ready_for_activation") is True
    chemistry_review["reference_chemistry_audit_path"] = stored_path
    chemistry_review["reference_chemistry_audit"] = {
        "audit_fingerprint": audit.get("audit_fingerprint", ""),
        "reference_reaction_ids": audit.get("reference_reaction_ids", {}),
        "reference_equations_balanced": audit.get(
            "reference_equations_balanced", False
        ),
        "ready_for_activation": ready,
        "preliminary_regressed_reactions": len(
            audit.get("preliminary_regressed_reaction_ids", [])
        ),
        "final_regressed_reactions": len(audit.get("regressed_reaction_ids", [])),
        "bookkeeping_changes": len(
            audit.get("microspecies_bookkeeping_changes", [])
        ),
        "bookkeeping_rejection_counts": rejection_counts,
        **expected_hashes,
    }
    chemistry_review["status"] = (
        "reference_ready_not_applied" if ready else "blocked_component_migration"
    )
    chemistry_review["model_sha256"] = expected_hashes["model_sha256"]

    temporary_path = dossier_path.with_suffix(dossier_path.suffix + ".tmp")
    temporary_path.write_text(
        json.dumps(dossier, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary_path.replace(dossier_path)
    return dossier


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit Rhea-backed essentiality reaction chemistry"
    )
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--case-id",
        action="append",
        dest="case_ids",
        help="case ID to audit; repeat as needed (default: all proposals)",
    )
    parser.add_argument(
        "--table", type=Path, default=DEFAULT_REACTION_CHEMISTRY_TABLE
    )
    parser.add_argument(
        "--microspecies-table", type=Path, default=DEFAULT_MICROSPECIES_TABLE
    )
    parser.add_argument(
        "--attach-dossiers",
        action="store_true",
        help="guardedly attach each audit summary to its durable case dossier",
    )
    parser.add_argument(
        "--dossier-dir", type=Path, default=DEFAULT_DOSSIER_DIR
    )
    parser.add_argument(
        "--global-component",
        action="store_true",
        help="also write global_component_migration.json for the full deferred batch",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    paths = write_reference_chemistry_audits(
        args.model,
        args.output_dir,
        case_ids=args.case_ids,
        table_path=args.table,
        microspecies_table_path=args.microspecies_table,
    )
    if args.global_component:
        global_path = write_global_reference_component_audit(
            args.model,
            args.output_dir / "global_component_migration.json",
            table_path=args.table,
            microspecies_table_path=args.microspecies_table,
        )
        print(global_path)
    for path in paths:
        if args.attach_dossiers:
            attach_reference_chemistry_audit(
                args.dossier_dir / path.name,
                path,
                model_path=args.model,
                table_path=args.table,
                microspecies_table_path=args.microspecies_table,
            )
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
