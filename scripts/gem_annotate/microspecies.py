"""Curated pH 7.3 microspecies convention for the iYali26 build pipeline.

The source workbook/model mixes fully protonated formulae with ionic charge
fields.  This module provides one auditable representation layer based on the
Rhea/ChEBI major microspecies at pH 7.3.  Formula and charge are always applied
as one atomic pair.  Rows that still require a connected-component migration
are deliberately visible in the curation table but are never applied here.

The module also implements the two acid/base bookkeeping rules used by the
pipeline:

* explicit hydroxide is normalised to ``H2O - H+``;
* H+/H2O may be added only when mass and charge can both become zero.

None of these functions writes SBML.  They mutate only the in-memory model
passed by the normal ``data/iyali26.xml -> model.xml`` build pipeline.
"""

from __future__ import annotations

import copy
import csv
import datetime as dt
import hashlib
import json
import logging
import math
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .config import CURATION_DATA_DIR

logger = logging.getLogger(__name__)

DEFAULT_MICROSPECIES_TABLE = CURATION_DATA_DIR / "metabolite_microspecies.csv"
REFERENCE_PH = 7.3

_ACTIVE_STATUS = "active"
_DEFERRED_STATUS = "component_review"
_VERIFIED_CURRENT_STATUS = "verified_current"
_NORMALIZE_STATUS = "normalize"
_ALLOWED_STATUSES = {
    _ACTIVE_STATUS,
    _DEFERRED_STATUS,
    _VERIFIED_CURRENT_STATUS,
    _NORMALIZE_STATUS,
}
_ALLOWED_SELECTORS = {"base_name", "bigg.metabolite", "metanetx.chemical"}
_FORMULA_RE = re.compile(r"^[A-Z][A-Za-z0-9]*$")
_CHEBI_RE = re.compile(r"^CHEBI:\d+$")
_FORMULA_SUFFIX_RE = re.compile(
    r"^(?:[A-Z][A-Za-z0-9\.]*|(?i:[pnq][+-]\d+|charge[+-]\d+))$",
)
_MNXM_CHEMISTRY_CONFLICT_NOTE = "metanetx_formula_charge_conflict"
_REQUIRED_COLUMNS = {
    "schema_version",
    "status",
    "family_id",
    "selector_type",
    "selector_value",
    "target_formula",
    "target_charge",
    "reference_ph",
    "chebi_id",
    "rhea_id",
    "source_url",
    "min_matches",
    "expected_metabolite_ids",
    "allowed_current_pairs",
    "rationale",
}


@dataclass(frozen=True)
class MicrospeciesRow:
    schema_version: int
    status: str
    family_id: str
    selector_type: str
    selector_value: str
    target_formula: str
    target_charge: int
    reference_ph: float
    chebi_id: str
    rhea_id: str
    source_url: str
    min_matches: int
    expected_metabolite_ids: frozenset[str]
    allowed_current_pairs: frozenset[tuple[str | None, int | None]]
    rationale: str

    @property
    def target_pair(self) -> tuple[str, int]:
        return self.target_formula, self.target_charge


def _parse_optional_charge(value: str) -> int | None:
    value = value.strip()
    if value in {"", "<missing>", "none", "None"}:
        return None
    return int(value)


def _parse_allowed_pairs(raw: str) -> frozenset[tuple[str | None, int | None]]:
    pairs: set[tuple[str | None, int | None]] = set()
    for token in raw.split(";"):
        token = token.strip()
        if not token:
            continue
        if "|" not in token:
            raise ValueError(
                f"invalid allowed_current_pairs token {token!r}; expected formula|charge"
            )
        formula_raw, charge_raw = token.split("|", 1)
        formula_raw = formula_raw.strip()
        formula = (
            None if formula_raw in {"", "<missing>", "none", "None"} else formula_raw
        )
        pairs.add((formula, _parse_optional_charge(charge_raw)))
    return frozenset(pairs)


def _parse_expected_ids(raw: str) -> frozenset[str]:
    values = [value.strip() for value in raw.split(";") if value.strip()]
    if len(values) != len(set(values)):
        raise ValueError("expected_metabolite_ids contains duplicates")
    return frozenset(values)


def load_curated_microspecies(
    table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> list[MicrospeciesRow]:
    """Load and fully validate the durable microspecies curation table."""

    path = Path(table_path)
    if not path.exists():
        raise FileNotFoundError(f"microspecies curation table not found: {path}")

    errors: list[str] = []
    rows: list[MicrospeciesRow] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        headers = set(reader.fieldnames or [])
        missing = sorted(_REQUIRED_COLUMNS - headers)
        if missing:
            raise ValueError(f"microspecies table lacks required columns: {missing}")

        seen_families: set[str] = set()
        for line_number, raw in enumerate(reader, start=2):
            prefix = f"row {line_number}"
            try:
                schema_version = int(raw["schema_version"].strip())
                status = raw["status"].strip()
                family_id = raw["family_id"].strip()
                selector_type = raw["selector_type"].strip()
                selector_value = raw["selector_value"].strip()
                target_formula = raw["target_formula"].strip()
                target_charge = int(raw["target_charge"].strip())
                reference_ph = float(raw["reference_ph"].strip())
                chebi_id = raw["chebi_id"].strip()
                rhea_id = raw["rhea_id"].strip()
                source_url = raw["source_url"].strip()
                min_matches = int(raw["min_matches"].strip())
                expected_ids = _parse_expected_ids(raw["expected_metabolite_ids"])
                allowed_pairs = _parse_allowed_pairs(raw["allowed_current_pairs"])
                rationale = raw["rationale"].strip()
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"{prefix}: cannot parse row ({exc})")
                continue

            if schema_version != 1:
                errors.append(f"{prefix}: unsupported schema_version {schema_version}")
            if status not in _ALLOWED_STATUSES:
                errors.append(f"{prefix}: unsupported status {status!r}")
            if not family_id:
                errors.append(f"{prefix}: family_id is required")
            elif family_id in seen_families:
                errors.append(f"{prefix}: duplicate family_id {family_id!r}")
            seen_families.add(family_id)
            if selector_type not in _ALLOWED_SELECTORS:
                errors.append(f"{prefix}: unsupported selector_type {selector_type!r}")
            if not selector_value:
                errors.append(f"{prefix}: selector_value is required")
            if not _FORMULA_RE.fullmatch(target_formula):
                errors.append(f"{prefix}: invalid exact formula {target_formula!r}")
            if not math.isclose(reference_ph, REFERENCE_PH, abs_tol=1e-12):
                errors.append(
                    f"{prefix}: reference_ph must be {REFERENCE_PH}, got {reference_ph}"
                )
            if not _CHEBI_RE.fullmatch(chebi_id):
                errors.append(f"{prefix}: invalid ChEBI identifier {chebi_id!r}")
            if rhea_id and not re.fullmatch(r"RHEA:\d+", rhea_id):
                errors.append(f"{prefix}: invalid Rhea identifier {rhea_id!r}")
            if not source_url.startswith(("https://", "http://")):
                errors.append(f"{prefix}: source_url must be HTTP(S)")
            if min_matches < 0:
                errors.append(f"{prefix}: min_matches cannot be negative")
            if min_matches != len(expected_ids):
                errors.append(
                    f"{prefix}: min_matches must equal the exact expected target "
                    f"count ({len(expected_ids)})"
                )
            if not allowed_pairs:
                errors.append(f"{prefix}: row requires allowed_current_pairs")
            if (target_formula, target_charge) not in allowed_pairs:
                errors.append(
                    f"{prefix}: row must allow its target pair for idempotence"
                )
            if not rationale:
                errors.append(f"{prefix}: rationale is required")

            rows.append(
                MicrospeciesRow(
                    schema_version=schema_version,
                    status=status,
                    family_id=family_id,
                    selector_type=selector_type,
                    selector_value=selector_value,
                    target_formula=target_formula,
                    target_charge=target_charge,
                    reference_ph=reference_ph,
                    chebi_id=chebi_id,
                    rhea_id=rhea_id,
                    source_url=source_url,
                    min_matches=min_matches,
                    expected_metabolite_ids=expected_ids,
                    allowed_current_pairs=allowed_pairs,
                    rationale=rationale,
                )
            )

    if errors:
        raise ValueError(
            "invalid microspecies curation table:\n- " + "\n- ".join(errors)
        )
    return rows


def _annotation_values(metabolite, key: str) -> set[str]:
    annotation = (
        metabolite.annotation if isinstance(metabolite.annotation, dict) else {}
    )
    value = annotation.get(key, [])
    if isinstance(value, (list, tuple, set)):
        return {str(item) for item in value}
    return {str(value)} if value not in (None, "") else set()


def metabolite_base_name(metabolite) -> str:
    """Return the compartment-independent chemical name used by table selectors."""

    name = (metabolite.name or "").strip().strip("_")
    if "_" in name:
        candidate, suffix = name.rsplit("_", 1)
        if _FORMULA_SUFFIX_RE.fullmatch(suffix.strip()):
            name = candidate.strip().strip("_")
    return name.casefold()


def resolve_microspecies_targets(model, row: MicrospeciesRow) -> list:
    """Resolve one exact family selector deterministically across compartments."""

    wanted = row.selector_value.casefold()
    matches = []
    for metabolite in model.metabolites:
        if row.selector_type == "base_name":
            hit = metabolite_base_name(metabolite) == wanted
        else:
            hit = row.selector_value in _annotation_values(
                metabolite, row.selector_type
            )
        if hit:
            matches.append(metabolite)
    return sorted(matches, key=lambda metabolite: metabolite.id)


def _resolve_pinned_targets(model, row: MicrospeciesRow) -> list:
    """Resolve a selector and require its full ID set to match the curated set."""

    selected = resolve_microspecies_targets(model, row)
    selected_ids = {metabolite.id for metabolite in selected}
    expected_ids = set(row.expected_metabolite_ids)
    missing_ids = sorted(expected_ids - {metabolite.id for metabolite in model.metabolites})
    if missing_ids:
        raise ValueError(
            f"{row.family_id}: expected metabolite IDs are absent: {missing_ids}"
        )
    if selected_ids != expected_ids:
        raise ValueError(
            f"{row.family_id}: selector target set changed; "
            f"expected={sorted(expected_ids)!r}, observed={sorted(selected_ids)!r}"
        )
    return [model.metabolites.get_by_id(mid) for mid in sorted(expected_ids)]


def _target_set_fingerprint(rows: Iterable[MicrospeciesRow]) -> str:
    payload = "\n".join(
        f"{row.family_id}\t{row.selector_type}\t{row.selector_value}\t"
        f"{';'.join(sorted(row.expected_metabolite_ids))}"
        for row in sorted(rows, key=lambda item: item.family_id)
    )
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _metabolite_pair(metabolite) -> tuple[str | None, int | None]:
    formula = metabolite.formula.strip() if metabolite.formula else None
    charge = int(metabolite.charge) if metabolite.charge is not None else None
    return formula, charge


def _is_allowed_current_pair(row: MicrospeciesRow, metabolite) -> bool:
    """Accept curated pairs or an audited MetaNetX conflict-preserved pair.

    Atomic MetaNetX enrichment deliberately leaves the name/source pair
    untouched when its formula conflicts with the complete MetaNetX pair.  A
    deferred family may recognize that state only when the structured note
    proves both the preserved pair and the exact curated target.  Active rows
    never receive this exception.
    """

    current = _metabolite_pair(metabolite)
    if current in row.allowed_current_pairs:
        return True
    if row.status != _DEFERRED_STATUS:
        return False
    raw_conflict = (metabolite.notes or {}).get(_MNXM_CHEMISTRY_CONFLICT_NOTE)
    if not isinstance(raw_conflict, str):
        return False
    try:
        conflict = json.loads(raw_conflict)
    except (TypeError, ValueError):
        return False
    preserved_pair = (
        conflict.get("existing_formula"),
        conflict.get("existing_charge"),
    )
    preserved_or_curated_legacy = preserved_pair == current or (
        current[1] == 0
        and current[0]
        in {
            formula
            for formula, _charge in row.allowed_current_pairs
            if formula is not None
        }
    )
    return (
        conflict.get("action") == "preserved_existing_pair"
        and preserved_or_curated_legacy
        and (
            conflict.get("proposed_formula"),
            conflict.get("proposed_charge"),
        )
        == row.target_pair
    )


def _is_fully_checkable(reaction) -> bool:
    return all(
        metabolite.formula and metabolite.charge is not None
        for metabolite in reaction.metabolites
    )


def _balance_without_noise(reaction, tolerance: float = 1e-9) -> dict[str, float]:
    try:
        balance = reaction.check_mass_balance()
    except Exception:
        return {"__error__": math.nan}
    return {
        key: float(value)
        for key, value in balance.items()
        if not math.isclose(float(value), 0.0, abs_tol=tolerance)
    }


def _fully_balanced_internal(reactions: Iterable) -> set[str]:
    balanced: set[str] = set()
    for reaction in reactions:
        if len(reaction.metabolites) <= 1 or not _is_fully_checkable(reaction):
            continue
        if not _balance_without_noise(reaction):
            balanced.add(reaction.id)
    return balanced


def apply_curated_microspecies(
    model,
    table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> dict:
    """Apply only ``status=active`` formula/charge pairs transactionally.

    Every row and selector is validated before mutation.  The batch is rolled
    back if it would make any previously fully balanced internal reaction
    unbalanced.  ``component_review`` rows are returned for audit but can never
    be activated through this function.  ``verified_current`` rows document a
    source-backed formula/charge pair that already matches the model; they are
    pinned during preflight but are never treated as a mutation candidate.
    """

    rows = load_curated_microspecies(table_path)
    active_rows = [row for row in rows if row.status == _ACTIVE_STATUS]
    deferred_rows = [row for row in rows if row.status == _DEFERRED_STATUS]
    verified_current_rows = [
        row for row in rows if row.status == _VERIFIED_CURRENT_STATUS
    ]

    resolved: dict[str, list] = {}
    errors: list[str] = []
    target_by_metabolite: dict[str, tuple[str, int, str]] = {}
    for row in rows:
        try:
            matches = _resolve_pinned_targets(model, row)
        except ValueError as exc:
            errors.append(str(exc))
            matches = []
        resolved[row.family_id] = matches
        for metabolite in matches:
            current = _metabolite_pair(metabolite)
            if not _is_allowed_current_pair(row, metabolite):
                errors.append(
                    f"{row.family_id}: {metabolite.id} has stale/unexpected pair "
                    f"{current!r}; allowed={sorted(row.allowed_current_pairs, key=str)!r}"
                )
            previous = target_by_metabolite.get(metabolite.id)
            proposed = (row.target_formula, row.target_charge, row.family_id)
            if previous and previous[:2] != proposed[:2]:
                errors.append(
                    f"{metabolite.id}: conflicting family targets {previous!r} and "
                    f"{proposed!r}"
                )
            target_by_metabolite[metabolite.id] = proposed

    if errors:
        raise ValueError("microspecies preflight failed:\n- " + "\n- ".join(errors))

    touched_reactions = {
        reaction
        for matches in resolved.values()
        for metabolite in matches
        for reaction in metabolite.reactions
    }
    balanced_before = _fully_balanced_internal(touched_reactions)
    snapshots = {
        metabolite.id: (
            metabolite.formula,
            metabolite.charge,
            copy.deepcopy(metabolite.annotation),
            copy.deepcopy(metabolite.notes),
        )
        for matches in resolved.values()
        for metabolite in matches
    }

    changed = 0
    already = 0
    audit_rows: list[dict] = []
    try:
        for row in active_rows:
            for metabolite in resolved[row.family_id]:
                before = _metabolite_pair(metabolite)
                if before == row.target_pair:
                    already += 1
                else:
                    # Formula and charge are intentionally assigned together.
                    metabolite.formula = row.target_formula
                    metabolite.charge = row.target_charge
                    changed += 1

                annotation = dict(metabolite.annotation or {})
                before_chebi = copy.deepcopy(annotation.get("chebi"))
                # ChEBI IDs identify exact chemical entities, so retaining a
                # neutral atom or a different protonation state would conflict
                # with the curated formula/charge pair.  Other namespaces are
                # preserved, but ChEBI is deliberately replaced by the exact
                # pH-convention identity.
                annotation["chebi"] = [row.chebi_id]
                metabolite.annotation = annotation
                notes = dict(metabolite.notes or {})
                notes["microspecies_convention"] = (
                    f"{row.family_id}: Rhea/ChEBI major microspecies at pH "
                    f"{row.reference_ph:g}; {row.source_url}"
                )
                metabolite.notes = notes
                audit_rows.append(
                    {
                        "family_id": row.family_id,
                        "metabolite_id": metabolite.id,
                        "before_formula": before[0],
                        "before_charge": before[1],
                        "before_chebi": before_chebi,
                        "after_formula": row.target_formula,
                        "after_charge": row.target_charge,
                        "after_chebi": row.chebi_id,
                        "source_url": row.source_url,
                    }
                )

        balanced_after = _fully_balanced_internal(touched_reactions)
        regressions = sorted(balanced_before - balanced_after)
        if regressions:
            raise ValueError(
                "microspecies batch would break previously balanced internal "
                f"reactions: {regressions}"
            )
    except Exception:
        for metabolite_id, snapshot in snapshots.items():
            metabolite = model.metabolites.get_by_id(metabolite_id)
            metabolite.formula = snapshot[0]
            metabolite.charge = snapshot[1]
            metabolite.annotation = snapshot[2]
            metabolite.notes = snapshot[3]
        raise

    deferred = []
    for row in deferred_rows:
        matches = resolved[row.family_id]
        deferred.append(
            {
                "family_id": row.family_id,
                "matched_metabolites": [metabolite.id for metabolite in matches],
                "target_formula": row.target_formula,
                "target_charge": row.target_charge,
                "would_change": [
                    metabolite.id
                    for metabolite in matches
                    if _metabolite_pair(metabolite) != row.target_pair
                ],
                "rationale": row.rationale,
            }
        )

    return {
        "reference_ph": REFERENCE_PH,
        "target_set_fingerprint": _target_set_fingerprint(rows),
        "changed_metabolites": changed,
        "already_canonical": already,
        "active_families": len(active_rows),
        "deferred_families": len(deferred_rows),
        "verified_current_families": len(verified_current_rows),
        "balanced_reaction_regressions": [],
        "changes": audit_rows,
        "deferred": deferred,
    }


def _reaction_balance_record(reaction, tolerance: float = 1e-9) -> dict:
    """Return one deterministic balance record without mutating a reaction."""

    metabolite_ids = sorted(metabolite.id for metabolite in reaction.metabolites)
    if len(reaction.metabolites) <= 1:
        return {
            "status": "boundary",
            "residual": {},
            "missing_metabolite_ids": [],
            "metabolite_ids": metabolite_ids,
        }

    missing = sorted(
        metabolite.id
        for metabolite in reaction.metabolites
        if not metabolite.formula or metabolite.charge is None
    )
    if missing:
        return {
            "status": "uncheckable",
            "residual": {},
            "missing_metabolite_ids": missing,
            "metabolite_ids": metabolite_ids,
        }

    residual = _balance_without_noise(reaction, tolerance)
    if "__error__" in residual:
        status = "error"
    elif residual:
        status = "imbalanced"
    else:
        status = "balanced"
    return {
        "status": status,
        "residual": {key: residual[key] for key in sorted(residual)},
        "missing_metabolite_ids": [],
        "metabolite_ids": metabolite_ids,
    }


def _component_closure(
    model,
    selected_rows: list[MicrospeciesRow],
    deferred_rows: list[MicrospeciesRow],
) -> tuple[list[MicrospeciesRow], list]:
    """Find deferred families joined through shared reactions.

    The closure deliberately traverses only curated ``component_review``
    families.  Traversing every ordinary metabolite would turn most of a GEM
    into one graph component and would not identify which additional
    microspecies decisions are actually available for review.
    """

    expected_id_to_row: dict[str, MicrospeciesRow] = {}
    for row in deferred_rows:
        for metabolite_id in row.expected_metabolite_ids:
            previous = expected_id_to_row.get(metabolite_id)
            if previous is not None and previous.family_id != row.family_id:
                raise ValueError(
                    f"{metabolite_id}: belongs to multiple component_review "
                    f"families ({previous.family_id!r}, {row.family_id!r})"
                )
            expected_id_to_row[metabolite_id] = row

    rows_by_id = {row.family_id: row for row in deferred_rows}
    closure_ids = {row.family_id for row in selected_rows}
    reaction_ids: set[str] = set()
    pending = sorted(closure_ids)
    while pending:
        family_id = pending.pop(0)
        row = rows_by_id[family_id]
        targets = _resolve_pinned_targets(model, row)
        for metabolite in targets:
            for reaction in metabolite.reactions:
                reaction_ids.add(reaction.id)
                for participant in reaction.metabolites:
                    neighbour = expected_id_to_row.get(participant.id)
                    if neighbour is None or neighbour.family_id in closure_ids:
                        continue
                    # Validate selector scope as soon as a neighbouring family
                    # enters the connected candidate component.
                    _resolve_pinned_targets(model, neighbour)
                    closure_ids.add(neighbour.family_id)
                    pending.append(neighbour.family_id)
        pending.sort()

    closure_rows = [rows_by_id[family_id] for family_id in sorted(closure_ids)]
    # Include every reaction touched by every family that entered the closure,
    # including reactions discovered by the final newly added family.
    reaction_ids = {
        reaction.id
        for row in closure_rows
        for metabolite in _resolve_pinned_targets(model, row)
        for reaction in metabolite.reactions
    }
    reactions = [model.reactions.get_by_id(rid) for rid in sorted(reaction_ids)]
    return closure_rows, reactions


def audit_component_migration(
    model,
    family_ids: Iterable[str] | None = None,
    *,
    tolerance: float = 1e-9,
    table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> dict:
    """Plan a deferred microspecies migration without changing ``model``.

    Only explicitly selected ``component_review`` rows are tentatively applied
    to an in-memory copy.  The report expands the audit scope to other deferred
    families connected through shared reactions, exposes the closure frontier,
    and fails readiness whenever an internal candidate reaction regresses,
    remains imbalanced, or cannot be checked.  The function never activates a
    curation row and never writes SBML.
    """

    rows = load_curated_microspecies(table_path)
    rows_by_id = {row.family_id: row for row in rows}
    deferred_rows = sorted(
        (row for row in rows if row.status == _DEFERRED_STATUS),
        key=lambda row: row.family_id,
    )

    if family_ids is None:
        selected_ids = [row.family_id for row in deferred_rows]
    else:
        selected_ids = sorted(set(family_ids))
        if not selected_ids:
            raise ValueError("component migration requires at least one family")

    selected_rows: list[MicrospeciesRow] = []
    for family_id in selected_ids:
        row = rows_by_id.get(family_id)
        if row is None:
            raise ValueError(f"unknown microspecies family {family_id!r}")
        if row.status != _DEFERRED_STATUS:
            raise ValueError(
                f"{family_id}: expected status component_review, got {row.status!r}"
            )
        selected_rows.append(row)

    selected_targets: dict[str, list] = {}
    current_pairs: dict[str, dict[str, object]] = {}
    for row in selected_rows:
        targets = _resolve_pinned_targets(model, row)
        selected_targets[row.family_id] = targets
        for metabolite in targets:
            current = _metabolite_pair(metabolite)
            if not _is_allowed_current_pair(row, metabolite):
                raise ValueError(
                    f"{row.family_id}: {metabolite.id} has stale/unexpected pair "
                    f"{current!r}; allowed="
                    f"{sorted(row.allowed_current_pairs, key=str)!r}"
                )
            current_pairs[metabolite.id] = {
                "family_id": row.family_id,
                "formula": current[0],
                "charge": current[1],
                "target_formula": row.target_formula,
                "target_charge": row.target_charge,
            }

    closure_rows, candidate_reactions = _component_closure(
        model, selected_rows, deferred_rows
    )
    closure_ids = [row.family_id for row in closure_rows]
    closure_metabolite_to_family = {
        metabolite_id: row.family_id
        for row in closure_rows
        for metabolite_id in row.expected_metabolite_ids
    }

    selected_metabolite_ids = {
        metabolite.id
        for targets in selected_targets.values()
        for metabolite in targets
    }
    initially_touched_ids = sorted(
        {
            reaction.id
            for targets in selected_targets.values()
            for metabolite in targets
            for reaction in metabolite.reactions
        }
    )

    before = {
        reaction.id: _reaction_balance_record(reaction, tolerance)
        for reaction in candidate_reactions
    }
    changed_metabolite_ids: list[str] = []
    snapshots = {
        metabolite.id: (metabolite.formula, metabolite.charge)
        for targets in selected_targets.values()
        for metabolite in targets
    }
    try:
        for row in selected_rows:
            for metabolite in selected_targets[row.family_id]:
                if _metabolite_pair(metabolite) != row.target_pair:
                    metabolite.formula = row.target_formula
                    metabolite.charge = row.target_charge
                    changed_metabolite_ids.append(metabolite.id)

        after = {
            reaction.id: _reaction_balance_record(reaction, tolerance)
            for reaction in candidate_reactions
        }
    finally:
        for metabolite_id, (formula, charge) in snapshots.items():
            metabolite = model.metabolites.get_by_id(metabolite_id)
            metabolite.formula = formula
            metabolite.charge = charge
    changed_metabolite_ids.sort()

    regressions = sorted(
        reaction_id
        for reaction_id in before
        if before[reaction_id]["status"] == "balanced"
        and after[reaction_id]["status"] != "balanced"
    )
    fixed = sorted(
        reaction_id
        for reaction_id in before
        if before[reaction_id]["status"] != "balanced"
        and after[reaction_id]["status"] == "balanced"
    )
    unresolved = sorted(
        reaction_id
        for reaction_id, record in after.items()
        if record["status"] in {"imbalanced", "uncheckable", "error"}
    )

    selected_id_set = set(selected_ids)
    problem_ids = sorted(set(regressions) | set(unresolved))
    frontier: list[dict] = []
    frontier_family_ids: set[str] = set()
    frontier_metabolite_ids: set[str] = set()
    for reaction_id in problem_ids:
        reaction = model.reactions.get_by_id(reaction_id)
        outside_families = sorted(
            {
                closure_metabolite_to_family[metabolite.id]
                for metabolite in reaction.metabolites
                if metabolite.id in closure_metabolite_to_family
                and closure_metabolite_to_family[metabolite.id]
                not in selected_id_set
            }
        )
        outside_metabolites = sorted(
            metabolite.id
            for metabolite in reaction.metabolites
            if metabolite.id not in selected_metabolite_ids
        )
        frontier_family_ids.update(outside_families)
        frontier_metabolite_ids.update(outside_metabolites)
        frontier.append(
            {
                "reaction_id": reaction_id,
                "unselected_component_family_ids": outside_families,
                "outside_selected_metabolite_ids": outside_metabolites,
            }
        )

    reaction_context = []
    for reaction in candidate_reactions:
        reaction_context.append(
            {
                "reaction_id": reaction.id,
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
                "stoichiometry": [
                    [metabolite.id, float(coefficient)]
                    for metabolite, coefficient in sorted(
                        reaction.metabolites.items(), key=lambda item: item[0].id
                    )
                ],
                "before": before[reaction.id],
            }
        )
    fingerprint_payload = {
        "schema": "component_migration_plan_v1",
        "reference_ph": REFERENCE_PH,
        "selected_family_ids": selected_ids,
        "closure_family_ids": closure_ids,
        "current_pairs": current_pairs,
        "reaction_context": reaction_context,
    }
    plan_fingerprint = "sha256:" + hashlib.sha256(
        json.dumps(
            fingerprint_payload,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()

    reaction_balances = [
        {
            "reaction_id": reaction_id,
            "before": before[reaction_id],
            "after": after[reaction_id],
        }
        for reaction_id in sorted(before)
    ]
    return {
        "schema_version": 1,
        "reference_ph": REFERENCE_PH,
        "plan_fingerprint": plan_fingerprint,
        "selected_family_ids": selected_ids,
        "closure_family_ids": closure_ids,
        "changed_metabolite_ids": changed_metabolite_ids,
        "initially_touched_reaction_ids": initially_touched_ids,
        "candidate_reaction_ids": sorted(before),
        "fixed_reaction_ids": fixed,
        "regressed_reaction_ids": regressions,
        "unresolved_reaction_ids": unresolved,
        "frontier_family_ids": sorted(frontier_family_ids),
        "frontier_metabolite_ids": sorted(frontier_metabolite_ids),
        "closure_frontier": frontier,
        "reaction_balances": reaction_balances,
        "ready_for_activation": not regressions and not unresolved,
    }


def _sha256_path(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_component_migration_audit(
    model_path: str | Path,
    output_path: str | Path,
    family_ids: Iterable[str],
    *,
    table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
    tolerance: float = 1e-9,
) -> Path:
    """Write one provenance-rich component audit with an atomic replace.

    The curation table is read-only in this operation.  A failed serialization,
    flush, or replacement leaves any existing audit byte-for-byte unchanged.
    """

    from cobra.io import read_sbml_model

    model_path = Path(model_path)
    output_path = Path(output_path)
    table_path = Path(table_path)
    selected_family_ids = sorted(set(family_ids))
    model = read_sbml_model(str(model_path))
    audit = audit_component_migration(
        model,
        selected_family_ids,
        tolerance=tolerance,
        table_path=table_path,
    )
    payload = {
        "schema_version": 1,
        "generated_at": dt.datetime.now(dt.UTC).isoformat(),
        "provenance": {
            "model_path": str(model_path.resolve()),
            "model_sha256": _sha256_path(model_path),
            "microspecies_table_path": str(table_path.resolve()),
            "microspecies_table_sha256": _sha256_path(table_path),
        },
        "selected_family_ids": selected_family_ids,
        "audit": audit,
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    mode = output_path.stat().st_mode & 0o777 if output_path.exists() else 0o644
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
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary_path, mode)
        os.replace(temporary_path, output_path)
    except Exception:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise
    return output_path


def _family_row(rows: list[MicrospeciesRow], family_id: str) -> MicrospeciesRow:
    matches = [row for row in rows if row.family_id == family_id]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one microspecies row for {family_id!r}")
    return matches[0]


def normalize_hydroxide_reactions(
    model,
    table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> dict:
    """Replace exact OH- with H2O - H+ in safe internal reactions.

    Boundary reactions and reactions spanning more than one compartment are
    reported but never modified because choosing a membrane side would invent
    proton coupling.
    """

    rows = load_curated_microspecies(table_path)
    hydroxide_row = _family_row(rows, "hydroxide")
    proton_row = _family_row(rows, "proton")
    water_row = _family_row(rows, "water")
    hydroxides = _resolve_pinned_targets(model, hydroxide_row)
    protons = {
        metabolite.compartment: metabolite
        for metabolite in _resolve_pinned_targets(model, proton_row)
        if _metabolite_pair(metabolite) == proton_row.target_pair
    }
    waters = {
        metabolite.compartment: metabolite
        for metabolite in _resolve_pinned_targets(model, water_row)
        if _metabolite_pair(metabolite) == water_row.target_pair
    }

    changed: list[str] = []
    rejected: list[dict] = []
    reaction_hydroxides: dict[object, list] = {}
    for hydroxide in hydroxides:
        if _metabolite_pair(hydroxide) != hydroxide_row.target_pair:
            rejected.append({"metabolite_id": hydroxide.id, "reason": "not_exact_OH-"})
            continue
        for reaction in sorted(hydroxide.reactions, key=lambda item: item.id):
            reaction_hydroxides.setdefault(reaction, []).append(hydroxide)

    for reaction in sorted(reaction_hydroxides, key=lambda item: item.id):
        compartments = {metabolite.compartment for metabolite in reaction.metabolites}
        if len(reaction.metabolites) <= 1:
            rejected.append({"reaction_id": reaction.id, "reason": "boundary"})
            continue
        if len(compartments) != 1:
            rejected.append(
                {"reaction_id": reaction.id, "reason": "multi_compartment"}
            )
            continue
        compartment = next(iter(compartments))
        proton = protons.get(compartment)
        water = waters.get(compartment)
        if proton is None or water is None:
            rejected.append(
                {"reaction_id": reaction.id, "reason": "missing_H+_or_H2O"}
            )
            continue

        additions: dict = {}
        for hydroxide in reaction_hydroxides[reaction]:
            coefficient = float(reaction.metabolites.get(hydroxide, 0.0))
            if math.isclose(coefficient, 0.0, abs_tol=1e-12):
                continue
            additions[hydroxide] = -coefficient
            additions[water] = additions.get(water, 0.0) + coefficient
            additions[proton] = additions.get(proton, 0.0) - coefficient
        if not additions:
            continue

        before = dict(reaction.metabolites)
        try:
            reaction.add_metabolites(additions)
            changed.append(reaction.id)
        except Exception:
            current = dict(reaction.metabolites)
            reaction.add_metabolites(
                {metabolite: -value for metabolite, value in current.items()}
            )
            reaction.add_metabolites(before)
            raise

    return {
        "changed_reactions": len(changed),
        "rejected_reactions": len(rejected),
        "changes": changed,
        "rejected": rejected,
    }


def balance_protons_and_water(
    model,
    reaction_ids: Iterable[str] | None = None,
    *,
    tolerance: float = 1e-9,
    table_path: str | Path = DEFAULT_MICROSPECIES_TABLE,
) -> dict:
    """Balance eligible reactions with H+/H2O while respecting charge.

    If the current residuals are ``H=h``, ``O=o`` and ``charge=q``, adding
    ``-o H2O`` and ``-q H+`` is valid exactly when ``h - 2o == q``.  Heavy
    element residuals, boundary reactions and multi-compartment reactions are
    never modified.
    """

    rows = load_curated_microspecies(table_path)
    proton_row = _family_row(rows, "proton")
    water_row = _family_row(rows, "water")
    protons = {
        metabolite.compartment: metabolite
        for metabolite in _resolve_pinned_targets(model, proton_row)
        if _metabolite_pair(metabolite) == proton_row.target_pair
    }
    waters = {
        metabolite.compartment: metabolite
        for metabolite in _resolve_pinned_targets(model, water_row)
        if _metabolite_pair(metabolite) == water_row.target_pair
    }

    if reaction_ids is None:
        reactions = sorted(model.reactions, key=lambda reaction: reaction.id)
    else:
        reactions = []
        for reaction_id in sorted(set(reaction_ids)):
            try:
                reactions.append(model.reactions.get_by_id(reaction_id))
            except KeyError as exc:
                raise ValueError(
                    f"unknown reaction for proton balance: {reaction_id}"
                ) from exc

    changes: list[dict] = []
    rejected: list[dict] = []
    skipped_missing_formula = 0
    skipped_curated_locks: list[str] = []
    for reaction in reactions:
        if len(reaction.metabolites) <= 1:
            continue
        curated_correction = (reaction.notes or {}).get(
            "curated_reaction_correction"
        )
        if (
            isinstance(curated_correction, dict)
            and curated_correction.get("lock_proton_water_stoichiometry") is True
        ):
            skipped_curated_locks.append(reaction.id)
            continue
        if not _is_fully_checkable(reaction):
            skipped_missing_formula += 1
            continue
        balance = _balance_without_noise(reaction, tolerance)
        if not balance:
            continue
        if "__error__" in balance:
            rejected.append({"reaction_id": reaction.id, "reason": "balance_error"})
            continue
        heavy = {
            key: value
            for key, value in balance.items()
            if key not in {"H", "O", "charge"}
        }
        if heavy:
            rejected.append(
                {
                    "reaction_id": reaction.id,
                    "reason": "heavy_elements",
                    "balance": balance,
                }
            )
            continue
        compartments = {metabolite.compartment for metabolite in reaction.metabolites}
        if len(compartments) != 1:
            rejected.append(
                {
                    "reaction_id": reaction.id,
                    "reason": "multi_compartment",
                    "balance": balance,
                }
            )
            continue

        h_residual = balance.get("H", 0.0)
        o_residual = balance.get("O", 0.0)
        charge_residual = balance.get("charge", 0.0)
        if not math.isclose(
            h_residual - 2.0 * o_residual,
            charge_residual,
            abs_tol=tolerance,
        ):
            rejected.append(
                {
                    "reaction_id": reaction.id,
                    "reason": "proton_gate",
                    "balance": balance,
                }
            )
            continue

        compartment = next(iter(compartments))
        proton = protons.get(compartment)
        water = waters.get(compartment)
        if (not math.isclose(o_residual, 0.0, abs_tol=tolerance) and water is None) or (
            not math.isclose(charge_residual, 0.0, abs_tol=tolerance) and proton is None
        ):
            rejected.append(
                {
                    "reaction_id": reaction.id,
                    "reason": "missing_H+_or_H2O",
                    "balance": balance,
                }
            )
            continue

        additions = {}
        if not math.isclose(o_residual, 0.0, abs_tol=tolerance):
            additions[water] = -o_residual
        if not math.isclose(charge_residual, 0.0, abs_tol=tolerance):
            additions[proton] = -charge_residual
        if not additions:
            continue

        reaction.add_metabolites(additions)
        after = _balance_without_noise(reaction, tolerance)
        if after:
            reaction.add_metabolites(
                {
                    metabolite: -coefficient
                    for metabolite, coefficient in additions.items()
                }
            )
            rejected.append(
                {
                    "reaction_id": reaction.id,
                    "reason": "postcondition_failed",
                    "balance": balance,
                    "after": after,
                }
            )
            continue
        changes.append(
            {
                "reaction_id": reaction.id,
                "before": balance,
                "additions": {
                    metabolite.id: coefficient
                    for metabolite, coefficient in additions.items()
                },
            }
        )

    return {
        "changed_reactions": len(changes),
        "rejected_reactions": len(rejected),
        "skipped_missing_formula": skipped_missing_formula,
        "skipped_curated_lock_reaction_ids": skipped_curated_locks,
        "changes": changes,
        "rejected": rejected,
    }


__all__ = [
    "DEFAULT_MICROSPECIES_TABLE",
    "REFERENCE_PH",
    "MicrospeciesRow",
    "apply_curated_microspecies",
    "audit_component_migration",
    "balance_protons_and_water",
    "load_curated_microspecies",
    "metabolite_base_name",
    "normalize_hydroxide_reactions",
    "resolve_microspecies_targets",
    "write_component_migration_audit",
]
