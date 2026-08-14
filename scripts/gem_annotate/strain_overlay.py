"""Strict, in-memory strain overlays for PO1f simulations.

This module is deliberately separate from the canonical model-build pipeline.
Callers should apply the experimental medium to a disposable model copy and
then call :func:`apply_strain_overlay`.  The supplied model is modified in
memory; no SBML file is read or written here.

The uracil concentration-ratio calculation exposed by this module is a
formulation-derived surrogate.  It is not a measured uptake rate or a measured
``Vmax``.  A profile can therefore retain the formulation calculation as batch
supply evidence while using a separate, explicitly labelled non-limiting bound
for static FBA.
"""

from __future__ import annotations

import copy
import json
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any


PROFILE_SCHEMA_VERSION = "1.0"
PROFILE_STRAIN = "PO1f"
PROFILE_REFERENCE_BACKGROUND = "W29/CLIB89"

URA3_REACTION_ID = "R612"
URA3_LOCUS = "URA3"
URA3_GENE_ID = "YALI1E31685g"
URA3_LEGACY_GENE_ID = "YALI0E26741g"
URA3_ALLELE = "ura3-302"
URA3_EXPECTED_BOUNDS = (0.0, 1000.0)
URA3_DISABLED_BOUNDS = (0.0, 0.0)

LEU2_REACTION_ID = "R45"
LEU2_LOCUS = "LEU2"
LEU2_GENE_ID = "YALI1C00464g"
LEU2_LEGACY_GENE_ID = "YALI0C00407g"
LEU2_ALLELE = "leu2-270"
LEU2_EXPECTED_GPR = "YALI1C00464g"
LEU2_PLASMID_PSEUDO_GENE = "PO1f_plasmid_LEU2"

LEUCINE_EXCHANGE_ID = "R1219"
URACIL_EXCHANGE_ID = "R1354"

CONCENTRATION_RATIO_STATUS = "concentration_ratio_surrogate_not_measured_vmax"
RUNTIME_OVERRIDE_STATUS = "static_fba_nonlimiting_not_measured_vmax"

_TOP_LEVEL_KEYS = {
    "schema_version",
    "profile_id",
    "strain",
    "operations",
    "medium",
    "assay_confounded_loci",
    "provenance_only_variants",
    "sources",
}
_STRAIN_KEYS = {
    "name",
    "reference_background",
    "genotype",
    "assay_background",
}
_DISABLE_OPERATION_KEYS = {
    "type",
    "reaction_id",
    "expected_before",
    "set_bounds",
    "locus",
    "gene_id",
    "legacy_gene_id",
    "allele",
    "protein_function",
    "evidence_status",
}
_PLASMID_OPERATION_KEYS = {
    "type",
    "reaction_id",
    "expected_before",
    "set_gene_reaction_rule",
    "pseudo_gene",
    "locus",
    "gene_id",
    "legacy_gene_id",
    "allele",
    "complement_source",
    "protein_function",
    "evidence_status",
}
_BOUNDS_KEYS = {"lower_bound", "upper_bound"}
_EXPECTED_GPR_KEYS = {"gene_reaction_rule"}
_MEDIUM_KEYS = {
    "uptake_assertions",
    "runtime_uptake_overrides",
    "formulation",
}
_FORMULATION_KEYS = {
    "uracil_mg_per_l",
    "uracil_millimolar",
    "uracil_molecular_weight_g_per_mol",
    "glucose_g_per_l",
    "glucose_molecular_weight_g_per_mol",
    "glucose_uptake",
    "supply_ratio_surrogate_uptake",
    "supply_ratio_surrogate_status",
    "runtime_override_status",
}


class StrainProfileError(ValueError):
    """Raised when a strain profile fails its closed schema or invariants."""


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise StrainProfileError(f"Duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_nonfinite_json_constant(value: str) -> None:
    raise StrainProfileError(f"Non-finite JSON number {value!r} is not allowed")


def _require_object(value: Any, where: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise StrainProfileError(f"{where} must be a JSON object")
    return value


def _require_exact_keys(
    value: dict[str, Any],
    expected: set[str],
    where: str,
) -> None:
    actual = set(value)
    missing = sorted(expected - actual)
    unknown = sorted(actual - expected)
    if missing or unknown:
        details: list[str] = []
        if missing:
            details.append(f"missing keys {missing}")
        if unknown:
            details.append(f"unknown keys {unknown}")
        raise StrainProfileError(f"{where} has " + " and ".join(details))


def _require_nonempty_string(value: Any, where: str) -> str:
    if type(value) is not str or not value.strip():
        raise StrainProfileError(f"{where} must be a non-empty string")
    if value != value.strip():
        raise StrainProfileError(f"{where} must not have surrounding whitespace")
    return value


def _require_number(
    value: Any,
    where: str,
    *,
    minimum: float | None = None,
    strictly_positive: bool = False,
) -> float:
    if isinstance(value, bool) or type(value) not in {int, float}:
        raise StrainProfileError(f"{where} must be a finite JSON number")
    result = float(value)
    if not math.isfinite(result):
        raise StrainProfileError(f"{where} must be finite")
    if strictly_positive and result <= 0.0:
        raise StrainProfileError(f"{where} must be greater than zero")
    if minimum is not None and result < minimum:
        raise StrainProfileError(f"{where} must be at least {minimum}")
    return result


def _require_string_list(value: Any, where: str) -> list[str]:
    if type(value) is not list:
        raise StrainProfileError(f"{where} must be a JSON array of strings")
    items = [
        _require_nonempty_string(item, f"{where}[{index}]")
        for index, item in enumerate(value)
    ]
    if len(items) != len(set(items)):
        raise StrainProfileError(f"{where} must not contain duplicates")
    return items


def _require_close(
    actual: float,
    expected: float,
    where: str,
    *,
    rel_tol: float = 1e-9,
    abs_tol: float = 1e-12,
) -> None:
    if not math.isclose(actual, expected, rel_tol=rel_tol, abs_tol=abs_tol):
        raise StrainProfileError(f"{where} is {actual!r}; expected {expected!r}")


def _validate_bounds(value: Any, where: str) -> tuple[float, float]:
    bounds = _require_object(value, where)
    _require_exact_keys(bounds, _BOUNDS_KEYS, where)
    lower = _require_number(bounds["lower_bound"], f"{where}.lower_bound")
    upper = _require_number(bounds["upper_bound"], f"{where}.upper_bound")
    if lower > upper:
        raise StrainProfileError(f"{where}.lower_bound must not exceed upper_bound")
    return lower, upper


def derive_concentration_ratio_uptake(
    uracil_mg_per_l: float,
    uracil_molecular_weight_g_per_mol: float,
    glucose_g_per_l: float,
    glucose_uptake: float,
    *,
    glucose_molecular_weight_g_per_mol: float = 180.156,
) -> float:
    """Return a formulation-derived uracil uptake surrogate.

    The calculation converts uracil and glucose concentrations to mmol/L and
    scales ``glucose_uptake`` by their concentration ratio::

        glucose_uptake * uracil_mM / glucose_mM

    This is a pure unit-and-ratio calculation.  The result is **not** an
    experimentally measured transport rate, kinetic limit, or ``Vmax``.
    """

    uracil_concentration = _require_number(
        uracil_mg_per_l,
        "uracil_mg_per_l",
        minimum=0.0,
    )
    uracil_mw = _require_number(
        uracil_molecular_weight_g_per_mol,
        "uracil_molecular_weight_g_per_mol",
        strictly_positive=True,
    )
    glucose_concentration = _require_number(
        glucose_g_per_l,
        "glucose_g_per_l",
        strictly_positive=True,
    )
    glucose_bound = _require_number(
        glucose_uptake,
        "glucose_uptake",
        minimum=0.0,
    )
    glucose_mw = _require_number(
        glucose_molecular_weight_g_per_mol,
        "glucose_molecular_weight_g_per_mol",
        strictly_positive=True,
    )

    uracil_millimolar = uracil_concentration / uracil_mw
    glucose_millimolar = 1000.0 * glucose_concentration / glucose_mw
    return float(glucose_bound * uracil_millimolar / glucose_millimolar)


def _validate_disable_operation(operation: dict[str, Any]) -> None:
    where = "profile.operations[disable_reaction]"
    _require_exact_keys(operation, _DISABLE_OPERATION_KEYS, where)
    if operation["type"] != "disable_reaction":
        raise StrainProfileError(f"{where}.type must be 'disable_reaction'")
    if operation["reaction_id"] != URA3_REACTION_ID:
        raise StrainProfileError(f"{where}.reaction_id must be {URA3_REACTION_ID!r}")
    expected = _validate_bounds(
        operation["expected_before"],
        f"{where}.expected_before",
    )
    target = _validate_bounds(operation["set_bounds"], f"{where}.set_bounds")
    if expected != URA3_EXPECTED_BOUNDS:
        raise StrainProfileError(
            f"{where}.expected_before must be R612 bounds {URA3_EXPECTED_BOUNDS}"
        )
    if target != URA3_DISABLED_BOUNDS:
        raise StrainProfileError(
            f"{where}.set_bounds must disable R612 with bounds {URA3_DISABLED_BOUNDS}"
        )
    if operation["locus"] != URA3_LOCUS:
        raise StrainProfileError(f"{where}.locus must be {URA3_LOCUS!r}")
    if operation["gene_id"] != URA3_GENE_ID:
        raise StrainProfileError(f"{where}.gene_id must be {URA3_GENE_ID!r}")
    if operation["legacy_gene_id"] != URA3_LEGACY_GENE_ID:
        raise StrainProfileError(
            f"{where}.legacy_gene_id must be {URA3_LEGACY_GENE_ID!r}"
        )
    if operation["allele"] != URA3_ALLELE:
        raise StrainProfileError(f"{where}.allele must be {URA3_ALLELE!r}")
    _require_nonempty_string(
        operation["protein_function"],
        f"{where}.protein_function",
    )
    _require_nonempty_string(
        operation["evidence_status"],
        f"{where}.evidence_status",
    )


def _validate_plasmid_operation(operation: dict[str, Any]) -> None:
    where = "profile.operations[plasmid_complement]"
    _require_exact_keys(operation, _PLASMID_OPERATION_KEYS, where)
    if operation["type"] != "plasmid_complement":
        raise StrainProfileError(f"{where}.type must be 'plasmid_complement'")
    if operation["reaction_id"] != LEU2_REACTION_ID:
        raise StrainProfileError(f"{where}.reaction_id must be {LEU2_REACTION_ID!r}")
    expected = _require_object(
        operation["expected_before"],
        f"{where}.expected_before",
    )
    _require_exact_keys(
        expected,
        _EXPECTED_GPR_KEYS,
        f"{where}.expected_before",
    )
    if expected["gene_reaction_rule"] != LEU2_EXPECTED_GPR:
        raise StrainProfileError(
            f"{where}.expected_before.gene_reaction_rule must be {LEU2_EXPECTED_GPR!r}"
        )
    if operation["set_gene_reaction_rule"] != LEU2_PLASMID_PSEUDO_GENE:
        raise StrainProfileError(
            f"{where}.set_gene_reaction_rule must be {LEU2_PLASMID_PSEUDO_GENE!r}"
        )
    if operation["pseudo_gene"] != LEU2_PLASMID_PSEUDO_GENE:
        raise StrainProfileError(
            f"{where}.pseudo_gene must be {LEU2_PLASMID_PSEUDO_GENE!r}"
        )
    if operation["locus"] != LEU2_LOCUS:
        raise StrainProfileError(f"{where}.locus must be {LEU2_LOCUS!r}")
    if operation["gene_id"] != LEU2_GENE_ID:
        raise StrainProfileError(f"{where}.gene_id must be {LEU2_GENE_ID!r}")
    if operation["legacy_gene_id"] != LEU2_LEGACY_GENE_ID:
        raise StrainProfileError(
            f"{where}.legacy_gene_id must be {LEU2_LEGACY_GENE_ID!r}"
        )
    if operation["allele"] != LEU2_ALLELE:
        raise StrainProfileError(f"{where}.allele must be {LEU2_ALLELE!r}")
    _require_nonempty_string(
        operation["complement_source"],
        f"{where}.complement_source",
    )
    _require_nonempty_string(
        operation["protein_function"],
        f"{where}.protein_function",
    )
    _require_nonempty_string(
        operation["evidence_status"],
        f"{where}.evidence_status",
    )


def _validate_medium(medium_value: Any) -> None:
    medium = _require_object(medium_value, "profile.medium")
    _require_exact_keys(medium, _MEDIUM_KEYS, "profile.medium")

    assertions = _require_object(
        medium["uptake_assertions"],
        "profile.medium.uptake_assertions",
    )
    _require_exact_keys(
        assertions,
        {LEUCINE_EXCHANGE_ID, URACIL_EXCHANGE_ID},
        "profile.medium.uptake_assertions",
    )
    leucine_uptake = _require_number(
        assertions[LEUCINE_EXCHANGE_ID],
        f"profile.medium.uptake_assertions.{LEUCINE_EXCHANGE_ID}",
        minimum=0.0,
    )
    uracil_assertion = _require_number(
        assertions[URACIL_EXCHANGE_ID],
        f"profile.medium.uptake_assertions.{URACIL_EXCHANGE_ID}",
        minimum=0.0,
    )
    _require_close(
        leucine_uptake,
        0.0,
        f"profile.medium.uptake_assertions.{LEUCINE_EXCHANGE_ID}",
    )

    overrides = _require_object(
        medium["runtime_uptake_overrides"],
        "profile.medium.runtime_uptake_overrides",
    )
    _require_exact_keys(
        overrides,
        {URACIL_EXCHANGE_ID},
        "profile.medium.runtime_uptake_overrides",
    )
    _require_number(
        overrides[URACIL_EXCHANGE_ID],
        f"profile.medium.runtime_uptake_overrides.{URACIL_EXCHANGE_ID}",
        minimum=0.0,
    )

    formulation = _require_object(
        medium["formulation"],
        "profile.medium.formulation",
    )
    _require_exact_keys(
        formulation,
        _FORMULATION_KEYS,
        "profile.medium.formulation",
    )
    uracil_mg_per_l = _require_number(
        formulation["uracil_mg_per_l"],
        "profile.medium.formulation.uracil_mg_per_l",
        strictly_positive=True,
    )
    _require_close(
        uracil_mg_per_l,
        20.0,
        "profile.medium.formulation.uracil_mg_per_l",
    )
    uracil_mw = _require_number(
        formulation["uracil_molecular_weight_g_per_mol"],
        "profile.medium.formulation.uracil_molecular_weight_g_per_mol",
        strictly_positive=True,
    )
    glucose_g_per_l = _require_number(
        formulation["glucose_g_per_l"],
        "profile.medium.formulation.glucose_g_per_l",
        strictly_positive=True,
    )
    glucose_mw = _require_number(
        formulation["glucose_molecular_weight_g_per_mol"],
        "profile.medium.formulation.glucose_molecular_weight_g_per_mol",
        strictly_positive=True,
    )
    glucose_uptake = _require_number(
        formulation["glucose_uptake"],
        "profile.medium.formulation.glucose_uptake",
        minimum=0.0,
    )
    stated_uracil_millimolar = _require_number(
        formulation["uracil_millimolar"],
        "profile.medium.formulation.uracil_millimolar",
        minimum=0.0,
    )
    calculated_uracil_millimolar = uracil_mg_per_l / uracil_mw
    _require_close(
        stated_uracil_millimolar,
        calculated_uracil_millimolar,
        "profile.medium.formulation.uracil_millimolar",
        rel_tol=5e-4,
        abs_tol=5e-6,
    )

    stated_surrogate = _require_number(
        formulation["supply_ratio_surrogate_uptake"],
        "profile.medium.formulation.supply_ratio_surrogate_uptake",
        minimum=0.0,
    )
    calculated_surrogate = derive_concentration_ratio_uptake(
        uracil_mg_per_l,
        uracil_mw,
        glucose_g_per_l,
        glucose_uptake,
        glucose_molecular_weight_g_per_mol=glucose_mw,
    )
    _require_close(
        stated_surrogate,
        calculated_surrogate,
        "profile.medium.formulation.supply_ratio_surrogate_uptake",
        rel_tol=5e-4,
        abs_tol=5e-6,
    )
    _require_close(
        uracil_assertion,
        stated_surrogate,
        f"profile.medium.uptake_assertions.{URACIL_EXCHANGE_ID}",
        rel_tol=5e-4,
        abs_tol=5e-6,
    )
    if formulation["supply_ratio_surrogate_status"] != CONCENTRATION_RATIO_STATUS:
        raise StrainProfileError(
            "profile.medium.formulation.supply_ratio_surrogate_status must be "
            f"{CONCENTRATION_RATIO_STATUS!r}"
        )
    if formulation["runtime_override_status"] != RUNTIME_OVERRIDE_STATUS:
        raise StrainProfileError(
            "profile.medium.formulation.runtime_override_status must be "
            f"{RUNTIME_OVERRIDE_STATUS!r}"
        )


def _validate_profile(profile_value: Any) -> dict[str, Any]:
    profile = _require_object(profile_value, "profile")
    _require_exact_keys(profile, _TOP_LEVEL_KEYS, "profile")
    if profile["schema_version"] != PROFILE_SCHEMA_VERSION:
        raise StrainProfileError(
            "profile.schema_version must be "
            f"{PROFILE_SCHEMA_VERSION!r}, not {profile['schema_version']!r}"
        )
    _require_nonempty_string(profile["profile_id"], "profile.profile_id")
    strain = _require_object(profile["strain"], "profile.strain")
    _require_exact_keys(strain, _STRAIN_KEYS, "profile.strain")
    if strain["name"] != PROFILE_STRAIN:
        raise StrainProfileError(
            f"profile.strain.name must be {PROFILE_STRAIN!r}"
        )
    if strain["reference_background"] != PROFILE_REFERENCE_BACKGROUND:
        raise StrainProfileError(
            "profile.strain.reference_background must be "
            f"{PROFILE_REFERENCE_BACKGROUND!r}"
        )
    _require_nonempty_string(strain["genotype"], "profile.strain.genotype")
    _require_nonempty_string(
        strain["assay_background"],
        "profile.strain.assay_background",
    )

    operations = profile["operations"]
    if type(operations) is not list:
        raise StrainProfileError("profile.operations must be a JSON array")
    if len(operations) != 2:
        raise StrainProfileError(
            "profile.operations must contain exactly one disable_reaction and "
            "one plasmid_complement operation"
        )
    by_type: dict[str, dict[str, Any]] = {}
    for index, operation_value in enumerate(operations):
        operation = _require_object(
            operation_value,
            f"profile.operations[{index}]",
        )
        operation_type = operation.get("type")
        if operation_type not in {"disable_reaction", "plasmid_complement"}:
            raise StrainProfileError(
                f"profile.operations[{index}].type is unsupported: {operation_type!r}"
            )
        if operation_type in by_type:
            raise StrainProfileError(
                f"profile.operations contains duplicate type {operation_type!r}"
            )
        by_type[operation_type] = operation
    if set(by_type) != {"disable_reaction", "plasmid_complement"}:
        raise StrainProfileError(
            "profile.operations must contain disable_reaction and plasmid_complement"
        )
    _validate_disable_operation(by_type["disable_reaction"])
    _validate_plasmid_operation(by_type["plasmid_complement"])
    _validate_medium(profile["medium"])

    confounded = _require_string_list(
        profile["assay_confounded_loci"],
        "profile.assay_confounded_loci",
    )
    if set(confounded) != {URA3_LOCUS, LEU2_LOCUS}:
        raise StrainProfileError(
            "profile.assay_confounded_loci must contain exactly URA3 and LEU2"
        )
    _require_string_list(
        profile["provenance_only_variants"],
        "profile.provenance_only_variants",
    )
    _require_string_list(profile["sources"], "profile.sources")
    return profile


def load_strain_profile(path: str | Path) -> dict[str, Any]:
    """Read and validate one strict PO1f JSON profile.

    Duplicate keys, non-standard non-finite numbers, unknown fields, missing
    fields, wrong JSON types, and inconsistent formulation calculations are
    rejected.  The returned dictionary contains only the profile data.
    """

    profile_path = Path(path)
    try:
        text = profile_path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise StrainProfileError(
            f"Strain profile is not valid UTF-8: {profile_path}"
        ) from exc
    try:
        value = json.loads(
            text,
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_nonfinite_json_constant,
        )
    except json.JSONDecodeError as exc:
        raise StrainProfileError(
            f"Invalid JSON strain profile {profile_path}: {exc.msg}"
        ) from exc
    profile = _validate_profile(value)
    return copy.deepcopy(profile)


def _get_reaction(model: Any, reaction_id: str) -> Any:
    try:
        return model.reactions.get_by_id(reaction_id)
    except KeyError as exc:
        raise ValueError(f"PO1f overlay requires reaction {reaction_id}") from exc


def _normalize_medium(
    medium: Mapping[str, Any],
    where: str,
) -> dict[str, float]:
    if not isinstance(medium, Mapping):
        raise ValueError(f"{where} must be a mapping of reaction IDs to uptake")
    normalized: dict[str, float] = {}
    for reaction_id, uptake_value in medium.items():
        if type(reaction_id) is not str or not reaction_id:
            raise ValueError(f"{where} contains an invalid reaction ID")
        try:
            uptake = _require_number(
                uptake_value,
                f"{where}[{reaction_id!r}]",
                minimum=0.0,
            )
        except StrainProfileError as exc:
            raise ValueError(str(exc)) from exc
        normalized[reaction_id] = uptake
    return normalized


def _medium_for_audit(
    medium: Mapping[str, float],
    *,
    explicit_reaction_ids: set[str],
) -> dict[str, float]:
    result = {str(key): float(value) for key, value in medium.items()}
    for reaction_id in explicit_reaction_ids:
        result.setdefault(reaction_id, 0.0)
    return dict(sorted(result.items()))


def _bounds_dict(reaction: Any) -> dict[str, float]:
    return {
        "lower_bound": float(reaction.lower_bound),
        "upper_bound": float(reaction.upper_bound),
    }


def _stoichiometry(reaction: Any) -> dict[str, float]:
    return {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in sorted(
            reaction.metabolites.items(),
            key=lambda item: item[0].id,
        )
    }


def apply_strain_overlay(
    model: Any,
    profile: dict[str, Any],
    *,
    active_medium: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Apply the strict PO1f profile to a disposable model copy in memory.

    The caller-owned ``model`` is mutated.  To preserve a canonical or baseline
    instance, pass ``baseline.copy()``.  All profile and model preconditions are
    checked before mutation.  ``active_medium`` should be the complete medium
    mapping just applied to the model; when omitted, ``model.medium`` is used.

    The returned audit contains only JSON-compatible primitive values.
    """

    validated = _validate_profile(profile)
    disable_operation = next(
        operation
        for operation in validated["operations"]
        if operation["type"] == "disable_reaction"
    )
    plasmid_operation = next(
        operation
        for operation in validated["operations"]
        if operation["type"] == "plasmid_complement"
    )

    ura3_reaction = _get_reaction(model, URA3_REACTION_ID)
    leu2_reaction = _get_reaction(model, LEU2_REACTION_ID)
    leucine_exchange = _get_reaction(model, LEUCINE_EXCHANGE_ID)
    uracil_exchange = _get_reaction(model, URACIL_EXCHANGE_ID)

    expected_ura3_bounds = _validate_bounds(
        disable_operation["expected_before"],
        "profile.operations[disable_reaction].expected_before",
    )
    target_ura3_bounds = _validate_bounds(
        disable_operation["set_bounds"],
        "profile.operations[disable_reaction].set_bounds",
    )
    actual_ura3_bounds = (
        float(ura3_reaction.lower_bound),
        float(ura3_reaction.upper_bound),
    )
    expected_leu2_gpr = plasmid_operation["expected_before"]["gene_reaction_rule"]
    target_leu2_gpr = plasmid_operation["set_gene_reaction_rule"]
    actual_leu2_gpr = str(leu2_reaction.gene_reaction_rule)
    baseline_model_state = (
        actual_ura3_bounds == expected_ura3_bounds
        and actual_leu2_gpr == expected_leu2_gpr
    )
    overlay_model_state = (
        actual_ura3_bounds == target_ura3_bounds
        and actual_leu2_gpr == target_leu2_gpr
    )
    if not baseline_model_state and not overlay_model_state:
        raise ValueError(
            "PO1f model preconditions are neither the expected state before "
            "PO1f overlay nor "
            "the complete applied overlay: "
            f"{URA3_REACTION_ID} bounds={actual_ura3_bounds}, "
            f"{LEU2_REACTION_ID} GPR={actual_leu2_gpr!r}"
        )

    exchange_ids = {reaction.id for reaction in model.exchanges}
    for reaction in (leucine_exchange, uracil_exchange):
        if reaction.id not in exchange_ids:
            raise ValueError(
                f"PO1f medium assertion target {reaction.id} is not a model "
                "exchange reaction"
            )

    actual_model_medium = _normalize_medium(model.medium, "model.medium")
    if active_medium is None:
        medium_before = dict(actual_model_medium)
    else:
        medium_before = _normalize_medium(active_medium, "active_medium")
        for reaction_id in {LEUCINE_EXCHANGE_ID, URACIL_EXCHANGE_ID}:
            supplied = medium_before.get(reaction_id, 0.0)
            actual = actual_model_medium.get(reaction_id, 0.0)
            if not math.isclose(supplied, actual, rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError(
                    f"active_medium reports {reaction_id} uptake {supplied}, "
                    f"but the model has {actual}"
                )

    assertions = validated["medium"]["uptake_assertions"]
    runtime_overrides = {
        reaction_id: float(value)
        for reaction_id, value in validated["medium"][
            "runtime_uptake_overrides"
        ].items()
    }
    baseline_medium_state = all(
        math.isclose(
            medium_before.get(reaction_id, 0.0),
            float(expected_value),
            rel_tol=5e-4,
            abs_tol=5e-6,
        )
        for reaction_id, expected_value in assertions.items()
    )
    overlay_medium_state = all(
        math.isclose(
            medium_before.get(reaction_id, 0.0),
            float(
                runtime_overrides.get(
                    reaction_id,
                    assertions.get(reaction_id, 0.0),
                )
            ),
            rel_tol=1e-9,
            abs_tol=1e-12,
        )
        for reaction_id in set(assertions) | set(runtime_overrides)
    )
    already_applied = overlay_model_state and overlay_medium_state
    for reaction_id, expected_value in assertions.items():
        actual_value = medium_before.get(reaction_id, 0.0)
        required_value = (
            runtime_overrides.get(reaction_id, expected_value)
            if already_applied
            else expected_value
        )
        if not math.isclose(
            actual_value,
            float(required_value),
            rel_tol=5e-4,
            abs_tol=5e-6,
        ):
            raise ValueError(
                f"Medium assertion failed for {reaction_id}: uptake is "
                f"{actual_value}, expected {float(required_value)}"
            )
    if not (
        (baseline_model_state and baseline_medium_state)
        or already_applied
    ):
        raise ValueError(
            "PO1f overlay found a partial or mismatched model/medium state; "
            "refusing a non-atomic reapplication"
        )

    protected_stoichiometry = {
        reaction.id: _stoichiometry(reaction)
        for reaction in (ura3_reaction, leu2_reaction)
    }
    operations_audit: list[dict[str, Any]] = []

    ura3_before = _bounds_dict(ura3_reaction)
    ura3_reaction.bounds = target_ura3_bounds
    operations_audit.append(
        {
            "type": "disable_reaction",
            "reaction_id": URA3_REACTION_ID,
            "locus": disable_operation["locus"],
            "gene_id": disable_operation["gene_id"],
            "legacy_gene_id": disable_operation["legacy_gene_id"],
            "allele": disable_operation["allele"],
            "protein_function": disable_operation["protein_function"],
            "evidence_status": disable_operation["evidence_status"],
            "before": ura3_before,
            "after": _bounds_dict(ura3_reaction),
            "changed": ura3_before != _bounds_dict(ura3_reaction),
        }
    )

    leu2_before_gpr = actual_leu2_gpr
    leu2_reaction.gene_reaction_rule = target_leu2_gpr
    pseudo_gene = model.genes.get_by_id(LEU2_PLASMID_PSEUDO_GENE)
    pseudo_gene.name = "LEU2 plasmid complementation"
    pseudo_notes = (
        dict(pseudo_gene.notes) if isinstance(pseudo_gene.notes, dict) else {}
    )
    pseudo_notes.update(
        {
            "strain": PROFILE_STRAIN,
            "overlay_role": "plasmid_complementation_pseudo_gene",
            "complemented_locus": LEU2_LOCUS,
            "protein_function": plasmid_operation["protein_function"],
            "evidence_status": plasmid_operation["evidence_status"],
        }
    )
    pseudo_gene.notes = pseudo_notes
    operations_audit.append(
        {
            "type": "plasmid_complement",
            "reaction_id": LEU2_REACTION_ID,
            "locus": plasmid_operation["locus"],
            "gene_id": plasmid_operation["gene_id"],
            "legacy_gene_id": plasmid_operation["legacy_gene_id"],
            "allele": plasmid_operation["allele"],
            "complement_source": plasmid_operation["complement_source"],
            "protein_function": plasmid_operation["protein_function"],
            "evidence_status": plasmid_operation["evidence_status"],
            "pseudo_gene": LEU2_PLASMID_PSEUDO_GENE,
            "before": {"gene_reaction_rule": leu2_before_gpr},
            "after": {"gene_reaction_rule": str(leu2_reaction.gene_reaction_rule)},
            "changed": leu2_before_gpr != str(leu2_reaction.gene_reaction_rule),
        }
    )

    runtime_medium = dict(actual_model_medium)
    for reaction_id, uptake in runtime_overrides.items():
        if uptake == 0.0:
            runtime_medium.pop(reaction_id, None)
        else:
            runtime_medium[reaction_id] = uptake
    model.medium = runtime_medium
    final_model_medium = _normalize_medium(model.medium, "model.medium")
    for reaction_id, expected_uptake in runtime_overrides.items():
        actual_uptake = final_model_medium.get(reaction_id, 0.0)
        if not math.isclose(
            actual_uptake,
            expected_uptake,
            rel_tol=1e-9,
            abs_tol=1e-12,
        ):
            raise RuntimeError(
                f"Runtime medium override failed for {reaction_id}: "
                f"{actual_uptake} != {expected_uptake}"
            )

    for reaction in (ura3_reaction, leu2_reaction):
        if _stoichiometry(reaction) != protected_stoichiometry[reaction.id]:
            raise RuntimeError(
                f"PO1f overlay changed protected stoichiometry for {reaction.id}"
            )

    explicit_medium_ids = set(assertions) | set(runtime_overrides)
    active_medium_audit = _medium_for_audit(
        final_model_medium,
        explicit_reaction_ids=explicit_medium_ids,
    )
    audit: dict[str, Any] = {
        "schema_version": PROFILE_SCHEMA_VERSION,
        "workflow": "po1f_runtime_strain_overlay",
        "profile_id": validated["profile_id"],
        "strain": copy.deepcopy(validated["strain"]),
        "operations": operations_audit,
        "medium": {
            "assertions_checked": {
                reaction_id: float(value)
                for reaction_id, value in sorted(assertions.items())
            },
            "runtime_uptake_overrides": dict(sorted(runtime_overrides.items())),
            "active_medium_before": _medium_for_audit(
                medium_before,
                explicit_reaction_ids=explicit_medium_ids,
            ),
            "active_medium": active_medium_audit,
            "formulation": copy.deepcopy(validated["medium"]["formulation"]),
        },
        "active_medium": active_medium_audit,
        "assay_confounded_loci": list(validated["assay_confounded_loci"]),
        "provenance_only_variants": list(validated["provenance_only_variants"]),
        "sources": list(validated["sources"]),
        "overlay_already_applied": already_applied,
        "model_layer_modified_in_memory": not already_applied,
        "canonical_sbml_written": False,
    }
    json.dumps(audit, sort_keys=True, allow_nan=False)
    return audit


__all__ = [
    "CONCENTRATION_RATIO_STATUS",
    "LEU2_EXPECTED_GPR",
    "LEU2_PLASMID_PSEUDO_GENE",
    "PROFILE_SCHEMA_VERSION",
    "RUNTIME_OVERRIDE_STATUS",
    "StrainProfileError",
    "apply_strain_overlay",
    "derive_concentration_ratio_uptake",
    "load_strain_profile",
]
