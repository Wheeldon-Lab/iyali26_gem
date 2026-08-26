"""Runtime-only CoQ9 dilution support for reproducible dFBA sensitivity runs.

The helper never writes SBML.  Call it only inside ``with model:`` so the
temporary reaction and linear constraint are removed after one simulation.
"""

from __future__ import annotations

import math
import random

from cobra import Reaction


Q9_METABOLITE_ID = "m468[C_mi]"
BIOMASS_REACTION_ID = "biomass_C"
COQ9_DILUTION_ID = "COQ9_DILUTION"
COQ9_DILUTION_CONSTRAINT_ID = "COQ9_DILUTION_GROWTH_COUPLING"
COQ9_POOL_SOURCE_ID = "DFBA_Q9_POOL_SOURCE"
ALPHA_SAMPLER_ID = "coq9_log10_triangular_v1"
ALPHA_LOW_MMOL_GDW = 1e-6
ALPHA_MODE_MMOL_GDW = 1e-4
ALPHA_HIGH_MMOL_GDW = 1e-3
CALIBRATION_STATUS = "sensitivity_only_not_calibrated"


def _integer(value: int, *, name: str, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _alpha(value: float) -> float:
    alpha = float(value)
    if not math.isfinite(alpha) or alpha < 0:
        raise ValueError("CoQ9 alpha must be finite and non-negative")
    return alpha


def sample_coq9_alpha(base_seed: int, replicate_id: int) -> float:
    """Draw one reproducible, log-space triangular CoQ9 coefficient."""
    seed = _integer(base_seed, name="base_seed")
    replicate = _integer(replicate_id, name="replicate_id", minimum=0)
    rng = random.Random(f"{ALPHA_SAMPLER_ID}:{seed}:{replicate}")
    log10_alpha = rng.triangular(-6.0, -3.0, -4.0)
    return float(10.0**log10_alpha)


def alpha_sampling_record(base_seed: int, replicate_id: int) -> dict[str, object]:
    """Return the complete, manifest-ready record for one sampled coefficient."""
    return {
        "sampler_id": ALPHA_SAMPLER_ID,
        "base_seed": _integer(base_seed, name="base_seed"),
        "replicate_id": _integer(replicate_id, name="replicate_id", minimum=0),
        "distribution": "log10_triangular",
        "low_mmol_gDW": ALPHA_LOW_MMOL_GDW,
        "mode_mmol_gDW": ALPHA_MODE_MMOL_GDW,
        "high_mmol_gDW": ALPHA_HIGH_MMOL_GDW,
        "alpha_mmol_gDW": sample_coq9_alpha(base_seed, replicate_id),
    }


def _require_absent(model, identifier: str, *, kind: str) -> None:
    collection = model.reactions if kind == "reaction" else model.constraints
    if identifier in collection:
        raise ValueError(f"runtime CoQ9 {kind} already exists: {identifier}")


def apply_runtime_coq9_dilution(model, alpha: float) -> Reaction:
    """Add ``CoQ9 -> ∅`` and enforce ``v_dilution = alpha * v_biomass``."""
    coefficient = _alpha(alpha)
    _require_absent(model, COQ9_DILUTION_ID, kind="reaction")
    _require_absent(model, COQ9_DILUTION_CONSTRAINT_ID, kind="constraint")
    try:
        q9 = model.metabolites.get_by_id(Q9_METABOLITE_ID)
        biomass = model.reactions.get_by_id(BIOMASS_REACTION_ID)
    except KeyError as error:
        raise ValueError(
            "runtime CoQ9 dilution requires "
            f"{Q9_METABOLITE_ID} and {BIOMASS_REACTION_ID}"
        ) from error

    dilution = Reaction(COQ9_DILUTION_ID, lower_bound=0.0, upper_bound=1000.0)
    dilution.add_metabolites({q9: -1.0})
    model.add_reactions([dilution])
    coupling = model.problem.Constraint(
        dilution.flux_expression - coefficient * biomass.flux_expression,
        lb=0.0,
        ub=0.0,
        name=COQ9_DILUTION_CONSTRAINT_ID,
    )
    model.add_cons_vars(coupling)
    return dilution


def add_runtime_coq9_pool_source(model) -> Reaction:
    """Add a disabled one-way CoQ9 source for explicitly modeled finite reserves."""
    _require_absent(model, COQ9_POOL_SOURCE_ID, kind="reaction")
    try:
        q9 = model.metabolites.get_by_id(Q9_METABOLITE_ID)
    except KeyError as error:
        raise ValueError(f"runtime CoQ9 reserve requires {Q9_METABOLITE_ID}") from error
    source = Reaction(COQ9_POOL_SOURCE_ID, lower_bound=0.0, upper_bound=0.0)
    source.add_metabolites({q9: 1.0})
    model.add_reactions([source])
    return source


def runtime_coq9_dilution_record(alpha: float, *, pool_source_enabled: bool) -> dict[str, object]:
    """Describe the disposable mechanism without claiming biological calibration."""
    return {
        "scope": "runtime_only",
        "calibration_status": CALIBRATION_STATUS,
        "reaction_id": COQ9_DILUTION_ID,
        "metabolite_id": Q9_METABOLITE_ID,
        "biomass_reaction_id": BIOMASS_REACTION_ID,
        "coupling": "v_COQ9_DILUTION = alpha_mmol_gDW * v_biomass_C",
        "alpha_mmol_gDW": _alpha(alpha),
        "pool_source": {
            "enabled": pool_source_enabled,
            "reaction_id": COQ9_POOL_SOURCE_ID if pool_source_enabled else None,
            "direction": "source_to_CoQ9" if pool_source_enabled else None,
        },
    }
