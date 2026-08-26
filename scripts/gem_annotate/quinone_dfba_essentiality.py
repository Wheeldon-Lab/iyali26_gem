"""Runtime-only CoQ9 reserve dFBA screen for the PO1f SD-Leu model.

This is a sensitivity experiment, not a calibrated biomass change.  It adds a
temporary CoQ9 demand proportional to biomass and a finite, one-way CoQ9
reserve.  The reserve is enabled only after source-free growth is zero, so it
is a complete-block control rather than a general partial-loss rescue model.
Neither the canonical SBML nor the essentiality dossier is changed.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from cobra import Metabolite, Reaction
from cobra.exceptions import OptimizationError
from cobra.flux_analysis import pfba
from cobra.manipulation.delete import knock_out_model_genes

from .config import load_project_paths
from .coq9_dilution import (
    CALIBRATION_STATUS,
    add_runtime_coq9_pool_source,
    alpha_sampling_record,
    apply_runtime_coq9_dilution,
    runtime_coq9_dilution_record,
)
from .essentiality_simulation_context import (
    load_effective_simulation_context,
    sha256_file,
    sha256_payload,
)
from .essentiality_evidence import target_fingerprint
from .validate_essential_genes import DEFAULT_CUTOFFS, load_experimental, reaction_case_context


WORKFLOW = "quinone_dfba_essentiality"
SCHEMA_VERSION = "1.6"
BIOMASS_ID = "biomass_C"
URACIL_MODES = ("finite_batch", "po1f_nonlimiting")

# Values are the concentration statements already recorded in the SD-Leu
# medium comments.  They are finite batch inventories, not uptake kinetics.
INITIAL_POOLS_MMOL_L = {
    "R1070": 111.0,
    "R1003": 0.054,
    "R1202": 0.287,
    "R1204": 0.601,
    "R1215": 0.095,
    "R1217": 0.381,
    "R1220": 0.274,
    "R1222": 0.134,
    "R1223": 0.303,
    "R1231": 0.840,
    "R1232": 0.245,
    "R1233": 0.276,
    "R1234": 1.195,
    "R1354": 0.178428,
}
DEFAULT_ALPHAS = (1e-6, 1e-4, 1e-3)
DEFAULT_POOL_MULTIPLIERS = (0.0, 0.5, 1.0, 2.0)
GUROBI_FEASIBILITY_TOL_MIN = 1e-9
GUROBI_FEASIBILITY_TOL_MAX = 1e-2
Q9_RESERVE_POLICY = "source_enabled_only_when_source_free_growth_is_zero"
R39_R19_RUNTIME_TOPOLOGY_ID = "r39_mito_r19_split_v1"
R39_R19_RUNTIME_REACTION_IDS = (
    "R763", "R407", "R969", "R39", "R808", "R715", "R40", "R19",
    "DFBA_R19_HYDROXYLATION", "DFBA_R19_FORMAL_OXIDATION", "R18", "R695", "R385",
)
R39_R19_RUNTIME_MAPPING = {
    "id": R39_R19_RUNTIME_TOPOLOGY_ID,
    "scope": "runtime_only",
    "calibration_status": "sensitivity_only_not_calibrated",
    "r39": "m641[C_mi] + 0.5 m64[C_mi] -> m28[C_mi] + m939[C_mi]",
    "closed_transport_reactions": ["R969", "R808"],
    "closed_legacy_reaction": "R19",
    "r19_hydroxylation": "m63[C_mi] + 0.5 m64[C_mi] -> DFBA_DDMQH2[C_mi]",
    "r19_formal_oxidation": "DFBA_DDMQH2[C_mi] + 0.5 m64[C_mi] -> m26[C_mi] + m59[C_mi]",
    "formal_oxidation_limit": "O2 is a balanced runtime bookkeeping convention, not a verified biological oxidant or enzyme.",
}
R39_R19_RUNTIME_MAPPING_SHA256 = sha256_payload(R39_R19_RUNTIME_MAPPING)
RUNTIME_GPR_TARGET_REACTIONS = (
    "R39", "R715", "R40", "DFBA_R19_HYDROXYLATION", "R18", "R695", "R385",
)
RUNTIME_GPR_BASELINE_RULES = {
    "R39": "",
    "R715": "YALI1B20835g",
    "R40": "YALI1F34625g",
    "DFBA_R19_HYDROXYLATION": "",
    "R18": "YALI1C25352g",
    "R695": "YALI1E18269g",
    "R385": "YALI1B20835g",
}
RUNTIME_GPR_SCENARIOS = {
    "baseline_topology": {
        "mapping_id": "baseline_topology_v1",
        "gpr_updates": {},
        "added_candidate_genes": [],
        "expected_added_gene_ko_closed_reactions": [],
        "hypothesis_label": "topology baseline with no added candidate GPR",
        "risk_labels": [],
    },
    "coq6_r39": {
        "mapping_id": "coq6_r39_v1",
        "gpr_updates": {"R39": "YALI1A08781g"},
        "added_candidate_genes": ["YALI1A08781g"],
        "expected_added_gene_ko_closed_reactions": ["R39"],
        "hypothesis_label": "constrained runtime direct-catalyst sensitivity mapping",
        "risk_labels": ["database_comparative_only", "native_yarrowia_unverified"],
    },
    "coq6_r19_hydroxylation": {
        "mapping_id": "coq6_r19_hydroxylation_v1",
        "gpr_updates": {"DFBA_R19_HYDROXYLATION": "YALI1A08781g"},
        "added_candidate_genes": ["YALI1A08781g"],
        "expected_added_gene_ko_closed_reactions": ["DFBA_R19_HYDROXYLATION"],
        "hypothesis_label": "runtime R19 hydroxylation direct-catalyst counterfactual",
        "risk_labels": ["chemistry_redox_risk", "r19_product_redox_chemistry_unresolved"],
    },
    "coq6_both_hydroxylations": {
        "mapping_id": "coq6_both_hydroxylations_v1",
        "gpr_updates": {
            "R39": "YALI1A08781g",
            "DFBA_R19_HYDROXYLATION": "YALI1A08781g",
        },
        "added_candidate_genes": ["YALI1A08781g"],
        "expected_added_gene_ko_closed_reactions": ["R39", "DFBA_R19_HYDROXYLATION"],
        "hypothesis_label": "comparative/database two-hydroxylation hypothesis only",
        "risk_labels": ["comparative_database_only", "chemistry_redox_risk"],
    },
    "coq8_routewide_absolute": {
        "mapping_id": "coq8_routewide_absolute_v1",
        "gpr_updates": {
            "R39": "YALI1B20527g",
            "R715": "YALI1B20527g and YALI1B20835g",
            "R40": "YALI1B20527g and YALI1F34625g",
            "DFBA_R19_HYDROXYLATION": "YALI1B20527g",
            "R18": "YALI1B20527g and YALI1C25352g",
            "R695": "YALI1B20527g and YALI1E18269g",
            "R385": "YALI1B20527g and YALI1B20835g",
        },
        "added_candidate_genes": ["YALI1B20527g"],
        "expected_added_gene_ko_closed_reactions": list(RUNTIME_GPR_TARGET_REACTIONS),
        "hypothesis_label": "routewide-absolute synthome-member counterfactual; not reaction evidence",
        "risk_labels": ["routewide_absolute_counterfactual", "not_reaction_evidence"],
    },
    "coq9_r695_absolute": {
        "mapping_id": "coq9_r695_absolute_v1",
        "gpr_updates": {"R695": "YALI1F34675g and YALI1E18269g"},
        "added_candidate_genes": ["YALI1F34675g"],
        "expected_added_gene_ko_closed_reactions": ["R695"],
        "hypothesis_label": "absolute R695 accessory-dependence counterfactual",
        "risk_labels": ["cross_species_conflicted"],
    },
    "coq9_routewide_absolute": {
        "mapping_id": "coq9_routewide_absolute_v1",
        "gpr_updates": {
            "R39": "YALI1F34675g",
            "R715": "YALI1F34675g and YALI1B20835g",
            "R40": "YALI1F34675g and YALI1F34625g",
            "DFBA_R19_HYDROXYLATION": "YALI1F34675g",
            "R18": "YALI1F34675g and YALI1C25352g",
            "R695": "YALI1F34675g and YALI1E18269g",
            "R385": "YALI1F34675g and YALI1B20835g",
        },
        "added_candidate_genes": ["YALI1F34675g"],
        "expected_added_gene_ko_closed_reactions": list(RUNTIME_GPR_TARGET_REACTIONS),
        "hypothesis_label": "routewide-absolute synthome-member counterfactual; not reaction evidence",
        "risk_labels": ["routewide_absolute_counterfactual", "not_reaction_evidence"],
    },
}


def _reaction_stoichiometry(reaction: Reaction) -> dict[str, float]:
    return {metabolite.id: float(coefficient) for metabolite, coefficient in reaction.metabolites.items()}


def _expect_runtime_reaction(
    model, reaction_id: str, stoichiometry: dict[str, float], bounds: tuple[float, float]
) -> None:
    reaction = model.reactions.get_by_id(reaction_id)
    if reaction.gene_reaction_rule or reaction.bounds != bounds or _reaction_stoichiometry(reaction) != stoichiometry:
        raise ValueError(f"R39/R19 runtime topology requires the current canonical {reaction_id} definition")


def _replace_stoichiometry(reaction: Reaction, stoichiometry: dict[Metabolite, float]) -> None:
    reaction.add_metabolites({metabolite: -coefficient for metabolite, coefficient in reaction.metabolites.items()})
    reaction.add_metabolites(stoichiometry)


def _apply_r39_r19_runtime_topology(model) -> dict[str, Any]:
    """Apply the reviewed R39/R19 topology convention to one disposable model copy."""
    if any(reaction_id in model.reactions for reaction_id in R39_R19_RUNTIME_REACTION_IDS[8:10]):
        raise ValueError("R39/R19 runtime topology is already present")
    _expect_runtime_reaction(
        model, "R39",
        {"m108[C_cy]": -1.0, "m109[C_cy]": -0.5, "m10[C_cy]": 1.0, "m110[C_cy]": 1.0},
        (0.0, 1000.0),
    )
    _expect_runtime_reaction(model, "R969", {"m108[C_cy]": -1.0, "m641[C_mi]": 1.0}, (-1000.0, 1000.0))
    _expect_runtime_reaction(model, "R808", {"m110[C_cy]": -1.0, "m939[C_mi]": 1.0}, (-1000.0, 1000.0))
    _expect_runtime_reaction(
        model, "R19",
        {"m63[C_mi]": -1.0, "m64[C_mi]": -1.0, "m26[C_mi]": 1.0, "m59[C_mi]": 1.0},
        (0.0, 1000.0),
    )
    required_metabolites = {
        "m641[C_mi]": ("C52H78O3", 0, "C_mi"),
        "m939[C_mi]": ("C52H77O4", -1, "C_mi"),
        "m63[C_mi]": ("C52H80O2", 0, "C_mi"),
        "m59[C_mi]": ("C52H78O3", 0, "C_mi"),
    }
    for metabolite_id, expected in required_metabolites.items():
        metabolite = model.metabolites.get_by_id(metabolite_id)
        if (metabolite.formula, metabolite.charge, metabolite.compartment) != expected:
            raise ValueError(f"R39/R19 runtime topology requires the current canonical {metabolite_id} definition")

    r39 = model.reactions.get_by_id("R39")
    _replace_stoichiometry(r39, {
        model.metabolites.get_by_id("m641[C_mi]"): -1.0,
        model.metabolites.get_by_id("m64[C_mi]"): -0.5,
        model.metabolites.get_by_id("m28[C_mi]"): 1.0,
        model.metabolites.get_by_id("m939[C_mi]"): 1.0,
    })
    model.reactions.get_by_id("R969").bounds = (0.0, 0.0)
    model.reactions.get_by_id("R808").bounds = (0.0, 0.0)
    model.reactions.get_by_id("R19").bounds = (0.0, 0.0)

    ddmqh2 = Metabolite(
        "DFBA_DDMQH2[C_mi]", name="2-methoxy-6-nonaprenylbenzene-1,4-diol (runtime-only)",
        formula="C52H80O3", charge=0, compartment="C_mi",
    )
    hydroxylation = Reaction("DFBA_R19_HYDROXYLATION", lower_bound=0.0, upper_bound=1000.0)
    hydroxylation.add_metabolites({
        model.metabolites.get_by_id("m63[C_mi]"): -1.0,
        model.metabolites.get_by_id("m64[C_mi]"): -0.5,
        ddmqh2: 1.0,
    })
    oxidation = Reaction("DFBA_R19_FORMAL_OXIDATION", lower_bound=0.0, upper_bound=1000.0)
    oxidation.add_metabolites({
        ddmqh2: -1.0,
        model.metabolites.get_by_id("m64[C_mi]"): -0.5,
        model.metabolites.get_by_id("m26[C_mi]"): 1.0,
        model.metabolites.get_by_id("m59[C_mi]"): 1.0,
    })
    model.add_reactions([hydroxylation, oxidation])
    balances = {
        reaction_id: model.reactions.get_by_id(reaction_id).check_mass_balance()
        for reaction_id in ("R39", "DFBA_R19_HYDROXYLATION", "DFBA_R19_FORMAL_OXIDATION")
    }
    if any(balances.values()):
        raise ValueError(f"R39/R19 runtime topology is not balanced: {balances}")
    return {
        "enabled": True,
        **R39_R19_RUNTIME_MAPPING,
        "mapping_sha256": R39_R19_RUNTIME_MAPPING_SHA256,
        "mass_balance": balances,
    }


def _runtime_topology_fluxes(solution, enabled: bool) -> dict[str, float]:
    if not enabled:
        return {}
    if solution.status != "optimal":
        return {f"pathway_flux_{reaction_id}_mmol_gDW_h": math.nan for reaction_id in R39_R19_RUNTIME_REACTION_IDS}
    return {
        f"pathway_flux_{reaction_id}_mmol_gDW_h": float(solution.fluxes[reaction_id])
        for reaction_id in R39_R19_RUNTIME_REACTION_IDS
    }


def _r39_r19_canonical_target_fingerprint(model) -> str:
    return target_fingerprint(
        reaction_case_context(model.reactions.get_by_id(reaction_id))
        for reaction_id in ("R39", "R969", "R808", "R19")
    )


def _runtime_gpr_mapping_payload(scenario: str) -> dict[str, Any]:
    if scenario not in RUNTIME_GPR_SCENARIOS:
        raise ValueError(f"Unknown runtime GPR scenario: {scenario}")
    return {
        "scenario_id": scenario,
        "scope": "runtime_only",
        "calibration_status": "sensitivity_only_not_calibrated",
        "required_topology_id": R39_R19_RUNTIME_TOPOLOGY_ID,
        "required_topology_mapping_sha256": R39_R19_RUNTIME_MAPPING_SHA256,
        "protected_unmapped_reactions": ["R19", "DFBA_R19_FORMAL_OXIDATION"],
        "binary_gpr_interpretation": (
            "Tests the imposed complete-block assumption, not the mapped protein's function."
        ),
        "expected_baseline_gprs": RUNTIME_GPR_BASELINE_RULES,
        **RUNTIME_GPR_SCENARIOS[scenario],
    }


def _apply_runtime_gpr_scenario(model, scenario: str) -> dict[str, Any]:
    """Apply one preregistered GPR hypothesis after the runtime topology transform."""
    mapping = _runtime_gpr_mapping_payload(scenario)
    required = {*RUNTIME_GPR_TARGET_REACTIONS, "DFBA_R19_FORMAL_OXIDATION", "R19"}
    if not required <= {reaction.id for reaction in model.reactions}:
        raise ValueError("Runtime GPR scenarios require the R39/R19 runtime topology")
    observed = {
        reaction_id: model.reactions.get_by_id(reaction_id).gene_reaction_rule
        for reaction_id in RUNTIME_GPR_TARGET_REACTIONS
    }
    if observed != RUNTIME_GPR_BASELINE_RULES:
        raise ValueError(f"Runtime GPR scenario baseline mismatch: {observed}")
    if (
        model.reactions.R19.bounds != (0.0, 0.0)
        or model.reactions.R19.gene_reaction_rule
        or model.reactions.get_by_id("DFBA_R19_FORMAL_OXIDATION").gene_reaction_rule
    ):
        raise ValueError("Runtime GPR scenarios require closed GPR-less R19 and GPR-less formal oxidation")
    for gene_id in mapping["added_candidate_genes"]:
        if not model.genes.has_id(gene_id):
            raise ValueError(f"Runtime candidate gene is absent from the base model: {gene_id}")
        gene = model.genes.get_by_id(gene_id)
        if gene.reactions:
            raise ValueError(f"Runtime candidate gene is not orphaned in the base model: {gene_id}")

    target_fingerprint_before = target_fingerprint(
        reaction_case_context(model.reactions.get_by_id(reaction_id))
        for reaction_id in (*RUNTIME_GPR_TARGET_REACTIONS, "DFBA_R19_FORMAL_OXIDATION", "R19")
    )
    for reaction_id, rule in mapping["gpr_updates"].items():
        model.reactions.get_by_id(reaction_id).gene_reaction_rule = rule
    expected = {**RUNTIME_GPR_BASELINE_RULES, **mapping["gpr_updates"]}
    applied = {
        reaction_id: model.reactions.get_by_id(reaction_id).gene_reaction_rule
        for reaction_id in RUNTIME_GPR_TARGET_REACTIONS
    }
    if applied != expected:
        raise RuntimeError(f"Runtime GPR scenario did not apply exactly: {applied}")
    return {
        "enabled": True,
        **mapping,
        "mapping_sha256": sha256_payload(mapping),
        "runtime_target_fingerprint_before_mapping": target_fingerprint_before,
    }


def _finite_medium(base_medium: dict[str, float], pools: dict[str, float], biomass: float, dt: float) -> dict[str, float]:
    medium = dict(base_medium)
    for reaction_id, amount in pools.items():
        if reaction_id not in medium:
            continue
        medium[reaction_id] = min(float(medium[reaction_id]), max(0.0, amount / (biomass * dt)))
    return medium


def _optimize_minimal_pool(model, source: Reaction) -> tuple[float, Any]:
    """Use reserve only after source-free growth is zero (complete-block control)."""
    model.objective = model.reactions.get_by_id(BIOMASS_ID)
    source_limit = source.upper_bound
    source.upper_bound = 0.0
    try:
        primary = pfba(model)
    except OptimizationError:
        primary = model.optimize()
    source.upper_bound = source_limit
    if primary.status != "optimal":
        return 0.0, primary
    growth = max(0.0, float(primary.fluxes[BIOMASS_ID]))
    if growth > 1e-9:
        return growth, primary
    try:
        primary = pfba(model)
    except OptimizationError:
        primary = model.optimize()
    if primary.status != "optimal":
        return 0.0, primary
    return max(0.0, float(primary.fluxes[BIOMASS_ID])), primary


def _software_versions(solver: str) -> dict[str, str]:
    versions = {"python": sys.version.split()[0], "cobra": __import__("cobra").__version__}
    if solver.lower() == "gurobi":
        import gurobipy

        versions["gurobipy"] = ".".join(map(str, gurobipy.gurobi.version()))
    return versions


def _gurobi_feasibility_tolerance(value: str) -> float:
    tolerance = float(value)
    if not math.isfinite(tolerance) or not GUROBI_FEASIBILITY_TOL_MIN <= tolerance <= GUROBI_FEASIBILITY_TOL_MAX:
        raise argparse.ArgumentTypeError(
            f"Gurobi feasibility tolerance must be in [{GUROBI_FEASIBILITY_TOL_MIN:g}, {GUROBI_FEASIBILITY_TOL_MAX:g}]"
        )
    return tolerance


def _configure_solver(model, solver: str, feasibility_tolerance: float | None) -> float:
    model.solver = solver
    if feasibility_tolerance is not None:
        if solver.lower() != "gurobi":
            raise ValueError("--feasibility-tol requires --solver gurobi")
        model.solver.configuration.tolerances.feasibility = feasibility_tolerance
    return float(model.solver.configuration.tolerances.feasibility)


def _alpha_samples(args: argparse.Namespace) -> list[dict[str, Any]]:
    """Resolve one fixed coefficient per replicate before any WT/KO simulation."""
    if args.alpha_seed is not None:
        if args.alphas is not None:
            raise ValueError("--alpha-seed cannot be combined with --alphas")
        if args.alpha_replicates < 1:
            raise ValueError("--alpha-replicates must be positive with --alpha-seed")
        return [
            alpha_sampling_record(args.alpha_seed, replicate_id)
            for replicate_id in range(args.alpha_replicates)
        ]
    if args.alpha_replicates != 1:
        raise ValueError("--alpha-replicates requires --alpha-seed")
    alphas = DEFAULT_ALPHAS if args.alphas is None else tuple(args.alphas)
    samples = []
    for alpha in alphas:
        value = float(alpha)
        if not math.isfinite(value) or value < 0:
            raise ValueError("--alphas values must be finite and non-negative")
        samples.append(
            {
                "sampler_id": "explicit_alpha_v1",
                "base_seed": None,
                "replicate_id": None,
                "distribution": "explicit",
                "low_mmol_gDW": None,
                "mode_mmol_gDW": None,
                "high_mmol_gDW": None,
                "alpha_mmol_gDW": value,
            }
        )
    if not samples:
        raise ValueError("at least one alpha is required")
    return samples


def _pool_multiplier(value: float) -> float:
    multiplier = float(value)
    if not math.isfinite(multiplier) or multiplier < 0:
        raise ValueError("pool multipliers must be finite and non-negative")
    return multiplier


def simulate_gene(
    base_model,
    *,
    gene_id: str | None,
    alpha: float,
    pool_multiplier: float,
    hours: float,
    dt: float,
    initial_biomass: float,
    uracil_mode: str = "finite_batch",
    r39_r19_runtime_topology: bool = False,
    runtime_gpr_scenario: str | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Run a disposable Euler batch simulation for WT or one single-gene KO."""
    if not all(math.isfinite(value) and value > 0 for value in (hours, dt, initial_biomass)):
        raise ValueError("hours, dt, and initial_biomass must be positive")
    pool_multiplier = _pool_multiplier(pool_multiplier)
    if uracil_mode not in URACIL_MODES:
        raise ValueError(f"uracil_mode must be one of {URACIL_MODES}")
    if not math.isclose(hours / dt, round(hours / dt), rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("hours must be an exact multiple of dt")
    if runtime_gpr_scenario is not None and not r39_r19_runtime_topology:
        raise ValueError("--runtime-gpr-scenario requires --r39-r19-runtime-topology")
    with base_model as model:
        if r39_r19_runtime_topology:
            _apply_r39_r19_runtime_topology(model)
        if runtime_gpr_scenario is not None:
            _apply_runtime_gpr_scenario(model, runtime_gpr_scenario)
        if gene_id is not None:
            knock_out_model_genes(model, [gene_id])
        source = add_runtime_coq9_pool_source(model)
        apply_runtime_coq9_dilution(model, alpha)
        base_medium = dict(model.medium)
        pools = dict(INITIAL_POOLS_MMOL_L)
        if uracil_mode == "po1f_nonlimiting":
            pools.pop("R1354")
        initial_q9_pool = alpha * initial_biomass * pool_multiplier
        q9_pool = initial_q9_pool
        q9_tolerance = max(1e-12, initial_q9_pool * 1e-9)
        biomass = initial_biomass
        depleted_at: float | None = 0.0 if q9_pool == 0 else None
        termination_status = "completed"
        termination_time_h = hours
        trajectory: list[dict[str, Any]] = []
        for step in range(int(round(hours / dt))):
            time = step * dt
            model.medium = _finite_medium(base_medium, pools, biomass, dt)
            source.upper_bound = q9_pool / (biomass * dt) if q9_pool > 0 else 0.0
            growth, solution = _optimize_minimal_pool(model, source)
            if solution.status != "optimal":
                trajectory.append(
                    {
                        "gene_id": gene_id or "WT",
                        "uracil_mode": uracil_mode,
                        "time_h": time,
                        "biomass_gDW_L": biomass,
                        "growth_h-1": math.nan,
                        "q9_pool_mmol_L": q9_pool,
                        "q9_source_flux_mmol_gDW_h": math.nan,
                        "uracil_mmol_L": pools.get("R1354", math.nan),
                        "glucose_mmol_L": pools["R1070"],
                        "status": solution.status,
                        **_runtime_topology_fluxes(solution, r39_r19_runtime_topology),
                    }
                )
                termination_status = solution.status
                termination_time_h = time
                break
            source_flux = max(0.0, float(solution.fluxes[source.id]))
            consumed_q9 = min(q9_pool, source_flux * biomass * dt)
            q9_pool = max(0.0, q9_pool - consumed_q9)
            if q9_pool <= q9_tolerance:
                q9_pool = 0.0
            if q9_pool == 0.0 and depleted_at is None:
                depleted_at = time + dt
            for reaction_id in pools:
                uptake = max(0.0, -float(solution.fluxes[reaction_id]))
                pools[reaction_id] = max(0.0, pools[reaction_id] - uptake * biomass * dt)
            trajectory.append(
                {
                    "gene_id": gene_id or "WT",
                    "uracil_mode": uracil_mode,
                    "time_h": time,
                    "biomass_gDW_L": biomass,
                    "growth_h-1": growth,
                    "q9_pool_mmol_L": q9_pool,
                    "q9_source_flux_mmol_gDW_h": source_flux,
                    "uracil_mmol_L": pools.get("R1354", math.nan),
                    "glucose_mmol_L": pools["R1070"],
                    "status": solution.status,
                    **_runtime_topology_fluxes(solution, r39_r19_runtime_topology),
                }
            )
            biomass *= 1.0 + growth * dt
        return (
            {
                "gene_id": gene_id or "WT",
                "uracil_mode": uracil_mode,
                "final_biomass_gDW_L": biomass,
                "dynamic_doublings": math.log2(biomass / initial_biomass),
                "initial_growth_h-1": trajectory[0]["growth_h-1"],
                "q9_pool_depleted_h": depleted_at,
                "q9_source_total_mmol_L": initial_q9_pool - q9_pool,
                "final_glucose_mmol_L": pools["R1070"],
                "final_uracil_mmol_L": pools.get("R1354", math.nan),
                "termination_status": termination_status,
                "termination_time_h": termination_time_h,
            },
            trajectory,
        )


def _gene_ids(model, requested: str | None, excluded: tuple[str, ...], chunk_index: int, chunk_count: int) -> list[str]:
    if chunk_count <= 0 or not 0 <= chunk_index < chunk_count:
        raise ValueError("chunk-index must be in [0, chunk-count)")
    available = sorted(gene.id for gene in model.genes if gene.id not in set(excluded))
    if requested is None:
        selected = available
    else:
        selected = [item.strip() for item in requested.split(",")]
        if not all(selected):
            raise ValueError("Requested genes must not contain empty gene IDs")
        if len(selected) != len(set(selected)):
            raise ValueError("Requested genes must be unique")
    missing = sorted(set(selected) - set(available))
    if missing:
        raise ValueError(f"Requested genes are not screenable: {missing}")
    return selected[chunk_index::chunk_count]


def _run_chunk(args: argparse.Namespace) -> Path:
    if args.runtime_gpr_scenario is not None and not args.r39_r19_runtime_topology:
        raise ValueError("--runtime-gpr-scenario requires --r39-r19-runtime-topology")
    pool_multipliers = [_pool_multiplier(value) for value in args.pool_multipliers]
    paths = load_project_paths(args.research_root, required=True)
    media_path = paths.media / "sd_leu.csv"
    profile_path = paths.strain_profiles / "po1f_sd_leu.json"
    experimental_path = paths.essentiality / "consensus_essential_genes.csv"
    paths.require(media_path, profile_path, experimental_path, paths.output_model)
    context = load_effective_simulation_context(
        model_path=paths.output_model, media_path=media_path, strain_profile_path=profile_path
    )
    solver_feasibility_tolerance = _configure_solver(
        context.model, args.solver, args.feasibility_tol
    )
    base_r1354_bound = float(context.model.medium["R1354"])
    genes = _gene_ids(context.model, args.genes, context.excluded_runtime_genes, args.chunk_index, args.chunk_count)
    if not genes:
        raise ValueError("This chunk has no assigned genes; reduce --chunk-count for a targeted run")
    out = Path(args.output_dir or paths.results / WORKFLOW / args.run_id).resolve()
    out.mkdir(parents=True, exist_ok=True)
    experimental = load_experimental(experimental_path, positive_only=True)
    alpha_samples = _alpha_samples(args)
    rows: list[dict[str, Any]] = []
    trajectories: list[dict[str, Any]] = []
    topology_audit_rows: list[dict[str, Any]] = []
    runtime_topology: dict[str, Any] | None = None
    runtime_gpr_scenario: dict[str, Any] = {"enabled": False}
    if args.r39_r19_runtime_topology:
        with context.model as model:
            canonical_target = _r39_r19_canonical_target_fingerprint(model)
            runtime_topology = _apply_r39_r19_runtime_topology(model)
            runtime_topology["canonical_target_fingerprint"] = canonical_target
            if args.runtime_gpr_scenario is not None:
                runtime_gpr_scenario = _apply_runtime_gpr_scenario(model, args.runtime_gpr_scenario)
    runtime_gpr_mapping_sha = runtime_gpr_scenario.get("mapping_sha256", "")
    for alpha_sample in alpha_samples:
        alpha = float(alpha_sample["alpha_mmol_gDW"])
        replicate_id = alpha_sample["replicate_id"]
        for multiplier in pool_multipliers:
            wt, wt_trace = simulate_gene(context.model, gene_id=None, alpha=alpha, pool_multiplier=multiplier, hours=args.hours, dt=args.dt, initial_biomass=args.initial_biomass, uracil_mode=args.uracil_mode, r39_r19_runtime_topology=args.r39_r19_runtime_topology, runtime_gpr_scenario=args.runtime_gpr_scenario)
            trajectories.extend([
                {
                    **item,
                    "alpha_mmol_gDW": alpha,
                    "alpha_replicate_id": replicate_id,
                    "pool_multiplier": multiplier,
                }
                for item in wt_trace
            ])
            if args.r39_r19_runtime_topology:
                topology_audit_rows.append({
                    "alpha_mmol_gDW": alpha,
                    "alpha_replicate_id": replicate_id,
                    "pool_multiplier": multiplier,
                    "runtime_gpr_scenario": args.runtime_gpr_scenario or "",
                    "runtime_gpr_mapping_sha256": runtime_gpr_mapping_sha,
                    **{key: value for key, value in wt_trace[0].items() if key.startswith("pathway_flux_")},
                })
            for gene_id in genes:
                result, trace = simulate_gene(context.model, gene_id=gene_id, alpha=alpha, pool_multiplier=multiplier, hours=args.hours, dt=args.dt, initial_biomass=args.initial_biomass, uracil_mode=args.uracil_mode, r39_r19_runtime_topology=args.r39_r19_runtime_topology, runtime_gpr_scenario=args.runtime_gpr_scenario)
                ratio = result["dynamic_doublings"] / wt["dynamic_doublings"] if wt["dynamic_doublings"] > 0 else math.nan
                rows.append({
                    **result,
                    "alpha_mmol_gDW": alpha,
                    "alpha_replicate_id": replicate_id,
                    "pool_multiplier": multiplier,
                    "runtime_gpr_scenario": args.runtime_gpr_scenario or "",
                    "runtime_gpr_mapping_sha256": runtime_gpr_mapping_sha,
                    "wt_dynamic_doublings": wt["dynamic_doublings"], "dynamic_growth_ratio": ratio,
                    "wt_termination_status": wt["termination_status"],
                    "wt_termination_time_h": wt["termination_time_h"],
                    "experimental_essential": gene_id in set(experimental["gene_id"]),
                    **{f"essential_at_{cutoff * 100:g}pct": bool(ratio < cutoff) for cutoff in DEFAULT_CUTOFFS},
                })
                trajectories.extend([
                    {
                        **item,
                        "alpha_mmol_gDW": alpha,
                        "alpha_replicate_id": replicate_id,
                        "pool_multiplier": multiplier,
                    }
                    for item in trace
                ])
    calls_path = out / f"chunk_{args.chunk_index:03d}_calls.tsv"
    trajectory_path = out / f"chunk_{args.chunk_index:03d}_trajectory.tsv"
    pd.DataFrame(rows).to_csv(calls_path, sep="\t", index=False)
    pd.DataFrame(trajectories).to_csv(trajectory_path, sep="\t", index=False)
    output_files = {
        "calls": {"filename": calls_path.name, "sha256": sha256_file(calls_path)},
        "trajectory": {"filename": trajectory_path.name, "sha256": sha256_file(trajectory_path)},
    }
    if topology_audit_rows:
        topology_audit_path = out / f"chunk_{args.chunk_index:03d}_runtime_topology_audit.tsv"
        pd.DataFrame(topology_audit_rows).to_csv(topology_audit_path, sep="\t", index=False)
        output_files["runtime_topology_audit"] = {
            "filename": topology_audit_path.name, "sha256": sha256_file(topology_audit_path),
        }
    manifest = {
        "workflow": WORKFLOW, "schema_version": SCHEMA_VERSION, "run_id": args.run_id,
        "chunk_index": args.chunk_index, "chunk_count": args.chunk_count, "genes": genes,
        "solver": args.solver, "runtime_versions": _software_versions(args.solver),
        "solver_feasibility_tolerance": solver_feasibility_tolerance,
        "script_sha256": sha256_file(Path(__file__)),
        "coq9_dilution_tool_sha256": sha256_file(Path(__file__).with_name("coq9_dilution.py")),
        "hours": args.hours,
        "dt_h": args.dt,
        "initial_biomass_gDW_L": args.initial_biomass,
        "alphas_mmol_gDW": [item["alpha_mmol_gDW"] for item in alpha_samples],
        "alpha_sampling": alpha_samples,
        "pool_multipliers": pool_multipliers,
        "uracil_mode": args.uracil_mode,
        "initial_pools_mmol_L": {**INITIAL_POOLS_MMOL_L, "R1354": None} if args.uracil_mode == "po1f_nonlimiting" else INITIAL_POOLS_MMOL_L,
        "base_r1354_bound_mmol_gDW_h": base_r1354_bound,
        "q9_reserve_definition": "alpha * initial_biomass * pool_multiplier mmol/L",
        "q9_reserve_policy": Q9_RESERVE_POLICY,
        "optimizer": "pfba",
        "nonoptimal_policy": "terminal_record_and_stop",
        "calibration_status": CALIBRATION_STATUS,
        "coq9_dilution": [
            runtime_coq9_dilution_record(
                float(item["alpha_mmol_gDW"]), pool_source_enabled=True
            )
            for item in alpha_samples
        ],
        "runtime_topology": runtime_topology or {"enabled": False},
        "runtime_gpr_scenario": runtime_gpr_scenario,
        "simulation_context": context.provenance(),
        "input_sha256": {"model": context.canonical_model_sha256, "medium": sha256_file(media_path), "profile": sha256_file(profile_path), "experimental": sha256_file(experimental_path)},
        "output_files": output_files,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest["fingerprint"] = sha256_payload(manifest)
    (out / f"chunk_{args.chunk_index:03d}_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def _merge(args: argparse.Namespace) -> Path:
    out = Path(args.output_dir).resolve()
    manifests = sorted(out.glob("chunk_*_manifest.json"))
    if not manifests:
        raise FileNotFoundError(f"No chunk manifests in {out}")
    payloads = sorted(
        (json.loads(path.read_text(encoding="utf-8")) for path in manifests),
        key=lambda payload: payload["chunk_index"],
    )
    for payload in payloads:
        observed_fingerprint = payload.get("fingerprint")
        expected_fingerprint = sha256_payload({
            key: value for key, value in payload.items() if key != "fingerprint"
        })
        if observed_fingerprint != expected_fingerprint:
            raise ValueError(f"Chunk {payload['chunk_index']} manifest fingerprint mismatch")
    fingerprints = {payload["fingerprint"] for payload in payloads}
    comparable = [{key: value for key, value in payload.items() if key not in {"chunk_index", "genes", "output_files", "created_at_utc", "fingerprint"}} for payload in payloads]
    if len({sha256_payload(item) for item in comparable}) != 1:
        raise ValueError("Chunk manifests do not share one simulation context")
    expected = set(range(payloads[0]["chunk_count"]))
    observed = {payload["chunk_index"] for payload in payloads}
    if len(payloads) != len(expected) or observed != expected:
        raise ValueError(f"Incomplete chunks: missing {sorted(expected - observed)}")
    calls_paths = []
    for payload in payloads:
        index = payload["chunk_index"]
        outputs = payload.get("output_files")
        if not isinstance(outputs, dict) or not {"calls", "trajectory"} <= set(outputs):
            raise ValueError(f"Chunk {index} is missing output file descriptors")
        if not set(outputs) <= {"calls", "trajectory", "runtime_topology_audit"}:
            raise ValueError(f"Chunk {index} has unknown output file descriptors")
        topology_enabled = bool(payload.get("runtime_topology", {}).get("enabled"))
        if topology_enabled != ("runtime_topology_audit" in outputs):
            raise ValueError(f"Chunk {index} runtime topology output descriptor mismatch")
        for kind, descriptor in outputs.items():
            expected_filename = f"chunk_{index:03d}_{kind}.tsv"
            if not isinstance(descriptor, dict) or set(descriptor) != {"filename", "sha256"}:
                raise ValueError(f"Chunk {index} has invalid {kind} output descriptor")
            if descriptor["filename"] != expected_filename:
                raise ValueError(f"Chunk {index} {kind} filename does not match its manifest")
            path = out / expected_filename
            if not path.is_file() or sha256_file(path) != descriptor["sha256"]:
                raise ValueError(f"Chunk {index} {kind} output hash mismatch")
            if kind == "calls":
                calls_paths.append(path)
    calls = pd.concat([pd.read_csv(path, sep="\t") for path in calls_paths], ignore_index=True)
    calls.to_csv(out / "essentiality_dynamic_calls.tsv", sep="\t", index=False)
    summary = {
        "workflow": WORKFLOW,
        "chunks": len(manifests),
        "chunk_fingerprints": sorted(fingerprints),
        "calls": len(calls),
        "calibration_status": CALIBRATION_STATUS,
        "alpha_sampling": payloads[0]["alpha_sampling"],
        "coq9_dilution": payloads[0]["coq9_dilution"],
        "coq9_dilution_tool_sha256": payloads[0]["coq9_dilution_tool_sha256"],
        "q9_reserve_policy": payloads[0]["q9_reserve_policy"],
        "runtime_gpr_scenario": payloads[0].get("runtime_gpr_scenario", {"enabled": False}),
    }
    (out / "merge_manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-root", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--run-id", default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    parser.add_argument("--solver", default="glpk")
    parser.add_argument(
        "--feasibility-tol", type=_gurobi_feasibility_tolerance,
        help="explicit Gurobi primal feasibility tolerance; recorded in the manifest",
    )
    parser.add_argument("--genes", help="comma-separated pilot genes; defaults to every model gene")
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--chunk-count", type=int, default=1)
    parser.add_argument("--hours", type=float, default=24.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--initial-biomass", type=float, default=0.01)
    parser.add_argument(
        "--alphas",
        type=float,
        nargs="+",
        help="explicit CoQ9 coefficients; defaults to the legacy 1e-6, 1e-4, 1e-3 grid",
    )
    parser.add_argument(
        "--alpha-seed",
        type=int,
        help="sample reproducible CoQ9 coefficients in log space; mutually exclusive with --alphas",
    )
    parser.add_argument(
        "--alpha-replicates",
        type=int,
        default=1,
        help="number of seeded alpha replicates; all WT and KO simulations share each replicate value",
    )
    parser.add_argument("--pool-multipliers", type=float, nargs="+", default=list(DEFAULT_POOL_MULTIPLIERS))
    parser.add_argument("--uracil-mode", choices=URACIL_MODES, default="finite_batch")
    parser.add_argument(
        "--r39-r19-runtime-topology", action="store_true",
        help="apply the reviewed R39 mitochondrial/R19 split convention only to disposable runtime copies",
    )
    parser.add_argument(
        "--runtime-gpr-scenario", choices=tuple(RUNTIME_GPR_SCENARIOS),
        help="apply one preregistered runtime-only candidate GPR hypothesis after the R39/R19 topology",
    )
    parser.add_argument("--merge", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.merge:
        if not args.output_dir:
            raise ValueError("--merge requires --output-dir")
        print(_merge(args))
    else:
        print(_run_chunk(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
