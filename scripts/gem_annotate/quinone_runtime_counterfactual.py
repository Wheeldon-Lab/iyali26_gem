"""Runtime-only quinone reachability and essentiality counterfactuals.

This experiment loads the released B-group model, applies SD-Leu and the PO1f
overlay, then changes only disposable in-memory copies.  It never writes SBML,
curated tables, the FN ledger, or Obsidian.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis

from .config import load_project_paths
from .essentiality_simulation_context import (
    load_effective_simulation_context,
    sha256_file,
    sha256_payload,
)
from .trna_biomass_pipeline import EXPECTED_PO1F_SCREEN
from .validate_essential_genes import (
    DEFAULT_CUTOFFS,
    PRIMARY_CUTOFF,
    build_per_gene_table,
    load_experimental,
    run_single_gene_deletions,
)


WORKFLOW = "quinone_runtime_counterfactual"
SCHEMA_VERSION = "1.0"
EXPERIMENT_DATE = "2026-08-06"
EXPECTED_INPUT_SHA256 = {
    "model": "0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415",
    "medium": "ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85",
    "profile": "35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383",
    "overlay_effect": "d15acbde9438f5d2391c4da23705a34a3585833062d616517d5af052088606c2",
    "experimental": "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b",
}

ROUTE_IDS = (
    "R763",
    "R407",
    "R969",
    "R39",
    "R808",
    "R715",
    "R40",
    "R19",
    "R18",
    "R695",
    "R385",
)
REPEATED_AND_REACTIONS = ("R715", "R19", "R18", "R695", "R385")
REPEATED_AND_GENES = {
    "YALI1F34625g",
    "YALI1B20527g",
    "YALI1A08781g",
    "YALI1F34675g",
    "YALI1C25352g",
    "YALI1B20835g",
    "YALI1E18269g",
}
STEP_SPECIFIC_GPRS = {
    "R763": "YALI1C26017g",
    "R407": "YALI1F08349g",
    "R39": "YALI1A08781g",
    "R715": "YALI1B20835g",
    "R40": "YALI1F34625g",
    "R19": "YALI1A08781g",
    "R18": "YALI1C25352g",
    "R695": "YALI1E18269g",
    "R385": "YALI1B20835g",
}

Q9_CHEMISTRY = {
    "m640[C_mi]": ("nonaprenyl diphosphate", "C45H76O7P2"),
    "m641[C_mi]": ("3-nonaprenyl-4-hydroxybenzoic acid", "C52H78O3"),
    "m108[C_cy]": ("3-nonaprenyl-4-hydroxybenzoic acid", "C52H78O3"),
    "m110[C_cy]": ("3-nonaprenyl-4,5-dihydroxybenzoic acid", "C52H77O4"),
    "m939[C_mi]": ("3-nonaprenyl-4,5-dihydroxybenzoic acid", "C52H77O4"),
    "m111[C_mi]": ("3-nonaprenyl-4-hydroxy-5-methoxybenzoic acid", "C53H79O4"),
    "m63[C_mi]": ("2-nonaprenyl-6-methoxyphenol", "C52H80O2"),
    "m59[C_mi]": ("2-nonaprenyl-6-methoxy-1,4-benzoquinone", "C52H78O3"),
    "m61[C_mi]": (
        "2-nonaprenyl-6-methoxy-3-methyl-1,4-benzoquinone",
        "C53H80O3",
    ),
    "m611[C_mi]": (
        "2-nonaprenyl-5-hydroxy-6-methoxy-3-methyl-1,4-benzoquinone",
        "C53H80O4",
    ),
    "m468[C_mi]": ("ubiquinone-9", "C54H82O4"),
    "m471[C_mi]": ("ubiquinol-9", "C54H84O4"),
}

DILUTION_COEFFICIENTS = (0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
SCREEN_DILUTION_COEFFICIENT = 1e-3

GENE_INFO = {
    "YALI1C26017g": {
        "symbol": "COQ1 candidate",
        "function": "CoQ9 side-chain nonaprenyl-diphosphate synthase candidate",
        "evidence_status": "experimentally verified in a heterologous catalytic-core assay; native localization unverified",
        "model_role": "R763",
    },
    "YALI1F08349g": {
        "symbol": "COQ2 candidate",
        "function": "4-hydroxybenzoate polyprenyltransferase candidate",
        "evidence_status": "uncharacterized in native Yarrowia; homology-supported annotation",
        "model_role": "R407",
    },
    "YALI1B20835g": {
        "symbol": "COQ3 candidate",
        "function": "ubiquinone-biosynthesis O-methyltransferase candidate",
        "evidence_status": "uncharacterized in native Yarrowia; comparative step-function evidence",
        "model_role": "R715 and R385 direct-catalyst candidate",
    },
    "YALI1F34625g": {
        "symbol": "COQ4 candidate",
        "function": "CoQ-ring oxidative decarboxylase/scaffold candidate",
        "evidence_status": "curated annotation; native locus function not experimentally verified",
        "model_role": "R40 direct-catalyst candidate",
    },
    "YALI1C25352g": {
        "symbol": "COQ5 candidate",
        "function": "CoQ-ring C-methyltransferase candidate",
        "evidence_status": "uncharacterized in native Yarrowia; comparative step-function evidence",
        "model_role": "R18 direct-catalyst candidate",
    },
    "YALI1A08781g": {
        "symbol": "COQ6 candidate",
        "function": "FAD monooxygenase candidate for two CoQ hydroxylations",
        "evidence_status": "uncharacterized in native Yarrowia; comparative step-function evidence",
        "model_role": "R39 and R19 direct-catalyst candidate",
    },
    "YALI1E18269g": {
        "symbol": "COQ7 candidate",
        "function": "demethoxyubiquinone hydroxylase candidate",
        "evidence_status": "uncharacterized in native Yarrowia; comparative step-function evidence",
        "model_role": "R695 direct-catalyst candidate",
    },
    "YALI1B20527g": {
        "symbol": "COQ8 candidate",
        "function": "ADCK-family regulatory ATPase/kinase candidate",
        "evidence_status": "uncharacterized in native Yarrowia; conserved-family evidence",
        "model_role": "member of current repeated AND; not a direct atom-transfer candidate",
    },
    "YALI1F34675g": {
        "symbol": "COQ9 candidate",
        "function": "CoQ-synthome-associated accessory candidate",
        "evidence_status": "uncharacterized in native Yarrowia; conserved-family evidence",
        "model_role": "member of current repeated AND; possible Coq6/Coq7 accessory",
    },
}


def _experiment_design() -> dict[str, Any]:
    return {
        "route_ids": list(ROUTE_IDS),
        "repeated_and_reactions": list(REPEATED_AND_REACTIONS),
        "repeated_and_genes": sorted(REPEATED_AND_GENES),
        "step_specific_gprs": STEP_SPECIFIC_GPRS,
        "q9_chemistry": Q9_CHEMISTRY,
        "q9_R763_stoichiometry": (
            "4 isopentenyl diphosphate + pentaprenyl diphosphate -> "
            "4 diphosphate + nonaprenyl diphosphate"
        ),
        "oxidized_terminal_R385_variant": (
            "SAM + 3-demethylubiquinone -> SAH + ubiquinone"
        ),
        "diagnostic_demand_bounds": [0.0, 1000.0],
        "dilution_coefficients_mmol_per_gDW": list(DILUTION_COEFFICIENTS),
        "screen_dilution_coefficient_mmol_per_gDW": (
            SCREEN_DILUTION_COEFFICIENT
        ),
        "screen_cutoffs_fraction_of_wt": list(DEFAULT_CUTOFFS),
        "screen_positive_only": True,
        "screen_scenarios": [
            "Q0_B_PO1f",
            "Q9_current_repeated_AND",
            "Q9_step_specific_GPR",
        ],
    }


def _software_versions(solver: str) -> dict[str, str]:
    versions = {
        "python": sys.version.split()[0],
        "pandas": pd.__version__,
        "cobra": __import__("cobra").__version__,
        "solver": solver,
    }
    if solver.lower() == "gurobi":
        import gurobipy

        versions["gurobi"] = ".".join(map(str, gurobipy.gurobi.version()))
    return versions


def _validate_output_dir(output_dir: Path) -> None:
    resolved = output_dir.resolve()
    project_experiments = Path(__file__).resolve().parents[2] / "docs" / "experiments"
    allowed_roots = (project_experiments.resolve(), Path("/tmp").resolve())
    if not any(resolved == root or root in resolved.parents for root in allowed_roots):
        raise ValueError(
            "Output must be inside docs/experiments or /tmp for this runtime-only workflow"
        )
    if not resolved.name.startswith("quinone_runtime_counterfactual_"):
        raise ValueError("Output directory name must start with quinone_runtime_counterfactual_")
    if resolved.exists() and any(resolved.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output directory: {resolved}")


def _stoichiometry(reaction: Any) -> dict[str, float]:
    return {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in sorted(
            reaction.metabolites.items(), key=lambda item: item[0].id
        )
    }


def _balance(reaction: Any) -> dict[str, float]:
    return {
        element: float(value)
        for element, value in sorted(reaction.check_mass_balance().items())
    }


def _stable_ratio(value: Any) -> float:
    """Remove solver-scale noise without changing threshold interpretation."""
    ratio = float(value)
    if abs(ratio) < 1e-12:
        return 0.0
    if abs(ratio - 1.0) < 1e-12:
        return 1.0
    return round(ratio, 12)


def _assert_frozen_base(model: Any) -> None:
    if model.objective.expression.as_coefficients_dict() == {}:
        raise ValueError("The effective model has no objective")
    if model.reactions.get_by_id("biomass_C").notes.get(
        "experimental_trna_biomass_mode"
    ) != "split_v1":
        raise ValueError("The runtime experiment requires the formal B-group model")
    split_count = sum(reaction.id.startswith("TRNA_BIOMASS_") for reaction in model.reactions)
    if split_count != 20:
        raise ValueError(f"Expected 20 B-group split reactions, found {split_count}")
    expected_r763 = {
        "m204[C_mi]": 1.0,
        "m640[C_mi]": 1.0,
        "m984[C_mi]": -1.0,
        "m985[C_mi]": -1.0,
    }
    expected_r385 = {
        "m28[C_mi]": -1.0,
        "m471[C_mi]": 1.0,
        "m60[C_mi]": -1.0,
        "m611[C_mi]": -1.0,
        "m62[C_mi]": 1.0,
    }
    if _stoichiometry(model.reactions.get_by_id("R763")) != expected_r763:
        raise ValueError("R763 no longer matches the frozen counterfactual input")
    if _stoichiometry(model.reactions.get_by_id("R385")) != expected_r385:
        raise ValueError("R385 no longer matches the frozen counterfactual input")
    for reaction_id in REPEATED_AND_REACTIONS:
        genes = {gene.id for gene in model.reactions.get_by_id(reaction_id).genes}
        if genes != REPEATED_AND_GENES:
            raise ValueError(f"{reaction_id} repeated-AND input changed: {sorted(genes)}")


def _apply_oxidized_terminal_q6_variant(model: Any) -> None:
    """Apply an atom-balanced oxidized endpoint; native redox state is unresolved."""
    reaction = model.reactions.get_by_id("R385")
    reaction.add_metabolites(
        {
            model.metabolites.get_by_id("m28[C_mi]"): 1.0,
            model.metabolites.get_by_id("m471[C_mi]"): -1.0,
            model.metabolites.get_by_id("m468[C_mi]"): 1.0,
        }
    )
    reaction.name = (
        "3-demethylubiquinone-6 O-methyltransferase "
        "(oxidized-terminal counterfactual)"
    )
    if _balance(reaction):
        raise ValueError(f"Balanced-Q6 R385 remains imbalanced: {_balance(reaction)}")


def _apply_q9_chemistry(model: Any) -> None:
    """Reinterpret the frozen Q6 objects as a balanced, runtime-only Q9 route."""
    impacted = sorted(
        {
            reaction.id
            for metabolite_id in Q9_CHEMISTRY
            for reaction in model.metabolites.get_by_id(metabolite_id).reactions
        }
    )
    before = {reaction_id: _balance(model.reactions.get_by_id(reaction_id)) for reaction_id in impacted}

    for metabolite_id, (name, formula) in Q9_CHEMISTRY.items():
        metabolite = model.metabolites.get_by_id(metabolite_id)
        metabolite.name = name
        metabolite.formula = formula

    model.reactions.get_by_id("R763").add_metabolites(
        {
            model.metabolites.get_by_id("m984[C_mi]"): -3.0,
            model.metabolites.get_by_id("m204[C_mi]"): 3.0,
        }
    )
    model.reactions.get_by_id("R763").name = (
        "nonaprenyl-diphosphate synthase, four-IPP lump (counterfactual)"
    )
    _apply_oxidized_terminal_q6_variant(model)
    model.reactions.get_by_id("R385").name = (
        "3-demethylubiquinone-9 O-methyltransferase (counterfactual)"
    )

    for reaction_id in ROUTE_IDS:
        imbalance = _balance(model.reactions.get_by_id(reaction_id))
        if imbalance:
            raise ValueError(f"Q9 route reaction {reaction_id} is imbalanced: {imbalance}")
    for reaction_id in impacted:
        if reaction_id in {"R763", "R385"}:
            continue
        after = _balance(model.reactions.get_by_id(reaction_id))
        if after != before[reaction_id]:
            raise ValueError(
                f"Q9 substitution introduced a new balance change in {reaction_id}: "
                f"{before[reaction_id]} -> {after}"
            )


def _apply_step_specific_gprs(model: Any) -> None:
    for reaction_id, rule in STEP_SPECIFIC_GPRS.items():
        model.reactions.get_by_id(reaction_id).gene_reaction_rule = rule


def _add_demand(model: Any, metabolite_id: str, reaction_id: str) -> Any:
    if model.reactions.has_id(reaction_id):
        raise ValueError(f"Counterfactual reaction already exists: {reaction_id}")
    reaction = Reaction(reaction_id)
    reaction.name = f"runtime-only demand for {metabolite_id}"
    reaction.bounds = (0.0, 1000.0)
    reaction.add_metabolites({model.metabolites.get_by_id(metabolite_id): -1.0})
    model.add_reactions([reaction])
    return reaction


def _couple_demand_to_growth(model: Any, demand: Any, coefficient: float) -> None:
    biomass = model.reactions.get_by_id("biomass_C")
    constraint = model.problem.Constraint(
        demand.flux_expression - float(coefficient) * biomass.flux_expression,
        lb=0.0,
        ub=0.0,
        name="CF_Q9_GROWTH_DILUTION",
    )
    model.add_cons_vars([constraint])


def _optimize(model: Any) -> Any:
    solution = model.optimize()
    if solution.status != "optimal" or solution.objective_value is None:
        raise RuntimeError(f"Counterfactual optimization failed: {solution.status}")
    return solution


def _maximum_demand(
    base_model: Any,
    *,
    scenario: str,
    transform: str,
    demand_metabolite: str,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    model = base_model.copy()
    if transform == "oxidized_q6_variant":
        _apply_oxidized_terminal_q6_variant(model)
    elif transform == "q9":
        _apply_q9_chemistry(model)
    elif transform != "legacy_q6":
        raise ValueError(f"Unknown chemistry transform: {transform}")
    demand = _add_demand(model, demand_metabolite, f"CF_DM_{scenario}")
    model.objective = demand
    solution = _optimize(model)
    route_fluxes = [
        {
            "scenario": scenario,
            "reaction_id": reaction_id,
            "flux": float(solution.fluxes[reaction_id]),
        }
        for reaction_id in ROUTE_IDS
    ]
    chemistry = [
        {
            "scenario": scenario,
            "reaction_id": reaction_id,
            "equation": model.reactions.get_by_id(reaction_id).reaction,
            "mass_balance": json.dumps(
                _balance(model.reactions.get_by_id(reaction_id)),
                sort_keys=True,
                separators=(",", ":"),
            ),
        }
        for reaction_id in ROUTE_IDS
    ]

    blocked = model.copy()
    blocked.reactions.get_by_id("R385").bounds = (0.0, 0.0)
    blocked.objective = blocked.reactions.get_by_id(demand.id)
    blocked_solution = _optimize(blocked)
    blocked_value = float(blocked_solution.objective_value)
    if abs(blocked_value) > 1e-8:
        raise ValueError(
            f"{scenario} demand remains producible with R385 closed: {blocked_value}"
        )
    summary = {
        "scenario": scenario,
        "chemistry": transform,
        "demand_metabolite": demand_metabolite,
        "maximum_type": "zero_growth_stoichiometric_reachability_maximum",
        "maximum_demand": float(solution.objective_value),
        "biomass_at_maximum_demand": float(solution.fluxes["biomass_C"]),
        "R385_flux": float(solution.fluxes["R385"]),
        "maximum_demand_with_R385_closed": blocked_value,
    }
    return summary, route_fluxes, chemistry


def _dilution_scan(base_model: Any) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for coefficient in DILUTION_COEFFICIENTS:
        model = base_model.copy()
        _apply_q9_chemistry(model)
        demand = _add_demand(model, "m468[C_mi]", "CF_DM_Q9_DILUTION")
        _couple_demand_to_growth(model, demand, coefficient)
        model.objective = model.reactions.get_by_id("biomass_C")
        solution = _optimize(model)
        growth = float(solution.objective_value)
        demand_flux = float(solution.fluxes[demand.id])
        expected = float(coefficient) * growth
        if abs(demand_flux - expected) > 1e-8:
            raise ValueError(
                f"Growth coupling failed at cQ={coefficient}: {demand_flux} != {expected}"
            )
        rows.append(
            {
                "cQ_mmol_per_gDW": coefficient,
                "growth_h-1": growth,
                "demand_flux_mmol_per_gDW_h": demand_flux,
                "R763_flux": float(solution.fluxes["R763"]),
                "R385_flux": float(solution.fluxes["R385"]),
                "status": solution.status,
                "interpretation": "sensitivity_only_not_calibrated",
            }
        )
    return pd.DataFrame(rows)


def _screen(
    model: Any,
    *,
    scenario: str,
    experimental: pd.DataFrame,
    excluded_gene_ids: tuple[str, ...],
    solver: str,
) -> tuple[pd.DataFrame, list[dict[str, Any]], float]:
    predictions, wt_growth = run_single_gene_deletions(
        model,
        solver,
        excluded_gene_ids=excluded_gene_ids,
    )
    per_gene = build_per_gene_table(
        experimental,
        predictions,
        DEFAULT_CUTOFFS,
        PRIMARY_CUTOFF,
    )
    positives = per_gene[
        per_gene["experimental_essential"].eq(True) & per_gene["in_model"].eq(True)
    ]
    curve: list[dict[str, Any]] = []
    stable_ratios = positives["ko_growth_ratio"].map(_stable_ratio)
    for cutoff in DEFAULT_CUTOFFS:
        tp = int((stable_ratios < cutoff).sum())
        fn = int(len(positives) - tp)
        curve.append(
            {
                "scenario": scenario,
                "cutoff_fraction_of_wt": cutoff,
                "TP": tp,
                "FN": fn,
                "recall": tp / len(positives),
                "wt_growth_h-1": wt_growth,
            }
        )
    return per_gene, curve, wt_growth


def _screen_model(base_model: Any, *, step_specific: bool, coefficient: float) -> Any:
    model = base_model.copy()
    _apply_q9_chemistry(model)
    if step_specific:
        _apply_step_specific_gprs(model)
    demand = _add_demand(model, "m468[C_mi]", "CF_DM_Q9_DILUTION")
    _couple_demand_to_growth(model, demand, coefficient)
    model.objective = model.reactions.get_by_id("biomass_C")
    return model


def _scenario_deltas(
    screens: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    contrasts = (
        (
            "baseline_to_q9_current_AND",
            "Q0_B_PO1f",
            "Q9_current_repeated_AND",
            "Q9_chemistry_plus_growth_demand",
        ),
        (
            "q9_current_AND_to_step_specific",
            "Q9_current_repeated_AND",
            "Q9_step_specific_GPR",
            "route_wide_step_specific_GPR_remapping",
        ),
    )
    rows: list[dict[str, Any]] = []
    for contrast, left_name, right_name, mechanism in contrasts:
        left = screens[left_name].set_index("gene_id")
        right = screens[right_name].set_index("gene_id")
        for gene_id in sorted(set(left.index) & set(right.index)):
            left_row = left.loc[gene_id]
            right_row = right.loc[gene_id]
            experimental_essential = left_row["experimental_essential"]
            in_model = left_row["in_model"]
            if (
                pd.isna(experimental_essential)
                or not bool(experimental_essential)
                or pd.isna(in_model)
                or not bool(in_model)
            ):
                continue
            before_ratio = _stable_ratio(left_row["ko_growth_ratio"])
            after_ratio = _stable_ratio(right_row["ko_growth_ratio"])
            for cutoff in DEFAULT_CUTOFFS:
                before = "TP" if before_ratio < cutoff else "FN"
                after = "TP" if after_ratio < cutoff else "FN"
                rows.append(
                    {
                        "contrast": contrast,
                        "mechanism": mechanism,
                        "cutoff_fraction_of_wt": cutoff,
                        "gene_id": gene_id,
                        "before_ko_growth_ratio": before_ratio,
                        "after_ko_growth_ratio": after_ratio,
                        "before_call": before,
                        "after_call": after,
                        "transition": f"{before}->{after}",
                        "changed": before != after,
                    }
                )
    return pd.DataFrame(rows)


def _quinone_gene_table(screens: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for scenario, table in screens.items():
        indexed = table.set_index("gene_id")
        for gene_id, info in GENE_INFO.items():
            row = indexed.loc[gene_id]
            experimental_value = row["experimental_essential"]
            experimental_essential = (
                None if pd.isna(experimental_value) else bool(experimental_value)
            )
            ko_growth_ratio = _stable_ratio(row["ko_growth_ratio"])
            record: dict[str, Any] = {
                "scenario": scenario,
                "gene_id": gene_id,
                **info,
                "experimental_essential": experimental_essential,
                "ko_growth_ratio": ko_growth_ratio,
            }
            for cutoff in DEFAULT_CUTOFFS:
                record[f"call_{cutoff * 100:g}pct"] = (
                    "TP"
                    if experimental_essential is True
                    and ko_growth_ratio < cutoff
                    else "FN"
                    if experimental_essential is True
                    else "outside_positive_reference"
                )
            rows.append(record)
    return pd.DataFrame(rows)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_report(
    output_path: Path,
    *,
    context: Any,
    baseline_growth: float,
    reachability: pd.DataFrame,
    dilution: pd.DataFrame,
    thresholds: pd.DataFrame,
    changed_counts: dict[str, int],
    positive_reference_total: int,
    in_model_positive_count: int,
) -> None:
    values = reachability.set_index("scenario")
    threshold_lines = [
        "| Scenario | Cutoff | TP | FN | Recall |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in thresholds.itertuples(index=False):
        threshold_lines.append(
            f"| {row.scenario} | {100 * row.cutoff_fraction_of_wt:g}% | "
            f"{row.TP} | {row.FN} | {100 * row.recall:.2f}% |"
        )
    scan_lines = [
        "| cQ (mmol/gDW) | Growth (h⁻¹) | Q9 dilution flux |",
        "|---:|---:|---:|",
    ]
    for row in dilution.to_dict(orient="records"):
        scan_lines.append(
            f"| {row['cQ_mmol_per_gDW']:.6g} | {row['growth_h-1']:.9g} | "
            f"{row['demand_flux_mmol_per_gDW_h']:.9g} |"
        )

    text = f"""## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: run
- Origin Date: {EXPERIMENT_DATE}
- Verification Status: UNVERIFIED (deterministic replay required)
- Version Label: exp_result_v1

# Runtime-only Quinone counterfactual experiment

## Experiment Result

- **ID**: `quinone-runtime-counterfactual-20260806`
- **Type**: simulation
- **Status**: completed
- **Base context fingerprint**: `{context.simulation_context_fingerprint}`
- **Canonical model SHA-256**: `{context.canonical_model_sha256}`
- **Baseline B–PO1f growth**: `{baseline_growth:.12g} h^-1`

## Main diagnosis

The frozen route is inactive because it has no effective terminal net demand:
all route FVA ranges are `[0,0]` without a demand.  An irreversible diagnostic
demand activates the complete route, including the previously omitted connector
`R969` in reverse.

This diagnosis survives two atom-balanced terminal variants:

- Legacy Q6 zero-growth stoichiometric reachability maximum:
  `{values.loc['Q1_legacy_Q6_demand', 'maximum_demand']:.9g}` mmol/gDW/h,
  but legacy `R385` is chemically imbalanced.
- Atom-balanced oxidized-terminal Q6 zero-growth reachability maximum:
  `{values.loc['Q1b_balanced_Q6_demand', 'maximum_demand']:.9g}` mmol/gDW/h
  after making `R385` produce ubiquinone-6.
- Atom-balanced oxidized-terminal Q9 zero-growth reachability maximum:
  `{values.loc['Q2_balanced_Q9_demand', 'maximum_demand']:.9g}` mmol/gDW/h
  after the explicit four-IPP nonaprenyl counterfactual.
- Closing `R385` reduces all three diagnostic demands to zero, so the result is
  not supplied by an alternate modeled source.

Therefore **missing net demand is the immediate cause of inactivity in this
simulation context**. It is not the only model defect: the legacy endpoint is
imbalanced, the biologically appropriate terminal redox microspecies is
unresolved, the chain length is CoQ6 rather than the evidence-supported CoQ9
counterfactual, and the native localization/GPR dependency remains unresolved.

## Growth-dilution sensitivity

The following coefficients are numerical sensitivity points, not fitted
physiological parameters:

{chr(10).join(scan_lines)}

## Essentiality attribution

The full positive-only B–PO1f screen was run for the frozen baseline and for a
Q9 dilution coefficient of `{SCREEN_DILUTION_COEFFICIENT:g} mmol/gDW` under two
GPR interpretations. Recall is TP / `{in_model_positive_count}` experimental
positive genes that map into the screened model. The source positive list has
`{positive_reference_total}` genes in total; genes absent from the model are not
counted as TP or FN.

{chr(10).join(threshold_lines)}

- Positive-list call changes caused by adding balanced Q9 chemistry and demand:
  `{changed_counts['baseline_to_q9_current_AND']}` gene-threshold rows.
- Positive-list call changes caused by route-wide step-specific GPR remapping,
  including adding candidate GPRs to previously ungated `R39/R40` and replacing
  the repeated seven-gene `AND` assignments:
  `{changed_counts['q9_current_AND_to_step_specific']}` gene-threshold rows.

These counterfactual call changes are causal model diagnostics, not evidence
that the proposed GPRs or coefficient are biologically correct.

## Chemistry boundary

- Q9 is built by reinterpreting the existing chain-specific IDs in memory and
  adding `C15H24` to each Q6 intermediate.
- Counterfactual `R763` is
  `4 IPP + pentaprenyl-PP -> 4 PPi + nonaprenyl-PP`.
- The oxidized-terminal counterfactual `R385` is
  `SAM + 3-demethylubiquinone-9 -> SAH + ubiquinone-9`.
- Every reaction in the Q9 route passes elemental and charge balance.
- [KEGG R08781](https://www.kegg.jp/entry/R08781) and
  [Rhea 81218](https://www.rhea-db.org/rhea/81218) represent generic oxidized
  quinone equations, whereas [Rhea 44381](https://www.rhea-db.org/rhea/44381)
  and Q9-specific [Rhea 17049](https://www.rhea-db.org/rhea/17049) represent a
  reduced ubiquinol equation with proton production. These are
  database-compatible alternatives; they do not resolve the native Yarrowia
  terminal redox state or supply a PO1f biomass coefficient.

## Anomalies and limitations

1. The current canonical `R385` has residual `H:+1, charge:-1`.
2. The native terminal substrate/product redox microspecies remains unresolved.
3. No W29/PO1f SD-Leu whole-cell Q9 mmol/gDW value is available.
4. `R39/R969/R808` localization and transport remain unverified in native
   *Yarrowia*.
5. The step-specific GPR is the minimal direct-catalyst interpretation; native
   COQ8/COQ9 accessory necessity remains unresolved.
6. The `0...1000` demand is used only as a reachability objective and is not a
   biological demand coefficient.

## No-write verification

The experiment wrote reports only.  It did not write SBML, curated tables,
formal GPRs/bounds, the durable FN ledger, or Obsidian.  See the frozen evidence
review in `docs/research/quinone_review_2026-08-06/review_report.md`.
"""
    output_path.write_text(text, encoding="utf-8")


def run(*, output_dir: Path, solver: str, research_root: Path | None = None) -> dict[str, Any]:
    _validate_output_dir(output_dir)
    paths = load_project_paths(research_root=research_root, required=True)
    model_path = paths.output_model
    medium_path = paths.media / "sd_leu.csv"
    profile_path = paths.strain_profiles / "po1f_sd_leu.json"
    experimental_path = paths.essentiality / "consensus_essential_genes.csv"

    context = load_effective_simulation_context(
        model_path=model_path,
        media_path=medium_path,
        strain_profile_path=profile_path,
    )
    actual_hashes = {
        "model": context.canonical_model_sha256,
        "medium": context.medium_sha256,
        "profile": context.strain_profile_sha256,
        "overlay_effect": context.strain_overlay_effect_sha256,
        "experimental": sha256_file(experimental_path),
    }
    if actual_hashes != EXPECTED_INPUT_SHA256:
        raise ValueError(
            "Frozen Quinone counterfactual inputs changed: "
            + json.dumps(actual_hashes, sort_keys=True)
        )
    base_model = context.model
    base_model.solver = solver
    script_sha256 = sha256_file(Path(__file__))
    experiment_design = _experiment_design()
    experiment_design_sha256 = sha256_payload(experiment_design)
    software_versions = _software_versions(solver)
    _assert_frozen_base(base_model)
    baseline_solution = _optimize(base_model)
    baseline_growth = float(baseline_solution.objective_value)

    fva_rows: list[dict[str, Any]] = []
    for fraction, label in ((0.0, "structural"), (1.0, "at_optimum")):
        fva = flux_variability_analysis(
            base_model,
            reaction_list=list(ROUTE_IDS),
            fraction_of_optimum=fraction,
            processes=1,
        )
        for reaction_id, row in fva.iterrows():
            fva_rows.append(
                {
                    "scenario": "Q0_B_PO1f_no_demand",
                    "objective_fraction": fraction,
                    "analysis": label,
                    "reaction_id": reaction_id,
                    "minimum": float(row["minimum"]),
                    "maximum": float(row["maximum"]),
                }
            )
    if any(abs(row[key]) > 1e-8 for row in fva_rows for key in ("minimum", "maximum")):
        raise ValueError("Frozen Q0 quinone route is no longer [0,0]")

    reachability_rows: list[dict[str, Any]] = []
    route_flux_rows: list[dict[str, Any]] = []
    chemistry_rows: list[dict[str, Any]] = []
    for scenario, transform, metabolite_id in (
        ("Q1_legacy_Q6_demand", "legacy_q6", "m471[C_mi]"),
        ("Q1b_balanced_Q6_demand", "oxidized_q6_variant", "m468[C_mi]"),
        ("Q2_balanced_Q9_demand", "q9", "m468[C_mi]"),
    ):
        summary, fluxes, chemistry = _maximum_demand(
            base_model,
            scenario=scenario,
            transform=transform,
            demand_metabolite=metabolite_id,
        )
        summary["scenario_fingerprint"] = sha256_payload(
            {
                "base_simulation_context_fingerprint": context.simulation_context_fingerprint,
                "script_sha256": script_sha256,
                "experiment_design_sha256": experiment_design_sha256,
                "scenario": scenario,
                "chemistry": transform,
                "demand_bounds": [0.0, 1000.0],
                "demand_metabolite": metabolite_id,
                "q9_chemistry": Q9_CHEMISTRY if transform == "q9" else None,
            }
        )
        reachability_rows.append(summary)
        route_flux_rows.extend(fluxes)
        chemistry_rows.extend(chemistry)

    dilution = _dilution_scan(base_model)
    experimental = load_experimental(experimental_path, positive_only=True)
    positive_reference_total = int(len(experimental))
    screen_models = {
        "Q0_B_PO1f": base_model.copy(),
        "Q9_current_repeated_AND": _screen_model(
            base_model,
            step_specific=False,
            coefficient=SCREEN_DILUTION_COEFFICIENT,
        ),
        "Q9_step_specific_GPR": _screen_model(
            base_model,
            step_specific=True,
            coefficient=SCREEN_DILUTION_COEFFICIENT,
        ),
    }
    screens: dict[str, pd.DataFrame] = {}
    threshold_rows: list[dict[str, Any]] = []
    screen_wt: dict[str, float] = {}
    for scenario, model in screen_models.items():
        per_gene, curve, wt_growth = _screen(
            model,
            scenario=scenario,
            experimental=experimental,
            excluded_gene_ids=context.excluded_runtime_genes,
            solver=solver,
        )
        screens[scenario] = per_gene
        threshold_rows.extend(curve)
        screen_wt[scenario] = wt_growth
    thresholds = pd.DataFrame(threshold_rows)
    baseline_per_gene = screens["Q0_B_PO1f"]
    in_model_positive_count = int(
        (
            baseline_per_gene["experimental_essential"].eq(True)
            & baseline_per_gene["in_model"].eq(True)
        ).sum()
    )
    baseline_curve = thresholds[thresholds["scenario"].eq("Q0_B_PO1f")]
    for row in baseline_curve.itertuples(index=False):
        expected = EXPECTED_PO1F_SCREEN[float(row.cutoff_fraction_of_wt)]
        if (int(row.TP), int(row.FN)) != expected:
            raise ValueError(
                f"Baseline essentiality changed at {row.cutoff_fraction_of_wt}: "
                f"{(row.TP, row.FN)} != {expected}"
            )

    deltas = _scenario_deltas(screens)
    gene_table = _quinone_gene_table(screens)
    changed_counts = {
        contrast: int(
            deltas[deltas["contrast"].eq(contrast)]["changed"].astype(bool).sum()
        )
        for contrast in deltas["contrast"].unique()
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        "reachability.tsv": pd.DataFrame(reachability_rows),
        "route_flux_at_max_demand.tsv": pd.DataFrame(route_flux_rows),
        "route_fva.tsv": pd.DataFrame(fva_rows),
        "chemistry_audit.tsv": pd.DataFrame(chemistry_rows),
        "dilution_scan.tsv": dilution,
        "essentiality_threshold_summary.tsv": thresholds,
        "essentiality_scenario_deltas.tsv": deltas,
        "quinone_gene_effects.tsv": gene_table,
    }
    for name, table in tables.items():
        table.to_csv(output_dir / name, sep="\t", index=False)

    result = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "experiment_date": EXPERIMENT_DATE,
        "status": "completed",
        "verification_status": "reproducibility_pending",
        "solver": solver,
        "script_sha256": script_sha256,
        "experiment_design": experiment_design,
        "experiment_design_sha256": experiment_design_sha256,
        "software_versions": software_versions,
        "input_hashes": actual_hashes,
        "base_context": {
            **context.provenance(),
            "canonical_model_sha256": context.canonical_model_sha256,
            "medium_sha256": context.medium_sha256,
            "experimental_sha256": actual_hashes["experimental"],
            "B_group_split_reactions": 20,
        },
        "baseline_growth_h-1": baseline_growth,
        "reachability": reachability_rows,
        "dilution_coefficients_mmol_per_gDW": list(DILUTION_COEFFICIENTS),
        "screen_dilution_coefficient_mmol_per_gDW": SCREEN_DILUTION_COEFFICIENT,
        "screen_wt_growth_h-1": screen_wt,
        "positive_reference_total": positive_reference_total,
        "in_model_positive_count": in_model_positive_count,
        "changed_gene_threshold_rows": changed_counts,
        "scientific_boundary": {
            "diagnosis": "missing_terminal_net_demand_is_immediate_modeled_cause",
            "legacy_R385": "mass_and_charge_imbalanced",
            "Q9_chemistry": (
                "atom_and_charge_balanced_oxidized_terminal_variant_not_curated_patch"
            ),
            "terminal_redox_state": "unresolved_oxidized_and_reduced_database_variants",
            "dilution_coefficients": "sensitivity_only_not_calibrated",
            "step_specific_GPR": "minimal_direct_catalyst_counterfactual_not_native_dependency_claim",
            "step_specific_GPR_scope": (
                "route_wide_remapping_including_previously_ungated_R39_and_R40"
            ),
            "maximum_demand_values": "zero_growth_stoichiometric_reachability_only",
        },
    }
    (output_dir / "experiment_result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    _write_report(
        output_dir / "experiment_result.md",
        context=context,
        baseline_growth=baseline_growth,
        reachability=tables["reachability.tsv"],
        dilution=dilution,
        thresholds=thresholds,
        changed_counts=changed_counts,
        positive_reference_total=positive_reference_total,
        in_model_positive_count=in_model_positive_count,
    )

    artifact_hashes = {
        path.name: _file_sha256(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path.name != "manifest.json"
    }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "experiment_date": EXPERIMENT_DATE,
        "status": "completed",
        "verification_status": "reproducibility_pending",
        "script_sha256": script_sha256,
        "experiment_design_sha256": experiment_design_sha256,
        "software_versions": software_versions,
        "input_hashes": actual_hashes,
        "base_context_fingerprint": context.simulation_context_fingerprint,
        "scenario_fingerprints": {
            row["scenario"]: row["scenario_fingerprint"] for row in reachability_rows
        }
        | {
            scenario: sha256_payload(
                {
                    "base_simulation_context_fingerprint": (
                        context.simulation_context_fingerprint
                    ),
                    "script_sha256": script_sha256,
                    "experiment_design_sha256": experiment_design_sha256,
                    "scenario": scenario,
                }
            )
            for scenario in [
                "Q0_B_PO1f",
                "Q9_current_repeated_AND",
                "Q9_step_specific_GPR",
                "Q9_dilution_scan",
            ]
        },
        "artifact_hashes": artifact_hashes,
        "writes_forbidden": [
            "model.xml",
            "data/iyali26.xml",
            "curated tables",
            "durable FN ledger",
            "Obsidian",
        ],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--solver", default="gurobi")
    parser.add_argument("--research-root", type=Path)
    args = parser.parse_args()
    result = run(
        output_dir=args.output_dir.resolve(),
        solver=args.solver,
        research_root=args.research_root,
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "diagnosis": result["scientific_boundary"]["diagnosis"],
                "baseline_growth_h-1": result["baseline_growth_h-1"],
                "changed_gene_threshold_rows": result["changed_gene_threshold_rows"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
