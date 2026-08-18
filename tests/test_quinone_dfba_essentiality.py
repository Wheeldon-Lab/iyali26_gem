import math
from pathlib import Path

import pandas as pd
import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model
from cobra.manipulation.delete import knock_out_model_genes

from scripts.gem_annotate import quinone_dfba_essentiality as dfba
from scripts.gem_annotate.summarize_quinone_dfba import FIVE_COQ_CONTROLS, summarize_calls
from scripts.gem_annotate.quinone_dfba_essentiality import _add_runtime_q9, _optimize_minimal_pool


def _toy_model():
    model = Model("q9_reserve")
    nutrient, uracil, q9 = (
        Metabolite(name, compartment="e") for name in ("a", "uracil", "m468[C_mi]")
    )
    uptake = Reaction("R1070", lower_bound=-2000, upper_bound=0)
    uptake.add_metabolites({nutrient: -1})
    uracil_uptake = Reaction("R1354", lower_bound=-1000, upper_bound=0)
    uracil_uptake.add_metabolites({uracil: -1})
    synthase = Reaction("Q9_SYNTHASE", lower_bound=0, upper_bound=1000)
    synthase.add_metabolites({nutrient: -1, q9: 1})
    synthase.gene_reaction_rule = "qgene"
    growth = Reaction("biomass_C", lower_bound=0, upper_bound=1000)
    growth.add_metabolites({nutrient: -1, uracil: -1})
    model.add_reactions([uptake, uracil_uptake, synthase, growth])
    return model


def _toy_model_with_maintenance():
    model = _toy_model()
    maintenance = Reaction("xMAINTENANCE", lower_bound=1, upper_bound=1000)
    maintenance.add_metabolites({model.metabolites.get_by_id("a"): -1})
    model.add_reactions([maintenance])
    return model


def test_runtime_q9_reserve_is_used_only_after_q9_synthesis_knockout():
    model = _toy_model()
    source, _ = _add_runtime_q9(model, 0.1)
    source.upper_bound = 1.0
    _, solution = _optimize_minimal_pool(model, source)
    assert solution.fluxes[source.id] == 0
    with model:
        knock_out_model_genes(model, ["qgene"])
        _, solution = _optimize_minimal_pool(model, source)
        assert solution.fluxes[source.id] > 0


def test_po1f_nonlimiting_uracil_keeps_base_bound_and_has_no_finite_pool(monkeypatch):
    monkeypatch.setattr(dfba, "INITIAL_POOLS_MMOL_L", {"R1070": 1e9, "R1354": 0.1})
    finite, _ = dfba.simulate_gene(
        _toy_model(), gene_id=None, alpha=0, pool_multiplier=0, hours=0.5, dt=0.5,
        initial_biomass=1, uracil_mode="finite_batch",
    )
    nonlimiting, trace = dfba.simulate_gene(
        _toy_model(), gene_id=None, alpha=0, pool_multiplier=0, hours=0.5, dt=0.5,
        initial_biomass=1, uracil_mode="po1f_nonlimiting",
    )

    assert nonlimiting["dynamic_doublings"] > finite["dynamic_doublings"]
    assert trace[0]["growth_h-1"] == 1000
    assert math.isnan(nonlimiting["final_uracil_mmol_L"])
    assert all(row["uracil_mode"] == "po1f_nonlimiting" for row in trace)
    assert all(math.isnan(row["uracil_mmol_L"]) for row in trace)


def test_nonoptimal_step_is_recorded_once_without_reading_invalid_fluxes(monkeypatch):
    monkeypatch.setattr(dfba, "INITIAL_POOLS_MMOL_L", {"R1070": 0.5, "R1354": 1e9})
    result, trace = dfba.simulate_gene(
        _toy_model_with_maintenance(), gene_id=None, alpha=0, pool_multiplier=0,
        hours=1, dt=0.5, initial_biomass=1,
    )

    assert [row["status"] for row in trace] == ["optimal", "infeasible"]
    assert math.isnan(trace[-1]["growth_h-1"])
    assert math.isnan(trace[-1]["q9_source_flux_mmol_gDW_h"])
    assert result["termination_status"] == "infeasible"
    assert result["termination_time_h"] == 0.5


def test_gurobi_feasibility_tolerance_is_bounded_and_explicit():
    parser = dfba._parser()
    args = parser.parse_args([
        "--research-root", "research", "--solver", "gurobi", "--feasibility-tol", "1e-9",
    ])

    assert args.feasibility_tol == 1e-9
    with pytest.raises(SystemExit):
        parser.parse_args(["--research-root", "research", "--feasibility-tol", "1e-10"])
    with pytest.raises(ValueError, match="requires --solver gurobi"):
        dfba._configure_solver(_toy_model(), "glpk", 1e-9)


def test_r39_r19_runtime_topology_is_balanced_and_leaves_canonical_model_unchanged():
    model = read_sbml_model(Path(__file__).parents[1] / "model.xml")
    reaction_ids = ("R39", "R969", "R808", "R19")
    before = {
        reaction_id: (
            dfba._reaction_stoichiometry(model.reactions.get_by_id(reaction_id)),
            model.reactions.get_by_id(reaction_id).bounds,
            model.reactions.get_by_id(reaction_id).gene_reaction_rule,
        )
        for reaction_id in reaction_ids
    }
    counts = (len(model.reactions), len(model.metabolites), len(model.genes))

    with model:
        audit = dfba._apply_r39_r19_runtime_topology(model)
        assert audit["mapping_sha256"] == dfba.R39_R19_RUNTIME_MAPPING_SHA256
        assert audit["mass_balance"] == {
            "R39": {}, "DFBA_R19_HYDROXYLATION": {}, "DFBA_R19_FORMAL_OXIDATION": {},
        }
        assert model.reactions.R969.bounds == (0.0, 0.0)
        assert model.reactions.R808.bounds == (0.0, 0.0)
        assert model.reactions.R19.bounds == (0.0, 0.0)
        assert model.reactions.get_by_id("DFBA_R19_HYDROXYLATION").gene_reaction_rule == ""
        assert model.reactions.get_by_id("DFBA_R19_FORMAL_OXIDATION").gene_reaction_rule == ""

    after = {
        reaction_id: (
            dfba._reaction_stoichiometry(model.reactions.get_by_id(reaction_id)),
            model.reactions.get_by_id(reaction_id).bounds,
            model.reactions.get_by_id(reaction_id).gene_reaction_rule,
        )
        for reaction_id in reaction_ids
    }
    assert counts == (len(model.reactions), len(model.metabolites), len(model.genes))
    assert before == after


def test_dfba_calls_summary_recomputes_grid_controls_and_monotonicity():
    rows = []
    for alpha in (1e-6, 1e-4, 1e-3):
        for pool in (0.0, 0.5, 1.0, 2.0):
            doublings = math.log2(1 + pool)
            for control in FIVE_COQ_CONTROLS:
                ratio = doublings / 10
                rows.append({
                    "gene_id": control["gene_id"], "alpha_mmol_gDW": alpha,
                    "pool_multiplier": pool, "dynamic_doublings": doublings,
                    "dynamic_growth_ratio": ratio,
                    "q9_source_total_mmol_L": alpha * 0.01 * pool,
                    "experimental_essential": False,
                    **{f"essential_at_{cutoff * 100:g}pct": ratio < cutoff for cutoff in (0.01, 0.05, 0.1, 0.15)},
                })

    tables, summary = summarize_calls(pd.DataFrame(rows), initial_biomass=0.01)

    assert len(tables["grid_summary"]) == 12
    assert len(tables["five_gene_pool_summary"]) == 60
    assert summary["pool_monotonicity_pass"]
    assert summary["q9_source_calls_bound_pass"]
    assert summary["five_control_max_abs_theory_error_doublings"] == 0
