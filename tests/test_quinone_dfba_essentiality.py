import math

from cobra import Metabolite, Model, Reaction
from cobra.manipulation.delete import knock_out_model_genes

from scripts.gem_annotate import quinone_dfba_essentiality as dfba
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
