import math

import pytest
from cobra import Metabolite, Model, Reaction
from cobra.flux_analysis import pfba
from cobra.manipulation.delete import knock_out_model_genes

from scripts.gem_annotate.coq9_dilution import (
    ALPHA_HIGH_MMOL_GDW,
    ALPHA_LOW_MMOL_GDW,
    COQ9_DILUTION_CONSTRAINT_ID,
    COQ9_DILUTION_ID,
    COQ9_POOL_SOURCE_ID,
    add_runtime_coq9_pool_source,
    apply_runtime_coq9_dilution,
    sample_coq9_alpha,
)


def _toy_model() -> Model:
    model = Model("coq9_dilution_toy")
    nutrient = Metabolite("nutrient", compartment="c")
    q9 = Metabolite("m468[C_mi]", compartment="C_mi")
    synthase = Reaction("Q9_SYNTHASE", lower_bound=0.0, upper_bound=1000.0)
    synthase.add_metabolites({nutrient: -1.0, q9: 1.0})
    synthase.gene_reaction_rule = "q9_gene"
    biomass = Reaction("biomass_C", lower_bound=0.0, upper_bound=1000.0)
    biomass.add_metabolites({nutrient: -1.0})
    supply = Reaction("NUTRIENT_SOURCE", lower_bound=0.0, upper_bound=1000.0)
    supply.add_metabolites({nutrient: 1.0})
    model.add_reactions([supply, synthase, biomass])
    model.objective = biomass
    return model


def test_seeded_alpha_is_reproducible_and_bounded():
    first = sample_coq9_alpha(20260826, 4)
    assert first == sample_coq9_alpha(20260826, 4)
    assert ALPHA_LOW_MMOL_GDW <= first <= ALPHA_HIGH_MMOL_GDW
    assert math.isfinite(first)
    with pytest.raises(ValueError, match="replicate_id"):
        sample_coq9_alpha(20260826, -1)


def test_runtime_dilution_is_coupled_and_restores_the_model():
    model = _toy_model()
    before = (tuple(model.reactions.list_attr("id")), tuple(model.metabolites.list_attr("id")))
    alpha = 0.1

    with model:
        dilution = apply_runtime_coq9_dilution(model, alpha)
        solution = pfba(model)
        assert solution.fluxes[dilution.id] == pytest.approx(
            alpha * solution.fluxes["biomass_C"]
        )
        assert solution.fluxes["Q9_SYNTHASE"] == pytest.approx(solution.fluxes[dilution.id])
        assert COQ9_POOL_SOURCE_ID not in model.reactions
        with pytest.raises(ValueError, match="already exists"):
            apply_runtime_coq9_dilution(model, alpha)

    assert before == (
        tuple(model.reactions.list_attr("id")), tuple(model.metabolites.list_attr("id"))
    )
    assert COQ9_DILUTION_ID not in model.reactions
    assert COQ9_DILUTION_CONSTRAINT_ID not in model.constraints


def test_pool_source_is_explicit_and_only_rescues_when_enabled():
    model = _toy_model()
    with model:
        dilution = apply_runtime_coq9_dilution(model, 0.1)
        with model:
            knock_out_model_genes(model, ["q9_gene"])
            assert pfba(model).fluxes["biomass_C"] == pytest.approx(0.0)

        source = add_runtime_coq9_pool_source(model)
        source.upper_bound = 1.0
        with model:
            knock_out_model_genes(model, ["q9_gene"])
            solution = pfba(model)
            assert solution.fluxes["biomass_C"] > 0.0
            assert solution.fluxes[source.id] > 0.0
            assert solution.fluxes[dilution.id] == pytest.approx(
                0.1 * solution.fluxes["biomass_C"]
            )
