from cobra import Metabolite, Model, Reaction
from cobra.manipulation.delete import knock_out_model_genes

from scripts.gem_annotate.quinone_dfba_essentiality import _add_runtime_q9, _optimize_minimal_pool


def _toy_model():
    model = Model("q9_reserve")
    nutrient, q9 = (Metabolite(name) for name in ("a", "m468[C_mi]"))
    uptake = Reaction("UPTAKE", lower_bound=0, upper_bound=10)
    uptake.add_metabolites({nutrient: 1})
    synthase = Reaction("Q9_SYNTHASE", lower_bound=0, upper_bound=1000)
    synthase.add_metabolites({nutrient: -1, q9: 1})
    synthase.gene_reaction_rule = "qgene"
    growth = Reaction("biomass_C", lower_bound=0, upper_bound=1000)
    growth.add_metabolites({nutrient: -1})
    model.add_reactions([uptake, synthase, growth])
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
