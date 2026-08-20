from pathlib import Path

import pandas as pd
import pytest
from cobra import Metabolite, Model, Reaction

from scripts.gem_annotate.essentiality_evidence import target_fingerprint
from scripts.gem_annotate.isozyme_capacity import (
    load_isozyme_capacity_scan,
    run_isozyme_capacity_scan,
    validate_isozyme_capacity_scan,
)


def _toy_isozyme_model() -> Model:
    model = Model("toy_isozyme")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")

    uptake = Reaction("UPTAKE")
    uptake.bounds = (0.0, 10.0)
    uptake.add_metabolites({a: 1.0})

    isozyme = Reaction("R_ISO")
    isozyme.bounds = (0.0, 1000.0)
    isozyme.add_metabolites({a: -1.0, b: 1.0})
    isozyme.gene_reaction_rule = "g_backup or g_main"

    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1000.0)
    biomass.add_metabolites({b: -1.0})

    model.add_reactions([uptake, isozyme, biomass])
    model.objective = biomass
    return model


def _fingerprint(model: Model, reaction_id: str) -> str:
    reaction = model.reactions.get_by_id(reaction_id)
    return target_fingerprint(
        [
            {
                "reaction_id": reaction.id,
                "stoichiometry": {
                    metabolite.id: float(coefficient)
                    for metabolite, coefficient in reaction.metabolites.items()
                },
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
                "gpr": reaction.gene_reaction_rule,
            }
        ]
    )


def _scan_frame(model: Model, fraction: float = 0.05) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "scenario_id": "toy_backup_cap",
                "case_id": "EGC-toy",
                "target_gene": "g_main",
                "reaction_id": "R_ISO",
                "capacity_fraction_of_wt_flux": fraction,
                "model_sha256": "model-sha",
                "media_sha256": "media-sha",
                "target_fingerprint": _fingerprint(model, "R_ISO"),
                "exploratory_only": True,
                "evidence_level": "hypothesis_sensitivity",
                "basis": "test",
                "rationale": "test only",
            }
        ]
    )


def test_ko_specific_capacity_scan_leaves_wt_and_model_unchanged() -> None:
    model = _toy_isozyme_model()
    model.solver = "glpk"
    wt_growth = float(model.slim_optimize())
    before_bounds = model.reactions.get_by_id("R_ISO").bounds
    scan = _scan_frame(model, fraction=0.05)
    validate_isozyme_capacity_scan(model, scan, "model-sha", "media-sha")

    baseline = pd.DataFrame(
        [
            {
                "gene_id": "g_main",
                "ko_growth": wt_growth,
                "ko_growth_ratio": 1.0,
            },
            {
                "gene_id": "g_backup",
                "ko_growth": wt_growth,
                "ko_growth_ratio": 1.0,
            },
        ]
    )
    experimental = pd.DataFrame(
        [{"gene_id": "g_main", "essential": True}]
    )

    result, summary = run_isozyme_capacity_scan(
        model,
        scan,
        baseline,
        experimental,
        wt_growth,
        (0.01, 0.05, 0.10, 0.15),
        "glpk",
    )

    assert float(result.loc[0, "scenario_ko_growth_ratio"]) == pytest.approx(0.05)
    assert bool(result.loc[0, "scenario_essential_at_10pct"]) is True
    assert bool(result.loc[0, "scenario_essential_at_5pct"]) is False
    assert bool(result.loc[0, "scenario_essential_at_1pct"]) is False
    assert summary["not_a_model_patch"] is True
    assert model.reactions.get_by_id("R_ISO").bounds == before_bounds
    assert float(model.slim_optimize()) == pytest.approx(wt_growth)


def test_capacity_scan_rejects_stale_sha_and_blocked_reaction() -> None:
    model = _toy_isozyme_model()
    scan = _scan_frame(model)
    with pytest.raises(ValueError, match="model SHA mismatch"):
        validate_isozyme_capacity_scan(model, scan, "different-sha", "media-sha")

    model.reactions.get_by_id("R_ISO").upper_bound = 0.0
    scan.loc[0, "target_fingerprint"] = _fingerprint(model, "R_ISO")
    with pytest.raises(ValueError, match="blocked"):
        validate_isozyme_capacity_scan(model, scan, "model-sha", "media-sha")


def test_capacity_scan_loader_requires_exploratory_only(tmp_path: Path) -> None:
    path = tmp_path / "scan.csv"
    path.write_text(
        "scenario_id,case_id,target_gene,reaction_id,"
        "capacity_fraction_of_wt_flux,model_sha256,media_sha256,"
        "target_fingerprint,exploratory_only,evidence_level,basis,rationale\n"
        "s,EGC-x,g1,R1,0.1,model,media,sha256:x,false,test,test,test\n"
    )
    with pytest.raises(ValueError, match="exploratory_only=true"):
        load_isozyme_capacity_scan(path)
