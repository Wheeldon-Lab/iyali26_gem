from pathlib import Path

import pytest
from cobra.io import read_sbml_model, write_sbml_model

from scripts.gem_annotate.patches import split_trna_charging_from_biomass
from scripts.gem_annotate.validate_essential_genes import (
    DEFAULT_MEDIA,
    apply_media,
    load_media,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
UNSPLIT_MODEL_PATH = REPO_ROOT / "data" / "iyali26.xml"
CANONICAL_MODEL_PATH = REPO_ROOT / "model.xml"


@pytest.mark.integration
def test_group_b_splits_all_twenty_trna_requirements_and_carries_flux(
    tmp_path: Path,
) -> None:
    model = read_sbml_model(str(UNSPLIT_MODEL_PATH))
    before_counts = (len(model.reactions), len(model.metabolites))

    audit = split_trna_charging_from_biomass(model)

    assert len(audit) == 20
    assert (len(model.reactions), len(model.metabolites)) == (
        before_counts[0] + 20,
        before_counts[1] + 20,
    )
    biomass = model.reactions.get_by_id("biomass_C")
    assert biomass.notes["experimental_trna_biomass_mode"] == "split_v1"
    for row in audit:
        amino_acid = model.metabolites.get_by_id(str(row["amino_acid_id"]))
        charged = model.metabolites.get_by_id(str(row["charged_trna_id"]))
        uncharged = model.metabolites.get_by_id(str(row["uncharged_trna_id"]))
        residue = model.metabolites.get_by_id(str(row["protein_residue_id"]))
        reaction = model.reactions.get_by_id(str(row["split_reaction_id"]))
        amount = float(row["coefficient"])

        assert amino_acid not in biomass.metabolites
        assert biomass.metabolites[residue] == pytest.approx(-amount)
        assert reaction.bounds == (0.0, 1000.0)
        assert dict(reaction.metabolites) == {
            charged: -1.0,
            uncharged: 1.0,
            residue: 1.0,
        }

    apply_media(model, load_media(DEFAULT_MEDIA))
    solution = model.optimize()
    assert solution.status == "optimal"
    assert solution.objective_value == pytest.approx(1.3284010586120676, rel=1e-9)
    for row in audit:
        expected_flux = float(row["coefficient"]) * solution.objective_value
        assert solution.fluxes[str(row["split_reaction_id"])] == pytest.approx(
            expected_flux,
            rel=1e-9,
        )
        assert solution.fluxes[str(row["charging_reaction_id"])] == pytest.approx(
            expected_flux,
            rel=1e-9,
        )

    # The exact split state survives SBML serialization and remains idempotent.
    roundtrip_path = tmp_path / "group_b_split.xml"
    write_sbml_model(model, str(roundtrip_path))
    roundtrip = read_sbml_model(str(roundtrip_path))
    roundtrip_counts = (len(roundtrip.reactions), len(roundtrip.metabolites))
    assert split_trna_charging_from_biomass(roundtrip) == []
    assert (len(roundtrip.reactions), len(roundtrip.metabolites)) == roundtrip_counts


def test_group_b_rejects_a_partial_split_state() -> None:
    model = read_sbml_model(str(UNSPLIT_MODEL_PATH))
    audit = split_trna_charging_from_biomass(model)
    model.remove_reactions([str(audit[0]["split_reaction_id"])])

    with pytest.raises(ValueError, match="Partial B-group tRNA biomass state"):
        split_trna_charging_from_biomass(model)


@pytest.mark.integration
def test_canonical_model_contains_the_promoted_b_group_translation_layer() -> None:
    model = read_sbml_model(str(CANONICAL_MODEL_PATH))
    biomass = model.reactions.get_by_id("biomass_C")
    split_reactions = [
        reaction
        for reaction in model.reactions
        if reaction.id.startswith("TRNA_BIOMASS_")
    ]

    assert biomass.notes["canonical_trna_biomass_mode"] == "split_v1"
    assert len(split_reactions) == 20
    assert split_trna_charging_from_biomass(model) == []

    apply_media(model, load_media(DEFAULT_MEDIA))
    solution = model.optimize()
    assert solution.status == "optimal"
    assert solution.objective_value == pytest.approx(1.4492618988553534, rel=1e-9)
    for reaction in split_reactions:
        amount = float(reaction.notes["biomass_coefficient"])
        assert solution.fluxes[reaction.id] == pytest.approx(
            amount * solution.objective_value,
            rel=1e-9,
        )
