from pathlib import Path

import pytest
from cobra.io import read_sbml_model

from scripts.gem_annotate.patches import split_trna_charging_from_biomass
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.trna_biomass_pipeline import (
    build_group_b_overlay,
    restore_free_amino_acid_biomass,
    validate_group_b_overlay,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
UNSPLIT_MODEL_PATH = REPO_ROOT / "data" / "iyali26.xml"
CANONICAL_MODEL_PATH = REPO_ROOT / "model.xml"


def _write_group_b(tmp_path: Path) -> Path:
    model = read_sbml_model(str(UNSPLIT_MODEL_PATH))
    split_trna_charging_from_biomass(model)
    path = tmp_path / "group_b.xml"
    write_deterministic_sbml_model(model, path)
    return path


def test_group_b_pipeline_overlay_audit_accepts_complete_split(tmp_path: Path) -> None:
    group_b_path = _write_group_b(tmp_path)

    audit = validate_group_b_overlay(UNSPLIT_MODEL_PATH, group_b_path)

    assert audit["carrier_conserving_split_reactions"] == 20
    assert audit["changed_preexisting_reactions"] == ["biomass_C"]
    assert audit["group_b_reactions"] == audit["canonical_reactions"] + 20
    assert audit["group_b_metabolites"] == audit["canonical_metabolites"] + 20
    assert audit["objective_reaction"] == "biomass_C"


def test_group_b_pipeline_builds_an_isolated_overlay_from_canonical_model(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "group_b.xml"

    split_audit = build_group_b_overlay(UNSPLIT_MODEL_PATH, output_path)
    audit = validate_group_b_overlay(UNSPLIT_MODEL_PATH, output_path)

    assert len(split_audit) == 20
    assert audit["changed_preexisting_reactions"] == ["biomass_C"]
    assert audit["group_b_reactions"] == audit["canonical_reactions"] + 20


def test_group_b_pipeline_copies_the_promoted_canonical_model(tmp_path: Path) -> None:
    output_path = tmp_path / "canonical_b_group.xml"

    split_audit = build_group_b_overlay(CANONICAL_MODEL_PATH, output_path)
    audit = validate_group_b_overlay(CANONICAL_MODEL_PATH, output_path)

    assert split_audit == []
    assert audit["canonical_b_group"] is True
    assert audit["changed_preexisting_reactions"] == []
    assert audit["group_b_reactions"] == audit["canonical_reactions"]


def test_free_amino_acid_counterfactual_exactly_reverses_promoted_b_group() -> None:
    model = read_sbml_model(str(CANONICAL_MODEL_PATH))
    reaction_count = len(model.reactions)
    metabolite_count = len(model.metabolites)

    inverse_audit = restore_free_amino_acid_biomass(model)

    assert len(inverse_audit) == 20
    assert not [
        reaction
        for reaction in model.reactions
        if reaction.id.startswith("TRNA_BIOMASS_")
    ]
    assert model.reactions.biomass_C.notes.get("canonical_trna_biomass_mode") is None
    assert model.reactions.biomass_C.notes.get("experimental_trna_biomass_mode") is None
    assert len(model.reactions) == reaction_count - 20
    assert len(model.metabolites) == metabolite_count - 20


def test_group_b_pipeline_overlay_audit_rejects_unrelated_change(tmp_path: Path) -> None:
    group_b_path = _write_group_b(tmp_path)
    model = read_sbml_model(str(group_b_path))
    model.reactions.get_by_id("R1").upper_bound = 999.0
    write_deterministic_sbml_model(model, group_b_path)

    with pytest.raises(
        ValueError,
        match="Only biomass_C may change among pre-existing reactions",
    ):
        validate_group_b_overlay(UNSPLIT_MODEL_PATH, group_b_path)
