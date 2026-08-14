from pathlib import Path

import pytest
from cobra import Metabolite, Model, Reaction

from scripts.gem_annotate.iron_uptake_pair1 import (
    FETC_GENE,
    FTRA_GENE,
    INHERITED_GPR,
    MODE_GPRS,
    apply_pair1_overlay,
    run_pair1_overlay,
)
from scripts.gem_annotate.sbml import write_deterministic_sbml_model


def _model() -> Model:
    model = Model("pair1-test")
    extracellular = Metabolite("iron_ex", compartment="ex")
    cytosolic = Metabolite("iron_cy", compartment="cy")
    reaction = Reaction("R2193")
    reaction.bounds = (0.0, 1000.0)
    reaction.add_metabolites({extracellular: -1.0, cytosolic: 1.0})
    reaction.gene_reaction_rule = INHERITED_GPR
    model.add_reactions([reaction])
    return model


def _protected_signature(model: Model):
    reaction = model.reactions.get_by_id("R2193")
    return (
        {
            metabolite.id: coefficient
            for metabolite, coefficient in reaction.metabolites.items()
        },
        reaction.bounds,
    )


def test_ftra_only_uses_directly_supported_gene_and_preserves_reaction() -> None:
    model = _model()
    protected = _protected_signature(model)

    audit = apply_pair1_overlay(model, mode="ftra_only")

    reaction = model.reactions.get_by_id("R2193")
    assert reaction.gene_reaction_rule == FTRA_GENE
    assert _protected_signature(model) == protected
    assert audit["evidence_status"] == "same_species_experimental_support"
    assert audit["inherited_gene_retained_as_orphan"] is True
    annotation = model.genes.get_by_id(FTRA_GENE).annotation
    assert annotation["ncbigene"] == "2910500"
    assert annotation["kegg.genes"] == "yli:2910500"
    assert annotation["uniprot"] == "Q6CA15"
    assert annotation["refseq"] == "XP_502497.1"
    assert annotation["sbo"] == "SBO:0000243"


def test_pair1_and_is_explicitly_labelled_as_inferred_sensitivity() -> None:
    model = _model()

    audit = apply_pair1_overlay(model, mode="pair1_and")

    reaction = model.reactions.get_by_id("R2193")
    assert reaction.gene_reaction_rule == MODE_GPRS["pair1_and"]
    assert {gene.id for gene in reaction.genes} == {FTRA_GENE, FETC_GENE}
    assert audit["evidence_status"] == "sensitivity_hypothesis"
    assert "AND relationship is inferred" in audit["limitation"]
    fetc = model.genes.get_by_id(FETC_GENE).annotation
    assert fetc["ncbigene"] == "2910503"
    assert fetc["uniprot"] == "Q6CA12"
    assert fetc["refseq"] == "XP_502500.2"


@pytest.mark.parametrize("mode", ["ftra_only", "pair1_and"])
def test_overlay_is_idempotent(mode: str) -> None:
    model = _model()

    first = apply_pair1_overlay(model, mode=mode)
    second = apply_pair1_overlay(model, mode=mode)

    assert first["changed"] is True
    assert second["changed"] is False
    assert model.reactions.get_by_id("R2193").gene_reaction_rule == MODE_GPRS[mode]


def test_overlay_refuses_unrelated_existing_gpr() -> None:
    model = _model()
    model.reactions.get_by_id("R2193").gene_reaction_rule = "UNRELATED_GENE"

    with pytest.raises(ValueError, match="unexpected GPR"):
        apply_pair1_overlay(model)


def test_runner_refuses_canonical_and_existing_outputs(tmp_path: Path) -> None:
    input_path = tmp_path / "input.xml"
    write_deterministic_sbml_model(_model(), input_path)

    with pytest.raises(ValueError, match="canonical model.xml"):
        run_pair1_overlay(
            input_model=input_path,
            output_model=Path(__file__).resolve().parents[1] / "model.xml",
        )

    existing = tmp_path / "existing.xml"
    existing.write_text("preserve", encoding="utf-8")
    with pytest.raises(FileExistsError, match="refuses to overwrite"):
        run_pair1_overlay(input_model=input_path, output_model=existing)


def test_runner_writes_roundtrippable_model_and_audit(tmp_path: Path) -> None:
    input_path = tmp_path / "input.xml"
    output_path = tmp_path / "output.xml"
    write_deterministic_sbml_model(_model(), input_path)

    audit = run_pair1_overlay(
        input_model=input_path,
        output_model=output_path,
        mode="ftra_only",
    )

    assert output_path.is_file()
    assert output_path.with_suffix(".xml.pair1-audit.json").is_file()
    assert audit["canonical_model"]["unchanged"] is True
    assert audit["overlay"]["after"]["gpr"] == FTRA_GENE
