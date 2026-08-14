import csv
from pathlib import Path

import pytest
from cobra.io import read_sbml_model

from scripts.gem_annotate.patch_runner import _patches
from scripts.gem_annotate.patches import (
    add_direct_enzyme_like_gprs,
    add_r612_ura3_gpr,
    correct_external_ndh2_gpr_and_remove_duplicate,
    remove_spurious_quinone_branches,
    replace_coq6_route_with_coq9,
)
from scripts.gem_annotate.sbml import write_deterministic_sbml_model


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_MODEL = REPO_ROOT / "data" / "iyali26.xml"


def _source_model():
    return read_sbml_model(str(SOURCE_MODEL))


def test_r612_gets_supported_ura3_gpr_and_is_idempotent() -> None:
    model = _source_model()

    assert add_r612_ura3_gpr(model) == 1
    reaction = model.reactions.get_by_id("R612")
    gene = model.genes.get_by_id("YALI1E31685g")
    assert reaction.gene_reaction_rule == "YALI1E31685g"
    assert gene.name == "URA3"
    assert gene.annotation["uniprot"] == "A0A1H6PUU4"
    assert gene.annotation["ec-code"] == "4.1.1.23"
    assert "ura3 genetic evidence" in reaction.notes["curated_gpr_correction"]
    assert "component_review" in reaction.notes["chemistry_status"]
    assert add_r612_ura3_gpr(model) == 0


def test_external_ndh2_correction_removes_only_verified_duplicate() -> None:
    model = _source_model()
    r1889_rule = model.reactions.get_by_id("R1889").gene_reaction_rule

    assert correct_external_ndh2_gpr_and_remove_duplicate(model) == 2
    r570 = model.reactions.get_by_id("R570")
    gene = model.genes.get_by_id("YALI1F32476g")
    assert r570.gene_reaction_rule == "YALI1F32476g"
    assert gene.name == "NDH2"
    assert gene.annotation["uniprot"] == "F2Z699"
    assert "R1889 has no GPR" in r570.notes["complex_i_scope"]
    assert "R2063" not in model.reactions
    assert model.reactions.get_by_id("R1889").gene_reaction_rule == r1889_rule == ""
    assert correct_external_ndh2_gpr_and_remove_duplicate(model) == 0


def test_external_ndh2_correction_refuses_an_unrecognised_source_rule() -> None:
    model = _source_model()
    model.reactions.get_by_id("R570").gene_reaction_rule = "unexpected_gene"

    with pytest.raises(ValueError, match="unexpected GPR"):
        correct_external_ndh2_gpr_and_remove_duplicate(model)


def test_direct_enzyme_like_gprs_cover_lip2_and_all_b6_kinase_vitamers() -> None:
    model = _source_model()

    assert add_direct_enzyme_like_gprs(model) == 3
    assert model.reactions.get_by_id("R2274").gene_reaction_rule == "YALI1A21372g"
    assert model.reactions.get_by_id("R1302").gene_reaction_rule == "YALI1A08512g"
    assert model.reactions.get_by_id("R1303").gene_reaction_rule == "YALI1A08512g"
    assert model.reactions.get_by_id("R1306").gene_reaction_rule == "YALI1A08512g"
    assert model.genes.get_by_id("YALI1A21372g").name == "LIP2"
    assert model.genes.get_by_id("YALI1A08512g").annotation["ec-code"] == "2.7.1.35"
    assert (
        model.reactions.get_by_id("R2274").notes["gpr_evidence_status"]
        == "experimentally_verified"
    )
    assert (
        model.reactions.get_by_id("R1302").notes["gpr_evidence_status"]
        == "curated_annotation"
    )
    assert add_direct_enzyme_like_gprs(model) == 0


def test_direct_enzyme_like_gprs_refuse_an_unrecognised_source_rule() -> None:
    model = _source_model()
    model.reactions.get_by_id("R1303").gene_reaction_rule = "unexpected_gene"

    with pytest.raises(ValueError, match="unexpected GPR"):
        add_direct_enzyme_like_gprs(model)


def test_spurious_quinone_cleanup_removes_only_reviewed_dead_branches() -> None:
    model = _source_model()
    before_growth = model.slim_optimize()

    assert remove_spurious_quinone_branches(model) == 6
    assert {
        "R189",
        "R2242",
        "R2247",
        "R2248",
        "R2249",
        "R2250",
    }.isdisjoint(model.reactions)
    assert {"R2243", "R2244", "R2245", "R2246"} <= {
        reaction.id for reaction in model.reactions
    }
    assert {
        "m367[C_cy]",
        "m368[C_cy]",
        "m1923[C_nu]",
        "m1924[C_nu]",
        "m1927[C_cy]",
        "m1928[C_cy]",
        "m1929[C_cy]",
        "m1930[C_cy]",
    }.isdisjoint(model.metabolites)
    unresolved_genes = {
        "YALI1B21088g",
        "YALI1D17983g",
        "YALI1E01159g",
        "YALI1E11415g",
        "YALI1E16694g",
        "YALI1E33302g",
    }
    assert unresolved_genes <= {gene.id for gene in model.genes}
    assert all(not model.genes.get_by_id(gene_id).reactions for gene_id in unresolved_genes)
    assert model.slim_optimize() == pytest.approx(before_growth, abs=1e-9)

    retained = model.reactions.get_by_id("R407")
    assert retained.gene_reaction_rule == "YALI1F08349g"
    assert {met.compartment for met in retained.metabolites} == {"C_mi"}
    assert "native Yarrowia CoQ9" in retained.notes["remaining_chain_length_gate"]

    coq2 = model.genes.get_by_id("YALI1F08349g")
    coq3 = model.genes.get_by_id("YALI1B20835g")
    assert coq2.name == "COQ2"
    assert coq2.annotation["uniprot"] == ["A0A1H6PM88", "Q6C2S2"]
    assert coq2.annotation["ec-code"] == "2.5.1.39"
    assert coq3.name == "COQ3"
    assert coq3.annotation["uniprot"] == ["A0A1D8N802", "Q6CEG2"]
    assert coq3.annotation["refseq"] == "XP_500950.3"
    assert coq3.annotation["ec-code"] == ["2.1.1.64", "2.1.1.114"]

    assert remove_spurious_quinone_branches(model) == 0


def test_spurious_quinone_cleanup_fails_closed_on_partial_or_changed_branch() -> None:
    partial = _source_model()
    partial.remove_reactions([partial.reactions.get_by_id("R189")])
    with pytest.raises(ValueError, match="only partially present"):
        remove_spurious_quinone_branches(partial)

    changed = _source_model()
    changed.reactions.get_by_id("R2242").gene_reaction_rule = "unexpected_gene"
    with pytest.raises(ValueError, match="unexpected GPR"):
        remove_spurious_quinone_branches(changed)

    changed_stoichiometry = _source_model()
    r2248 = changed_stoichiometry.reactions.get_by_id("R2248")
    oxygen = changed_stoichiometry.metabolites.get_by_id("m109[C_cy]")
    r2248.add_metabolites({oxygen: -0.5})
    with pytest.raises(ValueError, match="stoichiometry/bounds/compartment"):
        remove_spurious_quinone_branches(changed_stoichiometry)

    changed_bounds = _source_model()
    changed_bounds.reactions.get_by_id("R2250").upper_bound = 999.0
    with pytest.raises(ValueError, match="stoichiometry/bounds/compartment"):
        remove_spurious_quinone_branches(changed_bounds)


def test_spurious_quinone_cleanup_accepts_canonical_r189_balance_variant() -> None:
    model = _source_model()
    r189 = model.reactions.get_by_id("R189")
    proton = model.metabolites.get_by_id("m10[C_cy]")
    r189.add_metabolites({proton: -2.0})

    assert remove_spurious_quinone_branches(model) == 6


def test_spurious_quinone_cleanup_survives_sbml_roundtrip(tmp_path: Path) -> None:
    model = _source_model()
    remove_spurious_quinone_branches(model)
    output = tmp_path / "quinone-cleanup.xml"

    write_deterministic_sbml_model(model, output)
    reloaded = read_sbml_model(str(output))

    assert {
        "R189",
        "R2242",
        "R2247",
        "R2248",
        "R2249",
        "R2250",
    }.isdisjoint(reloaded.reactions)
    assert reloaded.reactions.get_by_id("R407").gene_reaction_rule == "YALI1F08349g"
    assert set(reloaded.genes.get_by_id("YALI1F08349g").annotation["uniprot"]) == {
        "A0A1H6PM88",
        "Q6C2S2",
    }
    assert not reloaded.genes.get_by_id("YALI1E01159g").reactions


_COQ_ROUTE_IDS = (
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
_LEGACY_COQ6_FORMULAS = {
    "m640[C_mi]": "C30H52O7P2",
    "m641[C_mi]": "C37H54O3",
    "m108[C_cy]": "C37H54O3",
    "m110[C_cy]": "C37H53O4",
    "m939[C_mi]": "C37H53O4",
    "m111[C_mi]": "C38H55O4",
    "m63[C_mi]": "C37H56O2",
    "m59[C_mi]": "C37H54O3",
    "m61[C_mi]": "C38H56O3",
    "m611[C_mi]": "C38H56O4",
    "m468[C_mi]": "C39H58O4",
    "m471[C_mi]": "C39H60O4",
}
_TARGET_COQ9_FORMULAS = {
    "m640[C_mi]": ("C45H76O7P2", 0),
    "m641[C_mi]": ("C52H78O3", 0),
    "m108[C_cy]": ("C52H78O3", 0),
    "m110[C_cy]": ("C52H77O4", -1),
    "m939[C_mi]": ("C52H77O4", -1),
    "m111[C_mi]": ("C53H79O4", -1),
    "m63[C_mi]": ("C52H80O2", 0),
    "m59[C_mi]": ("C52H78O3", 0),
    "m61[C_mi]": ("C53H80O3", 0),
    "m611[C_mi]": ("C53H80O4", 0),
    "m468[C_mi]": ("C54H82O4", 0),
    "m471[C_mi]": ("C54H84O4", 0),
}


def _set_stoichiometry(model, reaction_id: str, target: dict[str, float]) -> None:
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.add_metabolites(
        {metabolite: -coefficient for metabolite, coefficient in reaction.metabolites.items()}
    )
    reaction.add_metabolites(
        {
            model.metabolites.get_by_id(metabolite_id): coefficient
            for metabolite_id, coefficient in target.items()
        }
    )


def _restore_legacy_coq6_fixture(model) -> None:
    for metabolite_id, formula in _LEGACY_COQ6_FORMULAS.items():
        model.metabolites.get_by_id(metabolite_id).formula = formula
    _set_stoichiometry(
        model,
        "R763",
        {
            "m984[C_mi]": -1.0,
            "m985[C_mi]": -1.0,
            "m204[C_mi]": 1.0,
            "m640[C_mi]": 1.0,
        },
    )
    _set_stoichiometry(
        model,
        "R385",
        {
            "m28[C_mi]": -1.0,
            "m60[C_mi]": -1.0,
            "m611[C_mi]": -1.0,
            "m471[C_mi]": 1.0,
            "m62[C_mi]": 1.0,
        },
    )


def test_formal_coq9_patch_is_balanced_idempotent_and_structure_neutral() -> None:
    model = read_sbml_model(str(REPO_ROOT / "model.xml"))
    _restore_legacy_coq6_fixture(model)
    before_growth = model.slim_optimize()
    before_counts = (len(model.reactions), len(model.metabolites), len(model.genes))
    before_gprs = {reaction.id: reaction.gene_reaction_rule for reaction in model.reactions}
    before_bounds = {reaction.id: tuple(reaction.bounds) for reaction in model.reactions}
    before_biomass = {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in model.reactions.biomass_C.metabolites.items()
    }
    before_demands = {reaction.id for reaction in model.demands}
    before_sinks = {reaction.id for reaction in model.sinks}

    assert replace_coq6_route_with_coq9(model) > 0
    assert replace_coq6_route_with_coq9(model) == 0

    assert {
        metabolite_id: (
            model.metabolites.get_by_id(metabolite_id).formula,
            model.metabolites.get_by_id(metabolite_id).charge,
        )
        for metabolite_id in _TARGET_COQ9_FORMULAS
    } == _TARGET_COQ9_FORMULAS
    assert {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in model.reactions.R763.metabolites.items()
    } == {
        "m984[C_mi]": -4.0,
        "m985[C_mi]": -1.0,
        "m204[C_mi]": 4.0,
        "m640[C_mi]": 1.0,
    }
    assert {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in model.reactions.R385.metabolites.items()
    } == {
        "m60[C_mi]": -1.0,
        "m611[C_mi]": -1.0,
        "m468[C_mi]": 1.0,
        "m62[C_mi]": 1.0,
    }
    assert model.reactions.R385.annotation["ec-code"] == "2.1.1.64"
    assert all(
        model.reactions.get_by_id(reaction_id).check_mass_balance() == {}
        for reaction_id in _COQ_ROUTE_IDS
    )
    assert (len(model.reactions), len(model.metabolites), len(model.genes)) == before_counts
    assert {reaction.id: reaction.gene_reaction_rule for reaction in model.reactions} == before_gprs
    assert {reaction.id: tuple(reaction.bounds) for reaction in model.reactions} == before_bounds
    assert {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in model.reactions.biomass_C.metabolites.items()
    } == before_biomass
    assert {reaction.id for reaction in model.demands} == before_demands
    assert {reaction.id for reaction in model.sinks} == before_sinks
    assert model.slim_optimize() == pytest.approx(before_growth, abs=1e-9)


def test_formal_coq9_patch_fails_before_mutation_on_changed_r763() -> None:
    model = read_sbml_model(str(REPO_ROOT / "model.xml"))
    _restore_legacy_coq6_fixture(model)
    ipp = model.metabolites.get_by_id("m984[C_mi]")
    model.reactions.R763.add_metabolites({ipp: -1.0})
    formulas_before = {
        metabolite_id: model.metabolites.get_by_id(metabolite_id).formula
        for metabolite_id in _LEGACY_COQ6_FORMULAS
    }
    stoichiometry_before = {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in model.reactions.R763.metabolites.items()
    }

    with pytest.raises(ValueError, match="R763 no longer matches"):
        replace_coq6_route_with_coq9(model)

    assert {
        metabolite_id: model.metabolites.get_by_id(metabolite_id).formula
        for metabolite_id in _LEGACY_COQ6_FORMULAS
    } == formulas_before
    assert {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in model.reactions.R763.metabolites.items()
    } == stoichiometry_before


def test_individual_reviewed_corrections_are_exposed_by_patch_runner() -> None:
    available = _patches(allow_network=False)
    assert available["r612-ura3-gpr"] is add_r612_ura3_gpr
    assert available["external-ndh2-correction"] is correct_external_ndh2_gpr_and_remove_duplicate
    assert available["direct-enzyme-like-gprs"] is add_direct_enzyme_like_gprs
    assert available["quinone-branch-cleanup"] is remove_spurious_quinone_branches


def test_complex_i_evidence_table_is_explicitly_deferred() -> None:
    evidence_path = REPO_ROOT / "docs" / "curation" / "complex_i_gpr_evidence.csv"
    with evidence_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 42
    assert len({row["gene_id"] for row in rows}) == 42
    assert {row["gpr_decision"] for row in rows} == {"deferred"}
    assert {"YALI1D18037g", "YALI1D32594g", "YALI1M00472r"} <= {
        row["gene_id"] for row in rows if row["identity_status"] == "identity_conflict"
    }
