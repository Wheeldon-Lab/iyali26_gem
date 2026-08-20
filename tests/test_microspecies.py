from __future__ import annotations

import hashlib
from pathlib import Path

import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.microspecies import (
    DEFAULT_MICROSPECIES_TABLE,
    REFERENCE_PH,
    apply_curated_microspecies,
    audit_component_migration,
    balance_protons_and_water,
    load_curated_microspecies,
    normalize_hydroxide_reactions,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _met(mid: str, name: str, formula: str | None, charge: int, compartment="c"):
    return Metabolite(
        mid,
        name=name,
        formula=formula,
        charge=charge,
        compartment=compartment,
    )


def _model_with_common_ions() -> Model:
    model = Model("ions")
    metabolites = [
        _met("h_c", "H+_p+1", "H", 1),
        _met("h2o_c", "H2O_H2O", "H2O", 0),
        _met("hco3_c", "bicarbonate_CHO3", "CHO3", -1),
        _met("k_c", "potassium_K.H", "KH", 0),
        _met("k_e", "potassium_K.H", "KH", 0, "e"),
        _met("na_c", "sodium_Na.H", "NaH", 0),
        _met("na_e", "sodium_Na.H", "NaH", 0, "e"),
        _met("mg_c", "Magnesium", None, 0),
        _met("mg_e", "Magnesium", None, 0, "e"),
        _met("nh4_c", "ammonium_H3N", "H3N", 1),
    ]
    model.add_metabolites(metabolites)
    for rid, left, right in (
        ("TK", "k_e", "k_c"),
        ("TNA", "na_e", "na_c"),
        ("TMG", "mg_e", "mg_c"),
    ):
        reaction = Reaction(rid)
        reaction.add_metabolites(
            {
                model.metabolites.get_by_id(left): -1,
                model.metabolites.get_by_id(right): 1,
            }
        )
        model.add_reactions([reaction])
    return model


def _one_row_table(tmp_path: Path, row: str) -> Path:
    path = tmp_path / "microspecies.csv"
    path.write_text(
        "schema_version,status,family_id,selector_type,selector_value,"
        "target_formula,target_charge,reference_ph,chebi_id,rhea_id,source_url,"
        "min_matches,expected_metabolite_ids,allowed_current_pairs,rationale\n"
        + row
        + "\n",
        encoding="utf-8",
    )
    return path


def _multi_row_table(tmp_path: Path, rows: list[str]) -> Path:
    path = tmp_path / "microspecies.csv"
    path.write_text(
        "schema_version,status,family_id,selector_type,selector_value,"
        "target_formula,target_charge,reference_ph,chebi_id,rhea_id,source_url,"
        "min_matches,expected_metabolite_ids,allowed_current_pairs,rationale\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )
    return path


def _common_ion_table(tmp_path: Path) -> Path:
    return _multi_row_table(
        tmp_path,
        [
            "1,active,proton,base_name,H+,H,1,7.3,CHEBI:15378,,https://example.org,1,h_c,H|1,test",
            "1,active,water,base_name,H2O,H2O,0,7.3,CHEBI:15377,,https://example.org,1,h2o_c,H2O|0,test",
            "1,active,bicarbonate,base_name,bicarbonate,CHO3,-1,7.3,CHEBI:17544,,https://example.org,1,hco3_c,CHO3|-1,test",
            "1,active,potassium,base_name,potassium,K,1,7.3,CHEBI:29103,,https://example.org,2,k_c;k_e,KH|0;K|1,test",
            "1,active,sodium,base_name,sodium,Na,1,7.3,CHEBI:29101,,https://example.org,2,na_c;na_e,NaH|0;Na|1,test",
            "1,active,magnesium,base_name,Magnesium,Mg,2,7.3,CHEBI:18420,,https://example.org,2,mg_c;mg_e,<missing>|0;Mg|2,test",
        ],
    )


def _acid_base_table(tmp_path: Path, hydroxide_ids: str = "") -> Path:
    hydroxide_count = len([value for value in hydroxide_ids.split(";") if value])
    return _multi_row_table(
        tmp_path,
        [
            "1,active,proton,base_name,H+,H,1,7.3,CHEBI:15378,,https://example.org,1,h_c,H|1,test",
            "1,active,water,base_name,H2O,H2O,0,7.3,CHEBI:15377,,https://example.org,1,h2o_c,H2O|0,test",
            f"1,normalize,hydroxide,base_name,hydroxide,HO,-1,7.3,CHEBI:16234,,https://example.org,{hydroxide_count},{hydroxide_ids},HO|-1,test",
        ],
    )


def test_default_table_is_valid_and_defers_connected_families() -> None:
    rows = load_curated_microspecies()
    assert {row.status for row in rows} == {
        "active",
        "component_review",
        "verified_current",
        "normalize",
    }
    assert all(row.reference_ph == pytest.approx(REFERENCE_PH) for row in rows)
    assert {row.family_id for row in rows if row.status == "component_review"} >= {
        "ammonium",
        "phosphate",
        "diphosphate",
        "atp",
        "gtp",
        "mannose_1_phosphate",
    }
    assert {row.family_id for row in rows if row.status == "verified_current"} == {
        "glyceraldehyde_3_phosphate",
        "l_glutamine",
        "one_three_bisphospho_d_glycerate",
        "oxygen",
    }


def test_active_pairs_apply_atomically_and_are_idempotent(tmp_path: Path) -> None:
    model = _model_with_common_ions()
    table = _common_ion_table(tmp_path)
    model.metabolites.mg_c.annotation = {"chebi": "CHEBI:39128", "keep": "yes"}
    first = apply_curated_microspecies(model, table)

    assert first["changed_metabolites"] == 6
    assert first["balanced_reaction_regressions"] == []
    assert (model.metabolites.k_c.formula, model.metabolites.k_c.charge) == ("K", 1)
    assert (model.metabolites.na_c.formula, model.metabolites.na_c.charge) == ("Na", 1)
    assert (model.metabolites.mg_c.formula, model.metabolites.mg_c.charge) == ("Mg", 2)
    assert model.metabolites.mg_c.annotation["chebi"] == ["CHEBI:18420"]
    assert model.metabolites.mg_c.annotation["keep"] == "yes"
    # Verified chemistry is visible but cannot be applied as an isolated patch.
    assert (model.metabolites.nh4_c.formula, model.metabolites.nh4_c.charge) == (
        "H3N",
        1,
    )

    second = apply_curated_microspecies(model, table)
    assert second["changed_metabolites"] == 0
    assert second["already_canonical"] == first["already_canonical"] + 6


def test_stale_pair_is_rejected_before_any_mutation(tmp_path: Path) -> None:
    table = _one_row_table(
        tmp_path,
        "1,active,potassium,base_name,potassium,K,1,7.3,CHEBI:29103,,"
        "https://www.ebi.ac.uk/chebi/CHEBI:29103,1,k_c,KH|0;K|1,test",
    )
    model = Model("stale")
    potassium = _met("k_c", "potassium", "K2", 0)
    model.add_metabolites([potassium])

    with pytest.raises(ValueError, match="stale/unexpected pair"):
        apply_curated_microspecies(model, table)
    assert (potassium.formula, potassium.charge) == ("K2", 0)


def test_selector_scope_expansion_is_rejected(tmp_path: Path) -> None:
    table = _one_row_table(
        tmp_path,
        "1,active,potassium,base_name,potassium,K,1,7.3,CHEBI:29103,,"
        "https://www.ebi.ac.uk/chebi/CHEBI:29103,1,k_c,KH|0;K|1,test",
    )
    model = Model("scope")
    cytosolic = _met("k_c", "potassium", "KH", 0)
    unexpected = _met("k_e", "potassium", "KH", 0, "e")
    model.add_metabolites([cytosolic, unexpected])

    with pytest.raises(ValueError, match="selector target set changed"):
        apply_curated_microspecies(model, table)
    assert (cytosolic.formula, cytosolic.charge) == ("KH", 0)
    assert (unexpected.formula, unexpected.charge) == ("KH", 0)


def test_balanced_reaction_regression_rolls_back_batch(tmp_path: Path) -> None:
    table = _one_row_table(
        tmp_path,
        "1,active,potassium,base_name,potassium,K,1,7.3,CHEBI:29103,,"
        "https://www.ebi.ac.uk/chebi/CHEBI:29103,1,k_c,KH|0;K|1,test",
    )
    model = Model("rollback")
    potassium = _met("k_c", "potassium", "KH", 0)
    unchanged = _met("salt_c", "salt", "KH", 0)
    reaction = Reaction("R")
    reaction.add_metabolites({potassium: -1, unchanged: 1})
    model.add_reactions([reaction])
    assert reaction.check_mass_balance() == {}

    with pytest.raises(ValueError, match="break previously balanced"):
        apply_curated_microspecies(model, table)
    assert (potassium.formula, potassium.charge) == ("KH", 0)
    assert reaction.check_mass_balance() == {}


def _component_pair_model() -> Model:
    model = Model("component_pair")
    left = _met("a_c", "A", "H", 0)
    right = _met("b_c", "B", "H", 0)
    reaction = Reaction("R_COMPONENT")
    reaction.add_metabolites({left: -1, right: 1})
    model.add_reactions([reaction])
    return model


def _component_pair_table(tmp_path: Path) -> Path:
    return _multi_row_table(
        tmp_path,
        [
            "1,component_review,family_a,base_name,A,C,0,7.3,CHEBI:1,,"
            "https://example.org/a,1,a_c,H|0;C|0,test",
            "1,component_review,family_b,base_name,B,C,0,7.3,CHEBI:2,,"
            "https://example.org/b,1,b_c,H|0;C|0,test",
        ],
    )


def test_component_migration_audit_is_read_only_and_deterministic(
    tmp_path: Path,
) -> None:
    model = _component_pair_model()
    table = _component_pair_table(tmp_path)
    before_pairs = {
        metabolite.id: (metabolite.formula, metabolite.charge)
        for metabolite in model.metabolites
    }
    before_balance = model.reactions.R_COMPONENT.check_mass_balance()

    first = audit_component_migration(
        model, ["family_b", "family_a"], table_path=table
    )
    second = audit_component_migration(
        model, ["family_a", "family_b"], table_path=table
    )

    assert first == second
    assert first["plan_fingerprint"].startswith("sha256:")
    assert first["selected_family_ids"] == ["family_a", "family_b"]
    assert first["closure_family_ids"] == ["family_a", "family_b"]
    assert first["changed_metabolite_ids"] == ["a_c", "b_c"]
    assert first["ready_for_activation"] is True
    assert first["regressed_reaction_ids"] == []
    assert first["unresolved_reaction_ids"] == []
    assert {
        metabolite.id: (metabolite.formula, metabolite.charge)
        for metabolite in model.metabolites
    } == before_pairs
    assert model.reactions.R_COMPONENT.check_mass_balance() == before_balance


def test_incomplete_component_migration_reports_regression_and_frontier(
    tmp_path: Path,
) -> None:
    model = _component_pair_model()
    table = _component_pair_table(tmp_path)

    report = audit_component_migration(model, ["family_a"], table_path=table)

    assert report["ready_for_activation"] is False
    assert report["regressed_reaction_ids"] == ["R_COMPONENT"]
    assert report["unresolved_reaction_ids"] == ["R_COMPONENT"]
    assert report["frontier_family_ids"] == ["family_b"]
    assert report["frontier_metabolite_ids"] == ["b_c"]
    assert report["closure_frontier"] == [
        {
            "reaction_id": "R_COMPONENT",
            "unselected_component_family_ids": ["family_b"],
            "outside_selected_metabolite_ids": ["b_c"],
        }
    ]
    assert (model.metabolites.a_c.formula, model.metabolites.a_c.charge) == (
        "H",
        0,
    )


def test_component_migration_refuses_unknown_or_nondeferred_family(
    tmp_path: Path,
) -> None:
    model = _component_pair_model()
    table = _multi_row_table(
        tmp_path,
        [
            "1,active,family_a,base_name,A,C,0,7.3,CHEBI:1,,"
            "https://example.org/a,1,a_c,H|0;C|0,test",
            "1,component_review,family_b,base_name,B,C,0,7.3,CHEBI:2,,"
            "https://example.org/b,1,b_c,H|0;C|0,test",
        ],
    )

    with pytest.raises(ValueError, match="unknown microspecies family"):
        audit_component_migration(model, ["missing"], table_path=table)
    with pytest.raises(ValueError, match="expected status component_review"):
        audit_component_migration(model, ["family_a"], table_path=table)


def test_hydroxide_normalizes_only_safe_single_compartment_reactions(
    tmp_path: Path,
) -> None:
    model = _model_with_common_ions()
    hydroxide = _met("oh_c", "hydroxide", "HO", -1)
    second_hydroxide = _met("oh_alt_c", "hydroxide_HO", "HO", -1)
    same_species = _met("x_c", "x", "HO", -1)
    external_species = _met("x_e", "x external", "HO", -1, "e")
    model.add_metabolites(
        [hydroxide, second_hydroxide, same_species, external_species]
    )

    internal = Reaction("ROH")
    internal.add_metabolites(
        {hydroxide: -1, second_hydroxide: -2, same_species: 3}
    )
    boundary = Reaction("EX_OH")
    boundary.add_metabolites({hydroxide: -1})
    transport = Reaction("T_OH")
    transport.add_metabolites({hydroxide: -1, external_species: 1})
    model.add_reactions([internal, boundary, transport])

    table = _acid_base_table(tmp_path, "oh_alt_c;oh_c")
    apply_curated_microspecies(model, table)
    report = normalize_hydroxide_reactions(model, table)
    assert report["changed_reactions"] == 1
    assert report["rejected_reactions"] == 2
    assert hydroxide not in internal.metabolites
    assert second_hydroxide not in internal.metabolites
    assert internal.metabolites[model.metabolites.h2o_c] == -3
    assert internal.metabolites[model.metabolites.h_c] == 3
    assert internal.check_mass_balance() == {}
    assert hydroxide in boundary.metabolites
    assert hydroxide in transport.metabolites


def test_proton_water_gate_requires_simultaneous_mass_and_charge_balance(
    tmp_path: Path,
) -> None:
    model = _model_with_common_ions()
    water_like = _met("a_c", "A", "H2O", 0)
    oxide = _met("b_c", "B", "O", -2)
    neutral_oxide = _met("c_c", "C", "O", 0)
    model.add_metabolites([water_like, oxide, neutral_oxide])
    allowed = Reaction("R_ALLOWED")
    allowed.add_metabolites({water_like: -1, oxide: 1})
    rejected = Reaction("R_REJECTED")
    rejected.add_metabolites({water_like: -1, neutral_oxide: 1})
    model.add_reactions([allowed, rejected])

    table = _acid_base_table(tmp_path)
    apply_curated_microspecies(model, table)
    report = balance_protons_and_water(
        model,
        reaction_ids=["R_ALLOWED", "R_REJECTED"],
        table_path=table,
    )
    assert report["changed_reactions"] == 1
    assert allowed.metabolites[model.metabolites.h_c] == 2
    assert allowed.check_mass_balance() == {}
    assert model.metabolites.h_c not in rejected.metabolites
    assert any(
        row["reaction_id"] == "R_REJECTED" and row["reason"] == "proton_gate"
        for row in report["rejected"]
    )


def test_real_model_active_subset_is_safe_and_r540_deferred() -> None:
    model_path = REPO_ROOT / "model.xml"
    before_sha = hashlib.sha256(model_path.read_bytes()).hexdigest()
    model = read_sbml_model(str(model_path))
    r540_before = model.reactions.get_by_id("R540").reaction

    report = apply_curated_microspecies(model)
    assert report["balanced_reaction_regressions"] == []
    assert report["active_families"] == 7
    balance_report = balance_protons_and_water(model, reaction_ids=["R540"])
    assert model.reactions.get_by_id("R540").reaction == r540_before
    assert balance_report["changed_reactions"] == 0
    assert balance_report["rejected"][0]["reason"] == "proton_gate"
    assert model.metabolites.get_by_id("m1099[C_ex]").formula == "K"
    assert model.metabolites.get_by_id("m1108[C_ex]").charge == 1
    assert model.metabolites.get_by_id("m1893[C_cy]").charge == 2
    assert model.metabolites.get_by_id("m38[C_cy]").formula == "H3N"
    assert any(
        row["family_id"] == "gtp" and row["would_change"] for row in report["deferred"]
    )
    # Pipeline helpers operate in memory; model.xml is still generated only by main().
    assert hashlib.sha256(model_path.read_bytes()).hexdigest() == before_sha


def test_source_r540_is_not_rewritten_by_generic_balance() -> None:
    model = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))
    reaction = model.reactions.get_by_id("R540")
    before = {metabolite.id: coefficient for metabolite, coefficient in reaction.metabolites.items()}

    apply_curated_microspecies(model)
    report = balance_protons_and_water(model, reaction_ids=["R540"])

    after = {metabolite.id: coefficient for metabolite, coefficient in reaction.metabolites.items()}
    assert after == before
    assert report["changed_reactions"] == 0
    assert report["skipped_missing_formula"] == 1


def test_table_rejects_noncanonical_reference_ph(tmp_path: Path) -> None:
    table = _one_row_table(
        tmp_path,
        "1,active,potassium,base_name,potassium,K,1,7.0,CHEBI:29103,,"
        "https://www.ebi.ac.uk/chebi/CHEBI:29103,1,k_c,KH|0;K|1,test",
    )
    with pytest.raises(ValueError, match="reference_ph must be 7.3"):
        load_curated_microspecies(table)


def test_default_table_path_exists() -> None:
    assert DEFAULT_MICROSPECIES_TABLE.exists()
