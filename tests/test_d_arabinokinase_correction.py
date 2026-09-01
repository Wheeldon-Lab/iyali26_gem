from __future__ import annotations

from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate import patch_runner
from scripts.gem_annotate.patches import (
    fix_d_arabinokinase_direction_and_proton,
)
from scripts.gem_annotate.microspecies import balance_protons_and_water
from scripts.gem_annotate.sbml import write_deterministic_sbml_model


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_d_arabinokinase_is_forward_only_with_verified_product_proton() -> None:
    model = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))
    reaction = model.reactions.get_by_id("R2041")
    proton = model.metabolites.get_by_id("m10[C_cy]")

    assert reaction.bounds == (-1000.0, 1000.0)
    assert proton not in reaction.metabolites

    assert fix_d_arabinokinase_direction_and_proton(model) == 1
    assert reaction.bounds == (0.0, 1000.0)
    assert reaction.metabolites[proton] == 1.0
    assert reaction.notes["curated_reaction_correction"]["source"] == (
        "https://enzyme.expasy.org/EC/2.7.1.54"
    )
    assert reaction.notes["curated_reaction_correction"]["remaining_gate"] == (
        "ATP/ADP connected-component microspecies migration"
    )
    assert reaction.notes["curated_reaction_correction"][
        "lock_proton_water_stoichiometry"
    ] is True

    assert fix_d_arabinokinase_direction_and_proton(model) == 0

    assert reaction.metabolites[proton] == 1.0


def test_curated_proton_lock_survives_generic_balance(tmp_path: Path) -> None:
    model = Model("curated-lock")
    proton = Metabolite("h_c", name="H+_p+1", formula="H", charge=1, compartment="c")
    water = Metabolite(
        "h2o_c", name="H2O_H2O", formula="H2O", charge=0, compartment="c"
    )
    substrate = Metabolite("a_c", formula="CH2", charge=0, compartment="c")
    product = Metabolite("b_c", formula="C", charge=0, compartment="c")
    reaction = Reaction("LOCKED")
    reaction.add_metabolites({substrate: -1, product: 1, proton: 1})
    reaction.notes["curated_reaction_correction"] = {
        "lock_proton_water_stoichiometry": True
    }
    model.add_metabolites([water])
    model.add_reactions([reaction])

    table_path = tmp_path / "microspecies.csv"
    table_path.write_text(
        "schema_version,status,family_id,selector_type,selector_value,"
        "target_formula,target_charge,reference_ph,chebi_id,rhea_id,source_url,"
        "min_matches,expected_metabolite_ids,allowed_current_pairs,rationale\n"
        "1,active,proton,base_name,H+,H,1,7.3,CHEBI:15378,,"
        "https://example.org,1,h_c,H|1,test\n"
        "1,active,water,base_name,H2O,H2O,0,7.3,CHEBI:15377,,"
        "https://example.org,1,h2o_c,H2O|0,test\n",
        encoding="utf-8",
    )

    before = dict(reaction.metabolites)
    report = balance_protons_and_water(
        model, reaction_ids=[reaction.id], table_path=table_path
    )
    assert dict(reaction.metabolites) == before
    assert report["skipped_curated_lock_reaction_ids"] == [reaction.id]


def test_d_arabinokinase_patch_runner_writes_only_a_new_output(
    tmp_path: Path,
) -> None:
    source = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))
    source_path = tmp_path / "source.xml"
    output_path = tmp_path / "patched.xml"
    write_deterministic_sbml_model(source, source_path)

    audit = patch_runner.run_patch(
        "d-arabinokinase-direction",
        input_model=source_path,
        output_model=output_path,
    )

    assert audit["changes"] == 1
    patched = read_sbml_model(str(output_path))
    reaction = patched.reactions.get_by_id("R2041")
    assert reaction.bounds == (0.0, 1000.0)
    assert reaction.metabolites[
        patched.metabolites.get_by_id("m10[C_cy]")
    ] == 1.0
