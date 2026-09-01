from __future__ import annotations

from pathlib import Path

import pytest
from cobra.io import read_sbml_model

from scripts.gem_annotate.microspecies import (
    _is_allowed_current_pair,
    _resolve_pinned_targets,
    load_curated_microspecies,
)
from scripts.gem_annotate.reaction_chemistry import audit_reference_reaction_chemistry


REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = REPO_ROOT / "model.xml"


@pytest.mark.external_data
@pytest.mark.integration
def test_component_curation_explicitly_allows_the_current_model_baseline(
    external_data_file,
) -> None:
    model = read_sbml_model(str(MODEL_PATH))

    table_path = external_data_file("data/metabolite_microspecies.csv")
    for row in load_curated_microspecies(table_path):
        for metabolite in _resolve_pinned_targets(model, row):
            assert _is_allowed_current_pair(row, metabolite), (
                row.family_id,
                metabolite.id,
                metabolite.formula,
                metabolite.charge,
            )


@pytest.mark.external_data
@pytest.mark.integration
def test_r742_reference_requires_two_cytosolic_protons(external_data_file) -> None:
    report = audit_reference_reaction_chemistry(
        read_sbml_model(str(MODEL_PATH)),
        "EGC-r742-ssadh",
        table_path=external_data_file(
            "data/essentiality/reaction_chemistry_curation.csv"
        ),
        microspecies_table_path=external_data_file(
            "data/metabolite_microspecies.csv"
        ),
    )

    assert report["reference_reaction_ids"] == {"R742": "RHEA:13213"}
    assert report["target_reaction_balances"]["R742"]["after"]["status"] == "balanced"
    assert report["coefficient_changes"] == [
        {
            "reaction_id": "R742",
            "metabolite_id": "m10[C_cy]",
            "before_coefficient": 0.0,
            "after_coefficient": 2.0,
        }
    ]


@pytest.mark.external_data
@pytest.mark.integration
def test_r_g6pper_reference_has_no_er_proton_coefficient(
    external_data_file,
) -> None:
    report = audit_reference_reaction_chemistry(
        read_sbml_model(str(MODEL_PATH)),
        "EGC-g6pper-chemistry",
        table_path=external_data_file(
            "data/essentiality/reaction_chemistry_curation.csv"
        ),
        microspecies_table_path=external_data_file(
            "data/metabolite_microspecies.csv"
        ),
    )

    assert report["reference_reaction_ids"] == {"R_G6PPer": "RHEA:16689"}
    assert report["target_reaction_balances"]["R_G6PPer"]["after"]["status"] == "balanced"
    assert report["coefficient_changes"] == []
