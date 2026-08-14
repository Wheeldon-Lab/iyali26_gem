from __future__ import annotations

from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.microspecies import (
    _is_allowed_current_pair,
    _resolve_pinned_targets,
    load_curated_microspecies,
)
from scripts.gem_annotate.reaction_chemistry import audit_reference_reaction_chemistry


REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = REPO_ROOT / "model.xml"


def test_component_curation_explicitly_allows_the_current_model_baseline() -> None:
    model = read_sbml_model(str(MODEL_PATH))

    for row in load_curated_microspecies():
        for metabolite in _resolve_pinned_targets(model, row):
            assert _is_allowed_current_pair(row, metabolite), (
                row.family_id,
                metabolite.id,
                metabolite.formula,
                metabolite.charge,
            )


def test_r742_reference_requires_two_cytosolic_protons() -> None:
    report = audit_reference_reaction_chemistry(
        read_sbml_model(str(MODEL_PATH)), "EGC-r742-ssadh"
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


def test_r_g6pper_reference_has_no_er_proton_coefficient() -> None:
    report = audit_reference_reaction_chemistry(
        read_sbml_model(str(MODEL_PATH)), "EGC-g6pper-chemistry"
    )

    assert report["reference_reaction_ids"] == {"R_G6PPer": "RHEA:16689"}
    assert report["target_reaction_balances"]["R_G6PPer"]["after"]["status"] == "balanced"
    assert report["coefficient_changes"] == []
