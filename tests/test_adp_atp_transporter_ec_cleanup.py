from __future__ import annotations

from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate import patch_runner
from scripts.gem_annotate.patches import remove_stale_adp_atp_transporter_ec_codes
from scripts.gem_annotate.sbml import write_deterministic_sbml_model


REPO_ROOT = Path(__file__).resolve().parents[1]


def _seed_stale_transport_ec_annotations(model) -> None:
    """Model the post-xref state that the late pipeline cleanup receives."""

    model.reactions.R815.annotation = {
        "ec-code": ["2.7.4.6", "7.6.2.1"],
        "rhea": ["34999"],
    }
    model.reactions.R816.annotation = {
        "ec-code": "2.7.4.6",
        "rhea": ["35000"],
    }


def _reaction_state(reaction) -> tuple:
    """Everything the metadata-only cleanup must leave untouched."""

    return (
        reaction.name,
        reaction.lower_bound,
        reaction.upper_bound,
        reaction.gene_reaction_rule,
        tuple(
            sorted(
                (metabolite.id, coefficient)
                for metabolite, coefficient in reaction.metabolites.items()
            )
        ),
        {
            key: value
            for key, value in reaction.annotation.items()
            if key != "ec-code"
        },
    )


def test_adp_atp_transporter_ec_cleanup_is_precise_and_idempotent() -> None:
    model = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))
    _seed_stale_transport_ec_annotations(model)
    before = {
        reaction_id: _reaction_state(model.reactions.get_by_id(reaction_id))
        for reaction_id in ("R815", "R816")
    }

    assert remove_stale_adp_atp_transporter_ec_codes(model) == 2
    for reaction_id, expected_state in before.items():
        reaction = model.reactions.get_by_id(reaction_id)
        assert _reaction_state(reaction) == expected_state
    assert model.reactions.R815.annotation["ec-code"] == ["7.6.2.1"]
    assert "ec-code" not in model.reactions.R816.annotation

    assert remove_stale_adp_atp_transporter_ec_codes(model) == 0


def test_patch_runner_exposes_adp_atp_ec_cleanup(tmp_path: Path) -> None:
    source = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))
    _seed_stale_transport_ec_annotations(source)
    source_path = tmp_path / "source.xml"
    output_path = tmp_path / "patched.xml"
    write_deterministic_sbml_model(source, source_path)

    audit = patch_runner.run_patch(
        "adp-atp-transporter-ec-cleanup",
        input_model=source_path,
        output_model=output_path,
    )

    assert audit["changes"] == 2
    patched = read_sbml_model(str(output_path))
    assert patched.reactions.R815.annotation["ec-code"] == "7.6.2.1"
    assert "ec-code" not in patched.reactions.R816.annotation
