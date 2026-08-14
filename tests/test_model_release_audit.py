from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import write_sbml_model

from scripts.gem_annotate.model_release_audit import compare_models


def _write_model(path: Path, *, extra_change: bool = False) -> None:
    model = Model("release_audit")
    a = Metabolite("a_c", formula="C", charge=0, compartment="c")
    b = Metabolite("b_c", formula="C", charge=0, compartment="c")
    reaction = Reaction("R1")
    reaction.bounds = (0.0, 1000.0)
    reaction.add_metabolites({a: -1.0, b: 1.0})
    if extra_change:
        reaction.upper_bound = 2.0
    model.add_reactions([reaction])
    write_sbml_model(model, str(path))


def test_release_audit_rejects_nonwhitelisted_bound_change(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline.xml"
    candidate = tmp_path / "candidate.xml"
    _write_model(baseline)
    _write_model(candidate, extra_change=True)

    report = compare_models(baseline, candidate)

    assert report["reaction_diff"]["bounds_changed"] == ["R1"]
    assert report["release_gate"]["passed"] is False
    assert any(
        violation.startswith("reaction.bounds_changed")
        for violation in report["release_gate"]["violations"]
    )


def test_release_audit_records_repeat_sha_identity(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline.xml"
    candidate = tmp_path / "candidate.xml"
    repeat = tmp_path / "repeat.xml"
    _write_model(baseline)
    _write_model(candidate)
    repeat.write_bytes(candidate.read_bytes())

    report = compare_models(baseline, candidate, repeat)

    assert report["candidate_reproducibility"]["sha256_identical"] is True
    assert report["candidate_reproducibility"]["repeat_sha256"] == (
        report["candidate"]["sha256"]
    )
