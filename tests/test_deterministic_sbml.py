from __future__ import annotations

import stat
from pathlib import Path

import pytest
from cobra import Model, Reaction
from cobra.core import Group
from cobra.io import read_sbml_model

from scripts.gem_annotate import sbml


def _grouped_model(member_order: tuple[str, ...]) -> tuple[Model, Group]:
    model = Model("deterministic-groups")
    reactions = [Reaction("R1"), Reaction("R2"), Reaction("R3")]
    model.add_reactions(reactions)
    group = Group(
        "pathway",
        name="pathway",
        members=reactions,
        kind="collection",
    )
    model.add_groups([group])
    # Deliberately emulate two possible orders exposed by COBRApy's set-backed
    # Group.members property.  The production helper must erase this difference.
    group._members = tuple(model.reactions.get_by_id(rid) for rid in member_order)
    return model, group


def test_group_member_order_does_not_change_sbml_bytes(tmp_path: Path) -> None:
    first, _ = _grouped_model(("R3", "R1", "R2"))
    second, _ = _grouped_model(("R2", "R3", "R1"))
    first_path = tmp_path / "first.xml"
    second_path = tmp_path / "second.xml"

    sbml.write_deterministic_sbml_model(first, first_path)
    sbml.write_deterministic_sbml_model(second, second_path)

    assert first_path.read_bytes() == second_path.read_bytes()
    assert stat.S_IMODE(first_path.stat().st_mode) == 0o644
    roundtrip = read_sbml_model(str(first_path))
    assert {member.id for member in roundtrip.groups.pathway.members} == {
        "R1",
        "R2",
        "R3",
    }


def test_writer_restores_original_member_container(tmp_path: Path) -> None:
    model, group = _grouped_model(("R2", "R1", "R3"))
    original_members = group._members
    target = tmp_path / "model.xml"
    target.write_text("old model", encoding="utf-8")
    target.chmod(0o640)

    sbml.write_deterministic_sbml_model(model, target)

    assert group._members is original_members
    assert stat.S_IMODE(target.stat().st_mode) == 0o640


def test_writer_preserves_compression_suffix(tmp_path: Path) -> None:
    model, _ = _grouped_model(("R3", "R2", "R1"))
    target = tmp_path / "model.xml.gz"

    sbml.write_deterministic_sbml_model(model, target)

    assert target.read_bytes().startswith(b"\x1f\x8b")
    assert {reaction.id for reaction in read_sbml_model(str(target)).reactions} == {
        "R1",
        "R2",
        "R3",
    }


def test_failed_write_is_atomic_and_restores_groups(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    model, group = _grouped_model(("R2", "R1", "R3"))
    original_members = group._members
    target = tmp_path / "model.xml"
    target.write_text("existing model", encoding="utf-8")

    def fail_after_partial_write(_model, filename, **_kwargs) -> None:
        Path(filename).write_text("partial", encoding="utf-8")
        raise RuntimeError("simulated writer failure")

    monkeypatch.setattr(sbml, "write_sbml_model", fail_after_partial_write)

    with pytest.raises(RuntimeError, match="simulated writer failure"):
        sbml.write_deterministic_sbml_model(model, target)

    assert target.read_text(encoding="utf-8") == "existing model"
    assert group._members is original_members
    assert [
        path for path in tmp_path.iterdir() if path.name.startswith(".model.xml.")
    ] == []
