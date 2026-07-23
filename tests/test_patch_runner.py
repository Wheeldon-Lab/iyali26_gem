from __future__ import annotations

from pathlib import Path

import pytest
from cobra import Model, Reaction

from scripts.gem_annotate import patch_runner
from scripts.gem_annotate.sbml import write_deterministic_sbml_model


def _input_model(path: Path) -> Path:
    model = Model("patch-runner-test")
    model.add_reactions([Reaction("R_test")])
    write_deterministic_sbml_model(model, path)
    return path


def test_patch_runner_refuses_canonical_model_overwrite(tmp_path: Path) -> None:
    source = _input_model(tmp_path / "input.xml")

    with pytest.raises(ValueError, match="canonical model.xml"):
        patch_runner.run_patch(
            "c161-pool-extension",
            input_model=source,
            output_model=patch_runner.REPO_ROOT / "model.xml",
        )


def test_patch_runner_refuses_existing_output_without_touching_it(tmp_path: Path) -> None:
    source = _input_model(tmp_path / "input.xml")
    target = tmp_path / "already-there.xml"
    target.write_text("preserve me", encoding="utf-8")

    with pytest.raises(FileExistsError, match="refuses to overwrite"):
        patch_runner.run_patch(
            "c161-pool-extension", input_model=source, output_model=target
        )
    assert target.read_text(encoding="utf-8") == "preserve me"


def test_legacy_patch_entrypoints_only_delegate_to_shared_runner() -> None:
    root = Path(__file__).resolve().parents[1]
    for filename in (
        "apply_c161_pool_extension.py",
        "apply_ec_overload_cleanup.py",
        "apply_isozyme_gprs.py",
        "annotate_new_isozyme_genes.py",
    ):
        text = (root / "scripts" / filename).read_text(encoding="utf-8")
        assert "main_for_legacy" in text
        assert "write_sbml_model" not in text
        assert "model.xml\"" not in text


def test_no_script_directly_writes_the_canonical_model() -> None:
    root = Path(__file__).resolve().parents[1]
    violations = []
    for source in (root / "scripts").rglob("*.py"):
        text = source.read_text(encoding="utf-8")
        if "model.xml" in text and "write_sbml_model" in text:
            violations.append(source.relative_to(root))
    assert violations == []
