from pathlib import Path

import pytest

from scripts.gem_annotate.config import load_project_paths


def test_explicit_research_root_maps_legacy_paths(tmp_path: Path) -> None:
    root = tmp_path / "research"
    root.mkdir()
    paths = load_project_paths(root, required=True)

    assert paths.resolve_legacy_path(
        "data/essentiality/curation_cases.csv"
    ) == (root / "state/essentiality/repository/curation_cases.csv").resolve()
    assert paths.resolve_legacy_path("data/media/sd_leu.csv") == (
        root / "state/media/sd_leu.csv"
    ).resolve()
    assert paths.resolve_legacy_path("data/metabolite_microspecies.csv") == (
        root / "state/curation/data/metabolite_microspecies.csv"
    ).resolve()
    assert paths.resolve_legacy_path("model.xml") == (
        paths.repo_root / "model.xml"
    ).resolve()
    assert paths.resolve_legacy_path("data/iyali26.xml") == (
        paths.repo_root / "data/iyali26.xml"
    ).resolve()


def test_required_research_root_fails_closed(tmp_path: Path) -> None:
    missing = tmp_path / "missing"
    with pytest.raises(FileNotFoundError, match="Research workspace not found"):
        load_project_paths(missing, required=True)


def test_explicit_research_root_overrides_environment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    environment_root = tmp_path / "environment"
    explicit_root = tmp_path / "explicit"
    environment_root.mkdir()
    explicit_root.mkdir()
    monkeypatch.setenv("IYALI26_RESEARCH_ROOT", str(environment_root))

    assert load_project_paths(explicit_root, required=True).research_root == (
        explicit_root.resolve()
    )
