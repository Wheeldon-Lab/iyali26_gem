from pathlib import Path

import pytest

from scripts.gem_annotate.config import RESEARCH_ROOT_ENV, load_project_paths


@pytest.fixture
def external_data_file(request):
    """Resolve one required legacy path from the external research workspace."""

    if request.node.get_closest_marker("external_data") is None:
        pytest.fail("external_data_file requires @pytest.mark.external_data")

    def require(relative_legacy_path: str | Path) -> Path:
        relative_path = Path(relative_legacy_path)
        if relative_path.is_absolute():
            raise ValueError("external_data_file requires a relative legacy path")

        paths = load_project_paths()
        if not paths.configured:
            pytest.skip(
                f"{RESEARCH_ROOT_ENV} is not configured; required external data: "
                f"{relative_path.as_posix()}"
            )

        resolved = paths.resolve_legacy_path(relative_path)
        if not resolved.is_file():
            pytest.fail(
                f"Required external data file is missing: {resolved}", pytrace=False
            )
        return resolved

    return require
