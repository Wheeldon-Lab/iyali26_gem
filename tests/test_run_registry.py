from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.gem_annotate.run_registry import (
    DuplicateRunError,
    backfill,
    build_run_key,
    guard_duplicate_run,
    read_records,
    register_run,
)


def _details() -> tuple[dict[str, object], dict[str, str], dict[str, object]]:
    return ({"model": {"sha256": "a" * 64}}, {"runner.py": "b" * 64}, {"cutoff": 0.1})


def test_matching_successful_run_is_rejected_and_force_links_reproduction(tmp_path: Path) -> None:
    inputs, code, configuration = _details()
    key = build_run_key("essentiality", inputs=inputs, code_sources=code, configuration=configuration)
    original = tmp_path / "results" / "essentiality" / "original"
    reproduction = tmp_path / "results" / "essentiality" / "reproduction"
    register_run(tmp_path, workflow="essentiality", run_key=key, output_dir=original,
                 inputs=inputs, code_sources=code, configuration=configuration, status="complete")

    with pytest.raises(DuplicateRunError):
        guard_duplicate_run(tmp_path, workflow="essentiality", run_key=key, output_dir=reproduction)
    previous = guard_duplicate_run(
        tmp_path, workflow="essentiality", run_key=key, output_dir=reproduction,
        force_rerun=True, reproduction_reason="independent solver reproduction",
    )
    register_run(tmp_path, workflow="essentiality", run_key=key, output_dir=reproduction,
                 inputs=inputs, code_sources=code, configuration=configuration, status="complete",
                 previous=previous, reproduction_reason="independent solver reproduction")
    records = read_records(tmp_path)
    assert records[-1]["relationship"]["type"] == "reproduction_of"
    assert records[-1]["relationship"]["record_id"] == records[0]["record_id"]


def test_changed_code_sha_creates_a_distinct_run_key(tmp_path: Path) -> None:
    inputs, code, configuration = _details()
    first = build_run_key("flow", inputs=inputs, code_sources=code, configuration=configuration)
    changed = build_run_key("flow", inputs=inputs, code_sources={"runner.py": "c" * 64}, configuration=configuration)
    assert first != changed
    register_run(tmp_path, workflow="flow", run_key=first, output_dir=tmp_path / "first",
                 inputs=inputs, code_sources=code, configuration=configuration, status="complete")
    assert guard_duplicate_run(tmp_path, workflow="flow", run_key=changed, output_dir=tmp_path / "changed") is None


def test_backfill_appends_without_rewriting_historical_manifest(tmp_path: Path) -> None:
    result_dir = tmp_path / "artifacts" / "results" / "essentiality" / "historical"
    result_dir.mkdir(parents=True)
    manifest_path = result_dir / "run_manifest.json"
    manifest = {
        "inputs": {"model": {"sha256": "a" * 64}},
        "configuration": {"positive_only": True},
        "cutoffs": {"primary_fraction_of_wt": 0.1},
        "solver": {"name": "glpk"},
        "git": {"commit": "example"},
    }
    manifest_path.write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")
    before = manifest_path.read_bytes()
    (tmp_path / "relocation_manifest.csv").write_text(
        "source_path,destination_path,sha256,duplicate_of\nold.xlsx,new.xlsx,deadbeef,old.xlsx\n",
        encoding="utf-8",
    )

    assert backfill(tmp_path) == {"runs": 1, "duplicates": 1}
    assert manifest_path.read_bytes() == before
    assert backfill(tmp_path) == {"runs": 0, "duplicates": 0}
