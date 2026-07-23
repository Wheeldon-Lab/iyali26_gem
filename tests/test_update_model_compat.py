from __future__ import annotations

import importlib

import scripts.update_model as legacy_entrypoint

canonical_pipeline = importlib.import_module("scripts.gem_annotate.main")


def test_legacy_update_model_entrypoint_delegates_to_canonical_pipeline(
    monkeypatch,
):
    sentinel = object()
    monkeypatch.setattr(canonical_pipeline, "main", lambda: sentinel)

    assert legacy_entrypoint.main() is sentinel
