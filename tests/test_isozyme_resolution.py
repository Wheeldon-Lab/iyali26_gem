from __future__ import annotations

import json

import pandas as pd
import pytest

from scripts.gem_annotate.isozyme_resolution import (
    build_isozyme_resolution_ledger,
    classify_isozyme_counterfactual,
)


def _diagnostic(ko_ratio: float, closed_ratio: float) -> dict:
    return {
        "ko_growth_ratio": ko_ratio,
        "all_linked_reactions_closed_growth_ratio": closed_ratio,
    }


def _case(
    case_id: str,
    *,
    ko_ratio: float,
    closed_ratio: float,
    balance: dict | None = None,
    balance_error: str = "",
) -> dict:
    return {
        "case_id": case_id,
        "category": "isozyme_redundancy",
        "gene_ids": [f"g-{case_id}"],
        "reaction_ids": [f"R-{case_id}"],
        "model_sha256": "model-sha",
        "target_fingerprint": f"sha256:{case_id}",
        "model_context": {
            "diagnostics": [_diagnostic(ko_ratio, closed_ratio)],
            "reactions": [
                {
                    "reaction_id": f"R-{case_id}",
                    "gpr": f"g-{case_id} or backup-{case_id}",
                    "mass_balance": balance or {},
                    "mass_balance_error": balance_error,
                }
            ],
        },
    }


def test_counterfactual_requires_normal_gene_ko_and_lethal_reaction_ko() -> None:
    assert (
        classify_isozyme_counterfactual(
            [_diagnostic(1.0, 0.0)], primary_cutoff=0.10
        )
        == "model_causal_isozyme_candidate"
    )
    assert (
        classify_isozyme_counterfactual(
            [_diagnostic(0.90, 0.0)], primary_cutoff=0.10
        )
        == "noncausal_redundancy_signal"
    )
    assert (
        classify_isozyme_counterfactual(
            [_diagnostic(1.0, 0.10)], primary_cutoff=0.10
        )
        == "noncausal_redundancy_signal"
    )
    assert (
        classify_isozyme_counterfactual(
            [_diagnostic(1.0, 0.0), _diagnostic(1.0, 1.0)],
            primary_cutoff=0.10,
        )
        == "mixed_counterfactual_signal"
    )


def test_resolution_ledger_keeps_chemistry_identity_and_evidence_separate(
    tmp_path,
) -> None:
    cases = [
        _case("causal-balanced", ko_ratio=1.0, closed_ratio=0.0),
        _case(
            "causal-imbalanced",
            ko_ratio=1.0,
            closed_ratio=0.0,
            balance={"H": -1.0},
        ),
        _case(
            "noncausal-uncheckable",
            ko_ratio=1.0,
            closed_ratio=1.0,
            balance_error="missing polymer formula",
        ),
    ]
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    (evidence_dir / "causal-balanced.json").write_text(
        json.dumps(
            {
                "identity_review": {"status": "verified"},
                "verdict": "supported_patch_candidate",
                "proposed_operation": {"operation": "set_gpr"},
                "human_decision": {"decision": "pending"},
            }
        ),
        encoding="utf-8",
    )
    (evidence_dir / "causal-imbalanced.json").write_text(
        json.dumps(
            {
                "chemistry_review": {
                    "status": "blocked_component_migration",
                    "model_sha256": "model-sha",
                    "current_residual": {"H": -1.0},
                },
                "verdict": "needs_more_evidence",
            }
        ),
        encoding="utf-8",
    )

    table = build_isozyme_resolution_ledger(
        cases,
        primary_cutoff=0.10,
        ledger_path=tmp_path / "missing-ledger.csv",
        evidence_dir=evidence_dir,
        patch_table_path=tmp_path / "missing-patches.csv",
    ).set_index("case_id")

    assert table.loc["causal-balanced", "chemistry_status"] == "balanced"
    assert table.loc["causal-balanced", "identity_status"] == "verified"
    assert table.loc["causal-balanced", "resolution_status"] == "open"
    assert "workflow" in table.loc["causal-balanced", "blocker"]

    assert table.loc["causal-imbalanced", "chemistry_status"] == (
        "blocked_component_migration"
    )
    assert table.loc["causal-imbalanced", "blocker"] == (
        "reaction chemistry is blocked_component_migration"
    )

    assert table.loc["noncausal-uncheckable", "chemistry_status"] == (
        "uncheckable"
    )
    assert table.loc["noncausal-uncheckable", "resolution_status"] == (
        "reclassified_noncausal"
    )


@pytest.mark.external_data
@pytest.mark.integration
def test_current_sha_isozyme_inventory_has_33_groups_and_10_causal(
    external_data_file,
) -> None:
    table = pd.read_csv(
        external_data_file("data/essentiality/isozyme_resolution_ledger.csv")
    )
    assert len(table) == 33
    assert set(table["model_sha256"]) == {
        "39f4cae11c3f270400c8a227c78b6af3ed412e85b1ade6cb604b0f85c3d8b1d9"
    }
    assert (
        table["model_causal_class"]
        .eq("model_causal_isozyme_candidate")
        .sum()
        == 10
    )
    assert table["resolution_status"].eq("reclassified_noncausal").sum() == 23
