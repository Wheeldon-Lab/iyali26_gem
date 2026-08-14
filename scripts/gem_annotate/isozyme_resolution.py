"""Deterministic resolution ledger for essentiality isozyme signals.

The FN diagnostic calls a gene ``isozyme_redundancy`` whenever its linked
reactions retain FVA capacity after the gene knockout.  That signal alone does
not show that the annotated isozyme is the reason the knockout survives.  This
module groups the fresh cases by reaction/GPR and separates the high-confidence
counterfactual signal (gene KO grows normally, but closing all linked reactions
is lethal) from non-causal residual capacity.

The ledger is read-only with respect to the model.  It combines fresh simulation
packets with the durable evidence/state records so chemistry, gene identity,
literature, human approval and pipeline regression remain distinct gates.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .config import ESSENTIALITY_DIR
from .essentiality_evidence import (
    DEFAULT_EVIDENCE_DIR,
    DEFAULT_LEDGER,
    canonical_json,
    read_ledger,
)


DEFAULT_PATCH_TABLE = ESSENTIALITY_DIR / "curated_model_patches.csv"
ISOZYME_SURVIVAL_FLOOR = 0.90

ISOZYME_RESOLUTION_FIELDS = (
    "case_id",
    "model_sha256",
    "target_fingerprint",
    "gene_ids",
    "reaction_ids",
    "full_gpr",
    "n_fn_genes",
    "ko_ratio_min",
    "ko_ratio_max",
    "all_linked_closed_ratio_min",
    "all_linked_closed_ratio_max",
    "model_causal_class",
    "chemistry_status",
    "chemistry_imbalances",
    "identity_status",
    "evidence_verdict",
    "workflow_status",
    "proposed_operation",
    "human_decision",
    "pipeline_patch_status",
    "regression_status",
    "resolution_status",
    "blocker",
    "notes",
)


def classify_isozyme_counterfactual(
    diagnostic_rows: Iterable[dict[str, Any]],
    *,
    primary_cutoff: float,
    survival_floor: float = ISOZYME_SURVIVAL_FLOOR,
) -> str:
    """Classify whether residual capacity causally explains the simulated FN.

    ``model_causal_isozyme_candidate`` is deliberately a simulation statement,
    not a biological verdict.  Biology still requires balanced chemistry,
    verified gene identities and direct operation-specific evidence.
    """

    rows = list(diagnostic_rows)
    flags = [
        float(row.get("ko_growth_ratio", 0.0)) > survival_floor
        and float(row.get("all_linked_reactions_closed_growth_ratio", 1.0))
        < primary_cutoff
        for row in rows
    ]
    if flags and all(flags):
        return "model_causal_isozyme_candidate"
    if any(flags):
        return "mixed_counterfactual_signal"
    return "noncausal_redundancy_signal"


def _chemistry_summary(
    reaction_contexts: Iterable[dict[str, Any]],
) -> tuple[str, str]:
    imbalances: dict[str, Any] = {}
    uncheckable = False
    for reaction in reaction_contexts:
        reaction_id = str(reaction.get("reaction_id", ""))
        error = str(reaction.get("mass_balance_error", "") or "").strip()
        raw_balance = reaction.get("mass_balance", {})
        balance = {
            str(key): float(value)
            for key, value in (raw_balance.items() if isinstance(raw_balance, dict) else [])
            if abs(float(value)) > 1e-9
        }
        if error:
            uncheckable = True
            imbalances[reaction_id] = {"error": error}
        elif balance:
            imbalances[reaction_id] = balance

    if uncheckable:
        status = "uncheckable"
    elif imbalances:
        status = "imbalanced"
    else:
        status = "balanced"
    return status, canonical_json(imbalances)


def _identity_status(dossier: dict[str, Any]) -> str:
    explicit = dossier.get("identity_review", {})
    if isinstance(explicit, dict) and explicit.get("status"):
        return str(explicit["status"])

    crosschecks = dossier.get("identity_crosschecks", [])
    if not isinstance(crosschecks, list) or not crosschecks:
        return "pending"
    statuses = {
        str(row.get("status", "")).strip().casefold()
        for row in crosschecks
        if isinstance(row, dict)
    }
    if statuses & {"conflict", "mismatch", "wrong_mapping", "failed"}:
        return "conflict"
    # Database mappings are useful but cannot by themselves certify that every
    # GPR partner has the represented function in Y. lipolytica.
    return "partial_crosscheck"


def _reviewed_chemistry_summary(
    dossier: dict[str, Any],
    *,
    model_sha256: str,
    fallback_status: str,
    fallback_imbalances: str,
) -> tuple[str, str]:
    """Prefer an explicit current-SHA chemistry review over a raw residual."""

    review = dossier.get("chemistry_review", {})
    if not isinstance(review, dict) or not str(review.get("status", "")).strip():
        return fallback_status, fallback_imbalances
    if str(review.get("model_sha256", "")).strip() != model_sha256:
        return "stale_review", fallback_imbalances

    reviewed_imbalances = review.get(
        "current_residuals", review.get("current_residual", None)
    )
    if isinstance(reviewed_imbalances, dict):
        imbalance_text = canonical_json(reviewed_imbalances)
    else:
        imbalance_text = fallback_imbalances
    return str(review["status"]).strip(), imbalance_text


def _load_dossier(evidence_dir: Path, case_id: str) -> dict[str, Any]:
    path = evidence_dir / f"{case_id}.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _load_patch_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _pipeline_patch_status(
    case_id: str, patch_rows: Iterable[dict[str, str]]
) -> str:
    matching = [row for row in patch_rows if row.get("case_id") == case_id]
    if not matching:
        return "none"
    values = [
        f"{row.get('patch_id', '')}:{row.get('status', '')}" for row in matching
    ]
    return ";".join(sorted(values))


def _resolution(
    *,
    causal_class: str,
    chemistry_status: str,
    identity_status: str,
    evidence_verdict: str,
    workflow_status: str,
    patch_status: str,
) -> tuple[str, str, str]:
    if causal_class == "noncausal_redundancy_signal":
        return (
            "reclassified_noncausal",
            "",
            "Residual FVA capacity is not the cause of survival; hand the FN to an alternative-mechanism diagnosis.",
        )
    if causal_class == "mixed_counterfactual_signal":
        return (
            "open",
            "mixed gene-level counterfactual results",
            "Split or re-diagnose the grouped case before biological curation.",
        )
    if evidence_verdict in {
        "retain_current_model",
        "experimental_conflict",
        "outside_metabolic_scope",
    } and workflow_status in {
        "reviewed",
        "needs_more_evidence",
        "rejected",
        "regression_passed",
    }:
        return "resolved_no_patch", "", "Adversarially reviewed no-patch outcome."
    if chemistry_status not in {"balanced", "verified_balanced"}:
        return (
            "open",
            f"reaction chemistry is {chemistry_status}",
            "Resolve formulas, charges and stoichiometry before interpreting the GPR.",
        )
    if identity_status != "verified":
        return (
            "open",
            f"gene identity is {identity_status}",
            "Verify target and partner accessions/functions before proposing a patch.",
        )
    if evidence_verdict != "supported_patch_candidate":
        return (
            "open",
            "operation-specific direct evidence is incomplete",
            "Run literature review and independent skeptic review.",
        )
    if workflow_status == "regression_passed":
        return "resolved_regression_passed", "", "Accepted pipeline patch passed regression."
    if patch_status != "none":
        return (
            "open",
            "curated patch has not completed regression",
            "Rebuild through the pipeline and run the full Gurobi regression suite.",
        )
    if workflow_status == "awaiting_human":
        return "awaiting_human", "explicit human decision required", ""
    return "open", "evidence/approval workflow incomplete", ""


def build_isozyme_resolution_ledger(
    cases: Iterable[dict[str, Any]],
    *,
    primary_cutoff: float,
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
    patch_table_path: str | Path = DEFAULT_PATCH_TABLE,
) -> pd.DataFrame:
    """Combine fresh isozyme case packets with durable review state."""

    evidence_path = Path(evidence_dir)
    ledger_rows = {
        row["case_id"]: row for row in read_ledger(Path(ledger_path))
    }
    patch_rows = _load_patch_rows(Path(patch_table_path))
    output: list[dict[str, Any]] = []

    for case in sorted(cases, key=lambda row: str(row.get("case_id", ""))):
        if case.get("category") != "isozyme_redundancy":
            continue
        case_id = str(case["case_id"])
        model_context = case.get("model_context", {})
        diagnostic_rows = list(model_context.get("diagnostics", []))
        reactions = list(model_context.get("reactions", []))
        causal_class = classify_isozyme_counterfactual(
            diagnostic_rows,
            primary_cutoff=primary_cutoff,
        )
        dossier = _load_dossier(evidence_path, case_id)
        chemistry_status, chemistry_imbalances = _chemistry_summary(reactions)
        chemistry_status, chemistry_imbalances = _reviewed_chemistry_summary(
            dossier,
            model_sha256=str(case.get("model_sha256", "")),
            fallback_status=chemistry_status,
            fallback_imbalances=chemistry_imbalances,
        )
        durable = ledger_rows.get(case_id, {})
        identity_status = _identity_status(dossier)
        evidence_verdict = str(dossier.get("verdict", "") or "not_reviewed")
        workflow_status = str(
            durable.get("status", dossier.get("workflow_status", "detected"))
            or "detected"
        )
        proposed_operation = canonical_json(dossier.get("proposed_operation", {}))
        decision = dossier.get("human_decision", {})
        human_decision = (
            str(decision.get("decision", "pending"))
            if isinstance(decision, dict)
            else str(decision or "pending")
        )
        patch_status = _pipeline_patch_status(case_id, patch_rows)
        resolution_status, blocker, notes = _resolution(
            causal_class=causal_class,
            chemistry_status=chemistry_status,
            identity_status=identity_status,
            evidence_verdict=evidence_verdict,
            workflow_status=workflow_status,
            patch_status=patch_status,
        )

        ko_ratios = [float(row["ko_growth_ratio"]) for row in diagnostic_rows]
        closed_ratios = [
            float(row["all_linked_reactions_closed_growth_ratio"])
            for row in diagnostic_rows
        ]
        gprs = {
            str(reaction.get("reaction_id", "")): str(reaction.get("gpr", ""))
            for reaction in reactions
        }
        output.append(
            {
                "case_id": case_id,
                "model_sha256": case.get("model_sha256", ""),
                "target_fingerprint": case.get("target_fingerprint", ""),
                "gene_ids": ";".join(case.get("gene_ids", [])),
                "reaction_ids": ";".join(case.get("reaction_ids", [])),
                "full_gpr": canonical_json(gprs),
                "n_fn_genes": len(case.get("gene_ids", [])),
                "ko_ratio_min": min(ko_ratios),
                "ko_ratio_max": max(ko_ratios),
                "all_linked_closed_ratio_min": min(closed_ratios),
                "all_linked_closed_ratio_max": max(closed_ratios),
                "model_causal_class": causal_class,
                "chemistry_status": chemistry_status,
                "chemistry_imbalances": chemistry_imbalances,
                "identity_status": identity_status,
                "evidence_verdict": evidence_verdict,
                "workflow_status": workflow_status,
                "proposed_operation": proposed_operation,
                "human_decision": human_decision,
                "pipeline_patch_status": patch_status,
                "regression_status": (
                    "passed" if workflow_status == "regression_passed" else "not_run"
                ),
                "resolution_status": resolution_status,
                "blocker": blocker,
                "notes": notes,
            }
        )

    return pd.DataFrame(output, columns=ISOZYME_RESOLUTION_FIELDS)


__all__ = [
    "ISOZYME_RESOLUTION_FIELDS",
    "ISOZYME_SURVIVAL_FLOOR",
    "build_isozyme_resolution_ledger",
    "classify_isozyme_counterfactual",
]
