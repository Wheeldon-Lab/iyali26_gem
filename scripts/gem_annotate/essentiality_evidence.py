"""Durable evidence and human-approval controls for essentiality curation.

This module deliberately contains no model-editing code.  It validates evidence
dossiers, tracks the review state machine, and records explicit human decisions.
The build pipeline imports the same validators before applying any schema-v2
essentiality patch.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
import os
import re
import stat
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from .config import ESSENTIALITY_DIR, MEDIA_DIR, REPO_ROOT

DEFAULT_LEDGER = ESSENTIALITY_DIR / "curation_cases.csv"
DEFAULT_EVIDENCE_DIR = ESSENTIALITY_DIR / "evidence"
DEFAULT_MODEL = REPO_ROOT / "model.xml"
DEFAULT_EXPERIMENTAL = ESSENTIALITY_DIR / "consensus_essential_genes.csv"
DEFAULT_MEDIA = MEDIA_DIR / "sd_leu.csv"

CASE_STATUSES = (
    "detected",
    "queued",
    "researching",
    "reviewed",
    "awaiting_human",
    "needs_more_evidence",
    "rejected",
    "accepted",
    "implemented",
    "regression_passed",
)

VERDICTS = (
    "supported_patch_candidate",
    "retain_current_model",
    "experimental_conflict",
    "outside_metabolic_scope",
    "needs_more_evidence",
)

ALLOWED_TRANSITIONS = {
    "detected": {"queued"},
    "queued": {"researching"},
    "researching": {"reviewed"},
    "reviewed": {"awaiting_human", "needs_more_evidence", "rejected"},
    "awaiting_human": {"accepted", "needs_more_evidence", "rejected"},
    "needs_more_evidence": {"queued", "rejected"},
    "accepted": {"implemented"},
    "implemented": {"regression_passed"},
    "rejected": set(),
    "regression_passed": set(),
}

LEDGER_FIELDS = (
    "case_id",
    "status",
    "category",
    "gene_ids",
    "reaction_ids",
    "model_sha256",
    "experimental_sha256",
    "media_sha256",
    "evidence_schema_version",
    "simulation_context_fingerprint_version",
    "simulation_context_fingerprint",
    "strain_overlay_enabled",
    "strain_profile_id",
    "strain_profile_sha256",
    "strain_overlay_effect_fingerprint_version",
    "strain_overlay_effect_sha256",
    "case_packet_sha256",
    "target_fingerprint",
    "chemistry_fingerprint",
    "evidence_path",
    "detected_at",
    "updated_at",
    "previous_status",
    "stale_reason",
    "human_decision",
    "approved_by",
    "approved_at",
)

_LEGACY_OPTIONAL_LEDGER_FIELDS = {
    "chemistry_fingerprint",
    "evidence_schema_version",
    "simulation_context_fingerprint_version",
    "simulation_context_fingerprint",
    "strain_overlay_enabled",
    "strain_profile_id",
    "strain_profile_sha256",
    "strain_overlay_effect_fingerprint_version",
    "strain_overlay_effect_sha256",
    "case_packet_sha256",
}

SIMULATION_CONTEXT_FIELDS = (
    "simulation_context_fingerprint_version",
    "simulation_context_fingerprint",
    "strain_overlay_enabled",
    "strain_profile_id",
    "strain_profile_sha256",
    "strain_overlay_effect_fingerprint_version",
    "strain_overlay_effect_sha256",
)

SEARCH_AUDIT_FIELDS = {
    "searched_at",
    "databases",
    "queries",
    "inclusion_criteria",
    "exclusion_criteria",
    "screened_sources",
    "excluded_sources",
    "direct_evidence_found",
    "direct_evidence_absence_note",
}

SOURCE_RECORD_FIELDS = {
    "source_id",
    "title",
    "year",
    "species",
    "strain",
    "culture_conditions",
    "source_type",
    "evidence_type",
    "stance",
    "claim",
    "location",
    "evidence_tags",
    "genes",
    "reactions",
    "pathways",
    "methods",
    "result",
    "condition_match",
    "condition_mismatch_reason",
    "relevance",
    "confidence",
}

SOURCE_STANCES = {"supports", "contradicts", "contextual"}
CONDITION_MATCH_VALUES = {"exact", "partial", "mismatch", "unknown"}


def utc_now() -> str:
    """Return a stable UTC timestamp suitable for ledgers and dossiers."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: str | Path) -> str:
    """Hash a file without loading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_source_id(source: dict[str, Any]) -> str:
    """Return a deterministic citation ID from a DOI or permanent URL."""
    doi = str(source.get("doi", "") or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.startswith(prefix):
            doi = doi[len(prefix) :]
            break
    url = str(source.get("url", "") or "").strip().rstrip("/")
    if doi:
        locator = f"doi:{doi}"
    elif url:
        locator = f"url:{url}"
    else:
        raise ValueError("A literature source requires a DOI or permanent URL")
    return "SRC-" + hashlib.sha256(locator.encode("utf-8")).hexdigest()[:16]


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(value[key]) for key in sorted(value, key=str)}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float):
        return float(f"{value:.12g}")
    if value is None or isinstance(value, (str, int, bool)):
        return value
    return str(value)


def canonical_json(value: Any) -> str:
    """Serialize a value deterministically for IDs and fingerprints."""
    return json.dumps(
        _jsonable(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def _json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def _ledger_text(rows: Iterable[dict[str, Any]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=LEDGER_FIELDS, lineterminator="\r\n")
    writer.writeheader()
    for row in sorted(rows, key=lambda item: str(item.get("case_id", ""))):
        writer.writerow({field: row.get(field, "") for field in LEDGER_FIELDS})
    return buffer.getvalue()


def _stage_bytes(path: Path, payload: bytes, mode: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary_path, mode)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
    return temporary_path


def _atomic_write_bundle(payloads: dict[str | Path, str]) -> None:
    """Failure-safely replace a related ledger/dossier bundle.

    Every payload is staged and flushed before the first target changes.  If a
    later replacement fails, already-replaced targets are restored from their
    exact pre-transaction bytes.  This gives the multi-file workflow the same
    all-or-nothing behaviour as a single atomic ``os.replace`` under ordinary
    filesystem failures.
    """
    normalized = {Path(path): text for path, text in payloads.items()}
    originals: dict[Path, tuple[bool, bytes, int]] = {}
    staged: dict[Path, Path] = {}
    replaced: list[Path] = []
    try:
        for path, text in normalized.items():
            exists = path.exists()
            mode = stat.S_IMODE(path.stat().st_mode) if exists else 0o644
            originals[path] = (
                exists,
                path.read_bytes() if exists else b"",
                mode,
            )
            staged[path] = _stage_bytes(path, text.encode("utf-8"), mode)

        for path in normalized:
            os.replace(staged[path], path)
            replaced.append(path)

    except Exception as exc:
        rollback_errors: list[str] = []
        for path in reversed(replaced):
            existed, payload, mode = originals[path]
            try:
                if existed:
                    restore = _stage_bytes(path, payload, mode)
                    os.replace(restore, path)
                else:
                    path.unlink(missing_ok=True)
            except Exception as rollback_exc:  # pragma: no cover - catastrophic I/O
                rollback_errors.append(f"{path}: {rollback_exc}")
        if rollback_errors:
            raise RuntimeError(
                "Evidence transaction failed and rollback was incomplete: "
                + "; ".join(rollback_errors)
            ) from exc
        raise
    finally:
        for temporary_path in staged.values():
            temporary_path.unlink(missing_ok=True)


def stable_case_id(
    category: str, gene_ids: Iterable[str], reaction_ids: Iterable[str]
) -> str:
    """Build the required EGC-<12 hex> ID from biological grouping fields."""
    payload = {
        "category": str(category),
        "gene_ids": sorted({str(gene_id) for gene_id in gene_ids}),
        "reaction_ids": sorted({str(reaction_id) for reaction_id in reaction_ids}),
    }
    digest = hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()
    return f"EGC-{digest[:12]}"


def target_fingerprint(reaction_contexts: Iterable[dict[str, Any]]) -> str:
    """Fingerprint only model fields whose change invalidates simulation evidence."""
    targets = []
    for context in reaction_contexts:
        targets.append(
            {
                "reaction_id": context.get("reaction_id", ""),
                "stoichiometry": context.get("stoichiometry", {}),
                "lower_bound": context.get("lower_bound"),
                "upper_bound": context.get("upper_bound"),
                "gpr": context.get("gpr", ""),
            }
        )
    payload = sorted(targets, key=lambda item: str(item["reaction_id"]))
    digest = hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


def chemistry_fingerprint(reaction_contexts: Iterable[dict[str, Any]]) -> str:
    """Fingerprint target stoichiometry and every involved microspecies.

    Unlike :func:`target_fingerprint`, this fingerprint deliberately excludes
    bounds and GPRs.  It invalidates a chemistry audit whenever a target
    reaction's coefficients or any involved metabolite's formula, charge, or
    compartment changes.
    """
    targets: list[dict[str, Any]] = []
    for context in reaction_contexts:
        stoichiometry = context.get("stoichiometry", {})
        if not isinstance(stoichiometry, dict):
            stoichiometry = {}
        chemistry = context.get("metabolite_chemistry", {})
        if not isinstance(chemistry, dict):
            chemistry = {}
        metabolite_ids = sorted(str(item) for item in stoichiometry)
        targets.append(
            {
                "reaction_id": context.get("reaction_id", ""),
                "stoichiometry": {
                    metabolite_id: stoichiometry.get(metabolite_id)
                    for metabolite_id in metabolite_ids
                },
                "metabolite_chemistry": {
                    metabolite_id: {
                        "formula": (
                            chemistry.get(metabolite_id, {}).get("formula")
                            if isinstance(chemistry.get(metabolite_id), dict)
                            else None
                        ),
                        "charge": (
                            chemistry.get(metabolite_id, {}).get("charge")
                            if isinstance(chemistry.get(metabolite_id), dict)
                            else None
                        ),
                        "compartment": (
                            chemistry.get(metabolite_id, {}).get("compartment")
                            if isinstance(chemistry.get(metabolite_id), dict)
                            else None
                        ),
                    }
                    for metabolite_id in metabolite_ids
                },
            }
        )
    payload = sorted(targets, key=lambda item: str(item["reaction_id"]))
    digest = hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()
    return f"sha256:{digest}"


def assert_transition(current: str, target: str) -> None:
    if current not in CASE_STATUSES:
        raise ValueError(f"Unknown current essentiality state: {current!r}")
    if target not in CASE_STATUSES:
        raise ValueError(f"Unknown target essentiality state: {target!r}")
    if target not in ALLOWED_TRANSITIONS[current]:
        raise ValueError(
            f"Illegal essentiality state transition: {current} -> {target}"
        )


def transition_case_status(
    case_id: str,
    target_status: str,
    *,
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
    artifact_manifest: str | Path | None = None,
) -> dict[str, str]:
    """Advance a case without allowing the human-acceptance gate to be bypassed."""
    if target_status == "accepted":
        raise ValueError("Use record_human_decision for acceptance")
    rows = read_ledger(ledger_path)
    matches = [row for row in rows if row["case_id"] == case_id]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one ledger row for {case_id}; found {len(matches)}"
        )
    row = matches[0]
    assert_transition(row["status"], target_status)

    artifact: dict[str, Any] | None = None
    if target_status in {"implemented", "regression_passed"}:
        if artifact_manifest is None:
            raise ValueError(
                f"Cannot mark {target_status} without an --artifact-manifest"
            )
        manifest_path = Path(artifact_manifest)
        artifact = load_json_document(manifest_path)
        if not isinstance(artifact, dict):
            raise ValueError("Artifact manifest must be a JSON object")
        for field in (
            "case_id",
            *_record_provenance_fields(row),
            "target_fingerprint",
            "chemistry_fingerprint",
        ):
            expected = case_id if field == "case_id" else row.get(field, "")
            if str(artifact.get(field, "")) != str(expected):
                raise ValueError(
                    f"Artifact manifest {field} does not match the current case"
                )
        if target_status == "regression_passed" and artifact.get("regression_passed") is not True:
            raise ValueError("Regression manifest must set regression_passed=true")

    dossier_path = Path(evidence_dir) / f"{case_id}.json"
    if target_status in {"reviewed", "awaiting_human"}:
        if not dossier_path.exists():
            raise ValueError(f"Evidence dossier does not exist: {dossier_path}")
        dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
        for field in _record_provenance_fields(row):
            if dossier.get(field) != row.get(field):
                raise ValueError(
                    f"Evidence provenance {field} does not match the current ledger"
                )
        if dossier.get("target_fingerprint") != row.get("target_fingerprint"):
            raise ValueError(
                "Evidence target fingerprint does not match the current ledger"
            )
        if dossier.get("chemistry_fingerprint") != row.get(
            "chemistry_fingerprint"
        ):
            raise ValueError(
                "Evidence chemistry fingerprint does not match the current ledger"
            )
        adversarial = dossier.get("adversarial_review", {})
        if target_status == "reviewed" and (
            not isinstance(adversarial, dict)
            or str(adversarial.get("status", "")).strip().lower() != "complete"
        ):
            raise ValueError(
                "Cannot mark reviewed before the skeptic review is complete"
            )
        if target_status == "awaiting_human":
            require_valid_evidence_dossier(dossier, require_human_approval=False)

    now = utc_now()
    row.update({"status": target_status, "updated_at": now})
    payloads: dict[str | Path, str] = {}
    if dossier_path.exists():
        dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
        dossier["workflow_status"] = target_status
        dossier["workflow_updated_at"] = now
        if artifact is not None and artifact_manifest is not None:
            key = "implementation_manifest" if target_status == "implemented" else "regression_manifest"
            dossier[key] = {
                "path": str(Path(artifact_manifest).resolve()),
                "sha256": sha256_file(artifact_manifest),
                "recorded_at": now,
            }
        payloads[dossier_path] = _json_text(dossier)
    payloads[Path(ledger_path)] = _ledger_text(rows)
    _atomic_write_bundle(payloads)
    return row


def _nonempty(source: dict[str, Any], field: str) -> bool:
    value = source.get(field)
    return value is not None and str(value).strip() != ""


def _requires_simulation_context(record: dict[str, Any]) -> bool:
    """Schema-2.1 records are bound to an effective-model context."""
    return str(record.get("schema_version", "2.0")).strip() >= "2.1"


def _provenance_fields(record: dict[str, Any]) -> tuple[str, ...]:
    fields = (
        "model_sha256",
        "experimental_sha256",
        "media_sha256",
    )
    return fields + (SIMULATION_CONTEXT_FIELDS if _requires_simulation_context(record) else ())


def _nonempty_string_list(value: Any) -> bool:
    return (
        isinstance(value, list)
        and bool(value)
        and all(isinstance(item, str) and item.strip() for item in value)
    )


def validate_search_audit(search_audit: Any) -> list[str]:
    """Validate the reproducible search log stored with one reviewer result."""
    if not isinstance(search_audit, dict):
        return ["search_audit must be an object"]
    errors: list[str] = []
    missing = sorted(SEARCH_AUDIT_FIELDS - set(search_audit))
    if missing:
        errors.append(f"search_audit missing fields: {missing}")

    searched_at = str(search_audit.get("searched_at", "")).strip()
    try:
        parsed = datetime.fromisoformat(searched_at.replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            raise ValueError("timezone is required")
    except ValueError:
        errors.append("search_audit searched_at must be a timezone-aware ISO timestamp")

    for field in ("databases", "queries", "inclusion_criteria", "exclusion_criteria"):
        if not _nonempty_string_list(search_audit.get(field)):
            errors.append(f"search_audit {field} must be a non-empty string list")

    screened = search_audit.get("screened_sources", [])
    if not isinstance(screened, list) or any(not isinstance(item, dict) for item in screened):
        errors.append("search_audit screened_sources must be a list of objects")
        screened = []
    screened_ids: set[str] = set()
    for index, item in enumerate(screened):
        if not _nonempty(item, "source_id"):
            errors.append(f"search_audit screened_sources[{index}] has no source_id")
        else:
            screened_ids.add(str(item["source_id"]).strip())
        if str(item.get("disposition", "")).strip().lower() not in {
            "included",
            "excluded",
        }:
            errors.append(
                f"search_audit screened_sources[{index}] has invalid disposition"
            )

    excluded = search_audit.get("excluded_sources", [])
    if not isinstance(excluded, list) or any(not isinstance(item, dict) for item in excluded):
        errors.append("search_audit excluded_sources must be a list of objects")
        excluded = []
    for index, item in enumerate(excluded):
        if not _nonempty(item, "source_id") or not _nonempty(item, "reason"):
            errors.append(
                f"search_audit excluded_sources[{index}] requires source_id and reason"
            )
        elif str(item["source_id"]).strip() not in screened_ids:
            errors.append(
                f"search_audit excluded_sources[{index}] is absent from screened_sources"
            )

    direct_found = search_audit.get("direct_evidence_found")
    if not isinstance(direct_found, bool):
        errors.append("search_audit direct_evidence_found must be boolean")
    if direct_found is False and not _nonempty(
        search_audit, "direct_evidence_absence_note"
    ):
        errors.append(
            "search_audit must explain the absence of direct evidence"
        )
    return errors


def validate_source_records(sources: Any) -> list[str]:
    """Validate structured claim-level literature records and stable IDs."""
    if not isinstance(sources, list):
        return ["primary_sources must be a list"]
    errors: list[str] = []
    claim_keys: set[tuple[str, str, str]] = set()
    for index, source in enumerate(sources):
        if not isinstance(source, dict):
            errors.append(f"primary_sources[{index}] must be an object")
            continue
        missing = sorted(field for field in SOURCE_RECORD_FIELDS if field not in source)
        if not (_nonempty(source, "url") or _nonempty(source, "doi")):
            missing.append("url_or_doi")
        for field in SOURCE_RECORD_FIELDS - {
            "evidence_tags",
            "genes",
            "reactions",
            "pathways",
        }:
            if field in source and not _nonempty(source, field):
                missing.append(field)
        if missing:
            errors.append(
                f"primary_sources[{index}] missing: {sorted(set(missing))}"
            )

        for field in ("evidence_tags", "genes", "reactions", "pathways"):
            value = source.get(field)
            if not isinstance(value, list) or any(
                not isinstance(item, str) or not item.strip() for item in value
            ):
                errors.append(
                    f"primary_sources[{index}] {field} must be a string list"
                )

        stance = str(source.get("stance", "")).strip().lower()
        if stance not in SOURCE_STANCES:
            errors.append(f"primary_sources[{index}] has invalid stance {stance!r}")
        condition_match = str(source.get("condition_match", "")).strip().lower()
        if condition_match not in CONDITION_MATCH_VALUES:
            errors.append(
                f"primary_sources[{index}] has invalid condition_match {condition_match!r}"
            )

        if _nonempty(source, "url") or _nonempty(source, "doi"):
            expected_source_id = stable_source_id(source)
            if str(source.get("source_id", "")).strip() != expected_source_id:
                errors.append(
                    f"primary_sources[{index}] source_id must be {expected_source_id}"
                )
        source_id = str(source.get("source_id", "")).strip()
        claim_key = (
            source_id,
            str(source.get("claim", "")).strip(),
            str(source.get("location", "")).strip(),
        )
        if source_id and claim_key in claim_keys:
            errors.append(
                f"primary_sources[{index}] duplicates an existing source claim"
            )
        claim_keys.add(claim_key)
    return errors


def validate_contradiction_records(contradictions: Any) -> list[str]:
    if not isinstance(contradictions, list):
        return ["contradictions must be a list"]
    errors: list[str] = []
    for index, item in enumerate(contradictions):
        if not isinstance(item, dict):
            errors.append(f"contradictions[{index}] must be an object")
            continue
        missing = [
            field
            for field in ("source_id", "claim", "resolution_status", "rationale")
            if not _nonempty(item, field)
        ]
        if missing:
            errors.append(f"contradictions[{index}] missing: {missing}")
        if str(item.get("resolution_status", "")).strip().lower() not in {
            "resolved",
            "unresolved",
            "not_applicable",
        }:
            errors.append(
                f"contradictions[{index}] has invalid resolution_status"
            )
    return errors


def _is_yarrowia_lipolytica(value: Any) -> bool:
    normalized = str(value).strip().lower().replace("*", "")
    return normalized in {"yarrowia lipolytica", "y. lipolytica"}


_GPR_TOKEN = re.compile(r"[A-Za-z0-9_.:-]+")


def _gpr_gene_ids(rule: Any) -> set[str]:
    """Extract gene identifiers from a COBRA Boolean GPR expression."""
    return {
        token
        for token in _GPR_TOKEN.findall(str(rule or ""))
        if token.lower() not in {"and", "or"}
    }


def _target_reaction_contexts(dossier: dict[str, Any]) -> list[dict[str, Any]]:
    model_context = dossier.get("model_context", {})
    if not isinstance(model_context, dict):
        return []
    reactions = model_context.get("reactions", [])
    if not isinstance(reactions, list):
        return []
    return [context for context in reactions if isinstance(context, dict)]


def _target_reaction_ids(
    dossier: dict[str, Any], operation: dict[str, Any]
) -> set[str]:
    reaction_ids = {
        str(context.get("reaction_id", "")).strip()
        for context in _target_reaction_contexts(dossier)
        if str(context.get("reaction_id", "")).strip()
    }
    target_id = str(operation.get("target_id", "")).strip()
    if target_id:
        reaction_ids.add(target_id)
    return reaction_ids


def _required_identity_gene_ids(
    dossier: dict[str, Any], operation: dict[str, Any]
) -> tuple[set[str], list[str]]:
    required: set[str] = set()
    errors: list[str] = []
    for context in _target_reaction_contexts(dossier):
        declared = context.get("gpr_gene_ids")
        if isinstance(declared, list):
            required.update(str(gene_id).strip() for gene_id in declared if str(gene_id).strip())
        else:
            required.update(_gpr_gene_ids(context.get("gpr", "")))

    if str(operation.get("operation", "")).strip() == "set_gpr":
        proposed_gpr = str(operation.get("value", "")).strip()
        if not proposed_gpr:
            errors.append("proposed set_gpr operation has no value")
        else:
            required.update(_gpr_gene_ids(proposed_gpr))
    return required, errors


def _nested_activation_audit_errors(chemistry: dict[str, Any]) -> list[str]:
    """Reject any reference/component audit that is not explicitly ready."""
    errors: list[str] = []

    def visit(value: Any, path: str) -> None:
        if not isinstance(value, dict):
            return
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else str(key)
            key_text = str(key).strip().lower()
            if (
                "audit" in key_text
                and ("reference" in key_text or "component" in key_text)
                and (
                    not isinstance(child, dict)
                    or child.get("ready_for_activation") is not True
                )
            ):
                errors.append(f"{child_path} is not ready for activation")
            visit(child, child_path)

    visit(chemistry, "chemistry_review")
    return errors


def validate_evidence_dossier(
    dossier: dict[str, Any],
    *,
    require_supported_patch: bool = False,
    require_human_approval: bool = False,
) -> list[str]:
    """Return all evidence-policy violations found in a dossier.

    Database annotations and evidence from related yeasts may be useful context,
    but neither can satisfy the direct-evidence requirement for a patch.
    """
    errors: list[str] = []
    required_fields = {
        "case_id",
        "model_sha256",
        "experimental_sha256",
        "media_sha256",
        "target_fingerprint",
        "claim_under_review",
        "model_context",
        "primary_sources",
        "identity_crosschecks",
        "contradictions",
        "verdict",
        "proposed_operation",
        "confidence",
        "adversarial_review",
        "human_decision",
    }
    missing = sorted(required_fields - set(dossier))
    if missing:
        errors.append(f"missing dossier fields: {missing}")
    if _requires_simulation_context(dossier):
        context_missing = [
            field
            for field in SIMULATION_CONTEXT_FIELDS
            if field not in dossier
        ]
        if context_missing:
            errors.append(f"missing simulation-context fields: {context_missing}")
        elif bool(dossier.get("strain_overlay_enabled")):
            for field in (
                "simulation_context_fingerprint_version",
                "simulation_context_fingerprint",
                "strain_profile_id",
                "strain_profile_sha256",
                "strain_overlay_effect_fingerprint_version",
                "strain_overlay_effect_sha256",
            ):
                if not _nonempty(dossier, field):
                    errors.append(f"overlay dossier has empty {field}")

    verdict = str(dossier.get("verdict", ""))
    if verdict not in VERDICTS:
        errors.append(f"invalid verdict: {verdict!r}")
    if require_supported_patch and verdict != "supported_patch_candidate":
        errors.append("verdict is not supported_patch_candidate")

    sources = dossier.get("primary_sources", [])
    if not isinstance(sources, list):
        errors.append("primary_sources must be a list")
        sources = []
    errors.extend(validate_source_records(sources))
    direct_sources: list[dict[str, Any]] = []
    for source in sources:
        if not isinstance(source, dict):
            continue
        evidence_type = str(source.get("evidence_type", "")).strip().lower()
        source_type = str(source.get("source_type", "primary_research")).strip().lower()
        stance = str(source.get("stance", "supports")).strip().lower()
        if (
            _is_yarrowia_lipolytica(source.get("species"))
            and source_type in {"primary_research", "original_research"}
            and evidence_type in {"direct_experiment", "primary_experiment"}
            and stance == "supports"
        ):
            direct_sources.append(source)

    crosschecks = dossier.get("identity_crosschecks", [])
    if not isinstance(crosschecks, list):
        errors.append("identity_crosschecks must be a list")
        crosschecks = []
    valid_crosschecks = [
        item
        for item in crosschecks
        if isinstance(item, dict)
        and str(item.get("database", "")).strip().lower() in {"uniprot", "kegg"}
        and _nonempty(item, "gene_id")
        and (_nonempty(item, "url") or _nonempty(item, "accession"))
        and str(item.get("status", "match")).strip().lower() == "match"
    ]

    contradictions = dossier.get("contradictions", [])
    if not isinstance(contradictions, list):
        errors.append("contradictions must be a list")
        contradictions = []
    errors.extend(validate_contradiction_records(contradictions))
    unresolved = [
        item
        for item in contradictions
        if not isinstance(item, dict)
        or str(item.get("resolution_status", "unresolved")).strip().lower()
        not in {"resolved", "not_applicable"}
    ]

    adversarial = dossier.get("adversarial_review", {})
    skeptic_pass = (
        isinstance(adversarial, dict)
        and str(adversarial.get("status", "")).strip().lower() == "complete"
        and str(adversarial.get("verdict", "")).strip().lower() == "pass"
    )

    literature_review = dossier.get("literature_review")
    if isinstance(literature_review, dict):
        errors.extend(validate_search_audit(literature_review.get("search_audit")))
        if not _nonempty(literature_review, "reasoning"):
            errors.append("literature_review reasoning is missing")
    elif verdict == "supported_patch_candidate" or require_supported_patch:
        errors.append("supported patch dossier has no durable literature_review")

    if verdict == "supported_patch_candidate" or require_supported_patch:
        if not direct_sources:
            errors.append(
                "no supporting direct primary experiment in Yarrowia lipolytica"
            )
        if not valid_crosschecks:
            errors.append("no matching UniProt/KEGG identity cross-check")
        if unresolved:
            errors.append("unresolved direct contradiction is present")
        if not skeptic_pass:
            errors.append("adversarial review has not passed")

        operation = dossier.get("proposed_operation", {})
        if not isinstance(operation, dict) or not _nonempty(operation, "operation"):
            errors.append("proposed_operation must name an operation")
            operation_name = ""
        else:
            operation_name = str(operation["operation"]).strip()
            tags = {
                str(tag).strip().lower()
                for source in direct_sources
                for tag in source.get("evidence_tags", [])
            }
            required_tag_sets = {
                "set_gpr": {"direct_enzyme_function", "complex_membership"},
                "set_bounds": {"directionality", "localization", "transport"},
                "remove_reaction": {
                    "reaction_activity",
                    "directionality",
                    "localization",
                    "transport",
                },
            }
            permitted = required_tag_sets.get(operation_name)
            if permitted is not None and not (tags & permitted):
                errors.append(
                    f"direct sources lack operation-specific evidence for {operation_name}"
                )
            required_all_tags = {
                "partition_cpa_ura2": {
                    "complex_membership",
                    "metabolic_channeling",
                },
                "couple_trna_biomass": {
                    "biomass_composition",
                    "carrier_conservation",
                },
            }.get(operation_name, set())
            missing_tags = sorted(required_all_tags - tags)
            if missing_tags:
                errors.append(
                    f"direct sources lack required evidence tags for {operation_name}: "
                    f"{missing_tags}"
                )

        chemistry = dossier.get("chemistry_review", {})
        if not isinstance(chemistry, dict):
            errors.append("chemistry_review must be an object")
            chemistry = {}
        chemistry_status = str(chemistry.get("status", "")).strip().lower()
        if chemistry_status != "verified_balanced":
            errors.append(
                "reaction chemistry review is not verified_balanced for the proposed "
                f"operation: {chemistry_status or 'missing'}"
            )
        if str(chemistry.get("model_sha256", "")).strip() != str(
            dossier.get("model_sha256", "")
        ).strip():
            errors.append("chemistry review model SHA does not match dossier")
        dossier_chemistry_fingerprint = str(
            dossier.get("chemistry_fingerprint", "")
        ).strip()
        review_chemistry_fingerprint = str(
            chemistry.get("chemistry_fingerprint", "")
        ).strip()
        if not dossier_chemistry_fingerprint:
            errors.append("dossier chemistry_fingerprint is missing")
        reaction_contexts = _target_reaction_contexts(dossier)
        for context in reaction_contexts:
            reaction_id = str(context.get("reaction_id", "")).strip() or "<missing>"
            stoichiometry = context.get("stoichiometry")
            metabolite_chemistry = context.get("metabolite_chemistry")
            if not isinstance(stoichiometry, dict):
                errors.append(
                    f"reaction {reaction_id} stoichiometry must be an object"
                )
                continue
            if not isinstance(metabolite_chemistry, dict):
                errors.append(
                    f"reaction {reaction_id} metabolite_chemistry must be an object"
                )
                continue
            for metabolite_id in sorted(stoichiometry):
                entry = metabolite_chemistry.get(metabolite_id)
                if not isinstance(entry, dict):
                    errors.append(
                        f"reaction {reaction_id} lacks chemistry for {metabolite_id}"
                    )
                    continue
                if not str(entry.get("formula") or "").strip():
                    errors.append(
                        f"reaction {reaction_id} metabolite {metabolite_id} formula "
                        "is missing"
                    )
                if entry.get("charge") is None:
                    errors.append(
                        f"reaction {reaction_id} metabolite {metabolite_id} charge "
                        "is missing"
                    )
                if not str(entry.get("compartment") or "").strip():
                    errors.append(
                        f"reaction {reaction_id} metabolite {metabolite_id} "
                        "compartment is missing"
                    )
        recomputed_chemistry_fingerprint = chemistry_fingerprint(reaction_contexts)
        if (
            dossier_chemistry_fingerprint
            and dossier_chemistry_fingerprint != recomputed_chemistry_fingerprint
        ):
            errors.append(
                "dossier chemistry_fingerprint does not match model_context chemistry"
            )
        if review_chemistry_fingerprint != dossier_chemistry_fingerprint:
            errors.append("chemistry review fingerprint does not match dossier")
        if chemistry.get("ready_for_activation") is not True:
            errors.append("chemistry review is not ready for activation")

        target_reaction_ids = _target_reaction_ids(dossier, operation)
        audited_reaction_ids = chemistry.get("audited_reaction_ids", [])
        if not isinstance(audited_reaction_ids, list):
            errors.append("chemistry audited_reaction_ids must be a list")
            audited_reaction_id_set: set[str] = set()
        else:
            audited_reaction_id_set = {
                str(reaction_id).strip()
                for reaction_id in audited_reaction_ids
                if str(reaction_id).strip()
            }
        missing_audited_reactions = sorted(
            target_reaction_ids - audited_reaction_id_set
        )
        if missing_audited_reactions:
            errors.append(
                "chemistry audit does not cover target reactions: "
                f"{missing_audited_reactions}"
            )

        residuals = chemistry.get("residuals_by_reaction")
        if not isinstance(residuals, dict):
            errors.append("chemistry residuals_by_reaction must be an object")
        else:
            missing_residuals = sorted(target_reaction_ids - set(residuals))
            if missing_residuals:
                errors.append(
                    "chemistry residuals are missing target reactions: "
                    f"{missing_residuals}"
                )
            nonempty_residuals = sorted(
                str(reaction_id)
                for reaction_id, residual in residuals.items()
                if residual not in ({}, [], None, "")
            )
            if nonempty_residuals:
                errors.append(
                    "chemistry audit has non-empty reaction residuals: "
                    f"{nonempty_residuals}"
                )
        if not _nonempty(chemistry, "audit_path"):
            errors.append("chemistry audit_path is missing")
        if not _nonempty(chemistry, "audit_sha256"):
            errors.append("chemistry audit_sha256 is missing")
        errors.extend(_nested_activation_audit_errors(chemistry))

        identity = dossier.get("identity_review", {})
        if not isinstance(identity, dict):
            errors.append("identity_review must be an object")
            identity = {}
        identity_status = str(identity.get("status", "")).strip().lower()
        if identity_status != "verified":
            errors.append(
                "gene identity review is not verified: "
                f"{identity_status or 'missing'}"
            )
        if str(identity.get("model_sha256", "")).strip() != str(
            dossier.get("model_sha256", "")
        ).strip():
            errors.append("identity review model SHA does not match dossier")

        required_gene_ids, identity_errors = _required_identity_gene_ids(
            dossier, operation
        )
        errors.extend(identity_errors)
        reviewed_gene_ids = identity.get("reviewed_gene_ids")
        if not isinstance(reviewed_gene_ids, list):
            errors.append("identity reviewed_gene_ids must be a list")
            reviewed_gene_id_set: set[str] = set()
        else:
            reviewed_gene_id_set = {
                str(gene_id).strip()
                for gene_id in reviewed_gene_ids
                if str(gene_id).strip()
            }
        missing_reviewed = sorted(required_gene_ids - reviewed_gene_id_set)
        if missing_reviewed:
            errors.append(
                "identity review does not cover required genes: "
                f"{missing_reviewed}"
            )

        gene_entries = identity.get("genes")
        if not isinstance(gene_entries, list):
            errors.append("identity genes must be a list")
            gene_entries = []
        entries_by_id: dict[str, dict[str, Any]] = {}
        for index, entry in enumerate(gene_entries):
            if not isinstance(entry, dict):
                errors.append(f"identity genes[{index}] must be an object")
                continue
            gene_id = str(entry.get("gene_id", "")).strip()
            if not gene_id:
                errors.append(f"identity genes[{index}] has no gene_id")
                continue
            if gene_id in entries_by_id:
                errors.append(f"identity review has duplicate gene entry: {gene_id}")
                continue
            entries_by_id[gene_id] = entry
            for status_field in ("identity_status", "function_status"):
                status = str(entry.get(status_field, "")).strip().lower()
                if status != "verified":
                    errors.append(
                        f"identity gene {gene_id} {status_field} is not verified"
                    )
            if not _nonempty(entry, "functional_role"):
                errors.append(f"identity gene {gene_id} functional_role is missing")
            evidence_refs = entry.get("evidence_refs")
            if not isinstance(evidence_refs, list) or not any(
                str(reference).strip() for reference in evidence_refs
            ):
                errors.append(f"identity gene {gene_id} evidence_refs are missing")
        missing_gene_entries = sorted(required_gene_ids - set(entries_by_id))
        if missing_gene_entries:
            errors.append(
                "identity genes entries do not cover required genes: "
                f"{missing_gene_entries}"
            )

        crosschecked_gene_ids = {
            str(item.get("gene_id", "")).strip() for item in valid_crosschecks
        }
        missing_crosschecks = sorted(required_gene_ids - crosschecked_gene_ids)
        if missing_crosschecks:
            errors.append(
                "identity cross-checks do not cover required genes: "
                f"{missing_crosschecks}"
            )

    if require_human_approval:
        human = dossier.get("human_decision", {})
        if not isinstance(human, dict):
            errors.append("human_decision must be an object")
        else:
            if str(human.get("decision", "")).strip().lower() != "accepted":
                errors.append("human decision is not accepted")
            if str(human.get("approved_by", "")).strip() != "human_user":
                errors.append("approved_by must be human_user")
            if not _nonempty(human, "approved_at"):
                errors.append("human approval timestamp is missing")
    return errors


def require_valid_evidence_dossier(
    dossier: dict[str, Any],
    *,
    require_human_approval: bool = False,
) -> None:
    errors = validate_evidence_dossier(
        dossier,
        require_supported_patch=True,
        require_human_approval=require_human_approval,
    )
    if errors:
        raise ValueError("Evidence dossier rejected: " + "; ".join(errors))


def read_ledger(path: str | Path = DEFAULT_LEDGER) -> list[dict[str, str]]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = (
            set(LEDGER_FIELDS)
            - set(reader.fieldnames or [])
            - _LEGACY_OPTIONAL_LEDGER_FIELDS
        )
        if missing:
            raise ValueError(
                f"Essentiality case ledger missing columns: {sorted(missing)}"
            )
        return [
            {field: row.get(field, "") for field in LEDGER_FIELDS} for row in reader
        ]


def write_ledger(
    rows: Iterable[dict[str, Any]], path: str | Path = DEFAULT_LEDGER
) -> None:
    path = Path(path)
    _atomic_write_bundle({path: _ledger_text(rows)})


def dossier_skeleton(case: dict[str, Any]) -> dict[str, Any]:
    """Create the durable, unreviewed evidence record for a generated case."""
    dossier = {
        "schema_version": str(case.get("schema_version", "2.0")),
        "case_id": case["case_id"],
        "model_sha256": case["model_sha256"],
        "experimental_sha256": case["experimental_sha256"],
        "media_sha256": case["media_sha256"],
        "target_fingerprint": case["target_fingerprint"],
        "chemistry_fingerprint": case["chemistry_fingerprint"],
        "case_packet_sha256": case["case_packet_sha256"],
        "claim_under_review": case["claim_under_review"],
        "model_context": case["model_context"],
        "primary_sources": [],
        "identity_crosschecks": [],
        "contradictions": [],
        "verdict": "needs_more_evidence",
        "proposed_operation": {},
        "confidence": "none",
        "adversarial_review": {
            "status": "not_run",
            "verdict": "",
            "findings": [],
        },
        "human_decision": {
            "decision": "pending",
            "approved_by": "",
            "approved_at": "",
        },
    }
    if _requires_simulation_context(case):
        dossier.update({field: case.get(field) for field in SIMULATION_CONTEXT_FIELDS})
    return dossier


def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path.resolve())


_FRESH_CASE_REQUIRED_FIELDS = (
    "case_id",
    "model_sha256",
    "experimental_sha256",
    "media_sha256",
    "target_fingerprint",
    "chemistry_fingerprint",
)


def _case_reaction_contexts(case: dict[str, Any]) -> list[dict[str, Any]]:
    model_context = case.get("model_context")
    if not isinstance(model_context, dict):
        raise ValueError(f"{case.get('case_id', '<unknown>')}: model_context is invalid")
    reactions = model_context.get("reactions")
    if not isinstance(reactions, list):
        raise ValueError(
            f"{case.get('case_id', '<unknown>')}: model_context.reactions is invalid"
        )
    if not all(isinstance(context, dict) for context in reactions):
        raise ValueError(
            f"{case.get('case_id', '<unknown>')}: reaction context is invalid"
        )
    return reactions


def _contexts_by_reaction_id(
    contexts: Iterable[dict[str, Any]], *, case_id: str
) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for context in contexts:
        reaction_id = str(context.get("reaction_id", "")).strip()
        if not reaction_id:
            raise ValueError(f"{case_id}: reaction context has an empty reaction_id")
        if reaction_id in indexed:
            raise ValueError(f"{case_id}: duplicate reaction context {reaction_id}")
        indexed[reaction_id] = context
    return indexed


def _target_context_projection(context: dict[str, Any]) -> dict[str, Any]:
    return {
        "reaction_id": context.get("reaction_id", ""),
        "stoichiometry": context.get("stoichiometry", {}),
        "lower_bound": context.get("lower_bound"),
        "upper_bound": context.get("upper_bound"),
        "gpr": context.get("gpr", ""),
    }


def _fresh_case_context_index(case: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Validate a generated packet before any durable workflow file is touched."""
    case_id = str(case.get("case_id", "")).strip() or "<unknown>"
    missing = [field for field in _FRESH_CASE_REQUIRED_FIELDS if not _nonempty(case, field)]
    if missing:
        raise ValueError(f"{case_id}: fresh case has empty fields: {sorted(missing)}")
    if _requires_simulation_context(case):
        missing_context = [
            field
            for field in SIMULATION_CONTEXT_FIELDS
            if field not in case
            or (
                field not in {"strain_profile_id", "strain_profile_sha256", "strain_overlay_effect_fingerprint_version", "strain_overlay_effect_sha256"}
                and not _nonempty(case, field)
            )
        ]
        if missing_context:
            raise ValueError(
                f"{case_id}: fresh case has incomplete simulation context: {missing_context}"
            )
        if bool(case.get("strain_overlay_enabled")):
            for field in (
                "strain_profile_id",
                "strain_profile_sha256",
                "strain_overlay_effect_fingerprint_version",
                "strain_overlay_effect_sha256",
            ):
                if not _nonempty(case, field):
                    raise ValueError(f"{case_id}: overlay case has empty {field}")

    contexts = _case_reaction_contexts(case)
    indexed = _contexts_by_reaction_id(contexts, case_id=case_id)
    raw_reaction_ids = case.get("reaction_ids")
    if not isinstance(raw_reaction_ids, list):
        raise ValueError(f"{case_id}: reaction_ids must be a list")
    declared_reaction_ids = [
        str(reaction_id).strip() for reaction_id in raw_reaction_ids
    ]
    if (
        any(not reaction_id for reaction_id in declared_reaction_ids)
        or len(declared_reaction_ids) != len(set(declared_reaction_ids))
        or set(declared_reaction_ids) != set(indexed)
    ):
        raise ValueError(
            f"{case_id}: declared reaction_ids do not exactly match model_context"
        )

    for reaction_id, context in indexed.items():
        stoichiometry = context.get("stoichiometry")
        chemistry = context.get("metabolite_chemistry")
        if not isinstance(stoichiometry, dict) or not isinstance(chemistry, dict):
            raise ValueError(
                f"{case_id}: fresh reaction {reaction_id} lacks complete chemistry"
            )
        if set(map(str, stoichiometry)) != set(map(str, chemistry)):
            raise ValueError(
                f"{case_id}: fresh reaction {reaction_id} chemistry keys do not "
                "match stoichiometry"
            )
        for metabolite_id in stoichiometry:
            entry = chemistry.get(metabolite_id)
            if not isinstance(entry, dict) or not {
                "formula",
                "charge",
                "compartment",
            }.issubset(entry):
                raise ValueError(
                    f"{case_id}: fresh reaction {reaction_id} has incomplete "
                    f"chemistry for {metabolite_id}"
                )

    recomputed_target = target_fingerprint(contexts)
    if recomputed_target != str(case["target_fingerprint"]):
        raise ValueError(
            f"{case_id}: fresh target_fingerprint is internally inconsistent"
        )
    recomputed_chemistry = chemistry_fingerprint(contexts)
    if recomputed_chemistry != str(case["chemistry_fingerprint"]):
        raise ValueError(
            f"{case_id}: fresh chemistry_fingerprint is internally inconsistent"
        )
    return indexed


def _same_nonempty_provenance(
    source: dict[str, Any], case: dict[str, Any]
) -> bool:
    return all(
        _nonempty(source, field)
        and str(source.get(field)) == str(case.get(field))
        for field in ("case_id", *_provenance_fields(case), "target_fingerprint")
    )


def _migrate_same_sha_legacy_chemistry(
    existing_dossier: dict[str, Any],
    existing_row: dict[str, Any] | None,
    case: dict[str, Any],
    fresh_contexts: dict[str, dict[str, Any]],
    *,
    migrated_at: str,
) -> dict[str, Any] | None:
    """Backfill chemistry only for an unambiguous pre-chemistry schema-v2 case.

    Any ambiguity returns ``None`` so the caller follows the ordinary stale-case
    archive/reset path.  In particular, partially populated chemistry is never
    merged with a fresh packet.
    """
    if existing_row is None:
        return None
    if str(existing_dossier.get("schema_version", "2.0")).startswith("1"):
        return None
    if not _same_nonempty_provenance(existing_dossier, case):
        return None
    if not _same_nonempty_provenance(existing_row, case):
        return None
    if _nonempty(existing_dossier, "chemistry_fingerprint") or _nonempty(
        existing_row, "chemistry_fingerprint"
    ):
        return None

    try:
        legacy_context_list = _case_reaction_contexts(existing_dossier)
        legacy_contexts = _contexts_by_reaction_id(
            legacy_context_list, case_id=str(case["case_id"])
        )
    except ValueError:
        return None
    if set(legacy_contexts) != set(fresh_contexts):
        return None
    if any("metabolite_chemistry" in context for context in legacy_context_list):
        return None
    if target_fingerprint(legacy_context_list) != str(case["target_fingerprint"]):
        return None
    for reaction_id, legacy_context in legacy_contexts.items():
        if canonical_json(_target_context_projection(legacy_context)) != canonical_json(
            _target_context_projection(fresh_contexts[reaction_id])
        ):
            return None

    migrated = copy.deepcopy(existing_dossier)
    migrated_contexts = _case_reaction_contexts(migrated)
    for context in migrated_contexts:
        fresh = fresh_contexts[str(context["reaction_id"])]
        context["metabolite_chemistry"] = copy.deepcopy(
            fresh["metabolite_chemistry"]
        )
    recomputed_chemistry = chemistry_fingerprint(migrated_contexts)
    if recomputed_chemistry != str(case["chemistry_fingerprint"]):
        return None

    migrated["chemistry_fingerprint"] = recomputed_chemistry
    migrated["chemistry_fingerprint_migration"] = {
        "migration_type": "same_sha_legacy_backfill",
        "migrated_at": migrated_at,
        "case_id": str(case["case_id"]),
        "model_sha256": str(case["model_sha256"]),
        "experimental_sha256": str(case["experimental_sha256"]),
        "media_sha256": str(case["media_sha256"]),
        "target_fingerprint": str(case["target_fingerprint"]),
        "chemistry_fingerprint": recomputed_chemistry,
    }
    return migrated


def _archive_path_for_dossier(
    evidence_dir: Path, case_id: str, dossier: dict[str, Any]
) -> Path:
    old_target = str(dossier.get("target_fingerprint", "unknown") or "unknown")
    old_chemistry = str(
        dossier.get("chemistry_fingerprint", "unknown") or "unknown"
    )
    old_context = str(dossier.get("simulation_context_fingerprint", "baseline") or "baseline")
    old_packet = str(dossier.get("case_packet_sha256", "unknown") or "unknown")
    return evidence_dir / "archive" / (
        f"{case_id}.{old_target.replace(':', '-')}."
        f"{old_chemistry.replace(':', '-')}.{old_context}.{old_packet}.json"
    )


def merge_detected_cases(
    cases: Iterable[dict[str, Any]],
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
) -> list[dict[str, str]]:
    """Merge fresh cases with one failure-safe dossier/archive/ledger commit."""
    now = utc_now()
    fresh_cases = list(cases)
    for case in fresh_cases:
        supplied_packet_sha = str(case.get("case_packet_sha256", "")).strip()
        if not supplied_packet_sha:
            case["case_packet_sha256"] = case_packet_sha256(case)
        elif supplied_packet_sha != case_packet_sha256(case):
            raise ValueError(
                f"{case.get('case_id', '<unknown>')}: supplied case_packet_sha256 "
                "does not match packet contents"
            )
    fresh_contexts_by_case: dict[str, dict[str, dict[str, Any]]] = {}
    for case in fresh_cases:
        case_id = str(case.get("case_id", "")).strip()
        if case_id in fresh_contexts_by_case:
            raise ValueError(f"Duplicate fresh essentiality case: {case_id}")
        fresh_contexts_by_case[case_id] = _fresh_case_context_index(case)

    ledger_path = Path(ledger_path)
    evidence_dir = Path(evidence_dir)
    rows = {row["case_id"]: row for row in read_ledger(ledger_path)}
    payloads: dict[str | Path, str] = {}

    for case in fresh_cases:
        case_id = str(case["case_id"])
        fresh_contexts = fresh_contexts_by_case[case_id]
        dossier_path = evidence_dir / f"{case_id}.json"
        existing = rows.get(case_id)
        existing_dossier: dict[str, Any] | None = None
        dossier_state = "new"
        if dossier_path.exists():
            loaded = json.loads(dossier_path.read_text(encoding="utf-8"))
            if not isinstance(loaded, dict):
                raise ValueError(f"{case_id}: evidence dossier is not an object")
            existing_dossier = loaded

            current_contexts_match = False
            try:
                current_context_list = _case_reaction_contexts(existing_dossier)
                current_contexts = _contexts_by_reaction_id(
                    current_context_list, case_id=case_id
                )
                current_contexts_match = (
                    existing_dossier.get("target_fingerprint")
                    == case.get("target_fingerprint")
                    and existing_dossier.get("chemistry_fingerprint")
                    == case.get("chemistry_fingerprint")
                    and set(current_contexts) == set(fresh_contexts)
                    and target_fingerprint(current_context_list)
                    == case.get("target_fingerprint")
                    and chemistry_fingerprint(current_context_list)
                    == case.get("chemistry_fingerprint")
                    and all(
                        canonical_json(
                            _target_context_projection(current_contexts[reaction_id])
                        )
                        == canonical_json(
                            _target_context_projection(fresh_contexts[reaction_id])
                        )
                        for reaction_id in fresh_contexts
                    )
                )
            except ValueError:
                current_contexts_match = False

            if current_contexts_match:
                provenance_changed = any(
                    existing_dossier.get(field) != case.get(field)
                    for field in (*_provenance_fields(case), "case_packet_sha256")
                )
                if provenance_changed:
                    refreshed = dossier_skeleton(case)
                    for field in (
                        "primary_sources",
                        "identity_crosschecks",
                        "contradictions",
                    ):
                        if isinstance(existing_dossier.get(field), list):
                            refreshed[field] = existing_dossier[field]
                    if isinstance(existing_dossier.get("literature_review"), dict):
                        refreshed["literature_review"] = existing_dossier[
                            "literature_review"
                        ]
                    refreshed["reused_literature_provenance"] = {
                        "previous_model_sha256": existing_dossier.get(
                            "model_sha256", ""
                        ),
                        "previous_experimental_sha256": existing_dossier.get(
                            "experimental_sha256", ""
                        ),
                        "previous_media_sha256": existing_dossier.get("media_sha256", ""),
                        "previous_simulation_context_fingerprint": existing_dossier.get(
                            "simulation_context_fingerprint", ""
                        ),
                        "refreshed_at": now,
                    }
                    payloads[dossier_path] = _json_text(refreshed)
                    dossier_state = "provenance_refreshed"
                else:
                    dossier_state = "current"
            else:
                migrated = _migrate_same_sha_legacy_chemistry(
                    existing_dossier,
                    existing,
                    case,
                    fresh_contexts,
                    migrated_at=now,
                )
                if migrated is not None:
                    payloads[dossier_path] = _json_text(migrated)
                    dossier_state = "legacy_chemistry_migrated"
                else:
                    archive_path = _archive_path_for_dossier(
                        evidence_dir, case_id, existing_dossier
                    )
                    if not archive_path.exists():
                        payloads[archive_path] = _json_text(existing_dossier)
                    payloads[dossier_path] = _json_text(dossier_skeleton(case))
                    dossier_state = "stale_reset"
        else:
            payloads[dossier_path] = _json_text(dossier_skeleton(case))

        if existing and dossier_state == "legacy_chemistry_migrated":
            existing.update(
                {
                    "category": str(case["category"]),
                    "gene_ids": ";".join(case["gene_ids"]),
                    "reaction_ids": ";".join(case["reaction_ids"]),
                    "chemistry_fingerprint": str(case["chemistry_fingerprint"]),
                    "evidence_path": _display_path(dossier_path),
                    "updated_at": now,
                }
            )
            continue

        hashes_match = bool(
            existing
            and existing.get("model_sha256") == case.get("model_sha256")
            and existing.get("experimental_sha256") == case.get("experimental_sha256")
            and all(
                existing.get(field) == str(case.get(field, ""))
                for field in _provenance_fields(case)
            )
            and existing.get("target_fingerprint") == case.get("target_fingerprint")
            and existing.get("chemistry_fingerprint")
            == case.get("chemistry_fingerprint")
            and existing.get("case_packet_sha256") == case.get("case_packet_sha256")
        )
        if existing and hashes_match and dossier_state == "current":
            existing.update(
                {
                    "category": case["category"],
                    "gene_ids": ";".join(case["gene_ids"]),
                    "reaction_ids": ";".join(case["reaction_ids"]),
                    "evidence_path": _display_path(dossier_path),
                    "updated_at": now,
                }
            )
            continue

        previous_status = existing.get("status", "") if existing else "detected"
        stale_reasons: list[str] = []
        if existing:
            for field in (
                *_provenance_fields(case),
                "target_fingerprint",
                "chemistry_fingerprint",
                "case_packet_sha256",
            ):
                if existing.get(field) != str(case.get(field, "")):
                    stale_reasons.append(field)
            if dossier_state == "stale_reset" and not stale_reasons:
                stale_reasons.append("evidence_dossier")
        rows[case_id] = {
            "case_id": case_id,
            "status": "queued",
            "category": str(case["category"]),
            "gene_ids": ";".join(case["gene_ids"]),
            "reaction_ids": ";".join(case["reaction_ids"]),
            "model_sha256": str(case["model_sha256"]),
            "experimental_sha256": str(case["experimental_sha256"]),
            "media_sha256": str(case["media_sha256"]),
            "evidence_schema_version": str(case.get("schema_version", "2.0")),
            **{field: str(case.get(field, "") or "") for field in SIMULATION_CONTEXT_FIELDS},
            "case_packet_sha256": str(case["case_packet_sha256"]),
            "target_fingerprint": str(case["target_fingerprint"]),
            "chemistry_fingerprint": str(case["chemistry_fingerprint"]),
            "evidence_path": _display_path(dossier_path),
            "detected_at": existing.get("detected_at", now) if existing else now,
            "updated_at": now,
            "previous_status": previous_status,
            "stale_reason": ";".join(stale_reasons),
            "human_decision": "",
            "approved_by": "",
            "approved_at": "",
        }
    merged = list(rows.values())
    payloads[ledger_path] = _ledger_text(merged)
    _atomic_write_bundle(payloads)
    return merged


_CASE_PACKET_PROVENANCE_FIELDS = (
    "model_sha256",
    "experimental_sha256",
    "media_sha256",
    "target_fingerprint",
    "chemistry_fingerprint",
)

_CASE_PACKET_BINDING_FIELDS = (*_CASE_PACKET_PROVENANCE_FIELDS, "case_packet_sha256")


def _packet_provenance_fields(packet: dict[str, Any]) -> tuple[str, ...]:
    return _CASE_PACKET_PROVENANCE_FIELDS + (
        SIMULATION_CONTEXT_FIELDS if _requires_simulation_context(packet) else ()
    )


def _packet_binding_fields(packet: dict[str, Any]) -> tuple[str, ...]:
    return (*_packet_provenance_fields(packet), "case_packet_sha256")


def _record_provenance_fields(record: dict[str, Any]) -> tuple[str, ...]:
    return _CASE_PACKET_PROVENANCE_FIELDS + (
        SIMULATION_CONTEXT_FIELDS
        if str(record.get("evidence_schema_version", record.get("schema_version", "2.0"))) >= "2.1"
        else ()
    )


def case_packet_sha256(packet: dict[str, Any]) -> str:
    """Return the deterministic digest of a packet without its own digest."""

    payload = {
        key: value for key, value in packet.items() if key != "case_packet_sha256"
    }
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def reviewer_result_sha256(result: dict[str, Any]) -> str:
    """Return the deterministic digest of a reviewer result without its own digest."""

    payload = {
        key: value for key, value in result.items() if key != "reviewer_result_sha256"
    }
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def case_packet_batch_sha256(case_packets: Iterable[dict[str, Any]]) -> str:
    """Fingerprint the exact ordered packet set supplied to a skeptic."""

    payload = [
        {
            "case_id": str(packet.get("case_id", "")),
            "case_packet_sha256": str(packet.get("case_packet_sha256", "")),
        }
        for packet in case_packets
    ]
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()

_REVIEWER_RESULT_FIELDS = {
    "case_id",
    "claim_under_review",
    *_CASE_PACKET_BINDING_FIELDS,
    "search_audit",
    "primary_sources",
    "identity_crosschecks",
    "contradictions",
    "verdict",
    "proposed_operation",
    "confidence",
    "reasoning",
    "unresolved_questions",
    "reviewer_result_sha256",
}

_SKEPTIC_CASE_FIELDS = {
    "case_id",
    "case_packet_sha256",
    "reviewer_result_sha256",
    "status",
    "verdict",
    "findings",
    "unresolved_contradictions",
    "corrected_verdict",
    "confidence",
}


def load_json_document(path: str | Path) -> Any:
    path = Path(path)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValueError(f"Could not read JSON document {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON document {path}: {exc}") from exc


def validate_case_batch(
    case_packets: Iterable[dict[str, Any]],
    *,
    expected_batch_size: int | None = 3,
) -> list[dict[str, Any]]:
    """Validate the exact deterministic packets selected for one research batch."""
    packets = list(case_packets)
    if expected_batch_size is not None and len(packets) != expected_batch_size:
        raise ValueError(
            f"Expected exactly {expected_batch_size} selected cases; found {len(packets)}"
        )
    if not packets:
        raise ValueError("Research batch is empty")

    case_ids: list[str] = []
    for index, packet in enumerate(packets):
        if not isinstance(packet, dict):
            raise ValueError(f"Case packet {index} must be an object")
        required = {
            "case_id",
            "claim_under_review",
            *_packet_binding_fields(packet),
        }
        missing = sorted(required - set(packet))
        if missing:
            raise ValueError(f"Case packet {index} missing fields: {missing}")
        case_id = str(packet["case_id"]).strip()
        if not case_id:
            raise ValueError(f"Case packet {index} has an empty case_id")
        case_ids.append(case_id)
        for field in _packet_binding_fields(packet):
            if not str(packet.get(field, "")).strip():
                raise ValueError(f"Case packet {case_id} has an empty {field}")
        if str(packet["case_packet_sha256"]) != case_packet_sha256(packet):
            raise ValueError(f"Case packet {case_id} SHA does not match its contents")
    duplicates = sorted(
        {case_id for case_id in case_ids if case_ids.count(case_id) > 1}
    )
    if duplicates:
        raise ValueError(f"Duplicate selected case IDs: {duplicates}")
    return packets


def _assert_packet_matches_row(
    packet: dict[str, Any], row: dict[str, str]
) -> None:
    case_id = str(packet["case_id"])
    for field in _packet_provenance_fields(packet):
        if str(packet.get(field, "")) != row.get(field, ""):
            raise ValueError(
                f"Selected packet {case_id} {field} does not match the durable ledger"
            )


def _load_current_dossier(
    case_id: str,
    row: dict[str, str],
    evidence_dir: str | Path,
) -> tuple[Path, dict[str, Any]]:
    dossier_path = Path(evidence_dir) / f"{case_id}.json"
    if not dossier_path.exists():
        raise ValueError(f"Evidence dossier does not exist: {dossier_path}")
    dossier = load_json_document(dossier_path)
    if not isinstance(dossier, dict):
        raise ValueError(f"Evidence dossier must be an object: {dossier_path}")
    if dossier.get("case_id") != case_id:
        raise ValueError(f"Evidence dossier case_id does not match {case_id}")
    for field in _record_provenance_fields(row):
        if dossier.get(field) != row.get(field):
            raise ValueError(
                f"Evidence dossier {case_id} {field} does not match the durable ledger"
            )
    return dossier_path, dossier


def transition_selected_batch_to_researching(
    case_packets: Iterable[dict[str, Any]],
    *,
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
    expected_batch_size: int | None = 3,
) -> list[dict[str, str]]:
    """Atomically move only the selected, current queued cases to researching."""
    packets = validate_case_batch(
        case_packets, expected_batch_size=expected_batch_size
    )
    rows = read_ledger(ledger_path)
    rows_by_id = {row["case_id"]: row for row in rows}
    if len(rows_by_id) != len(rows):
        raise ValueError("Durable ledger contains duplicate case IDs")

    dossier_records: list[tuple[Path, dict[str, Any]]] = []
    selected_rows: list[dict[str, str]] = []
    for packet in packets:
        case_id = str(packet["case_id"])
        row = rows_by_id.get(case_id)
        if row is None:
            raise ValueError(f"Selected case is absent from the durable ledger: {case_id}")
        if row["status"] != "queued":
            raise ValueError(
                f"Selected case {case_id} must be queued, not {row['status']}"
            )
        _assert_packet_matches_row(packet, row)
        dossier_record = _load_current_dossier(case_id, row, evidence_dir)
        if dossier_record[1].get("claim_under_review") != packet.get(
            "claim_under_review"
        ):
            raise ValueError(
                f"Evidence dossier {case_id} claim does not match the selected packet"
            )
        if dossier_record[1].get("case_packet_sha256") != packet.get(
            "case_packet_sha256"
        ):
            raise ValueError(
                f"Evidence dossier {case_id} case-packet SHA does not match the selected packet"
            )
        dossier_records.append(dossier_record)
        selected_rows.append(row)

    now = utc_now()
    payloads: dict[str | Path, str] = {}
    for row, (dossier_path, dossier) in zip(selected_rows, dossier_records):
        assert_transition(row["status"], "researching")
        row.update({"status": "researching", "updated_at": now})
        dossier["workflow_status"] = "researching"
        dossier["workflow_updated_at"] = now
        payloads[dossier_path] = _json_text(dossier)
    payloads[Path(ledger_path)] = _ledger_text(rows)
    _atomic_write_bundle(payloads)
    return selected_rows


def validate_reviewer_results(
    reviewer_results: Iterable[dict[str, Any]],
    case_packets: Iterable[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Require exactly one well-formed literature-review result per packet."""
    packets = list(case_packets)
    packet_by_id = {str(packet["case_id"]): packet for packet in packets}
    results = list(reviewer_results)
    if any(not isinstance(result, dict) for result in results):
        raise ValueError("Every reviewer result must be an object")
    result_by_id = {str(result.get("case_id", "")): result for result in results}
    if len(result_by_id) != len(results):
        raise ValueError("Reviewer results contain duplicate case IDs")
    if set(result_by_id) != set(packet_by_id):
        raise ValueError(
            "Reviewer result case IDs do not exactly cover the selected batch: "
            f"expected={sorted(packet_by_id)}, found={sorted(result_by_id)}"
        )

    for case_id, result in result_by_id.items():
        packet = packet_by_id[case_id]
        required_result_fields = _REVIEWER_RESULT_FIELDS | set(
            _packet_provenance_fields(packet)
        )
        missing = sorted(required_result_fields - set(result))
        if missing:
            raise ValueError(f"Reviewer result {case_id} missing fields: {missing}")
        if result["claim_under_review"] != packet_by_id[case_id]["claim_under_review"]:
            raise ValueError(
                f"Reviewer result {case_id} claim does not match the selected packet"
            )
        for field in _packet_binding_fields(packet):
            if str(result.get(field, "")) != str(packet[field]):
                raise ValueError(
                    f"Reviewer result {case_id} {field} does not match the selected packet"
                )
        if str(result["reviewer_result_sha256"]) != reviewer_result_sha256(result):
            raise ValueError(
                f"Reviewer result {case_id} SHA does not match its contents"
            )
        if result["verdict"] not in VERDICTS:
            raise ValueError(
                f"Reviewer result {case_id} has invalid verdict {result['verdict']!r}"
            )
        search_errors = validate_search_audit(result.get("search_audit"))
        if search_errors:
            raise ValueError(
                f"Reviewer result {case_id} invalid search audit: "
                + "; ".join(search_errors)
            )
        for field in (
            "primary_sources",
            "identity_crosschecks",
            "contradictions",
            "unresolved_questions",
        ):
            if not isinstance(result[field], list):
                raise ValueError(f"Reviewer result {case_id} {field} must be a list")
        source_errors = validate_source_records(result["primary_sources"])
        if source_errors:
            raise ValueError(
                f"Reviewer result {case_id} invalid primary sources: "
                + "; ".join(source_errors)
            )
        audit_included_ids = {
            str(item.get("source_id", "")).strip()
            for item in result["search_audit"]["screened_sources"]
            if isinstance(item, dict)
            and str(item.get("disposition", "")).strip().lower() == "included"
        }
        source_ids = {
            str(item.get("source_id", "")).strip()
            for item in result["primary_sources"]
            if isinstance(item, dict)
        }
        missing_audit_sources = sorted(source_ids - audit_included_ids)
        if missing_audit_sources:
            raise ValueError(
                f"Reviewer result {case_id} sources absent from search audit: "
                f"{missing_audit_sources}"
            )
        has_direct_source = any(
            isinstance(item, dict)
            and _is_yarrowia_lipolytica(item.get("species"))
            and str(item.get("source_type", "")).strip().lower()
            in {"primary_research", "original_research"}
            and str(item.get("evidence_type", "")).strip().lower()
            in {"direct_experiment", "primary_experiment"}
            and str(item.get("stance", "")).strip().lower() == "supports"
            for item in result["primary_sources"]
        )
        if result["search_audit"]["direct_evidence_found"] is not has_direct_source:
            raise ValueError(
                f"Reviewer result {case_id} direct_evidence_found does not match "
                "the structured source records"
            )
        if not all(isinstance(item, dict) for item in result["identity_crosschecks"]):
            raise ValueError(f"Reviewer result {case_id} has a non-object cross-check")
        contradiction_errors = validate_contradiction_records(
            result["contradictions"]
        )
        if contradiction_errors:
            raise ValueError(
                f"Reviewer result {case_id} invalid contradictions: "
                + "; ".join(contradiction_errors)
            )
        contradicting_source_ids = {
            str(item.get("source_id", "")).strip()
            for item in result["primary_sources"]
            if isinstance(item, dict)
            and str(item.get("stance", "")).strip().lower() == "contradicts"
        }
        recorded_contradiction_ids = {
            str(item.get("source_id", "")).strip()
            for item in result["contradictions"]
            if isinstance(item, dict)
        }
        unrecorded_contradictions = sorted(
            contradicting_source_ids - recorded_contradiction_ids
        )
        if unrecorded_contradictions:
            raise ValueError(
                f"Reviewer result {case_id} contradicting sources are not recorded "
                f"in contradictions: {unrecorded_contradictions}"
            )
        if not isinstance(result["proposed_operation"], dict):
            raise ValueError(
                f"Reviewer result {case_id} proposed_operation must be an object"
            )
        if not str(result["confidence"]).strip():
            raise ValueError(f"Reviewer result {case_id} has empty confidence")
        if not str(result["reasoning"]).strip():
            raise ValueError(f"Reviewer result {case_id} has empty reasoning")
    return result_by_id


def validate_skeptic_batch_result(
    skeptic_batch_result: dict[str, Any],
    case_packets: Iterable[dict[str, Any]],
    reviewer_results: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Validate one skeptic document that exactly covers the selected batch."""
    if not isinstance(skeptic_batch_result, dict):
        raise ValueError("Skeptic batch result must be one object")
    if not {
        "batch_case_ids",
        "case_packet_batch_sha256",
        "reviewer_result_sha256_by_case",
        "results",
    }.issubset(skeptic_batch_result):
        raise ValueError(
            "Skeptic batch result requires batch_case_ids, case_packet_batch_sha256, "
            "reviewer_result_sha256_by_case and results"
        )
    expected_ids = [str(packet["case_id"]) for packet in case_packets]
    declared_ids = skeptic_batch_result["batch_case_ids"]
    results = skeptic_batch_result["results"]
    if not isinstance(declared_ids, list) or not all(
        isinstance(case_id, str) and case_id.strip() for case_id in declared_ids
    ):
        raise ValueError("Skeptic batch_case_ids must be a list of non-empty strings")
    if len(set(declared_ids)) != len(declared_ids):
        raise ValueError("Skeptic batch_case_ids contain duplicates")
    if set(declared_ids) != set(expected_ids):
        raise ValueError(
            "Skeptic document does not exactly cover the selected batch: "
            f"expected={sorted(expected_ids)}, found={sorted(declared_ids)}"
        )
    if str(skeptic_batch_result["case_packet_batch_sha256"]) != case_packet_batch_sha256(
        case_packets
    ):
        raise ValueError("Skeptic batch case-packet SHA does not match selected packets")
    reviewer_hashes = skeptic_batch_result["reviewer_result_sha256_by_case"]
    if not isinstance(reviewer_hashes, dict) or set(reviewer_hashes) != set(expected_ids):
        raise ValueError(
            "Skeptic reviewer-result SHA map does not exactly cover selected cases"
        )
    for case_id in expected_ids:
        if case_id not in reviewer_results:
            raise ValueError(f"Skeptic result references unknown reviewer case {case_id}")
        expected_hash = reviewer_results[case_id]["reviewer_result_sha256"]
        if str(reviewer_hashes[case_id]) != str(expected_hash):
            raise ValueError(
                f"Skeptic reviewer-result SHA for {case_id} does not match reviewer output"
            )
    if not isinstance(results, list) or any(
        not isinstance(result, dict) for result in results
    ):
        raise ValueError("Skeptic results must be a list of objects")
    result_by_id = {str(result.get("case_id", "")): result for result in results}
    if len(result_by_id) != len(results):
        raise ValueError("Skeptic results contain duplicate case IDs")
    if set(result_by_id) != set(expected_ids):
        raise ValueError("Skeptic results do not exactly cover batch_case_ids")

    for case_id, result in result_by_id.items():
        missing = sorted(_SKEPTIC_CASE_FIELDS - set(result))
        if missing:
            raise ValueError(f"Skeptic result {case_id} missing fields: {missing}")
        packet = next(packet for packet in case_packets if packet["case_id"] == case_id)
        if str(result["case_packet_sha256"]) != str(packet["case_packet_sha256"]):
            raise ValueError(
                f"Skeptic result {case_id} case-packet SHA does not match selected packet"
            )
        if str(result["reviewer_result_sha256"]) != str(
            reviewer_results[case_id]["reviewer_result_sha256"]
        ):
            raise ValueError(
                f"Skeptic result {case_id} reviewer-result SHA does not match reviewer output"
            )
        if str(result["status"]).strip().lower() != "complete":
            raise ValueError(f"Skeptic result {case_id} is not complete")
        verdict = str(result["verdict"]).strip().lower()
        if verdict not in {"pass", "block", "downgrade"}:
            raise ValueError(
                f"Skeptic result {case_id} has invalid verdict {verdict!r}"
            )
        if not isinstance(result["findings"], list) or not isinstance(
            result["unresolved_contradictions"], list
        ):
            raise ValueError(
                f"Skeptic result {case_id} findings and unresolved_contradictions "
                "must be lists"
            )
        corrected = str(result["corrected_verdict"]).strip()
        if corrected and corrected not in VERDICTS:
            raise ValueError(
                f"Skeptic result {case_id} has invalid corrected_verdict {corrected!r}"
            )
        if verdict in {"block", "downgrade"} and not corrected:
            raise ValueError(
                f"Skeptic result {case_id} must supply a corrected_verdict"
            )
        if verdict == "pass" and result["unresolved_contradictions"]:
            raise ValueError(
                f"Skeptic result {case_id} cannot pass with unresolved contradictions"
            )
        if not str(result["confidence"]).strip():
            raise ValueError(f"Skeptic result {case_id} has empty confidence")
    return result_by_id


def import_research_batch_results(
    case_packets: Iterable[dict[str, Any]],
    reviewer_results: Iterable[dict[str, Any]],
    skeptic_batch_result: dict[str, Any],
    *,
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
    expected_batch_size: int | None = 3,
) -> dict[str, Any]:
    """Atomically import one complete reviewer+skeptic batch without acceptance."""
    packets = validate_case_batch(
        case_packets, expected_batch_size=expected_batch_size
    )
    reviewer_by_id = validate_reviewer_results(reviewer_results, packets)
    skeptic_by_id = validate_skeptic_batch_result(
        skeptic_batch_result, packets, reviewer_by_id
    )
    rows = read_ledger(ledger_path)
    rows_by_id = {row["case_id"]: row for row in rows}
    if len(rows_by_id) != len(rows):
        raise ValueError("Durable ledger contains duplicate case IDs")

    prepared: list[tuple[dict[str, str], Path, dict[str, Any], str]] = []
    for packet in packets:
        case_id = str(packet["case_id"])
        row = rows_by_id.get(case_id)
        if row is None:
            raise ValueError(f"Selected case is absent from the durable ledger: {case_id}")
        if row["status"] != "researching":
            raise ValueError(
                f"Selected case {case_id} must be researching, not {row['status']}"
            )
        _assert_packet_matches_row(packet, row)
        dossier_path, dossier = _load_current_dossier(case_id, row, evidence_dir)
        if dossier.get("claim_under_review") != packet.get("claim_under_review"):
            raise ValueError(
                f"Evidence dossier {case_id} claim does not match the selected packet"
            )
        if dossier.get("case_packet_sha256") != packet.get("case_packet_sha256"):
            raise ValueError(
                f"Evidence dossier {case_id} case-packet SHA does not match the selected packet"
            )
        reviewer = reviewer_by_id[case_id]
        skeptic = skeptic_by_id[case_id]

        for field in (
            "claim_under_review",
            "primary_sources",
            "identity_crosschecks",
            "contradictions",
            "verdict",
            "proposed_operation",
            "confidence",
            "unresolved_questions",
        ):
            dossier[field] = reviewer[field]
        dossier["literature_review"] = reviewer
        dossier["reviewer_verdict"] = reviewer["verdict"]

        skeptic_verdict = str(skeptic["verdict"]).strip().lower()
        corrected_verdict = str(skeptic["corrected_verdict"]).strip()
        if corrected_verdict:
            dossier["verdict"] = corrected_verdict
        if skeptic_verdict in {"block", "downgrade"}:
            if dossier["verdict"] == "supported_patch_candidate":
                raise ValueError(
                    f"Skeptic {skeptic_verdict} for {case_id} cannot retain a "
                    "supported_patch_candidate verdict"
                )
        dossier["adversarial_review"] = {
            "status": "complete",
            "verdict": skeptic_verdict,
            "findings": skeptic["findings"],
            "unresolved_contradictions": skeptic["unresolved_contradictions"],
            "corrected_verdict": corrected_verdict,
            "confidence": skeptic["confidence"],
        }

        target_status = "reviewed"
        assert_transition(row["status"], target_status)
        if (
            dossier["verdict"] == "supported_patch_candidate"
            and skeptic_verdict == "pass"
        ):
            require_valid_evidence_dossier(dossier, require_human_approval=False)
            assert_transition("reviewed", "awaiting_human")
            target_status = "awaiting_human"
        else:
            errors = validate_evidence_dossier(dossier)
            if errors:
                raise ValueError(
                    f"Evidence dossier {case_id} rejected: " + "; ".join(errors)
                )
        prepared.append((row, dossier_path, dossier, target_status))

    now = utc_now()
    payloads: dict[str | Path, str] = {}
    awaiting_human: list[str] = []
    reviewed: list[str] = []
    for row, dossier_path, dossier, target_status in prepared:
        row.update({"status": target_status, "updated_at": now})
        dossier["workflow_status"] = target_status
        dossier["workflow_updated_at"] = now
        dossier["review_imported_at"] = now
        payloads[dossier_path] = _json_text(dossier)
        if target_status == "awaiting_human":
            awaiting_human.append(row["case_id"])
        else:
            reviewed.append(row["case_id"])
    payloads[Path(ledger_path)] = _ledger_text(rows)
    _atomic_write_bundle(payloads)
    return {
        "batch_case_ids": [str(packet["case_id"]) for packet in packets],
        "reviewed_case_ids": reviewed,
        "awaiting_human_case_ids": awaiting_human,
        "accepted_case_ids": [],
    }


def _load_current_case_for_review(
    case_id: str,
    *,
    ledger_path: str | Path,
    evidence_dir: str | Path,
) -> tuple[list[dict[str, str]], dict[str, str], Path, dict[str, Any]]:
    """Load one current case while rejecting stale or duplicate durable state."""

    rows = read_ledger(ledger_path)
    matches = [row for row in rows if row["case_id"] == case_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one ledger row for {case_id}; found {len(matches)}")
    row = matches[0]
    dossier_path, dossier = _load_current_dossier(case_id, row, evidence_dir)
    return rows, row, dossier_path, dossier


def _require_review_provenance(review: dict[str, Any], dossier: dict[str, Any]) -> None:
    for field in (*_provenance_fields(dossier), "case_id"):
        if str(review.get(field, "")) != str(dossier.get(field, "")):
            raise ValueError(f"Review {field} does not match the current evidence dossier")


def _required_case_identity_gene_ids(dossier: dict[str, Any]) -> set[str]:
    required = {
        str(gene_id).strip()
        for gene_id in dossier.get("model_context", {}).get("genes", [])
        if isinstance(gene_id, dict) and str(gene_id.get("gene_id", "")).strip()
    }
    required.update(
        str(gene_id).strip()
        for gene_id in dossier.get("gene_ids", [])
        if str(gene_id).strip()
    )
    for reaction in dossier.get("model_context", {}).get("reactions", []):
        if not isinstance(reaction, dict):
            continue
        required.update(
            str(gene_id).strip()
            for gene_id in reaction.get("gpr_gene_ids", [])
            if str(gene_id).strip()
        )
    return required


def import_identity_review(
    case_id: str,
    review: dict[str, Any],
    *,
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
) -> dict[str, str]:
    """Atomically attach a current-SHA identity review to one durable dossier."""

    if not isinstance(review, dict):
        raise ValueError("Identity review must be a JSON object")
    rows, row, dossier_path, dossier = _load_current_case_for_review(
        case_id, ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    _require_review_provenance(review, dossier)
    identity = review.get("identity_review", review)
    if not isinstance(identity, dict):
        raise ValueError("identity_review must be an object")
    if str(identity.get("model_sha256", "")) != str(dossier["model_sha256"]):
        raise ValueError("Identity review model SHA does not match the current dossier")
    status = str(identity.get("status", "")).strip()
    if not status:
        raise ValueError("Identity review status is required")
    if status == "verified":
        required_gene_ids = _required_case_identity_gene_ids(dossier)
        reviewed = {
            str(gene_id).strip()
            for gene_id in identity.get("reviewed_gene_ids", [])
            if str(gene_id).strip()
        }
        if not required_gene_ids <= reviewed:
            raise ValueError(
                "Verified identity review does not cover target/GPR genes: "
                f"{sorted(required_gene_ids - reviewed)}"
            )
        genes = identity.get("genes")
        if not isinstance(genes, list):
            raise ValueError("Verified identity review requires a genes list")
        by_id = {
            str(entry.get("gene_id", "")).strip(): entry
            for entry in genes
            if isinstance(entry, dict) and str(entry.get("gene_id", "")).strip()
        }
        for gene_id in required_gene_ids:
            entry = by_id.get(gene_id, {})
            if (
                str(entry.get("identity_status", "")).strip() != "verified"
                or str(entry.get("function_status", "")).strip() != "verified"
                or not _nonempty(entry, "functional_role")
                or not _nonempty_string_list(entry.get("evidence_refs"))
            ):
                raise ValueError(f"Verified identity review is incomplete for {gene_id}")
    dossier["identity_review"] = identity
    dossier["identity_review_imported_at"] = utc_now()
    row["updated_at"] = dossier["identity_review_imported_at"]
    _atomic_write_bundle({dossier_path: _json_text(dossier), Path(ledger_path): _ledger_text(rows)})
    return row


def _resolve_chemistry_audit_path(value: str, evidence_dir: str | Path) -> Path:
    raw = Path(value)
    candidate = raw if raw.is_absolute() else (REPO_ROOT / raw)
    candidate = candidate.resolve()
    root = Path(evidence_dir).resolve()
    if root not in (candidate, *candidate.parents):
        raise ValueError("Chemistry audit must be stored under the evidence directory")
    return candidate


def import_chemistry_review(
    case_id: str,
    review: dict[str, Any],
    *,
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
) -> dict[str, str]:
    """Atomically attach a hash-verified current chemistry review to one dossier."""

    if not isinstance(review, dict):
        raise ValueError("Chemistry review must be a JSON object")
    rows, row, dossier_path, dossier = _load_current_case_for_review(
        case_id, ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    _require_review_provenance(review, dossier)
    chemistry = review.get("chemistry_review", review)
    if not isinstance(chemistry, dict):
        raise ValueError("chemistry_review must be an object")
    required = {
        "status", "model_sha256", "chemistry_fingerprint", "audit_path", "audit_sha256"
    }
    missing = sorted(field for field in required if not _nonempty(chemistry, field))
    if missing:
        raise ValueError(f"Chemistry review missing fields: {missing}")
    if str(chemistry["model_sha256"]) != str(dossier["model_sha256"]):
        raise ValueError("Chemistry review model SHA does not match the current dossier")
    if str(chemistry["chemistry_fingerprint"]) != str(dossier["chemistry_fingerprint"]):
        raise ValueError("Chemistry review fingerprint does not match the current dossier")
    audit_path = _resolve_chemistry_audit_path(str(chemistry["audit_path"]), evidence_dir)
    if not audit_path.is_file() or sha256_file(audit_path) != str(chemistry["audit_sha256"]):
        raise ValueError("Chemistry audit path or SHA does not match the supplied review")
    if str(chemistry["status"]) == "verified_balanced":
        target_ids = _target_reaction_ids(dossier, dossier.get("proposed_operation", {}))
        audited = {
            str(reaction_id).strip()
            for reaction_id in chemistry.get("audited_reaction_ids", [])
            if str(reaction_id).strip()
        }
        residuals = chemistry.get("residuals_by_reaction")
        if (
            chemistry.get("ready_for_activation") is not True
            or not target_ids <= audited
            or not isinstance(residuals, dict)
            or any(residuals.get(reaction_id) not in ({}, [], None, "") for reaction_id in target_ids)
        ):
            raise ValueError("Verified chemistry review lacks balanced target audit evidence")
    dossier["chemistry_review"] = chemistry
    dossier["chemistry_review_imported_at"] = utc_now()
    row["updated_at"] = dossier["chemistry_review_imported_at"]
    _atomic_write_bundle({dossier_path: _json_text(dossier), Path(ledger_path): _ledger_text(rows)})
    return row


def verify_live_acceptance_inputs(
    case_id: str,
    row: dict[str, str],
    dossier: dict[str, Any],
    *,
    model_path: str | Path = DEFAULT_MODEL,
    experimental_path: str | Path = DEFAULT_EXPERIMENTAL,
    media_path: str | Path = DEFAULT_MEDIA,
    strain_profile_path: str | Path | None = None,
) -> dict[str, str]:
    """Recompute every simulation provenance gate immediately before acceptance."""
    from .essentiality_simulation_context import load_effective_simulation_context

    paths = {
        "model_sha256": Path(model_path),
        "experimental_sha256": Path(experimental_path),
        "media_sha256": Path(media_path),
    }
    for field, path in paths.items():
        if not path.exists():
            raise ValueError(f"Current acceptance input does not exist: {path}")
        current_sha = sha256_file(path)
        if row.get(field) != current_sha:
            raise ValueError(
                f"Current {field} does not match the durable ledger for {case_id}: "
                f"{row.get(field, '')} != {current_sha}"
            )
        if dossier.get(field) != current_sha:
            raise ValueError(
                f"Current {field} does not match the evidence dossier for {case_id}"
            )

    overlay_expected = bool(dossier.get("strain_overlay_enabled", False))
    if overlay_expected and strain_profile_path is None:
        raise ValueError("Overlay dossier acceptance requires --strain-profile")
    if not overlay_expected and strain_profile_path is not None:
        raise ValueError("Baseline dossier acceptance must not receive a strain profile")
    simulation = load_effective_simulation_context(
        model_path=model_path,
        media_path=media_path,
        strain_profile_path=strain_profile_path,
    )
    if _requires_simulation_context(dossier):
        for field, expected in simulation.provenance().items():
            row_value = row.get(field, "")
            dossier_value = dossier.get(field)
            if field == "strain_overlay_enabled":
                row_matches = str(row_value).strip().lower() == str(bool(expected)).lower()
            elif expected is None:
                row_matches = str(row_value).strip() == ""
            else:
                row_matches = str(row_value) == str(expected)
            if not row_matches:
                raise ValueError(
                    f"Current {field} does not match the durable ledger for {case_id}"
                )
            if dossier_value != expected:
                raise ValueError(
                    f"Current {field} does not match the evidence dossier for {case_id}"
                )

    if dossier.get("case_id") != case_id:
        raise ValueError(f"Evidence dossier case_id does not match {case_id}")
    contexts = dossier.get("model_context", {}).get("reactions", [])
    reaction_ids = sorted(
        {
            str(context.get("reaction_id", "")).strip()
            for context in contexts
            if isinstance(context, dict)
            and str(context.get("reaction_id", "")).strip()
        }
    )
    if not reaction_ids:
        reaction_ids = sorted(
            reaction_id.strip()
            for reaction_id in row.get("reaction_ids", "").split(";")
            if reaction_id.strip()
        )
    if not reaction_ids:
        raise ValueError(f"No live target reactions are recorded for {case_id}")

    model = simulation.model
    live_contexts: list[dict[str, Any]] = []
    for reaction_id in reaction_ids:
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError as exc:
            raise ValueError(
                f"Current model is missing target reaction {reaction_id} for {case_id}"
            ) from exc
        live_contexts.append(
            {
                "reaction_id": reaction.id,
                "stoichiometry": {
                    metabolite.id: float(coefficient)
                    for metabolite, coefficient in sorted(
                        reaction.metabolites.items(), key=lambda item: item[0].id
                    )
                },
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
                "gpr": reaction.gene_reaction_rule,
                "gpr_gene_ids": sorted(gene.id for gene in reaction.genes),
                "metabolite_chemistry": {
                    metabolite.id: {
                        "formula": metabolite.formula,
                        "charge": metabolite.charge,
                        "compartment": metabolite.compartment,
                    }
                    for metabolite in sorted(
                        reaction.metabolites, key=lambda item: item.id
                    )
                },
            }
        )
    live_fingerprint = target_fingerprint(live_contexts)
    if row.get("target_fingerprint") != live_fingerprint:
        raise ValueError(
            f"Current target fingerprint does not match the durable ledger for "
            f"{case_id}: {row.get('target_fingerprint', '')} != {live_fingerprint}"
        )
    if dossier.get("target_fingerprint") != live_fingerprint:
        raise ValueError(
            f"Current target fingerprint does not match the evidence dossier for {case_id}"
        )
    live_chemistry_fingerprint = chemistry_fingerprint(live_contexts)
    if row.get("chemistry_fingerprint") != live_chemistry_fingerprint:
        raise ValueError(
            f"Current chemistry fingerprint does not match the durable ledger for "
            f"{case_id}: {row.get('chemistry_fingerprint', '')} != "
            f"{live_chemistry_fingerprint}"
        )
    if dossier.get("chemistry_fingerprint") != live_chemistry_fingerprint:
        raise ValueError(
            "Current chemistry fingerprint does not match the evidence dossier for "
            f"{case_id}"
        )
    chemistry_review = dossier.get("chemistry_review", {})
    if not isinstance(chemistry_review, dict) or chemistry_review.get(
        "chemistry_fingerprint"
    ) != live_chemistry_fingerprint:
        raise ValueError(
            "Current chemistry fingerprint does not match the chemistry review for "
            f"{case_id}"
        )
    return {
        "model_sha256": row["model_sha256"],
        "experimental_sha256": row["experimental_sha256"],
        "media_sha256": row["media_sha256"],
        "target_fingerprint": live_fingerprint,
        "chemistry_fingerprint": live_chemistry_fingerprint,
        **simulation.provenance(),
    }


def record_human_decision(
    case_id: str,
    decision: str,
    *,
    ledger_path: str | Path = DEFAULT_LEDGER,
    evidence_dir: str | Path = DEFAULT_EVIDENCE_DIR,
    model_path: str | Path = DEFAULT_MODEL,
    experimental_path: str | Path = DEFAULT_EXPERIMENTAL,
    media_path: str | Path = DEFAULT_MEDIA,
    strain_profile_path: str | Path | None = None,
) -> dict[str, str]:
    """Record one explicit accept/reject/defer decision in ledger and dossier."""
    normalized = decision.strip().lower()
    target_status = {
        "accept": "accepted",
        "accepted": "accepted",
        "reject": "rejected",
        "rejected": "rejected",
        "defer": "needs_more_evidence",
        "deferred": "needs_more_evidence",
    }.get(normalized)
    if target_status is None:
        raise ValueError("Decision must be accept, reject, or defer")

    rows = read_ledger(ledger_path)
    matches = [row for row in rows if row["case_id"] == case_id]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one ledger row for {case_id}; found {len(matches)}"
        )
    row = matches[0]
    dossier_path = Path(evidence_dir) / f"{case_id}.json"
    if not dossier_path.exists():
        raise ValueError(f"Evidence dossier does not exist: {dossier_path}")
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    if dossier.get("case_id") != case_id:
        raise ValueError("Evidence case ID does not match the current ledger")
    if dossier.get("target_fingerprint") != row.get("target_fingerprint"):
        raise ValueError("Evidence target fingerprint does not match the current ledger")
    if dossier.get("chemistry_fingerprint") != row.get("chemistry_fingerprint"):
        raise ValueError(
            "Evidence chemistry fingerprint does not match the current ledger"
        )
    for field, label in (
        ("model_sha256", "model"),
        ("experimental_sha256", "experimental"),
        ("media_sha256", "media"),
    ):
        if dossier.get(field) != row.get(field):
            raise ValueError(f"Evidence {label} SHA does not match the current ledger")
    for field in _record_provenance_fields(row):
        if dossier.get(field) != row.get(field):
            raise ValueError(
                f"Evidence provenance {field} does not match the current ledger"
            )

    if target_status == "accepted":
        if row["status"] != "awaiting_human":
            raise ValueError(
                f"Cannot accept {case_id} from state {row['status']}; review must finish first"
            )
        require_valid_evidence_dossier(dossier, require_human_approval=False)
        live_provenance = verify_live_acceptance_inputs(
            case_id,
            row,
            dossier,
            model_path=model_path,
            experimental_path=experimental_path,
            media_path=media_path,
            strain_profile_path=strain_profile_path,
        )
    elif row["status"] not in {"reviewed", "awaiting_human", "needs_more_evidence"}:
        raise ValueError(f"Cannot decide {case_id} from state {row['status']}")

    if row["status"] != target_status:
        assert_transition(row["status"], target_status)
    now = utc_now()
    human_decision = {
        "decision": (
            "accepted"
            if target_status == "accepted"
            else "rejected"
            if target_status == "rejected"
            else "deferred"
        ),
        "approved_by": "human_user" if target_status == "accepted" else "",
        "approved_at": now if target_status == "accepted" else "",
        "decided_at": now,
    }
    dossier["human_decision"] = human_decision
    dossier["workflow_status"] = target_status
    dossier["workflow_updated_at"] = now
    if target_status == "accepted":
        dossier["acceptance_live_provenance"] = {
            **live_provenance,
            "verified_at": now,
        }
    row.update(
        {
            "status": target_status,
            "updated_at": now,
            "human_decision": human_decision["decision"],
            "approved_by": human_decision["approved_by"],
            "approved_at": human_decision["approved_at"],
        }
    )
    _atomic_write_bundle(
        {
            dossier_path: _json_text(dossier),
            Path(ledger_path): _ledger_text(rows),
        }
    )
    return row


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Advance an essentiality case through guarded review states"
    )
    parser.add_argument(
        "action",
        choices=(
            "accept",
            "reject",
            "defer",
            "queued",
            "researching",
            "reviewed",
            "awaiting_human",
            "implemented",
            "regression_passed",
            "start-batch",
            "import-batch-review",
            "import-chemistry-review",
            "import-identity-review",
        ),
    )
    parser.add_argument("case_id", nargs="?")
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--evidence-dir", type=Path, default=DEFAULT_EVIDENCE_DIR)
    parser.add_argument("--batch-file", type=Path)
    parser.add_argument("--reviewer-results", type=Path)
    parser.add_argument("--skeptic-result", type=Path)
    parser.add_argument("--review-file", type=Path)
    parser.add_argument("--artifact-manifest", type=Path)
    parser.add_argument("--expected-batch-size", type=int, default=3)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--experimental", type=Path, default=DEFAULT_EXPERIMENTAL)
    parser.add_argument("--media", type=Path, default=DEFAULT_MEDIA)
    parser.add_argument(
        "--strain-profile",
        type=Path,
        help="Exact runtime strain profile required to replay an overlay-bound case",
    )
    args = parser.parse_args()
    try:
        if args.action == "start-batch":
            if args.case_id is not None:
                parser.error("start-batch does not accept a case_id")
            if args.batch_file is None:
                parser.error("start-batch requires --batch-file")
            batch = load_json_document(args.batch_file)
            if not isinstance(batch, list):
                parser.error("--batch-file must contain a JSON list of case packets")
            updated = transition_selected_batch_to_researching(
                batch,
                ledger_path=args.ledger,
                evidence_dir=args.evidence_dir,
                expected_batch_size=args.expected_batch_size,
            )
        elif args.action == "import-batch-review":
            if args.case_id is not None:
                parser.error("import-batch-review does not accept a case_id")
            missing_paths = [
                flag
                for flag, path in (
                    ("--batch-file", args.batch_file),
                    ("--reviewer-results", args.reviewer_results),
                    ("--skeptic-result", args.skeptic_result),
                )
                if path is None
            ]
            if missing_paths:
                parser.error(
                    "import-batch-review requires " + ", ".join(missing_paths)
                )
            batch = load_json_document(args.batch_file)
            reviewer_document = load_json_document(args.reviewer_results)
            skeptic_document = load_json_document(args.skeptic_result)
            if not isinstance(batch, list):
                parser.error("--batch-file must contain a JSON list of case packets")
            if isinstance(reviewer_document, dict):
                reviewer_document = reviewer_document.get("results")
            if not isinstance(reviewer_document, list):
                parser.error(
                    "--reviewer-results must contain a JSON list or {results: [...]}"
                )
            if not isinstance(skeptic_document, dict):
                parser.error("--skeptic-result must contain one JSON object")
            updated = import_research_batch_results(
                batch,
                reviewer_document,
                skeptic_document,
                ledger_path=args.ledger,
                evidence_dir=args.evidence_dir,
                expected_batch_size=args.expected_batch_size,
            )
        elif args.action in {"import-chemistry-review", "import-identity-review"}:
            if args.case_id is None or args.review_file is None:
                parser.error(f"{args.action} requires case_id and --review-file")
            review_document = load_json_document(args.review_file)
            if args.action == "import-chemistry-review":
                updated = import_chemistry_review(
                    args.case_id,
                    review_document,
                    ledger_path=args.ledger,
                    evidence_dir=args.evidence_dir,
                )
            else:
                updated = import_identity_review(
                    args.case_id,
                    review_document,
                    ledger_path=args.ledger,
                    evidence_dir=args.evidence_dir,
                )
        else:
            if args.case_id is None:
                parser.error(f"{args.action} requires case_id")
            if args.action in {"accept", "reject", "defer"}:
                updated = record_human_decision(
                    args.case_id,
                    args.action,
                    ledger_path=args.ledger,
                    evidence_dir=args.evidence_dir,
                    model_path=args.model,
                    experimental_path=args.experimental,
                    media_path=args.media,
                    strain_profile_path=args.strain_profile,
                )
            else:
                updated = transition_case_status(
                    args.case_id,
                    args.action,
                    ledger_path=args.ledger,
                    evidence_dir=args.evidence_dir,
                    artifact_manifest=args.artifact_manifest,
                )
    except ValueError as exc:
        parser.error(str(exc))
    print(json.dumps(updated, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
