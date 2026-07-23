"""Append-only registry for external iYali26 research runs.

The registry is a JSONL ledger in the external research workspace.  Existing
run manifests are treated as immutable evidence: backfill only adds ledger
records and never rewrites a result directory or manifest.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .essentiality_evidence import canonical_json


REGISTRY_RELATIVE_PATH = Path("artifacts") / "run_registry.jsonl"


class DuplicateRunError(RuntimeError):
    """Raised when a complete matching run already exists."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_json(payload: object) -> str:
    import hashlib

    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def registry_path(research_root: Path) -> Path:
    return Path(research_root).resolve() / REGISTRY_RELATIVE_PATH


def read_records(research_root: Path) -> list[dict[str, Any]]:
    path = registry_path(research_root)
    if not path.is_file():
        return []
    records: list[dict[str, Any]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid run registry JSON at {path}:{line_number}") from exc
        if not isinstance(record, dict) or "record_id" not in record:
            raise ValueError(f"Invalid run registry record at {path}:{line_number}")
        records.append(record)
    return records


def append_record(research_root: Path, record: Mapping[str, Any]) -> bool:
    """Append one immutable record, returning False when already registered."""

    path = registry_path(research_root)
    records = read_records(research_root)
    record_id = str(record["record_id"])
    if any(str(item["record_id"]) == record_id for item in records):
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(canonical_json(dict(record)) + "\n")
    return True


def build_run_key(
    workflow: str,
    *,
    inputs: Mapping[str, Any],
    code_sources: Mapping[str, Any],
    configuration: Mapping[str, Any],
) -> str:
    """Return the canonical fingerprint used to reject accidental reruns."""

    return sha256_json(
        {
            "schema_version": 1,
            "workflow": workflow,
            "inputs": dict(inputs),
            "code_sources": dict(code_sources),
            "configuration": dict(configuration),
        }
    )


def _successful_matching_record(
    records: Iterable[Mapping[str, Any]], workflow: str, run_key: str
) -> dict[str, Any] | None:
    for record in reversed(list(records)):
        if (
            record.get("record_type") == "run"
            and record.get("workflow") == workflow
            and record.get("run_key") == run_key
            and record.get("status") == "complete"
        ):
            return dict(record)
    return None


def guard_duplicate_run(
    research_root: Path,
    *,
    workflow: str,
    run_key: str,
    output_dir: Path,
    force_rerun: bool = False,
    reproduction_reason: str | None = None,
) -> dict[str, Any] | None:
    """Reject a matching complete result unless an explicit reproduction is justified."""

    previous = _successful_matching_record(read_records(research_root), workflow, run_key)
    if previous is None or Path(str(previous.get("output_dir", ""))).resolve() == Path(output_dir).resolve():
        return previous
    if not force_rerun:
        raise DuplicateRunError(
            "A successful matching run already exists at "
            f"{previous.get('output_dir')}; reuse it or pass --force-rerun with "
            "--reproduction-reason."
        )
    if not reproduction_reason or not reproduction_reason.strip():
        raise ValueError("--force-rerun requires a non-empty --reproduction-reason")
    return previous


def register_run(
    research_root: Path,
    *,
    workflow: str,
    run_key: str,
    output_dir: Path,
    inputs: Mapping[str, Any],
    code_sources: Mapping[str, Any],
    configuration: Mapping[str, Any],
    status: str,
    manifest_path: Path | None = None,
    previous: Mapping[str, Any] | None = None,
    reproduction_reason: str | None = None,
) -> bool:
    relationship = (
        {"type": "reproduction_of", "record_id": previous["record_id"], "reason": reproduction_reason}
        if previous is not None and Path(str(previous.get("output_dir", ""))).resolve() != Path(output_dir).resolve()
        else {"type": "original"}
    )
    record = {
        "schema_version": 1,
        "record_type": "run",
        "record_id": sha256_json(
            {"workflow": workflow, "run_key": run_key, "output_dir": str(Path(output_dir).resolve())}
        ),
        "recorded_at": utc_now(),
        "workflow": workflow,
        "run_key": run_key,
        "status": status,
        "output_dir": str(Path(output_dir).resolve()),
        "manifest_path": str(Path(manifest_path).resolve()) if manifest_path else None,
        "inputs": dict(inputs),
        "code_sources": dict(code_sources),
        "configuration": dict(configuration),
        "relationship": relationship,
    }
    return append_record(research_root, record)


def _essentiality_record(manifest_path: Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    inputs = dict(manifest.get("inputs", {}))
    configuration = {
        **dict(manifest.get("configuration", {})),
        **dict(manifest.get("cutoffs", {})),
        "solver": dict(manifest.get("solver", {})),
    }
    code_sources = {"git": dict(manifest.get("git", {}))}
    run_key = build_run_key(
        "essentiality",
        inputs=inputs,
        code_sources=code_sources,
        configuration=configuration,
    )
    return {
        "workflow": "essentiality",
        "run_key": run_key,
        "inputs": inputs,
        "code_sources": code_sources,
        "configuration": configuration,
        "status": "complete",
    }


def backfill(research_root: Path) -> dict[str, int]:
    """Append registry records for immutable historical runs and duplicate files."""

    root = Path(research_root).resolve()
    added_runs = 0
    added_duplicates = 0
    added_corrections = 0
    for manifest_path in sorted((root / "artifacts" / "results").rglob("run_manifest.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if "run_key" in manifest and "experiment_id" in manifest:
            workflow = "flow"
            details = {
                "workflow": workflow,
                "run_key": str(manifest["run_key"]),
                "inputs": dict(manifest.get("inputs", {})),
                "code_sources": dict(manifest.get("code_sources", {})),
                "configuration": dict(manifest.get("config", {})),
                "status": str(manifest.get("status", "complete")),
            }
        else:
            details = _essentiality_record(manifest_path, manifest)
        if details["status"] != "complete":
            continue
        if register_run(
            root,
            workflow=details["workflow"],
            run_key=details["run_key"],
            output_dir=manifest_path.parent,
            manifest_path=manifest_path,
            inputs=details["inputs"],
            code_sources=details["code_sources"],
            configuration=details["configuration"],
            status="complete",
        ):
            added_runs += 1
        for previous in read_records(root):
            if (
                previous.get("record_type") == "run"
                and previous.get("output_dir") == str(manifest_path.parent.resolve())
                and previous.get("workflow") != details["workflow"]
            ):
                correction = {
                    "schema_version": 1,
                    "record_type": "classification_correction",
                    "record_id": sha256_json(
                        {"supersedes": previous["record_id"], "workflow": details["workflow"]}
                    ),
                    "recorded_at": utc_now(),
                    "output_dir": str(manifest_path.parent.resolve()),
                    "relationship": {
                        "type": "supersedes_classification",
                        "record_id": previous["record_id"],
                        "workflow": details["workflow"],
                    },
                }
                if append_record(root, correction):
                    added_corrections += 1

    relocation = root / "relocation_manifest.csv"
    if relocation.is_file():
        import csv

        with relocation.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                duplicate_of = (row.get("duplicate_of") or "").strip()
                if not duplicate_of:
                    continue
                source = row.get("source_path", "")
                destination = row.get("destination_path", "")
                record = {
                    "schema_version": 1,
                    "record_type": "artifact_duplicate",
                    "record_id": sha256_json({"source": source, "destination": destination}),
                    "recorded_at": utc_now(),
                    "source_path": source,
                    "destination_path": destination,
                    "sha256": row.get("sha256"),
                    "relationship": {"type": "duplicate_of", "path": duplicate_of},
                }
                if append_record(root, record):
                    added_duplicates += 1
    return {
        "runs": added_runs,
        "duplicates": added_duplicates,
        "classification_corrections": added_corrections,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Manage append-only iYali26 run registry")
    subparsers = parser.add_subparsers(dest="command", required=True)
    backfill_parser = subparsers.add_parser("backfill")
    backfill_parser.add_argument("--research-root", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "backfill":
        print(json.dumps(backfill(args.research_root), sort_keys=True))
        return 0
    raise AssertionError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
