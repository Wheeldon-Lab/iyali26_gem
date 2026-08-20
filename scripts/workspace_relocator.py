#!/usr/bin/env python3
"""Relocate workspace artifacts without deleting or overwriting file content.

This utility is intentionally conservative:

* every regular file is hashed before and after an atomic rename;
* existing destinations are never overwritten;
* every moved file is recorded in an append-only CSV manifest;
* compatibility symlinks are created only for paths used by repository tools;
* no unlink, remove, rmtree, or replacement operation is used.

The repository and research root must be on the same filesystem so ``os.rename``
is an atomic move rather than a copy-and-delete operation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


MANIFEST_FIELDS = (
    "source_path",
    "destination_path",
    "size_bytes",
    "sha256",
    "duplicate_of",
    "moved_at_utc",
    "compatibility_symlink",
)


@dataclass(frozen=True)
class MoveSpec:
    source: Path
    destination: Path
    compatibility_symlink: bool = False


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def files_beneath(path: Path) -> list[Path]:
    if path.is_symlink():
        return []
    if path.is_file():
        return [path]
    return sorted(item for item in path.rglob("*") if item.is_file() and not item.is_symlink())


def directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in files_beneath(path))


def unique_destination(path: Path, source: Path) -> Path:
    if not path.exists() and not path.is_symlink():
        return path
    if path.is_dir() and not any(path.iterdir()):
        raise FileExistsError(
            f"Refusing to replace an existing empty directory: {path}. "
            "Choose a child destination instead."
        )
    fingerprint = (
        sha256_file(source)[:12]
        if source.is_file()
        else hashlib.sha256(str(source).encode("utf-8")).hexdigest()[:12]
    )
    candidate = path.with_name(f"{path.name}.source-{fingerprint}")
    counter = 2
    while candidate.exists() or candidate.is_symlink():
        candidate = path.with_name(f"{path.name}.source-{fingerprint}-{counter}")
        counter += 1
    return candidate


def assert_same_filesystem(repo_root: Path, research_root: Path) -> None:
    if repo_root.stat().st_dev != research_root.stat().st_dev:
        raise RuntimeError(
            "Zero-delete relocation requires repository and research root on the "
            "same filesystem so os.rename is atomic."
        )


def append_rows(manifest_handle, rows: Iterable[dict[str, object]]) -> None:
    writer = csv.DictWriter(manifest_handle, fieldnames=MANIFEST_FIELDS)
    for row in rows:
        writer.writerow(row)
    manifest_handle.flush()
    os.fsync(manifest_handle.fileno())


def move_one(
    spec: MoveSpec,
    *,
    repo_root: Path,
    manifest_handle,
    first_source_by_hash: dict[str, str],
) -> Path:
    source = spec.source
    if not source.exists() and not source.is_symlink():
        raise FileNotFoundError(f"Relocation source does not exist: {source}")
    if source.is_symlink():
        raise RuntimeError(f"Refusing to relocate an existing symlink: {source}")

    destination = unique_destination(spec.destination, source)
    destination.parent.mkdir(parents=True, exist_ok=True)

    source_files = files_beneath(source)
    source_records: list[tuple[Path, int, str]] = []
    for source_file in source_files:
        source_records.append(
            (
                source_file.relative_to(source) if source.is_dir() else Path(source.name),
                source_file.stat().st_size,
                sha256_file(source_file),
            )
        )

    os.rename(source, destination)

    rows: list[dict[str, object]] = []
    for relative_path, size_bytes, expected_hash in source_records:
        destination_file = (
            destination / relative_path if destination.is_dir() else destination
        )
        if not destination_file.is_file():
            raise RuntimeError(f"Moved file is missing: {destination_file}")
        actual_hash = sha256_file(destination_file)
        if actual_hash != expected_hash:
            raise RuntimeError(
                f"SHA-256 mismatch after moving {source}: "
                f"{expected_hash} != {actual_hash}"
            )
        duplicate_of = first_source_by_hash.get(actual_hash, "")
        first_source_by_hash.setdefault(actual_hash, str(source / relative_path))
        original_file = source / relative_path if source.is_dir() else source
        rows.append(
            {
                "source_path": str(original_file.relative_to(repo_root)),
                "destination_path": str(destination_file),
                "size_bytes": size_bytes,
                "sha256": actual_hash,
                "duplicate_of": duplicate_of,
                "moved_at_utc": utc_now(),
                "compatibility_symlink": str(spec.compatibility_symlink).lower(),
            }
        )
    append_rows(manifest_handle, rows)

    if spec.compatibility_symlink:
        os.symlink(destination, source, target_is_directory=destination.is_dir())
    return destination


def data_file_destination(research_root: Path, source: Path) -> Path:
    name = source.name
    diagnostic_tokens = (
        "audit",
        "blocked",
        "candidate",
        "classification",
        "critique",
        "diagnos",
        "factpack",
        "fullscan",
        "lacking",
        "memote",
        "offender",
        "result",
        "strategy",
        "summary",
        "transport",
        "unannotated",
        "vote",
        "landscape",
    )
    if source.suffix == ".py":
        return research_root / "artifacts" / "research_scripts" / "data" / name
    if (
        name.startswith("_")
        or source.suffix == ".log"
        or any(token in name.lower() for token in diagnostic_tokens)
    ):
        return research_root / "artifacts" / "diagnostics" / "data" / name
    return research_root / "state" / "curation" / "data" / name


def build_move_specs(repo_root: Path, research_root: Path) -> list[MoveSpec]:
    specs: list[MoveSpec] = []
    moved_directory_sources: list[Path] = []

    directory_moves = (
        ("results", "artifacts/results", True),
        ("explorer", "artifacts/explorer", False),
        ("external_models", "reference/external_models", True),
        ("experiments", "experiments/repository", True),
        ("data/metanetx", "reference/metanetx", True),
        ("data/kegg", "reference/kegg", True),
        ("data/ncbi", "reference/ncbi", True),
        ("data/expasy", "reference/expasy", True),
        ("data/cache", "cache/data", True),
        ("data/essentiality", "state/essentiality/repository", True),
        ("data/media", "state/media", True),
        ("data/yali1_yali0_map", "reference/locus_map", True),
        (
            "data/essential_gene_metabolic",
            "state/curation/essential_gene_metabolic",
            True,
        ),
        ("data/growth", "state/growth", True),
        ("docs/weekly_briefing", "artifacts/weekly_briefing", True),
        ("local", "archive/local", False),
        (".pytest_cache", "cache/tooling/pytest", False),
        (".ruff_cache", "cache/tooling/ruff", False),
    )
    for source_text, destination_text, compatibility in directory_moves:
        source = repo_root / source_text
        if source.exists() and not source.is_symlink():
            moved_directory_sources.append(source)
            specs.append(
                MoveSpec(
                    source=source,
                    destination=research_root / destination_text,
                    compatibility_symlink=compatibility,
                )
            )

    root_file_moves = (
        ("42003_2023_4996_MOESM10_ESM.xlsx", "raw/assays"),
        ("file.pdf", "raw/literature"),
        ("model_baseline.xml", "artifacts/legacy_models"),
        ("model_before.xml", "artifacts/legacy_models"),
        ("model_ensembl.xml", "artifacts/legacy_models"),
        ("model_with_esm_genes.xml", "artifacts/legacy_models"),
        ("change_gene_diff.html", "artifacts/reports"),
        ("change_metabolite.html", "artifacts/reports"),
        ("change_reaction_diff.html", "artifacts/reports"),
        ("with_esm_gene.html", "artifacts/reports"),
        ("diagnosis.txt", "artifacts/diagnostics/root"),
        ("yhong075@ucr.edu", "archive/local_config"),
        ("scripts/update_model.ipynb", "experiments/notebooks"),
    )
    for source_text, destination_dir in root_file_moves:
        source = repo_root / source_text
        if source.exists() and not source.is_symlink():
            specs.append(
                MoveSpec(
                    source=source,
                    destination=research_root / destination_dir / source.name,
                )
            )

    data_root = repo_root / "data"
    reserved_data_names = {
        "iyali26.xml",
        "metanetx",
        "kegg",
        "ncbi",
        "expasy",
        "cache",
        "essentiality",
        "media",
        "yali1_yali0_map",
        "essential_gene_metabolic",
        "growth",
        "__pycache__",
        ".DS_Store",
    }
    if data_root.is_dir():
        for source in sorted(data_root.iterdir()):
            if source.name in reserved_data_names or source.is_symlink():
                continue
            specs.append(
                MoveSpec(
                    source=source,
                    destination=data_file_destination(research_root, source),
                    compatibility_symlink=True,
                )
            )

    system_files = [
        repo_root / ".DS_Store",
        repo_root / ".Rhistory",
        repo_root / "data" / ".DS_Store",
        repo_root / "docs" / ".DS_Store",
        repo_root / "docs" / ".Rhistory",
        repo_root / "scripts" / ".DS_Store",
        repo_root / "scripts" / "gem_annotate" / ".DS_Store",
    ]
    for source in system_files:
        if source.exists() and not source.is_symlink():
            relative = source.relative_to(repo_root)
            specs.append(
                MoveSpec(
                    source=source,
                    destination=research_root / "archive" / "system" / relative,
                )
            )

    cache_directories: list[Path] = []
    for parent, directory_names, _ in os.walk(repo_root):
        parent_path = Path(parent)
        if parent_path == repo_root / ".git" or repo_root / ".git" in parent_path.parents:
            directory_names[:] = []
            continue
        if parent_path == repo_root / ".venv" or repo_root / ".venv" in parent_path.parents:
            directory_names[:] = []
            continue
        for directory_name in tuple(directory_names):
            if directory_name == "__pycache__":
                candidate = parent_path / directory_name
                nested_in_moved_directory = any(
                    moved_source == candidate or moved_source in candidate.parents
                    for moved_source in moved_directory_sources
                )
                if candidate != data_root / "__pycache__" and not nested_in_moved_directory:
                    cache_directories.append(candidate)
                directory_names.remove(directory_name)
    data_pycache = data_root / "__pycache__"
    if data_pycache.exists() and not data_pycache.is_symlink():
        cache_directories.append(data_pycache)
    for source in sorted(set(cache_directories)):
        if source.exists() and not source.is_symlink():
            relative = source.relative_to(repo_root)
            specs.append(
                MoveSpec(
                    source=source,
                    destination=research_root / "cache" / "python" / relative,
                )
            )

    return specs


def relocate(repo_root: Path, research_root: Path) -> Path:
    repo_root = repo_root.resolve()
    research_root = research_root.resolve()
    assert_same_filesystem(repo_root, research_root)

    manifest_path = research_root / "relocation_manifest.csv"
    partial_path = research_root / "relocation_manifest.partial.csv"
    if manifest_path.exists() or partial_path.exists():
        raise FileExistsError(
            "A relocation manifest already exists; refusing to overwrite it."
        )

    specs = build_move_specs(repo_root, research_root)
    first_source_by_hash: dict[str, str] = {}
    with partial_path.open("x", newline="", encoding="utf-8") as manifest_handle:
        csv.DictWriter(manifest_handle, fieldnames=MANIFEST_FIELDS).writeheader()
        manifest_handle.flush()
        os.fsync(manifest_handle.fileno())
        for index, spec in enumerate(specs, start=1):
            destination = move_one(
                spec,
                repo_root=repo_root,
                manifest_handle=manifest_handle,
                first_source_by_hash=first_source_by_hash,
            )
            print(f"[{index}/{len(specs)}] {spec.source} -> {destination}", flush=True)

    os.rename(partial_path, manifest_path)
    return manifest_path


def append_relocation(
    repo_root: Path,
    research_root: Path,
    source: Path,
    destination: Path,
    *,
    compatibility_symlink: bool = False,
) -> Path:
    """Append one late-discovered source to an existing relocation manifest."""

    repo_root = repo_root.resolve()
    research_root = research_root.resolve()
    source = source.resolve()
    destination = destination.resolve()
    assert_same_filesystem(repo_root, research_root)
    try:
        source.relative_to(repo_root)
        destination.relative_to(research_root)
    except ValueError as exc:
        raise ValueError(
            "Append relocation source must be inside the repository and "
            "destination must be inside the research workspace."
        ) from exc

    manifest_path = research_root / "relocation_manifest.csv"
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Relocation manifest not found: {manifest_path}"
        )

    first_source_by_hash: dict[str, str] = {}
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise ValueError(
                f"Unexpected relocation manifest columns: {reader.fieldnames}"
            )
        for row in reader:
            first_source_by_hash.setdefault(row["sha256"], row["source_path"])

    with manifest_path.open("a", newline="", encoding="utf-8") as manifest_handle:
        return move_one(
            MoveSpec(
                source=source,
                destination=destination,
                compatibility_symlink=compatibility_symlink,
            ),
            repo_root=repo_root,
            manifest_handle=manifest_handle,
            first_source_by_hash=first_source_by_hash,
        )


def verify_manifest(repo_root: Path, manifest_path: Path) -> int:
    failures: list[str] = []
    row_count = 0
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            row_count += 1
            destination = Path(row["destination_path"])
            if not destination.is_file():
                failures.append(f"missing destination: {destination}")
                continue
            actual_size = destination.stat().st_size
            expected_size = int(row["size_bytes"])
            if actual_size != expected_size:
                failures.append(
                    f"size mismatch: {destination}: {actual_size} != {expected_size}"
                )
                continue
            actual_hash = sha256_file(destination)
            if actual_hash != row["sha256"]:
                failures.append(
                    f"SHA-256 mismatch: {destination}: "
                    f"{actual_hash} != {row['sha256']}"
                )
            if row["compatibility_symlink"] == "true":
                source = repo_root / row["source_path"]
                if not source.exists():
                    failures.append(f"missing compatibility path: {source}")
    if failures:
        for failure in failures:
            print(failure, file=sys.stderr)
        return 1
    print(f"Verified {row_count} relocated files against {manifest_path}")
    return 0


def write_content_manifest(research_root: Path, output_path: Path) -> Path:
    research_root = research_root.resolve()
    output_path = output_path.resolve()
    partial_path = output_path.with_name(f"{output_path.name}.partial")
    if output_path.exists() or partial_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite an existing content manifest: {output_path}"
        )
    relocation_manifest = research_root / "relocation_manifest.csv"
    if not relocation_manifest.is_file():
        raise FileNotFoundError(
            f"Relocation manifest not found: {relocation_manifest}"
        )
    excluded = {output_path, partial_path}
    files = []
    for path in sorted(
        item
        for item in research_root.rglob("*")
        if item.is_file() and not item.is_symlink() and item.resolve() not in excluded
    ):
        files.append(
            {
                "path": str(path.relative_to(research_root)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    payload = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "research_root_name": research_root.name,
        "relocation_manifest_sha256": sha256_file(relocation_manifest),
        "files": files,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with partial_path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.rename(partial_path, output_path)
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    relocate_parser = subparsers.add_parser("relocate")
    relocate_parser.add_argument("--repo-root", type=Path, required=True)
    relocate_parser.add_argument("--research-root", type=Path, required=True)

    append_parser = subparsers.add_parser("append")
    append_parser.add_argument("--repo-root", type=Path, required=True)
    append_parser.add_argument("--research-root", type=Path, required=True)
    append_parser.add_argument("--source", type=Path, required=True)
    append_parser.add_argument("--destination", type=Path, required=True)
    append_parser.add_argument("--compatibility-symlink", action="store_true")

    plan_parser = subparsers.add_parser("plan")
    plan_parser.add_argument("--repo-root", type=Path, required=True)
    plan_parser.add_argument("--research-root", type=Path, required=True)

    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--repo-root", type=Path, required=True)
    verify_parser.add_argument("--manifest", type=Path, required=True)

    manifest_parser = subparsers.add_parser("manifest")
    manifest_parser.add_argument("--research-root", type=Path, required=True)
    manifest_parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "plan":
        specs = build_move_specs(args.repo_root.resolve(), args.research_root.resolve())
        total_size = 0
        for spec in specs:
            size = (
                spec.source.stat().st_size
                if spec.source.is_file()
                else directory_size(spec.source)
            )
            total_size += size
            print(
                f"{size}\t{spec.compatibility_symlink}\t"
                f"{spec.source}\t{spec.destination}"
            )
        print(f"Planned {len(specs)} moves totaling {total_size} bytes")
        return 0
    if args.command == "relocate":
        manifest = relocate(args.repo_root, args.research_root)
        print(f"Relocation manifest: {manifest}")
        return 0
    if args.command == "append":
        destination = append_relocation(
            args.repo_root,
            args.research_root,
            args.source,
            args.destination,
            compatibility_symlink=args.compatibility_symlink,
        )
        print(f"Appended relocation: {args.source} -> {destination}")
        return 0
    if args.command == "manifest":
        output = write_content_manifest(args.research_root, args.output)
        print(f"Research content manifest: {output}")
        return 0
    return verify_manifest(args.repo_root.resolve(), args.manifest.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
