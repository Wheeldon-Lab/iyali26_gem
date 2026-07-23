"""Safe locus-tag resolution across the YALI0 and YALI1 assemblies.

Spelling normalisation is deliberately assembly preserving: case and the
optional underscore after ``YALI0``/``YALI1`` may vary, but a YALI0 tag is
never converted to a same-suffix YALI1 tag.  Cross-assembly aliases come only
from the two version-controlled mapping tables and are admitted only when the
combined mapping is one-to-one in both directions.
"""

from __future__ import annotations

import csv
import hashlib
import re
from collections import defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from .config import ESSENTIALITY_DIR, LOCUS_MAP_DIR, load_project_paths


_YALI_LOCUS_RE = re.compile(
    r"^(YALI[01])_?([A-Za-z]\d+[A-Za-z])$", re.IGNORECASE
)

_DEFAULT_IYLI21_MAP = LOCUS_MAP_DIR / "iyli21_genes_vs_S2.csv"
_DEFAULT_S2_MAP = LOCUS_MAP_DIR / "s2_metabolic_genes.csv"
_DEFAULT_IDENTITY_EXCLUSIONS = (
    ESSENTIALITY_DIR / "curated_locus_identity_exclusions.csv"
)


class AmbiguousLocusMappingError(ValueError):
    """Raised when a requested cross-assembly alias is not one-to-one."""


class LocusLookupConflictError(ValueError):
    """Raised when multiple available genes share one canonical spelling."""


def _parsed_yali_locus(raw: str) -> tuple[str, str] | None:
    match = _YALI_LOCUS_RE.fullmatch(str(raw).strip())
    if match is None:
        return None
    prefix = match.group(1).upper()
    raw_suffix = match.group(2)
    suffix = raw_suffix[0].upper() + raw_suffix[1:-1] + raw_suffix[-1].lower()
    return prefix, suffix


def canonical_locus_key(raw: str) -> str:
    """Return a lowercase lookup key without changing genome assembly.

    Non-YALI identifiers are treated as opaque case-insensitive identifiers.
    This also accommodates explicitly mapped legacy mitochondrial names such
    as ``YALIfMp03`` without inventing any relationship for them.
    """
    value = str(raw).strip()
    parsed = _parsed_yali_locus(value)
    if parsed is None:
        return value.casefold()
    prefix, suffix = parsed
    return f"{prefix}{suffix}".casefold()


def locus_spelling_variants(raw: str) -> set[str]:
    """Return exact-assembly spellings differing only by case/underscore."""
    value = str(raw).strip()
    parsed = _parsed_yali_locus(value)
    if parsed is None:
        return {value}
    prefix, suffix = parsed
    return {f"{prefix}{suffix}", f"{prefix}_{suffix}"}


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_pairs(path: Path, yali1_column: str, yali0_column: str) -> list[tuple[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Locus crosswalk source is missing: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = {yali1_column, yali0_column} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Locus crosswalk {path} is missing columns: {sorted(missing)}"
            )
        pairs: list[tuple[str, str]] = []
        for row in reader:
            yali1 = row[yali1_column].strip()
            yali0 = row[yali0_column].strip()
            if yali1 and yali0:
                pairs.append((canonical_locus_key(yali1), canonical_locus_key(yali0)))
    return pairs


def _read_identity_exclusions(path: Path) -> set[tuple[str, str]]:
    """Load evidence-backed cross-assembly pairs that must never be aliases."""

    if not path.exists():
        raise FileNotFoundError(f"Locus identity exclusion source is missing: {path}")
    required = {
        "schema_version",
        "status",
        "yali1",
        "yali0",
        "case_id",
        "evidence_path",
        "source_url",
        "reason",
    }
    exclusions: set[tuple[str, str]] = set()
    errors: list[str] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Locus identity exclusions {path} are missing columns: "
                f"{sorted(missing)}"
            )
        for line_number, row in enumerate(reader, start=2):
            prefix = f"{path.name} row {line_number}"
            yali1 = canonical_locus_key(row["yali1"])
            yali0 = canonical_locus_key(row["yali0"])
            if row["schema_version"].strip() != "1":
                errors.append(f"{prefix}: unsupported schema_version")
            if row["status"].strip() != "active":
                errors.append(f"{prefix}: status must be active")
            if not yali1.startswith("yali1") or not yali0.startswith("yali0"):
                errors.append(f"{prefix}: expected one YALI1 and one YALI0 locus")
            if not row["case_id"].strip().startswith("EGC-"):
                errors.append(f"{prefix}: evidence case_id is required")
            evidence_path = Path(row["evidence_path"].strip())
            if not evidence_path.is_absolute():
                evidence_path = load_project_paths().resolve_legacy_path(evidence_path)
            if not evidence_path.is_file():
                errors.append(f"{prefix}: evidence_path does not exist")
            if not row["source_url"].strip().startswith(("https://", "http://")):
                errors.append(f"{prefix}: source_url must be HTTP(S)")
            if not row["reason"].strip():
                errors.append(f"{prefix}: reason is required")
            pair = (yali1, yali0)
            if pair in exclusions:
                errors.append(f"{prefix}: duplicate excluded pair {pair!r}")
            exclusions.add(pair)
    if errors:
        raise ValueError("Invalid locus identity exclusions:\n- " + "\n- ".join(errors))
    return exclusions


@dataclass(frozen=True)
class LocusCrosswalk:
    """Conflict-aware, one-to-one aliases backed by curated CSV sources."""

    aliases: Mapping[str, str]
    ambiguous_components: Mapping[str, tuple[str, ...]]
    excluded_pairs: frozenset[tuple[str, str]]
    fingerprint: str
    source_sha256: Mapping[str, str]

    @classmethod
    def from_csvs(
        cls,
        iyli21_path: str | Path = _DEFAULT_IYLI21_MAP,
        s2_metabolic_path: str | Path = _DEFAULT_S2_MAP,
        identity_exclusions_path: str | Path | None = None,
    ) -> "LocusCrosswalk":
        iyli21 = Path(iyli21_path)
        s2_metabolic = Path(s2_metabolic_path)
        identity_exclusions = (
            Path(identity_exclusions_path)
            if identity_exclusions_path is not None
            else None
        )
        excluded_pairs = (
            _read_identity_exclusions(identity_exclusions)
            if identity_exclusions is not None
            else set()
        )
        sources = (
            (iyli21, "yali1_s2", "yali0"),
            (s2_metabolic, "yali1", "yali0"),
        )

        graph: dict[str, set[str]] = defaultdict(set)
        source_hashes: dict[str, str] = {}
        if identity_exclusions is not None:
            source_hashes[identity_exclusions.name] = _file_sha256(
                identity_exclusions
            )
        observed_pairs: set[tuple[str, str]] = set()
        for path, yali1_column, yali0_column in sources:
            source_hashes[path.name] = _file_sha256(path)
            for yali1, yali0 in _read_pairs(path, yali1_column, yali0_column):
                observed_pairs.add((yali1, yali0))
                if (yali1, yali0) in excluded_pairs:
                    continue
                if yali1 == yali0:
                    continue
                graph[yali1].add(yali0)
                graph[yali0].add(yali1)

        stale_exclusions = sorted(excluded_pairs - observed_pairs)
        if stale_exclusions:
            raise ValueError(
                "Locus identity exclusions do not match the current crosswalk "
                f"sources: {stale_exclusions}"
            )

        aliases: dict[str, str] = {}
        unsafe: set[str] = set()
        for locus, neighbours in graph.items():
            if len(neighbours) != 1:
                unsafe.add(locus)
                continue
            counterpart = next(iter(neighbours))
            if graph.get(counterpart) != {locus}:
                unsafe.add(locus)

        # Mark the whole connected component unsafe, not just the branch node.
        # Otherwise a leaf in a one-to-many mapping could still look valid.
        queue: deque[str] = deque(unsafe)
        while queue:
            locus = queue.popleft()
            for neighbour in graph.get(locus, ()):
                if neighbour not in unsafe:
                    unsafe.add(neighbour)
                    queue.append(neighbour)

        ambiguous_components: dict[str, tuple[str, ...]] = {}
        for locus in sorted(unsafe):
            component: set[str] = set()
            pending = [locus]
            while pending:
                current = pending.pop()
                if current in component:
                    continue
                component.add(current)
                pending.extend(graph.get(current, ()))
            ambiguous_components[locus] = tuple(sorted(component - {locus}))

        for locus, neighbours in graph.items():
            if locus in unsafe or len(neighbours) != 1:
                continue
            aliases[locus] = next(iter(neighbours))

        digest = hashlib.sha256()
        digest.update(b"safe-yali-crosswalk-v1\0")
        for name, sha256 in sorted(source_hashes.items()):
            digest.update(f"{name}\0{sha256}\0".encode())
        for locus, counterpart in sorted(aliases.items()):
            digest.update(f"alias\0{locus}\0{counterpart}\0".encode())
        for locus, component in sorted(ambiguous_components.items()):
            digest.update(f"ambiguous\0{locus}\0{','.join(component)}\0".encode())
        for yali1, yali0 in sorted(excluded_pairs):
            digest.update(f"excluded\0{yali1}\0{yali0}\0".encode())

        return cls(
            aliases=aliases,
            ambiguous_components=ambiguous_components,
            excluded_pairs=frozenset(excluded_pairs),
            fingerprint=f"sha256:{digest.hexdigest()}",
            source_sha256=source_hashes,
        )

    def counterpart(self, raw: str) -> str | None:
        """Return the explicit one-to-one counterpart or refuse ambiguity."""
        key = canonical_locus_key(raw)
        if key in self.ambiguous_components:
            alternatives = ", ".join(self.ambiguous_components[key])
            raise AmbiguousLocusMappingError(
                f"Ambiguous cross-assembly mapping for {raw!r}: {alternatives}"
            )
        return self.aliases.get(key)

    def lookup_keys(self, raw: str, *, fail_on_ambiguity: bool = False) -> set[str]:
        """Return the exact canonical key plus a safe explicit counterpart."""
        key = canonical_locus_key(raw)
        keys = {key}
        try:
            counterpart = self.counterpart(raw)
        except AmbiguousLocusMappingError:
            if fail_on_ambiguity:
                raise
            return keys
        if counterpart:
            keys.add(counterpart)
        return keys

    def query_spellings(self, raw: str) -> set[str]:
        """Return safe exact/crosswalk spellings for external database queries."""
        spellings = set(locus_spelling_variants(raw))
        try:
            counterpart = self.counterpart(raw)
        except AmbiguousLocusMappingError:
            return spellings
        if counterpart:
            spellings.update(locus_spelling_variants(counterpart))
        return spellings

    def build_lookup(self, gene_ids: Iterable[str]) -> dict[str, str]:
        """Build canonical exact/explicit-alias keys for available model genes.

        Exact spellings take precedence.  Any collision is omitted or raised,
        never resolved by iteration order.
        """
        exact_candidates: dict[str, set[str]] = defaultdict(set)
        ids = [str(gene_id) for gene_id in gene_ids]
        for gene_id in ids:
            exact_candidates[canonical_locus_key(gene_id)].add(gene_id)

        exact_conflicts = {
            key: values for key, values in exact_candidates.items() if len(values) > 1
        }
        if exact_conflicts:
            details = "; ".join(
                f"{key}={sorted(values)}" for key, values in sorted(exact_conflicts.items())
            )
            raise LocusLookupConflictError(
                f"Multiple model genes share a canonical locus spelling: {details}"
            )

        lookup = {key: next(iter(values)) for key, values in exact_candidates.items()}
        alias_candidates: dict[str, set[str]] = defaultdict(set)
        for gene_id in ids:
            exact_key = canonical_locus_key(gene_id)
            try:
                counterpart = self.counterpart(gene_id)
            except AmbiguousLocusMappingError:
                continue
            if counterpart and counterpart != exact_key:
                alias_candidates[counterpart].add(gene_id)

        for key, candidates in alias_candidates.items():
            if key in lookup:
                # Preserve a real exact model gene over an alias.
                continue
            if len(candidates) == 1:
                lookup[key] = next(iter(candidates))
        return lookup


def load_default_locus_crosswalk() -> LocusCrosswalk:
    """Load the repository's two authoritative YALI1/YALI0 mapping tables."""
    return LocusCrosswalk.from_csvs(
        identity_exclusions_path=_DEFAULT_IDENTITY_EXCLUSIONS
    )


def model_gene_fingerprint(gene_ids: Iterable[str]) -> str:
    """Fingerprint the exact model-gene namespace used by a derived cache."""
    digest = hashlib.sha256()
    for gene_id in sorted(str(gene_id) for gene_id in gene_ids):
        digest.update(gene_id.encode())
        digest.update(b"\0")
    return f"sha256:{digest.hexdigest()}"


__all__ = [
    "AmbiguousLocusMappingError",
    "LocusCrosswalk",
    "LocusLookupConflictError",
    "canonical_locus_key",
    "load_default_locus_crosswalk",
    "locus_spelling_variants",
    "model_gene_fingerprint",
]
