"""Normalize the Cas9/Cas12a essentiality workbook with row provenance.

The source workbook is a fixed experimental artifact.  This module deliberately
fails on layout drift, missing values, duplicate genes, or unexpected row
counts so a partially read workbook cannot silently become calibration data.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import re
from pathlib import Path
from typing import Mapping

import pandas as pd


ASSAY_SHEETS = ("Cas9", "Cas12a")
SOURCE_COLUMNS = (
    "Gene ID",
    "FS",
    "Raw p-value",
    "Corrected p-value",
    "Essentiality",
)
NORMALIZED_COLUMNS = (
    "gene_id",
    "source_gene_id",
    "assay",
    "fitness_score",
    "raw_p_value",
    "q_value",
    "experimental_call",
    "source_sheet",
    "source_row",
    "source_sha256",
)
DEFAULT_EXPECTED_ROWS = {"Cas9": 7_854, "Cas12a": 7_795}
CALL_MAP = {"Essential": "essential", "Non-essential": "nonessential"}
_SOURCE_GENE_PATTERN = re.compile(r"^YALI[012]_.+\S$")


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_gene_id(value: object) -> str:
    """Convert source IDs such as ``YALI1_E10500g`` to GEM form."""
    return re.sub(r"^(YALI[012])_", r"\1", str(value).strip())


def _require_complete(raw: pd.DataFrame, sheet: str) -> None:
    missing = raw[list(SOURCE_COLUMNS)].isna()
    if missing.any().any():
        row_positions, column_positions = missing.to_numpy().nonzero()
        locations = [
            f"{SOURCE_COLUMNS[column]} row {int(row) + 2}"
            for row, column in zip(row_positions[:10], column_positions[:10])
        ]
        raise ValueError(
            f"{sheet} contains missing values: {', '.join(locations[:10])}"
        )


def _numeric_column(raw: pd.DataFrame, column: str, sheet: str) -> pd.Series:
    numeric = pd.to_numeric(raw[column], errors="coerce")
    invalid = numeric.isna() | ~numeric.map(
        lambda value: bool(pd.notna(value) and math.isfinite(float(value)))
    )
    if invalid.any():
        rows = (invalid[invalid].index + 2).tolist()
        raise ValueError(
            f"{sheet} column {column!r} contains non-numeric values at "
            f"source rows {rows[:10]}"
        )
    return numeric.astype(float)


def normalize_assay_sheet(
    raw: pd.DataFrame,
    *,
    sheet: str,
    expected_rows: int,
    source_sha256: str,
) -> pd.DataFrame:
    """Validate and normalize one already-loaded workbook sheet."""
    if list(raw.columns) != list(SOURCE_COLUMNS):
        raise ValueError(
            f"{sheet} columns do not match the expected source layout. "
            f"Expected {list(SOURCE_COLUMNS)}, found {list(raw.columns)}"
        )
    if len(raw) != expected_rows:
        raise ValueError(
            f"{sheet} row count is {len(raw)}; expected exactly {expected_rows}"
        )
    _require_complete(raw, sheet)

    source_gene_ids = raw["Gene ID"].astype(str).str.strip()
    malformed = ~source_gene_ids.map(
        lambda value: bool(_SOURCE_GENE_PATTERN.fullmatch(value))
    )
    if malformed.any():
        examples = source_gene_ids[malformed].head(10).tolist()
        raise ValueError(f"{sheet} contains malformed source gene IDs: {examples}")

    gene_ids = source_gene_ids.map(normalize_gene_id)
    duplicates = gene_ids[gene_ids.duplicated(keep=False)]
    if not duplicates.empty:
        raise ValueError(
            f"{sheet} contains duplicate gene IDs after normalization: "
            f"{sorted(duplicates.unique())[:10]}"
        )

    calls = raw["Essentiality"].astype(str).str.strip()
    unknown_calls = sorted(set(calls) - set(CALL_MAP))
    if unknown_calls:
        raise ValueError(
            f"{sheet} contains unrecognized Essentiality values: {unknown_calls}"
        )

    fitness = _numeric_column(raw, "FS", sheet)
    raw_p = _numeric_column(raw, "Raw p-value", sheet)
    q_value = _numeric_column(raw, "Corrected p-value", sheet)
    for label, values in (("Raw p-value", raw_p), ("Corrected p-value", q_value)):
        # The published workbook contains floating-point round-off a few parts
        # in 10^12 above one.  Accept only that numerical noise, then clamp.
        invalid = (values < -1e-12) | (values > 1 + 1e-9)
        if invalid.any():
            rows = (invalid[invalid].index + 2).tolist()
            raise ValueError(
                f"{sheet} column {label!r} has values outside [0, 1] at "
                f"source rows {rows[:10]}"
            )

    result = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "source_gene_id": source_gene_ids,
            "assay": sheet,
            "fitness_score": fitness,
            "raw_p_value": raw_p.clip(0.0, 1.0),
            "q_value": q_value.clip(0.0, 1.0),
            "experimental_call": calls.map(CALL_MAP),
            "source_sheet": sheet,
            "source_row": raw.index.to_series().astype(int) + 2,
            "source_sha256": source_sha256,
        }
    )
    if result.isna().any().any():
        raise ValueError(f"{sheet} normalization unexpectedly produced missing values")
    return result[list(NORMALIZED_COLUMNS)]


def normalize_workbook(
    input_workbook: Path,
    *,
    expected_rows: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """Normalize both assay sheets and return one provenance-rich table."""
    input_workbook = Path(input_workbook)
    if not input_workbook.is_file():
        raise ValueError(f"Input workbook not found: {input_workbook}")
    counts = dict(DEFAULT_EXPECTED_ROWS if expected_rows is None else expected_rows)
    if set(counts) != set(ASSAY_SHEETS):
        raise ValueError(
            f"Expected row counts must be provided for {list(ASSAY_SHEETS)}"
        )
    if any(not isinstance(value, int) or value <= 0 for value in counts.values()):
        raise ValueError("Expected row counts must be positive integers")

    workbook = pd.ExcelFile(input_workbook)
    if tuple(workbook.sheet_names) != ASSAY_SHEETS:
        raise ValueError(
            f"Workbook sheets must be exactly {list(ASSAY_SHEETS)} in that order; "
            f"found {workbook.sheet_names}"
        )
    source_sha = sha256_file(input_workbook)
    tables = []
    for sheet in ASSAY_SHEETS:
        raw = pd.read_excel(workbook, sheet_name=sheet)
        tables.append(
            normalize_assay_sheet(
                raw,
                sheet=sheet,
                expected_rows=counts[sheet],
                source_sha256=source_sha,
            )
        )

    result = pd.concat(tables, ignore_index=True)
    duplicates = result.duplicated(["assay", "gene_id"], keep=False)
    if duplicates.any():
        examples = (
            result.loc[duplicates, ["assay", "gene_id"]].head(10).to_dict("records")
        )
        raise ValueError(f"Duplicate assay/gene rows after normalization: {examples}")
    expected_total = sum(counts.values())
    if len(result) != expected_total or result.isna().any().any():
        raise ValueError(
            "Normalized workbook failed row-count or completeness validation"
        )
    return result[list(NORMALIZED_COLUMNS)]


def write_normalized_assays(
    input_workbook: Path,
    output: Path,
    *,
    expected_rows: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """Normalize a workbook and write its stable CSV representation."""
    result = normalize_workbook(input_workbook, expected_rows=expected_rows)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False, float_format="%.17g", lineterminator="\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize Cas9/Cas12a screen fitness data with row provenance."
    )
    parser.add_argument("--input-workbook", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-cas9-rows", type=int, default=7_854)
    parser.add_argument("--expected-cas12a-rows", type=int, default=7_795)
    args = parser.parse_args()
    try:
        result = write_normalized_assays(
            args.input_workbook,
            args.output,
            expected_rows={
                "Cas9": args.expected_cas9_rows,
                "Cas12a": args.expected_cas12a_rows,
            },
        )
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(
        f"Wrote {len(result)} normalized rows "
        f"(Cas9={args.expected_cas9_rows}, Cas12a={args.expected_cas12a_rows}) "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
