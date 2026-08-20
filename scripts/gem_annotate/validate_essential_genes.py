"""Evidence-first gene-essentiality validation for the iYali26 GEM.

The primary reference bundled with this repository is a positive-only list of
1,612 consensus-essential genes from Ramesh et al. (2023).  Because it has no
confirmed non-essential controls, this module reports coverage, TP, FN and
recall without inventing TN/FP labels.

The experimental medium is applied at runtime and never written into model.xml.
All model exchanges are closed first; only exchanges listed in the medium CSV
are opened.  Optional diagnostics classify false negatives and prepare a
curation queue, but they never mutate the model.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
import json
import logging
import math
import re
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion
from cobra.io import read_sbml_model

from .config import (
    ESSENTIALITY_DIR,
    MEDIA_DIR,
    REPO_ROOT,
    RESULTS_DIR,
    load_project_paths,
)
from .essentiality_evidence import (
    DEFAULT_EVIDENCE_DIR,
    DEFAULT_LEDGER,
    canonical_json,
    chemistry_fingerprint,
    merge_detected_cases,
    read_ledger,
    sha256_file,
    stable_case_id,
    target_fingerprint,
)
from .isozyme_capacity import (
    load_isozyme_capacity_scan,
    run_isozyme_capacity_scan,
    validate_isozyme_capacity_scan,
)
from .isozyme_resolution import (
    ISOZYME_SURVIVAL_FLOOR,
    build_isozyme_resolution_ledger,
    classify_isozyme_counterfactual,
)
from .run_registry import (
    build_run_key,
    guard_duplicate_run,
    register_run,
    utc_now as registry_utc_now,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DEFAULT_EXPERIMENTAL = ESSENTIALITY_DIR / "consensus_essential_genes.csv"
DEFAULT_MEDIA = MEDIA_DIR / "sd_leu.csv"
DEFAULT_OUTPUT_DIR = RESULTS_DIR / "essentiality"
DEFAULT_CUTOFFS = (0.01, 0.05, 0.10, 0.15)
PRIMARY_CUTOFF = 0.10
FLUX_EPS = 1e-9
LEUCINE_EXCHANGE_ID = "R1219"
ASSAY_RATIO_METRIC_DECIMALS = 9
ASSAY_NAMES = ("Cas9", "Cas12a")
ASSAY_EXPECTED_ROWS = {"Cas9": 7_854, "Cas12a": 7_795}
ASSAY_FITNESS_COLUMNS = (
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
ASSAY_CALLS = {"essential", "nonessential"}

GENE_COLUMN_ALIASES = {
    "gene",
    "gene_id",
    "clib89 gene id",
    "source_gene_id",
    "yali1_id",
}
FUNCTION_COLUMN_ALIASES = {
    "function",
    "putative function",
    "putative function (obtained from patterson et al. metabolic engineering, 2018)",
    "description",
    "comment",
}


def normalize_gene_id(value: object) -> str:
    """Normalize source IDs such as YALI1_A00309g to SBML form."""
    text = str(value).strip()
    return re.sub(r"^(YALI[012])_", r"\1", text)


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, dtype=str)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, dtype=str)


def _find_column(columns: Iterable[object], aliases: set[str]) -> str | None:
    by_lower = {str(column).strip().lower(): str(column) for column in columns}
    for alias in aliases:
        if alias in by_lower:
            return by_lower[alias]
    return None


def _parse_essential(value: object) -> bool:
    truthy = {"1", "true", "yes", "essential", "y"}
    falsy = {"0", "false", "no", "non-essential", "nonessential", "n"}
    normalized = str(value).strip().lower()
    if normalized in truthy:
        return True
    if normalized in falsy:
        return False
    raise ValueError(f"Unrecognised essentiality value: {value!r}")


def load_experimental(path: Path, positive_only: bool = False) -> pd.DataFrame:
    """Load CSV/TSV/XLSX essentiality data and normalize gene identifiers.

    Positive-only input may omit an ``essential`` column when ``positive_only``
    is true.  The returned table preserves source/function/confidence metadata.
    """
    raw = _read_table(path).rename(columns=lambda c: str(c).strip())
    gene_column = _find_column(raw.columns, GENE_COLUMN_ALIASES)
    if gene_column is None:
        raise ValueError(
            "Experimental table has no recognized gene column. "
            f"Found: {list(raw.columns)}"
        )

    lower_columns = {str(column).strip().lower(): str(column) for column in raw.columns}
    essential_column = lower_columns.get("essential")
    if essential_column is None and not positive_only:
        raise ValueError(
            "Experimental table has no 'essential' column. Pass --positive-only "
            "only when every row is an experimentally essential positive."
        )

    function_column = _find_column(raw.columns, FUNCTION_COLUMN_ALIASES)
    source_column = lower_columns.get("source")
    confidence_column = lower_columns.get("confidence")
    source_gene_column = lower_columns.get("source_gene_id")

    result = pd.DataFrame()
    result["source_gene_id"] = raw[source_gene_column or gene_column].fillna("").str.strip()
    result["gene_id"] = raw[gene_column].map(normalize_gene_id)
    result["function"] = raw[function_column].fillna("") if function_column else ""
    result["source"] = raw[source_column].fillna("") if source_column else str(path)
    result["confidence"] = (
        raw[confidence_column].fillna("") if confidence_column else "unspecified"
    )
    result["essential"] = (
        True if essential_column is None else raw[essential_column].map(_parse_essential)
    )

    result = result[result["gene_id"].ne("")].copy()
    duplicates = result[result["gene_id"].duplicated(keep=False)]
    if not duplicates.empty:
        duplicate_ids = sorted(duplicates["gene_id"].unique())
        raise ValueError(f"Duplicate experimental gene IDs after normalization: {duplicate_ids[:10]}")
    return result.reset_index(drop=True)


def _finite_numeric(values: pd.Series, *, column: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.map(lambda value: bool(pd.notna(value) and math.isfinite(float(value))))
    if not finite.all():
        rows = (finite[~finite].index + 2).tolist()
        raise ValueError(
            f"Assay fitness column {column!r} contains missing or non-finite "
            f"values at file rows {rows[:10]}"
        )
    return numeric.astype(float)


def load_assay_fitness(
    path: Path,
    expected_rows: dict[str, int] | None = ASSAY_EXPECTED_ROWS,
) -> pd.DataFrame:
    """Load the normalized Cas9/Cas12a fitness table and verify provenance.

    Passing ``expected_rows=None`` is intended only for small synthetic unit
    fixtures.  Production callers use the exact published row-count contract.
    """
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Assay fitness file not found: {path}")
    raw = pd.read_csv(path, dtype=str, keep_default_na=False)
    if list(raw.columns) != list(ASSAY_FITNESS_COLUMNS):
        raise ValueError(
            "Assay fitness columns do not match the normalized schema. "
            f"Expected {list(ASSAY_FITNESS_COLUMNS)}, found {list(raw.columns)}"
        )
    empty = raw[list(ASSAY_FITNESS_COLUMNS)].apply(
        lambda column: column.astype(str).str.strip().eq("")
    )
    if empty.any().any():
        row_positions, column_positions = empty.to_numpy().nonzero()
        locations = [
            f"{ASSAY_FITNESS_COLUMNS[column]} row {int(row) + 2}"
            for row, column in zip(row_positions[:10], column_positions[:10])
        ]
        raise ValueError(
            f"Assay fitness file contains missing values: {', '.join(locations)}"
        )

    result = raw.copy()
    result["source_gene_id"] = result["source_gene_id"].str.strip()
    result["gene_id"] = result["gene_id"].map(normalize_gene_id)
    source_normalized = result["source_gene_id"].map(normalize_gene_id)
    id_mismatch = result["gene_id"].ne(source_normalized)
    if id_mismatch.any():
        rows = (id_mismatch[id_mismatch].index + 2).tolist()
        raise ValueError(
            f"gene_id does not match normalized source_gene_id at file rows {rows[:10]}"
        )
    assay_map = {name.casefold(): name for name in ASSAY_NAMES}
    normalized_assays = result["assay"].str.strip().str.casefold().map(assay_map)
    if normalized_assays.isna().any():
        unknown = sorted(result.loc[normalized_assays.isna(), "assay"].unique())
        raise ValueError(f"Unknown assay names in normalized fitness file: {unknown}")
    result["assay"] = normalized_assays
    result["source_sheet"] = result["source_sheet"].str.strip()
    sheet_mismatch = result["source_sheet"].ne(result["assay"])
    if sheet_mismatch.any():
        rows = (sheet_mismatch[sheet_mismatch].index + 2).tolist()
        raise ValueError(f"Assay/source_sheet mismatch at file rows {rows[:10]}")

    result["experimental_call"] = result["experimental_call"].str.strip().str.casefold()
    unknown_calls = sorted(set(result["experimental_call"]) - ASSAY_CALLS)
    if unknown_calls:
        raise ValueError(f"Unknown experimental assay calls: {unknown_calls}")
    for column in ("fitness_score", "raw_p_value", "q_value"):
        result[column] = _finite_numeric(result[column], column=column)
    for column in ("raw_p_value", "q_value"):
        invalid = (result[column] < 0.0) | (result[column] > 1.0)
        if invalid.any():
            rows = (invalid[invalid].index + 2).tolist()
            raise ValueError(
                f"Assay fitness column {column!r} is outside [0, 1] at "
                f"file rows {rows[:10]}"
            )

    source_rows = pd.to_numeric(result["source_row"], errors="coerce")
    invalid_source_rows = (
        source_rows.isna()
        | (source_rows < 2)
        | source_rows.map(lambda value: not float(value).is_integer())
    )
    if invalid_source_rows.any():
        rows = (invalid_source_rows[invalid_source_rows].index + 2).tolist()
        raise ValueError(f"Invalid source_row values at file rows {rows[:10]}")
    result["source_row"] = source_rows.astype(int)

    result["source_sha256"] = result["source_sha256"].str.strip()
    bad_sha = ~result["source_sha256"].str.fullmatch(r"[0-9a-f]{64}")
    if bad_sha.any():
        rows = (bad_sha[bad_sha].index + 2).tolist()
        raise ValueError(f"Invalid source_sha256 values at file rows {rows[:10]}")
    if result["source_sha256"].nunique() != 1:
        raise ValueError("Normalized assay rows must share one source workbook SHA-256")

    duplicates = result.duplicated(["assay", "gene_id"], keep=False)
    if duplicates.any():
        examples = result.loc[duplicates, ["assay", "gene_id"]].head(10).to_dict("records")
        raise ValueError(f"Duplicate assay/gene rows in fitness file: {examples}")
    if set(result["assay"]) != set(ASSAY_NAMES):
        raise ValueError(f"Assay fitness file must contain both {list(ASSAY_NAMES)}")
    if expected_rows is not None:
        if set(expected_rows) != set(ASSAY_NAMES):
            raise ValueError(f"Expected assay row counts must cover {list(ASSAY_NAMES)}")
        counts = result["assay"].value_counts().to_dict()
        for assay in ASSAY_NAMES:
            if counts.get(assay, 0) != expected_rows[assay]:
                raise ValueError(
                    f"{assay} normalized row count is {counts.get(assay, 0)}; "
                    f"expected exactly {expected_rows[assay]}"
                )
            expected_source_rows = set(range(2, expected_rows[assay] + 2))
            observed_source_rows = set(
                result.loc[result["assay"].eq(assay), "source_row"]
            )
            if observed_source_rows != expected_source_rows:
                raise ValueError(
                    f"{assay} source_row provenance is incomplete or duplicated"
                )
    return result[list(ASSAY_FITNESS_COLUMNS)].reset_index(drop=True)


def build_assay_fitness_table(
    assay_fitness: pd.DataFrame,
    predictions: pd.DataFrame,
    primary_cutoff: float,
) -> pd.DataFrame:
    """Join continuous assay phenotypes to model KO growth ratios."""
    required_predictions = {"gene_id", "ko_growth_ratio"}
    missing = required_predictions - set(predictions.columns)
    if missing:
        raise ValueError(f"Prediction table missing columns: {sorted(missing)}")
    if predictions["gene_id"].duplicated().any():
        raise ValueError("Prediction table contains duplicate gene IDs")
    prediction_columns = [
        column
        for column in ("gene_id", "ko_status", "ko_growth", "ko_growth_ratio")
        if column in predictions.columns
    ]
    joined = assay_fitness.merge(
        predictions[prediction_columns],
        how="left",
        on="gene_id",
        validate="many_to_one",
    )
    joined["in_model"] = joined["ko_growth_ratio"].notna()
    calls = pd.Series(pd.NA, index=joined.index, dtype="boolean")
    calls.loc[joined["in_model"]] = (
        joined.loc[joined["in_model"], "ko_growth_ratio"] < primary_cutoff
    )
    joined["predicted_essential_primary"] = calls
    return joined


def _finite_metric_pairs(group: pd.DataFrame) -> pd.DataFrame:
    pairs = group[["fitness_score", "ko_growth_ratio", "experimental_call"]].copy()
    valid = pairs[["fitness_score", "ko_growth_ratio"]].apply(
        lambda column: column.map(
            lambda value: bool(pd.notna(value) and math.isfinite(float(value)))
        )
    ).all(axis=1)
    return pairs.loc[valid]


def _correlation_or_none(left: pd.Series, right: pd.Series) -> float | None:
    if len(left) < 2 or left.nunique() < 2 or right.nunique() < 2:
        return None
    value = float(left.corr(right))
    return value if math.isfinite(value) else None


def _roc_auc(labels: pd.Series, scores: pd.Series) -> float | None:
    """Return tie-aware binary ROC AUC without adding a sklearn dependency."""
    positives = labels.astype(bool)
    n_positive = int(positives.sum())
    n_negative = int(len(positives) - n_positive)
    if n_positive == 0 or n_negative == 0:
        return None
    ranks = scores.astype(float).rank(method="average", ascending=True)
    positive_rank_sum = float(ranks.loc[positives].sum())
    auc = (
        positive_rank_sum - n_positive * (n_positive + 1) / 2
    ) / (n_positive * n_negative)
    return float(auc)


def make_assay_fitness_summary(
    joined: pd.DataFrame,
    model_gene_ids: Iterable[str],
    primary_cutoff: float,
    growth_cutoffs: tuple[float, ...] = DEFAULT_CUTOFFS,
) -> dict[str, object]:
    """Summarize continuous assay/model agreement without altering TP/FN recall."""
    required = {
        "gene_id",
        "assay",
        "fitness_score",
        "experimental_call",
        "ko_growth_ratio",
        "in_model",
    }
    missing = required - set(joined.columns)
    if missing:
        raise ValueError(f"Joined assay table missing columns: {sorted(missing)}")
    model_gene_ids = set(model_gene_ids)
    expected_in_model = joined["gene_id"].isin(model_gene_ids)
    if not expected_in_model.eq(joined["in_model"].astype(bool)).all():
        raise ValueError("Joined assay in_model flags do not match the supplied model genes")

    cutoffs = tuple(sorted(set(growth_cutoffs) | {primary_cutoff}))
    per_assay: dict[str, dict[str, object]] = {}
    for assay in ASSAY_NAMES:
        group = joined[joined["assay"].eq(assay)]
        metric_pairs = _finite_metric_pairs(group)
        # LP solvers can return numerically different values at about 1e-12 for
        # biologically identical KO ratios (especially the large ratio=1 tie).
        # Quantize only the continuous-assay diagnostics so Spearman ranks and
        # AUC ties are reproducible; essentiality calls retain the raw ratios.
        metric_ratios = metric_pairs["ko_growth_ratio"].round(
            ASSAY_RATIO_METRIC_DECIMALS
        )
        spearman = _correlation_or_none(
            metric_pairs["fitness_score"].rank(method="average"),
            metric_ratios.rank(method="average"),
        )
        pearson = _correlation_or_none(
            metric_pairs["fitness_score"], metric_ratios
        )
        auc = _roc_auc(
            metric_pairs["experimental_call"].eq("essential"),
            -metric_ratios,
        )
        in_model = group["in_model"].astype(bool)
        per_assay[assay] = {
            "assay_rows": int(len(group)),
            "assay_unique_genes": int(group["gene_id"].nunique()),
            "model_overlap": int(in_model.sum()),
            "fraction_assay_genes_in_model": (
                float(in_model.sum() / len(group)) if len(group) else None
            ),
            "essential_calls": int(group["experimental_call"].eq("essential").sum()),
            "nonessential_calls": int(
                group["experimental_call"].eq("nonessential").sum()
            ),
            "metric_pairs": int(len(metric_pairs)),
            "spearman_fitness_vs_ko_growth_ratio": spearman,
            "linear_r_squared_fitness_vs_ko_growth_ratio": (
                pearson * pearson if pearson is not None else None
            ),
            "assay_call_roc_auc": auc,
        }

    call_matrix = joined.pivot(
        index="gene_id", columns="assay", values="experimental_call"
    )
    concordant_mask = pd.Series(False, index=call_matrix.index)
    if set(ASSAY_NAMES).issubset(call_matrix.columns):
        concordant_mask = call_matrix[list(ASSAY_NAMES)].eq("nonessential").all(axis=1)
    proxy_ids = set(call_matrix.index[concordant_mask])
    proxy = (
        joined[joined["gene_id"].isin(proxy_ids) & joined["in_model"].astype(bool)]
        .drop_duplicates("gene_id")
        .copy()
    )
    curve: list[dict[str, object]] = []
    for cutoff in cutoffs:
        predicted_essential = int((proxy["ko_growth_ratio"] < cutoff).sum())
        safe = int(len(proxy) - predicted_essential)
        curve.append(
            {
                "cutoff_fraction_of_wt": cutoff,
                "model_proxy_genes": int(len(proxy)),
                "predicted_essential_count": predicted_essential,
                "safe_nonessential_count": safe,
                "safety_rate": float(safe / len(proxy)) if len(proxy) else None,
            }
        )
    primary = next(
        item for item in curve if item["cutoff_fraction_of_wt"] == primary_cutoff
    )
    return {
        "schema_version": "1.0",
        "semantics": {
            "fitness_correlation": "fitness_score versus model KO/WT growth ratio",
            "linear_r_squared": "ordinary least-squares R-squared with an intercept",
            "ko_ratio_metric_quantization": (
                f"rounded to {ASSAY_RATIO_METRIC_DECIMALS} decimal places for "
                "solver-stable continuous metrics; threshold calls use raw ratios"
            ),
            "roc_positive_class": "experimental essential",
            "roc_model_score": "negative KO/WT growth ratio",
            "proxy_definition": (
                "genes called nonessential by both Cas9 and Cas12a; safety means "
                "the model does not call the gene essential at the stated cutoff"
            ),
            "positive_only_recall_unchanged": True,
        },
        "per_assay": per_assay,
        "concordant_nonessential_proxy_safety": {
            "source_proxy_genes": int(len(proxy_ids)),
            "model_proxy_genes": int(len(proxy)),
            "primary": primary,
            "cutoff_curve": curve,
        },
    }


def load_media(path: Path) -> pd.DataFrame:
    media = pd.read_csv(path, dtype={"exchange": str})
    required = {"exchange", "uptake"}
    missing = required - set(media.columns)
    if missing:
        raise ValueError(f"Media file missing columns {sorted(missing)}: {path}")
    media["exchange"] = media["exchange"].str.strip()
    media["uptake"] = pd.to_numeric(media["uptake"], errors="raise")
    if media["exchange"].duplicated().any():
        repeated = sorted(media.loc[media["exchange"].duplicated(), "exchange"].unique())
        raise ValueError(f"Duplicate exchanges in media file: {repeated}")
    invalid = media[(media["uptake"] < 0) | (media["uptake"] > 1000)]
    if not invalid.empty:
        raise ValueError("Media uptake values must be between 0 and 1000")
    return media


def apply_media(model, media: pd.DataFrame) -> dict[str, float]:
    """Replace the model medium with the exact listed uptake definition."""
    model_reactions = {reaction.id for reaction in model.reactions}
    missing = sorted(set(media["exchange"]) - model_reactions)
    if missing:
        raise ValueError(f"Media exchange reactions not found in model: {missing}")

    allowed = {
        row.exchange: float(row.uptake)
        for row in media.itertuples(index=False)
        if float(row.uptake) > 0
    }
    model.medium = allowed
    if LEUCINE_EXCHANGE_ID in model.medium:
        raise ValueError(
            f"SD-Leu invariant violated: {LEUCINE_EXCHANGE_ID} is open for uptake"
        )
    return dict(model.medium)


def _gene_id_from_deletion_row(index: object, row: pd.Series) -> str:
    ids = row.get("ids")
    if isinstance(ids, (set, frozenset)) and ids:
        return str(next(iter(ids)))
    if isinstance(index, (set, frozenset)) and index:
        return str(next(iter(index)))
    return str(ids if ids is not None else index)


def run_single_gene_deletions(model, solver: str) -> tuple[pd.DataFrame, float]:
    model.solver = solver
    solution = model.optimize()
    if solution.status != "optimal" or solution.objective_value is None:
        raise RuntimeError(f"Wild-type FBA is not optimal: {solution.status}")
    wt_growth = float(solution.objective_value)
    if not (0.1 <= wt_growth <= 2.0):
        raise RuntimeError(
            f"Wild-type growth {wt_growth:.6g} h^-1 is outside the accepted "
            "0.1-2.0 h^-1 range under the experimental medium"
        )

    gene_ids = sorted(gene.id for gene in model.genes)
    logger.info("Running deterministic single-gene deletion for %d model genes", len(gene_ids))
    deletion = single_gene_deletion(model, gene_list=gene_ids, processes=1)
    rows: list[dict[str, object]] = []
    for index, row in deletion.iterrows():
        gene_id = _gene_id_from_deletion_row(index, row)
        status = str(row.get("status", "optimal"))
        growth_value = row.get("growth")
        if status != "optimal" or growth_value is None or pd.isna(growth_value):
            growth = 0.0
        else:
            growth = max(0.0, float(growth_value))
        rows.append(
            {
                "gene_id": gene_id,
                "ko_status": status,
                "ko_growth": growth,
                "ko_growth_ratio": growth / wt_growth,
            }
        )
    predictions = pd.DataFrame(rows).sort_values("gene_id").reset_index(drop=True)
    return predictions, wt_growth


def _cutoff_label(cutoff: float) -> str:
    percent = cutoff * 100
    return f"essential_at_{percent:g}pct".replace(".", "p")


def build_per_gene_table(
    experimental: pd.DataFrame,
    predictions: pd.DataFrame,
    cutoffs: tuple[float, ...],
    primary_cutoff: float,
) -> pd.DataFrame:
    prediction_map = predictions.set_index("gene_id")
    experimental_map = experimental.set_index("gene_id")
    union_ids = sorted(set(experimental["gene_id"]) | set(predictions["gene_id"]))
    rows: list[dict[str, object]] = []

    for gene_id in union_ids:
        in_experiment = gene_id in experimental_map.index
        in_model = gene_id in prediction_map.index
        exp = experimental_map.loc[gene_id] if in_experiment else None
        pred = prediction_map.loc[gene_id] if in_model else None
        ratio = float(pred["ko_growth_ratio"]) if in_model else math.nan
        experimental_essential = bool(exp["essential"]) if in_experiment else None
        predicted_primary = bool(in_model and ratio < primary_cutoff)

        if not in_model:
            classification = "outside_gem_scope"
        elif in_experiment and experimental_essential:
            classification = "TP" if predicted_primary else "FN"
        elif in_experiment:
            classification = "FP" if predicted_primary else "TN"
        elif predicted_primary:
            classification = "unverified_prediction"
        else:
            classification = "not_in_positive_reference"

        record: dict[str, object] = {
            "gene_id": gene_id,
            "source_gene_id": exp["source_gene_id"] if in_experiment else "",
            "function": exp["function"] if in_experiment else "",
            "source": exp["source"] if in_experiment else "",
            "confidence": exp["confidence"] if in_experiment else "",
            "experimental_essential": experimental_essential,
            "in_model": in_model,
            "ko_status": pred["ko_status"] if in_model else "",
            "ko_growth": float(pred["ko_growth"]) if in_model else math.nan,
            "ko_growth_ratio": ratio,
            "predicted_essential_primary": predicted_primary if in_model else None,
            "classification": classification,
        }
        for cutoff in cutoffs:
            record[_cutoff_label(cutoff)] = bool(in_model and ratio < cutoff) if in_model else None
        rows.append(record)

    return pd.DataFrame(rows)


def make_summary(
    model,
    experimental: pd.DataFrame,
    per_gene: pd.DataFrame,
    wt_growth: float,
    cutoffs: tuple[float, ...],
    primary_cutoff: float,
    experimental_path: Path,
    media_path: Path,
    active_medium: dict[str, float],
    positive_only: bool,
) -> dict[str, object]:
    experimental_ids = set(experimental["gene_id"])
    model_ids = {gene.id for gene in model.genes}
    overlap = experimental_ids & model_ids
    experimental_positive = per_gene[
        per_gene["gene_id"].isin(overlap) & per_gene["experimental_essential"].eq(True)
    ]

    curve = []
    for cutoff in cutoffs:
        calls = experimental_positive["ko_growth_ratio"] < cutoff
        tp = int(calls.sum())
        fn = int(len(calls) - tp)
        curve.append(
            {
                "cutoff_fraction_of_wt": cutoff,
                "TP": tp,
                "FN": fn,
                "recall": tp / len(calls) if len(calls) else None,
            }
        )

    primary = next(item for item in curve if item["cutoff_fraction_of_wt"] == primary_cutoff)
    return {
        "schema_version": "2.0",
        "model": {
            "id": model.id,
            "num_genes": len(model.genes),
            "num_reactions": len(model.reactions),
            "num_metabolites": len(model.metabolites),
        },
        "experimental": {
            "path": str(experimental_path.resolve()),
            "positive_only": positive_only,
            "num_rows": len(experimental),
            "num_positive": int(experimental["essential"].sum()),
        },
        "coverage": {
            "experimental_genes": len(experimental_ids),
            "model_genes": len(model_ids),
            "intersection": len(overlap),
            "outside_gem_scope": len(experimental_ids - model_ids),
            "model_only_unlabelled": len(model_ids - experimental_ids),
        },
        "medium": {
            "path": str(media_path.resolve()),
            "num_open_uptakes": len(active_medium),
            "leucine_exchange": LEUCINE_EXCHANGE_ID,
            "leucine_uptake_open": LEUCINE_EXCHANGE_ID in active_medium,
        },
        "wt_growth": wt_growth,
        "primary_cutoff_fraction_of_wt": primary_cutoff,
        "primary": primary,
        "cutoff_curve": curve,
        "unverified_predictions": int(
            per_gene["classification"].eq("unverified_prediction").sum()
        ),
        "notes": [
            "Positive-only reference: accuracy, MCC, precision and specificity are not reported.",
            "Genes outside the GEM are outside scope, not experimentally non-essential.",
        ],
    }


def _fva_capacity(model, reactions: list, fraction: float = 1.0) -> float:
    if not reactions:
        return 0.0
    fva = flux_variability_analysis(
        model,
        reaction_list=reactions,
        fraction_of_optimum=fraction,
        processes=1,
    )
    if fva.empty:
        return 0.0
    return float(fva[["minimum", "maximum"]].abs().to_numpy().max())


def _open_uptake_metabolites(model) -> set[str]:
    active = set(model.medium)
    metabolite_ids: set[str] = set()
    for reaction_id in active:
        reaction = model.reactions.get_by_id(reaction_id)
        metabolite_ids.update(met.id for met in reaction.metabolites)
    return metabolite_ids


def _closed_reaction_growth(model, reactions: Iterable[object], wt_growth: float) -> tuple[float, float]:
    """Return growth and WT ratio after closing all supplied reactions."""
    reactions = list(reactions)
    if not reactions:
        return wt_growth, 1.0
    with model:
        for reaction in reactions:
            reaction.bounds = (0.0, 0.0)
        value = model.slim_optimize()
        growth = 0.0 if value is None else max(0.0, float(value))
    return growth, growth / wt_growth


def _verify_bypass_candidates(
    model,
    gene,
    wt_fluxes: pd.Series,
    lethal_growth: float,
    max_candidates: int,
) -> list[str]:
    """Verify high-flux-delta single-reaction bypasses with a combined KO.

    Candidate generation is heuristic because alternate optima can change the
    chosen FBA vertex.  A reported bypass is nevertheless exact: blocking it
    together with the gene must reduce growth below the primary threshold.
    """
    with model:
        gene.knock_out()
        ko_solution = model.optimize()
        if ko_solution.status != "optimal":
            return []
        delta = (ko_solution.fluxes - wt_fluxes).abs().sort_values(ascending=False)
        gene_reactions = {reaction.id for reaction in gene.reactions}
        candidate_ids = [
            rid
            for rid, value in delta.items()
            if value > 1e-6 and rid not in gene_reactions
        ][:max_candidates]

    verified: list[str] = []
    for reaction_id in candidate_ids:
        with model:
            gene.knock_out()
            model.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)
            growth = model.slim_optimize()
            growth = 0.0 if growth is None else max(0.0, float(growth))
        if growth < lethal_growth:
            verified.append(reaction_id)
    return verified


def diagnose_false_negatives(
    model,
    per_gene: pd.DataFrame,
    wt_growth: float,
    primary_cutoff: float,
    max_bypass_candidates: int,
) -> pd.DataFrame:
    fn_ids = per_gene.loc[per_gene["classification"].eq("FN"), "gene_id"].tolist()
    logger.info("Diagnosing %d false-negative genes with FVA", len(fn_ids))
    uptake_metabolites = _open_uptake_metabolites(model)
    wt_solution = model.optimize()
    wt_fluxes = wt_solution.fluxes
    lethal_growth = primary_cutoff * wt_growth
    rows: list[dict[str, object]] = []

    for index, gene_id in enumerate(fn_ids, start=1):
        gene = model.genes.get_by_id(gene_id)
        reactions = sorted(gene.reactions, key=lambda reaction: reaction.id)
        wt_capacity = _fva_capacity(model, reactions)
        with model:
            gene.knock_out()
            ko_growth_value = model.slim_optimize()
            ko_growth = 0.0 if ko_growth_value is None else max(0.0, float(ko_growth_value))
            ko_capacity = _fva_capacity(model, reactions)
        reaction_ko_growth, reaction_ko_ratio = _closed_reaction_growth(
            model, reactions, wt_growth
        )

        if wt_capacity <= FLUX_EPS:
            category = "inactive_reaction"
        elif ko_capacity > FLUX_EPS:
            category = "isozyme_redundancy"
        else:
            reaction_metabolites = {met.id for reaction in reactions for met in reaction.metabolites}
            category = (
                "nutrient_bypass"
                if reaction_metabolites & uptake_metabolites
                else "metabolic_bypass"
            )

        verified_bypasses: list[str] = []
        if category == "metabolic_bypass" and max_bypass_candidates > 0:
            verified_bypasses = _verify_bypass_candidates(
                model,
                gene,
                wt_fluxes,
                lethal_growth,
                max_bypass_candidates,
            )

        rows.append(
            {
                "gene_id": gene_id,
                "category": category,
                "n_reactions": len(reactions),
                "reaction_ids": ";".join(reaction.id for reaction in reactions),
                "reaction_names": ";".join((reaction.name or "") for reaction in reactions),
                "gpr": " | ".join(reaction.gene_reaction_rule for reaction in reactions),
                "wt_capacity": wt_capacity,
                "ko_capacity": ko_capacity,
                "ko_growth": ko_growth,
                "ko_growth_ratio": ko_growth / wt_growth,
                "all_linked_reactions_closed_growth": reaction_ko_growth,
                "all_linked_reactions_closed_growth_ratio": reaction_ko_ratio,
                "model_causal_isozyme": (
                    category == "isozyme_redundancy"
                    and ko_growth / wt_growth > ISOZYME_SURVIVAL_FLOOR
                    and reaction_ko_ratio < primary_cutoff
                ),
                "isozyme_interpretation": (
                    classify_isozyme_counterfactual(
                        [
                            {
                                "ko_growth_ratio": ko_growth / wt_growth,
                                "all_linked_reactions_closed_growth_ratio": reaction_ko_ratio,
                            }
                        ],
                        primary_cutoff=primary_cutoff,
                    )
                    if category == "isozyme_redundancy"
                    else "not_applicable"
                ),
                "verified_single_bypasses": ";".join(verified_bypasses),
                "n_verified_single_bypasses": len(verified_bypasses),
            }
        )
        if index % 25 == 0 or index == len(fn_ids):
            logger.info("  diagnosed %d/%d", index, len(fn_ids))

    return pd.DataFrame(rows)


CURATION_RULES = {
    "isozyme_redundancy": (
        1,
        "Review OR/AND logic, complex membership, compartment and enzyme identity; change GPR only with external evidence.",
    ),
    "inactive_reaction": (
        2,
        "Audit zero bounds, dead ends, directionality and evidence-backed biomass requirements; do not force flux.",
    ),
    "nutrient_bypass": (
        3,
        "Confirm the nutrient is present in SD-Leu and that uptake/transport direction is biologically valid.",
    ),
    "metabolic_bypass": (
        4,
        "Review verified combined-KO bypasses; remove or constrain only with thermodynamic/localisation evidence.",
    ),
}


def build_curation_queue(diagnostics: pd.DataFrame) -> pd.DataFrame:
    queue = diagnostics.copy()
    queue["priority"] = queue["category"].map(lambda category: CURATION_RULES[category][0])
    queue["recommended_action"] = queue["category"].map(
        lambda category: CURATION_RULES[category][1]
    )
    queue["evidence_status"] = "requires_manual_review"
    queue["proposed_operation"] = ""
    queue["evidence_url"] = ""
    queue["rationale"] = ""
    columns = [
        "priority",
        "gene_id",
        "category",
        "reaction_ids",
        "reaction_names",
        "gpr",
        "ko_growth",
        "ko_growth_ratio",
        "all_linked_reactions_closed_growth",
        "all_linked_reactions_closed_growth_ratio",
        "model_causal_isozyme",
        "isozyme_interpretation",
        "wt_capacity",
        "ko_capacity",
        "verified_single_bypasses",
        "n_verified_single_bypasses",
        "recommended_action",
        "evidence_status",
        "proposed_operation",
        "evidence_url",
        "rationale",
    ]
    return queue[columns].sort_values(["priority", "ko_growth_ratio", "gene_id"])


def _split_semicolon_ids(value: object) -> list[str]:
    return sorted({item.strip() for item in str(value or "").split(";") if item.strip()})


def _json_safe(value: Any) -> Any:
    """Normalize COBRA annotations and numeric scalars for JSON packets."""
    return json.loads(canonical_json(value))


def reaction_case_context(reaction) -> dict[str, Any]:
    """Extract the complete local reaction context used by literature agents."""
    stoichiometry = {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in sorted(
            reaction.metabolites.items(), key=lambda item: item[0].id
        )
    }
    compartments = sorted({metabolite.compartment for metabolite in reaction.metabolites})
    metabolite_chemistry = {
        metabolite.id: {
            "formula": metabolite.formula,
            "charge": metabolite.charge,
            "compartment": metabolite.compartment,
        }
        for metabolite in sorted(reaction.metabolites, key=lambda item: item.id)
    }
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mass_balance = _json_safe(reaction.check_mass_balance())
        mass_balance_error = ""
    except (TypeError, ValueError) as exc:
        mass_balance = {}
        mass_balance_error = str(exc)
    return {
        "reaction_id": reaction.id,
        "name": reaction.name or "",
        "equation": reaction.reaction,
        "stoichiometry": stoichiometry,
        "metabolite_chemistry": _json_safe(metabolite_chemistry),
        "compartments": compartments,
        "lower_bound": float(reaction.lower_bound),
        "upper_bound": float(reaction.upper_bound),
        "gpr": reaction.gene_reaction_rule,
        "gpr_gene_ids": sorted(gene.id for gene in reaction.genes),
        "ec_codes": _json_safe((reaction.annotation or {}).get("ec-code", [])),
        "annotations": _json_safe(reaction.annotation or {}),
        "subsystem": reaction.subsystem or "",
        "mass_balance": mass_balance,
        "mass_balance_error": mass_balance_error,
    }


def _gene_case_context(gene, case_gene_ids: set[str]) -> dict[str, Any]:
    reactions = sorted(gene.reactions, key=lambda reaction: reaction.id)
    partners = sorted(
        {
            partner.id
            for reaction in reactions
            for partner in reaction.genes
            if partner.id not in case_gene_ids
        }
    )
    return {
        "gene_id": gene.id,
        "name": gene.name or "",
        "annotations": _json_safe(gene.annotation or {}),
        "reaction_ids": [reaction.id for reaction in reactions],
        "gpr_partners": partners,
    }


CASE_QUESTIONS = {
    "isozyme_redundancy": [
        "Are the proteins encoded by the case genes true isozymes, or obligatory subunits of one enzyme complex in Y. lipolytica?",
        "Does direct Y. lipolytica evidence support the complete OR/AND GPR and the represented compartment?",
    ],
    "metabolic_bypass": [
        "Is each simulated bypass reaction directly demonstrated in Y. lipolytica under a compatible condition?",
        "Are the bypass direction, compartment and any required transport step biologically supported?",
    ],
    "inactive_reaction": [
        "Why can the associated reaction carry no flux: bounds, direction, transport, dead-end connectivity, or absent biomass demand?",
        "Is there direct Y. lipolytica evidence for the reaction, its compartment and any proposed biomass or translation coupling?",
    ],
    "nutrient_bypass": [
        "Is the bypass nutrient actually present and bioavailable in the defined SD-Leu experiment?",
        "Is the modelled uptake or intracellular transport direction supported in Y. lipolytica?",
    ],
}


CASE_CLAIMS = {
    "isozyme_redundancy": "The current GPR may encode isozyme redundancy where the proteins are instead non-redundant complex subunits or differently localized enzymes.",
    "metabolic_bypass": "A modelled metabolic bypass may explain survival even though the experimental gene is essential.",
    "inactive_reaction": "The gene-associated reaction is inactive because of a connectivity, directionality, transport, biomass, or translation representation issue.",
    "nutrient_bypass": "An allowed nutrient uptake or transport route may bypass the experimentally essential function.",
}


def _case_priority(
    category: str,
    diagnostic_rows: list[dict[str, Any]],
    primary_cutoff: float,
) -> tuple[int, str]:
    causal_isozyme_signal = any(
        float(row.get("ko_growth_ratio", 0.0)) > ISOZYME_SURVIVAL_FLOOR
        and float(row.get("all_linked_reactions_closed_growth_ratio", 1.0))
        < primary_cutoff
        for row in diagnostic_rows
    )
    verified = any(
        int(row.get("n_verified_single_bypasses", 0)) > 0 for row in diagnostic_rows
    )
    if category == "isozyme_redundancy" and causal_isozyme_signal:
        return 1, "gene KO survives but closing every linked reaction is lethal"
    if category == "metabolic_bypass" and verified:
        return 2, "gene plus reaction KO verified at least one single-reaction bypass"
    if category == "inactive_reaction":
        return 3, "inactive reaction requires connectivity/biomass/translation review"
    if category == "nutrient_bypass":
        return 4, "possible nutrient uptake or transport bypass"
    return 5, "lower-priority unverified redundancy or bypass"


def build_agent_cases(
    model,
    diagnostics: pd.DataFrame,
    per_gene: pd.DataFrame,
    summary: dict[str, object],
    primary_cutoff: float,
) -> list[dict[str, Any]]:
    """Group fresh FN diagnostics into deterministic, self-contained packets."""
    per_gene_map = per_gene.set_index("gene_id")
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in diagnostics.to_dict(orient="records"):
        reaction_ids = tuple(_split_semicolon_ids(row.get("reaction_ids", "")))
        gpr_signature = tuple(
            (
                reaction_id,
                model.reactions.get_by_id(reaction_id).gene_reaction_rule,
            )
            for reaction_id in reaction_ids
        )
        key = (str(row["category"]), reaction_ids, gpr_signature)
        groups.setdefault(key, []).append(row)

    cases: list[dict[str, Any]] = []
    cutoff_labels = [_cutoff_label(float(item["cutoff_fraction_of_wt"])) for item in summary["cutoff_curve"]]
    for (category, reaction_ids, _gpr_signature), diagnostic_rows in groups.items():
        gene_ids = sorted({str(row["gene_id"]) for row in diagnostic_rows})
        case_gene_ids = set(gene_ids)
        reactions = [model.reactions.get_by_id(reaction_id) for reaction_id in reaction_ids]
        reaction_contexts = [reaction_case_context(reaction) for reaction in reactions]
        verified_bypasses = sorted(
            {
                bypass_id
                for row in diagnostic_rows
                for bypass_id in _split_semicolon_ids(row.get("verified_single_bypasses", ""))
            }
        )
        bypass_contexts = [
            reaction_case_context(model.reactions.get_by_id(reaction_id))
            for reaction_id in verified_bypasses
        ]
        priority, ranking_reason = _case_priority(
            category, diagnostic_rows, primary_cutoff
        )
        case_id = stable_case_id(category, gene_ids, reaction_ids)
        threshold_results = {}
        for gene_id in gene_ids:
            per_gene_row = per_gene_map.loc[gene_id]
            threshold_results[gene_id] = {
                label: bool(per_gene_row[label]) for label in cutoff_labels
            }
        packet: dict[str, Any] = {
            "schema_version": "2.0",
            "case_id": case_id,
            "category": category,
            "priority": priority,
            "ranking_reason": ranking_reason,
            "gene_ids": gene_ids,
            "reaction_ids": list(reaction_ids),
            "model_sha256": summary["model"]["sha256"],
            "experimental_sha256": summary["experimental"]["sha256"],
            "media_sha256": summary["medium"]["sha256"],
            "target_fingerprint": target_fingerprint(reaction_contexts),
            "chemistry_fingerprint": chemistry_fingerprint(reaction_contexts),
            "claim_under_review": CASE_CLAIMS[category],
            "falsifiable_questions": CASE_QUESTIONS[category],
            "threshold_results": threshold_results,
            "model_context": {
                "model_id": model.id,
                "wt_growth": summary["wt_growth"],
                "primary_cutoff_fraction_of_wt": primary_cutoff,
                "medium": summary["medium"],
                "genes": [
                    _gene_case_context(model.genes.get_by_id(gene_id), case_gene_ids)
                    for gene_id in gene_ids
                ],
                "reactions": reaction_contexts,
                "diagnostics": _json_safe(diagnostic_rows),
                "verified_bypasses": bypass_contexts,
            },
        }
        packet["case_packet_sha256"] = sha256_file_payload(packet)
        cases.append(packet)
    return sorted(cases, key=lambda case: (case["priority"], case["case_id"]))


def sha256_file_payload(packet: dict[str, Any]) -> str:
    """Hash an in-memory packet using the same SHA-256 convention as files."""
    import hashlib

    return hashlib.sha256(canonical_json(packet).encode("utf-8")).hexdigest()


_CASE_PROVENANCE_FIELDS = (
    "model_sha256",
    "experimental_sha256",
    "media_sha256",
    "target_fingerprint",
    "chemistry_fingerprint",
)


def _case_max_ko_growth_ratio(case: dict[str, Any]) -> float:
    diagnostics = case.get("model_context", {}).get("diagnostics", [])
    ratios = [
        float(row.get("ko_growth_ratio", 0.0))
        for row in diagnostics
        if isinstance(row, dict)
    ]
    return max(ratios, default=0.0)


def _case_batch_sort_key(case: dict[str, Any]) -> tuple[int, float, str]:
    """Prefer verified bypasses, then higher KO survival, deterministically."""
    return (
        int(case.get("priority", 99)),
        -_case_max_ko_growth_ratio(case),
        str(case.get("case_id", "")),
    )


def _case_will_be_queued(
    case: dict[str, Any], row: dict[str, str] | None
) -> bool:
    if row is None:
        return True
    hashes_match = all(
        str(row.get(field, "")) == str(case.get(field, ""))
        for field in _CASE_PROVENANCE_FIELDS
    )
    return not hashes_match or row.get("status") == "queued"


def prepare_agent_case_files(
    cases: list[dict[str, Any]],
    output_dir: Path,
    batch_size: int,
    ledger_path: Path = DEFAULT_LEDGER,
    evidence_dir: Path = DEFAULT_EVIDENCE_DIR,
    case_category: str | None = None,
    require_full_batch: bool = False,
) -> dict[str, Any]:
    """Persist durable case state and write fresh ignored agent packets."""
    if case_category is not None and case_category not in CURATION_RULES:
        raise ValueError(f"Unknown essentiality case category: {case_category}")
    cases = [
        case
        for case in cases
        if case_category is None or case.get("category") == case_category
    ]
    cases = sorted(cases, key=_case_batch_sort_key)

    current_rows = {row["case_id"]: row for row in read_ledger(ledger_path)}
    preflight_eligible = [
        case
        for case in cases
        if _case_will_be_queued(case, current_rows.get(str(case["case_id"])))
    ]
    if require_full_batch and len(preflight_eligible) < batch_size:
        raise ValueError(
            f"Only {len(preflight_eligible)} queued {case_category or 'essentiality'} "
            f"cases are available; {batch_size} are required"
        )

    ledger = merge_detected_cases(cases, ledger_path=ledger_path, evidence_dir=evidence_dir)
    ledger_by_id = {row["case_id"]: row for row in ledger}
    eligible = [
        case
        for case in cases
        if ledger_by_id[case["case_id"]]["status"] == "queued"
    ]
    selected_ids = {case["case_id"] for case in eligible[:batch_size]}

    cases_path = output_dir / "essentiality_agent_cases.jsonl"
    cases_path.write_text(
        "".join(json.dumps(case, sort_keys=True) + "\n" for case in cases),
        encoding="utf-8",
    )
    batch_path = output_dir / "essentiality_agent_batch.json"
    batch = [case for case in cases if case["case_id"] in selected_ids]
    batch_path.write_text(
        json.dumps(batch, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    queue_rows = []
    for case in cases:
        row = ledger_by_id[case["case_id"]]
        queue_rows.append(
            {
                "priority": case["priority"],
                "case_id": case["case_id"],
                "status": row["status"],
                "selected_for_batch": case["case_id"] in selected_ids,
                "category": case["category"],
                "gene_ids": ";".join(case["gene_ids"]),
                "reaction_ids": ";".join(case["reaction_ids"]),
                "ranking_reason": case["ranking_reason"],
                "target_fingerprint": case["target_fingerprint"],
                "model_sha256": case["model_sha256"],
                "evidence_path": row["evidence_path"],
            }
        )
    queue = pd.DataFrame(queue_rows).sort_values(["priority", "case_id"])
    queue_path = output_dir / "essentiality_agent_queue.tsv"
    _write_tsv(queue, queue_path)
    return {
        "num_cases": len(cases),
        "num_queued": len(eligible),
        "batch_size": batch_size,
        "case_category": case_category or "all",
        "selected_case_ids": sorted(selected_ids),
        "cases_path": str(cases_path.resolve()),
        "batch_path": str(batch_path.resolve()),
        "queue_path": str(queue_path.resolve()),
        "ledger_path": str(ledger_path.resolve()),
        "evidence_dir": str(evidence_dir.resolve()),
    }


def audit_translation_module(model) -> tuple[pd.DataFrame, dict[str, object]]:
    """Audit the locked tRNA-aware biomass reactions without enabling them."""
    candidate_ids = [rid for rid in ("R1387", "R1710") if rid in model.reactions]
    if not candidate_ids:
        return pd.DataFrame(), {
            "candidate_reactions": [],
            "charged_trna_count": 0,
            "ready_to_connect": False,
            "blocking_reason": "No tRNA-aware biomass candidate reaction exists",
        }

    candidate = model.reactions.get_by_id(candidate_ids[0])
    candidate_fluxes: dict[str, float] = {}
    for candidate_id in candidate_ids:
        with model:
            reaction = model.reactions.get_by_id(candidate_id)
            reaction.bounds = (0.0, 1000.0)
            model.objective = reaction
            flux = model.slim_optimize()
            candidate_fluxes[candidate_id] = (
                0.0 if flux is None else max(0.0, float(flux))
            )
    charged = sorted(
        [
            metabolite
            for metabolite, coefficient in candidate.metabolites.items()
            if coefficient < 0 and "trna" in (metabolite.name or "").lower()
        ],
        key=lambda metabolite: metabolite.id,
    )
    candidate_set = set(candidate_ids)
    rows: list[dict[str, object]] = []
    complete = True
    carriers_balanced = True

    for metabolite in charged:
        charging = sorted(
            [reaction for reaction in metabolite.reactions if reaction.id not in candidate_set],
            key=lambda reaction: reaction.id,
        )
        uncharged = sorted(
            {
                other.id
                for reaction in charging
                for other in reaction.metabolites
                if other.id != metabolite.id and "trna" in (other.name or "").lower()
            }
        )
        open_charging = [
            reaction.id
            for reaction in charging
            if reaction.lower_bound != 0 or reaction.upper_bound != 0
        ]
        uncharged_in_candidate = [
            metabolite_id
            for metabolite_id in uncharged
            if model.metabolites.get_by_id(metabolite_id) in candidate.metabolites
        ]
        carrier_balanced = False
        if len(uncharged_in_candidate) == 1:
            uncharged_metabolite = model.metabolites.get_by_id(uncharged_in_candidate[0])
            carrier_balanced = math.isclose(
                abs(float(candidate.metabolites[metabolite])),
                abs(float(candidate.metabolites[uncharged_metabolite])),
                rel_tol=1e-9,
                abs_tol=1e-12,
            )
        pair_complete = bool(uncharged and charging and open_charging)
        complete = complete and pair_complete
        carriers_balanced = carriers_balanced and carrier_balanced
        rows.append(
            {
                "charged_trna_id": metabolite.id,
                "charged_trna_name": metabolite.name or "",
                "uncharged_trna_ids": ";".join(uncharged),
                "charging_reaction_ids": ";".join(reaction.id for reaction in charging),
                "open_charging_reaction_ids": ";".join(open_charging),
                "charging_gprs": " | ".join(reaction.gene_reaction_rule for reaction in charging),
                "pair_complete": pair_complete,
                "carrier_balanced_in_candidate": carrier_balanced,
            }
        )

    expected_count = 20
    candidate_feasible = any(value > FLUX_EPS for value in candidate_fluxes.values())
    ready = (
        complete
        and carriers_balanced
        and len(charged) == expected_count
        and candidate_feasible
    )
    summary = {
        "candidate_reactions": candidate_ids,
        "charged_trna_count": len(charged),
        "expected_charged_trna_count": expected_count,
        "all_pairs_complete": complete,
        "all_carriers_balanced": carriers_balanced,
        "candidate_max_flux": candidate_fluxes,
        "candidate_feasible": candidate_feasible,
        "ready_to_connect": ready,
        "blocking_reason": (
            ""
            if ready
            else "Do not connect to biomass until all 20 carrier pairs are balanced and a tRNA-aware biomass candidate carries flux"
        ),
    }
    return pd.DataFrame(rows), summary


def _write_tsv(dataframe: pd.DataFrame, path: Path) -> None:
    dataframe.to_csv(path, sep="\t", index=False, na_rep="")


def _parse_cutoffs(value: str) -> tuple[float, ...]:
    cutoffs = tuple(sorted({float(item.strip()) for item in value.split(",") if item.strip()}))
    if not cutoffs or any(cutoff <= 0 or cutoff >= 1 for cutoff in cutoffs):
        raise argparse.ArgumentTypeError("growth cutoffs must be fractions strictly between 0 and 1")
    return cutoffs


def _distribution_version(distribution: str) -> str | None:
    try:
        return importlib_metadata.version(distribution)
    except importlib_metadata.PackageNotFoundError:
        return None


def _git_value(repo_root: Path, *arguments: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=repo_root,
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def _git_provenance(repo_root: Path) -> dict[str, object]:
    status = _git_value(repo_root, "status", "--porcelain")
    return {
        "commit": _git_value(repo_root, "rev-parse", "HEAD"),
        "branch": _git_value(repo_root, "rev-parse", "--abbrev-ref", "HEAD"),
        # Record only the boolean state.  File names and repository remotes are
        # intentionally excluded from the portable manifest.
        "dirty": bool(status) if status is not None else None,
    }


def _solver_provenance(model, requested_solver: str) -> dict[str, object]:
    interface = getattr(getattr(model, "solver", None), "interface", None)
    interface_name = getattr(interface, "__name__", None)
    if not interface_name:
        interface_name = type(getattr(model, "solver", None)).__module__
    final_component = str(interface_name).rsplit(".", 1)[-1]
    actual_solver = final_component.removesuffix("_interface")
    backend_distributions = {
        "gurobi": "gurobipy",
        "glpk": "swiglpk",
        "glpk_exact": "swiglpk",
        "cplex": "cplex",
        "scipy": "scipy",
    }
    distribution = backend_distributions.get(actual_solver)
    return {
        "requested": str(requested_solver),
        "actual": actual_solver,
        "optlang_interface": str(interface_name),
        "backend_distribution": distribution,
        "backend_version": _distribution_version(distribution) if distribution else None,
    }


def build_run_manifest(
    *,
    model,
    model_path: Path,
    experimental_path: Path,
    media_path: Path,
    requested_solver: str,
    primary_cutoff: float,
    growth_cutoffs: tuple[float, ...],
    positive_only: bool,
    assay_fitness_path: Path | None = None,
    assay_source_sha256s: Iterable[str] = (),
    repo_root: Path = REPO_ROOT,
    generated_at: str | None = None,
    run_key: str | None = None,
    code_sources: dict[str, str] | None = None,
) -> dict[str, object]:
    """Build a credential-free provenance manifest for one validation run."""
    inputs: dict[str, object] = {
        "model": {
            "path": str(Path(model_path).resolve()),
            "sha256": sha256_file(Path(model_path)),
        },
        "experimental": {
            "path": str(Path(experimental_path).resolve()),
            "sha256": sha256_file(Path(experimental_path)),
        },
        "medium": {
            "path": str(Path(media_path).resolve()),
            "sha256": sha256_file(Path(media_path)),
        },
    }
    if assay_fitness_path is not None:
        inputs["assay_fitness"] = {
            "path": str(Path(assay_fitness_path).resolve()),
            "sha256": sha256_file(Path(assay_fitness_path)),
            "source_workbook_sha256": sorted(set(assay_source_sha256s)),
        }
    manifest = {
        "schema_version": "1.0",
        "generated_at": generated_at
        or datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "inputs": inputs,
        "solver": _solver_provenance(model, requested_solver),
        "software": {
            "python": sys.version.split()[0],
            "cobra": _distribution_version("cobra"),
            "pandas": _distribution_version("pandas"),
            "optlang": _distribution_version("optlang"),
        },
        "git": _git_provenance(Path(repo_root)),
        "cutoffs": {
            "primary_fraction_of_wt": float(primary_cutoff),
            "growth_fractions_of_wt": [float(value) for value in growth_cutoffs],
        },
        "configuration": {
            "positive_only": bool(positive_only),
            "assay_fitness_enabled": assay_fitness_path is not None,
            "credentials_included": False,
        },
    }
    if run_key is not None:
        manifest["run_key"] = run_key
    if code_sources is not None:
        manifest["code_sources"] = dict(code_sources)
    return manifest


def validate_essential_genes(
    experimental_path: Path,
    model_path: Path,
    media_path: Path,
    output_dir: Path,
    primary_cutoff: float = PRIMARY_CUTOFF,
    growth_cutoffs: tuple[float, ...] = DEFAULT_CUTOFFS,
    positive_only: bool = False,
    diagnose: bool = False,
    prepare_agent_cases: bool = False,
    batch_size: int = 3,
    case_ledger_path: Path = DEFAULT_LEDGER,
    evidence_dir: Path = DEFAULT_EVIDENCE_DIR,
    solver: str = "gurobi",
    max_bypass_candidates: int = 25,
    isozyme_capacity_scan_path: Path | None = None,
    assay_fitness_path: Path | None = None,
    case_category: str | None = None,
    run_key: str | None = None,
    code_sources: dict[str, str] | None = None,
) -> dict[str, object]:
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if case_category is not None and case_category not in CURATION_RULES:
        raise ValueError(f"Unknown essentiality case category: {case_category}")
    if case_category is not None and not prepare_agent_cases:
        raise ValueError("case_category requires prepare_agent_cases")
    if isozyme_capacity_scan_path is not None:
        if diagnose or prepare_agent_cases:
            raise ValueError(
                "Exploratory isozyme capacity scans cannot be combined with "
                "--diagnose or --prepare-agent-cases; those outputs must use the "
                "unmodified model signal"
            )
        if output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve():
            raise ValueError(
                "Exploratory isozyme capacity scans require a dedicated "
                "--output-dir so official baseline outputs are not overwritten"
            )
    growth_cutoffs = tuple(sorted(set(growth_cutoffs) | {primary_cutoff}))
    experimental = load_experimental(experimental_path, positive_only=positive_only)
    assay_fitness = (
        load_assay_fitness(assay_fitness_path)
        if assay_fitness_path is not None
        else None
    )
    media = load_media(media_path)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = read_sbml_model(str(model_path))
    active_medium = apply_media(model, media)
    predictions, wt_growth = run_single_gene_deletions(model, solver)
    per_gene = build_per_gene_table(experimental, predictions, growth_cutoffs, primary_cutoff)
    summary = make_summary(
        model,
        experimental,
        per_gene,
        wt_growth,
        growth_cutoffs,
        primary_cutoff,
        experimental_path,
        media_path,
        active_medium,
        positive_only,
    )
    model_sha256 = sha256_file(model_path)
    summary["model"].update(
        {
            "path": str(model_path.resolve()),
            "sha256": model_sha256,
        }
    )
    summary["experimental"]["sha256"] = sha256_file(experimental_path)
    media_sha256 = sha256_file(media_path)
    summary["medium"]["sha256"] = media_sha256

    output_dir.mkdir(parents=True, exist_ok=True)
    per_gene_path = output_dir / "essentiality_per_gene.tsv"
    summary_path = output_dir / "essentiality_summary.json"
    _write_tsv(per_gene, per_gene_path)

    if assay_fitness is not None and assay_fitness_path is not None:
        joined_assays = build_assay_fitness_table(
            assay_fitness,
            predictions,
            primary_cutoff,
        )
        assay_report_path = output_dir / "essentiality_assay_fitness.tsv"
        _write_tsv(joined_assays, assay_report_path)
        assay_summary = make_assay_fitness_summary(
            joined_assays,
            (gene.id for gene in model.genes),
            primary_cutoff,
            growth_cutoffs,
        )
        assay_summary["input"] = {
            "path": str(assay_fitness_path.resolve()),
            "sha256": sha256_file(assay_fitness_path),
            "source_workbook_sha256": sorted(
                assay_fitness["source_sha256"].unique().tolist()
            ),
        }
        assay_summary["report"] = str(assay_report_path.resolve())
        summary["assay_fitness"] = assay_summary

    if isozyme_capacity_scan_path is not None:
        capacity_scan = load_isozyme_capacity_scan(isozyme_capacity_scan_path)
        validate_isozyme_capacity_scan(
            model,
            capacity_scan,
            model_sha256,
            media_sha256,
        )
        capacity_results, capacity_summary = run_isozyme_capacity_scan(
            model,
            capacity_scan,
            predictions,
            experimental,
            wt_growth,
            growth_cutoffs,
            solver,
        )
        capacity_tsv = output_dir / "isozyme_capacity_sensitivity.tsv"
        capacity_json = output_dir / "isozyme_capacity_sensitivity_summary.json"
        _write_tsv(capacity_results, capacity_tsv)
        capacity_summary["input"] = {
            "path": str(isozyme_capacity_scan_path.resolve()),
            "sha256": sha256_file(isozyme_capacity_scan_path),
        }
        capacity_summary["model_sha256"] = model_sha256
        capacity_json.write_text(
            json.dumps(capacity_summary, indent=2, sort_keys=True) + "\n"
        )
        summary["exploratory_isozyme_capacity"] = {
            "enabled": True,
            "exploratory_only": True,
            "official_baseline_unchanged": True,
            "table": str(capacity_tsv),
            "summary": str(capacity_json),
        }

    if diagnose or prepare_agent_cases:
        diagnostics = diagnose_false_negatives(
            model,
            per_gene,
            wt_growth,
            primary_cutoff,
            max_bypass_candidates,
        )
        queue = build_curation_queue(diagnostics)
        if case_category is not None:
            queue = queue.loc[queue["category"].eq(case_category)].copy()
        _write_tsv(diagnostics, output_dir / "essentiality_fn_diagnosis.tsv")
        _write_tsv(queue, output_dir / "essentiality_curation_queue.tsv")
        summary["diagnosis_counts"] = {
            str(category): int(count)
            for category, count in diagnostics["category"].value_counts().items()
        }
        if case_category is None:
            translation, translation_summary = audit_translation_module(model)
            _write_tsv(translation, output_dir / "translation_module_audit.tsv")
            summary["translation_module"] = translation_summary
        cases = build_agent_cases(
            model,
            diagnostics,
            per_gene,
            summary,
            primary_cutoff,
        )
        if prepare_agent_cases:
            summary["agent_cases"] = prepare_agent_case_files(
                cases,
                output_dir,
                batch_size,
                ledger_path=case_ledger_path,
                evidence_dir=evidence_dir,
                case_category=case_category,
                require_full_batch=case_category is not None,
            )
        if case_category in {None, "isozyme_redundancy"}:
            resolution_cases = [
                case
                for case in cases
                if case_category is None or case["category"] == case_category
            ]
            isozyme_resolution = build_isozyme_resolution_ledger(
                resolution_cases,
                primary_cutoff=primary_cutoff,
                ledger_path=case_ledger_path,
                evidence_dir=evidence_dir,
            )
            isozyme_resolution_path = output_dir / "essentiality_isozyme_resolution.tsv"
            _write_tsv(isozyme_resolution, isozyme_resolution_path)
            causal_counts = isozyme_resolution["model_causal_class"].value_counts()
            resolution_counts = isozyme_resolution["resolution_status"].value_counts()
            summary["isozyme_resolution"] = {
                "table": str(isozyme_resolution_path.resolve()),
                "case_groups": int(len(isozyme_resolution)),
                "model_causal_case_groups": int(
                    causal_counts.get("model_causal_isozyme_candidate", 0)
                ),
                "noncausal_case_groups": int(
                    causal_counts.get("noncausal_redundancy_signal", 0)
                ),
                "causal_class_counts": {
                    str(label): int(count) for label, count in causal_counts.items()
                },
                "resolution_status_counts": {
                    str(label): int(count) for label, count in resolution_counts.items()
                },
            }

    manifest = build_run_manifest(
        model=model,
        model_path=model_path,
        experimental_path=experimental_path,
        media_path=media_path,
        requested_solver=solver,
        primary_cutoff=primary_cutoff,
        growth_cutoffs=growth_cutoffs,
        positive_only=positive_only,
        assay_fitness_path=assay_fitness_path,
        assay_source_sha256s=(
            assay_fitness["source_sha256"].unique().tolist()
            if assay_fitness is not None
            else ()
        ),
        run_key=run_key,
        code_sources=code_sources,
    )
    manifest_path = output_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    summary["run_manifest"] = str(manifest_path.resolve())
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    logger.info("Per-gene report: %s", per_gene_path)
    logger.info("Summary: %s", summary_path)
    logger.info("Run manifest: %s", manifest_path)
    return summary


def print_summary(summary: dict[str, object]) -> None:
    coverage = summary["coverage"]
    primary = summary["primary"]
    print("\nGene Essentiality Validation (positive-only reference)")
    print("=" * 62)
    print(f"WT growth:             {summary['wt_growth']:.6f} h^-1")
    print(f"Experimental positives:{coverage['experimental_genes']:>8}")
    print(f"In-model intersection: {coverage['intersection']:>8}")
    print(f"Outside GEM scope:      {coverage['outside_gem_scope']:>8}")
    print(f"Primary cutoff:         {summary['primary_cutoff_fraction_of_wt']:.0%} of WT")
    print(f"TP / FN:                {primary['TP']:>4} / {primary['FN']:<4}")
    print(f"Recall:                 {primary['recall']:.3f}")
    print("Cutoff curve:")
    for item in summary["cutoff_curve"]:
        print(
            f"  {item['cutoff_fraction_of_wt']:>5.0%}: "
            f"TP={item['TP']:>3} FN={item['FN']:>3} recall={item['recall']:.3f}"
        )
    print("Accuracy/MCC/precision omitted: the reference contains no negative controls.")
    assay_summary = summary.get("assay_fitness")
    if assay_summary:
        print("\nSeparate assay-fitness diagnostics:")
        for assay, metrics in assay_summary["per_assay"].items():
            spearman = metrics["spearman_fitness_vs_ko_growth_ratio"]
            r_squared = metrics["linear_r_squared_fitness_vs_ko_growth_ratio"]
            auc = metrics["assay_call_roc_auc"]
            spearman_text = f"{spearman:.3f}" if spearman is not None else "n/a"
            r_squared_text = f"{r_squared:.3f}" if r_squared is not None else "n/a"
            auc_text = f"{auc:.3f}" if auc is not None else "n/a"
            print(
                f"  {assay}: overlap={metrics['model_overlap']}/"
                f"{metrics['assay_unique_genes']} Spearman={spearman_text} "
                f"R^2={r_squared_text} ROC_AUC={auc_text}"
            )
        proxy = assay_summary["concordant_nonessential_proxy_safety"]["primary"]
        print(
            "  Concordant-nonessential proxy safety at primary cutoff: "
            f"{proxy['safe_nonessential_count']}/{proxy['model_proxy_genes']} "
            f"({proxy['safety_rate']:.3f})"
            if proxy["safety_rate"] is not None
            else "  Concordant-nonessential proxy safety: n/a"
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate iYali26 gene essentiality against positive-only or labelled experiments."
    )
    parser.add_argument(
        "--research-root",
        type=Path,
        help="External research workspace; overrides IYALI26_RESEARCH_ROOT",
    )
    parser.add_argument("--experimental", "-e", type=Path)
    parser.add_argument("--model", "-m", type=Path, default=REPO_ROOT / "model.xml")
    parser.add_argument("--media", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--force-rerun",
        action="store_true",
        help="reproduce an identical successful run; requires --reproduction-reason",
    )
    parser.add_argument(
        "--reproduction-reason",
        help="why a matching successful result must be reproduced",
    )
    parser.add_argument("--positive-only", action="store_true")
    parser.add_argument(
        "--assay-fitness",
        type=Path,
        help=(
            "Normalized Cas9/Cas12a fitness CSV. Metrics are reported separately "
            "and do not change positive-only TP/FN/recall."
        ),
    )
    parser.add_argument("--primary-cutoff", type=float, default=PRIMARY_CUTOFF)
    parser.add_argument(
        "--growth-cutoffs",
        type=_parse_cutoffs,
        default=DEFAULT_CUTOFFS,
        help="Comma-separated fractions of WT (default: 0.01,0.05,0.10,0.15)",
    )
    parser.add_argument("--diagnose", action="store_true")
    parser.add_argument(
        "--prepare-agent-cases",
        action="store_true",
        help="Recompute FN diagnostics and prepare deterministic literature-review cases",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=3,
        help="Number of queued cases placed in the next agent batch (default: 3)",
    )
    parser.add_argument(
        "--case-category",
        choices=sorted(CURATION_RULES),
        help=(
            "Restrict durable case preparation to one diagnosis category. "
            "Requires --prepare-agent-cases."
        ),
    )
    parser.add_argument("--solver", default="gurobi")
    parser.add_argument(
        "--max-bypass-candidates",
        type=int,
        default=25,
        help="Maximum high-flux-delta candidates tested per metabolic-bypass FN",
    )
    parser.add_argument(
        "--isozyme-capacity-scan",
        type=Path,
        help=(
            "Run an exploratory KO-specific residual-isozyme capacity scan. "
            "Requires a dedicated --output-dir and cannot be combined with "
            "diagnostics or agent cases."
        ),
    )
    args = parser.parse_args()

    if args.force_rerun and not args.reproduction_reason:
        parser.error("--force-rerun requires --reproduction-reason")

    try:
        project_paths = load_project_paths(args.research_root, required=True)
        project_paths.require(project_paths.essentiality, project_paths.media)
    except (FileNotFoundError, RuntimeError) as exc:
        parser.error(str(exc))
    if args.experimental is None:
        args.experimental = (
            project_paths.essentiality / "consensus_essential_genes.csv"
        )
    if args.media is None:
        args.media = project_paths.media / "sd_leu.csv"
    for path, label in (
        (args.experimental, "experimental file"),
        (args.model, "model"),
        (args.media, "media file"),
    ):
        if not path.exists():
            parser.error(f"{label} not found: {path}")
    if not (0 < args.primary_cutoff < 1):
        parser.error("--primary-cutoff must be strictly between 0 and 1")
    if args.max_bypass_candidates < 0:
        parser.error("--max-bypass-candidates cannot be negative")
    if args.batch_size <= 0:
        parser.error("--batch-size must be positive")
    if args.case_category is not None and not args.prepare_agent_cases:
        parser.error("--case-category requires --prepare-agent-cases")
    if args.isozyme_capacity_scan is not None and not args.isozyme_capacity_scan.exists():
        parser.error(f"isozyme capacity scan not found: {args.isozyme_capacity_scan}")
    if args.assay_fitness is not None and not args.assay_fitness.exists():
        parser.error(f"assay fitness file not found: {args.assay_fitness}")

    inputs: dict[str, dict[str, str]] = {
        "model": {"path": str(args.model.resolve()), "sha256": sha256_file(args.model)},
        "experimental": {
            "path": str(args.experimental.resolve()),
            "sha256": sha256_file(args.experimental),
        },
        "medium": {"path": str(args.media.resolve()), "sha256": sha256_file(args.media)},
    }
    if args.assay_fitness is not None:
        inputs["assay_fitness"] = {
            "path": str(args.assay_fitness.resolve()),
            "sha256": sha256_file(args.assay_fitness),
        }
    code_sources = {
        str(path.resolve()): sha256_file(path)
        for path in (
            Path(__file__),
            Path(__file__).with_name("essentiality_evidence.py"),
            Path(__file__).with_name("patches.py"),
        )
        if path.is_file()
    }
    configuration = {
        "positive_only": bool(args.positive_only),
        "primary_cutoff": float(args.primary_cutoff),
        "growth_cutoffs": [float(value) for value in args.growth_cutoffs],
        "solver": args.solver,
        "diagnose": bool(args.diagnose),
        "prepare_agent_cases": bool(args.prepare_agent_cases),
        "case_category": args.case_category,
    }
    run_key = build_run_key(
        "essentiality",
        inputs=inputs,
        code_sources=code_sources,
        configuration=configuration,
    )
    if args.output_dir is None:
        args.output_dir = (
            project_paths.results
            / "essentiality"
            / f"{run_key[:12]}-{registry_utc_now().replace(':', '').replace('+00:00', 'Z')}"
        )
    if args.output_dir.exists():
        parser.error(
            "--output-dir already exists; use a new directory to preserve every prior result"
        )
    try:
        previous_run = guard_duplicate_run(
            project_paths.research_root,
            workflow="essentiality",
            run_key=run_key,
            output_dir=args.output_dir,
            force_rerun=args.force_rerun,
            reproduction_reason=args.reproduction_reason,
        )
    except (ValueError, RuntimeError) as exc:
        parser.error(str(exc))

    try:
        summary = validate_essential_genes(
            experimental_path=args.experimental,
            model_path=args.model,
            media_path=args.media,
            output_dir=args.output_dir,
            primary_cutoff=args.primary_cutoff,
            growth_cutoffs=args.growth_cutoffs,
            positive_only=args.positive_only,
            diagnose=args.diagnose,
            prepare_agent_cases=args.prepare_agent_cases,
            batch_size=args.batch_size,
            solver=args.solver,
            max_bypass_candidates=args.max_bypass_candidates,
            isozyme_capacity_scan_path=args.isozyme_capacity_scan,
            assay_fitness_path=args.assay_fitness,
            case_category=args.case_category,
            case_ledger_path=project_paths.essentiality / "curation_cases.csv",
            evidence_dir=project_paths.essentiality / "evidence",
            run_key=run_key,
            code_sources=code_sources,
        )
    except (ValueError, RuntimeError) as exc:
        logger.error("%s", exc)
        sys.exit(2)
    register_run(
        project_paths.research_root,
        workflow="essentiality",
        run_key=run_key,
        output_dir=args.output_dir,
        manifest_path=Path(str(summary["run_manifest"])),
        inputs=inputs,
        code_sources=code_sources,
        configuration=configuration,
        status="complete",
        previous=previous_run,
        reproduction_reason=args.reproduction_reason,
    )
    print_summary(summary)


if __name__ == "__main__":
    main()
