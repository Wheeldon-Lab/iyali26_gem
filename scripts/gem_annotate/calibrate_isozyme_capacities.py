"""Deterministic joint calibration of the provisional R4/R1846 capacities.

This command evaluates a small, pre-approved Cartesian grid.  Each scenario
globally partitions both OR-isozyme reaction families, runs every single-gene
knockout with one worker, and applies conservative collateral-phenotype gates.
It never writes an SBML model and never changes the canonical capacity profile.

The capacity fractions in the grid are relative to the corresponding reaction
flux in the unmodified SD-Leu wild type.  They are exploratory parameters, not
measured enzyme capacities.  A selected profile therefore retains the existing
proteomics/kcat replacement warnings and is written only as a separate output.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import tempfile
import warnings
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from cobra.core.gene import GPR
from cobra.io import read_sbml_model

from .config import ESSENTIALITY_DIR, MEDIA_DIR, REPO_ROOT, RESULTS_DIR
from .essentiality_evidence import canonical_json, sha256_file, utc_now
from .provisional_capacity import (
    REQUIRED_COLUMNS as PROFILE_REQUIRED_COLUMNS,
    apply_provisional_isozyme_capacities,
    load_provisional_capacity_table,
)
from .validate_essential_genes import (
    DEFAULT_CUTOFFS,
    apply_media,
    build_assay_fitness_table,
    load_experimental,
    load_media,
    make_assay_fitness_summary,
    run_single_gene_deletions,
)


DEFAULT_MODEL = REPO_ROOT / "model.xml"
DEFAULT_PROFILE = ESSENTIALITY_DIR / "provisional_isozyme_capacities.csv"
DEFAULT_GRID = (
    ESSENTIALITY_DIR / "scenarios" / "provisional_capacity_joint_grid.csv"
)
DEFAULT_EXPERIMENTAL = ESSENTIALITY_DIR / "consensus_essential_genes.csv"
DEFAULT_ASSAY_FITNESS = ESSENTIALITY_DIR / "cas9_cas12a_fitness.csv"
DEFAULT_MEDIA = MEDIA_DIR / "sd_leu.csv"
DEFAULT_OUTPUT_DIR = RESULTS_DIR / "essentiality" / "capacity_joint_grid"

GRID_REQUIRED_COLUMNS = {
    "scenario_id",
    "r4_capacity_fraction_of_wt_flux",
    "r1846_capacity_fraction_of_wt_flux",
    "exploratory_only",
    "basis",
    "rationale",
}
EXPECTED_R4_FRACTIONS = (0.025, 0.075, 0.15)
EXPECTED_R1846_FRACTIONS = (0.010, 0.025, 0.050)
TARGET_REACTIONS = ("R4", "R1846")
BLOCKED_CONTROL_REACTION = "R1843"
TARGET_LOWER_BOUND = 0.05
TARGET_UPPER_BOUND = 0.10
WT_DELTA_TOLERANCE = 1e-8
BACKUP_KO_RATIO_MINIMUM = 0.90
WORKERS = 1
RATIO_CHANGE_EPS = 1e-9
RESULT_DECIMALS = 9
SCHEMA_VERSION = 1

GRID_SUMMARY_NAME = "joint_capacity_grid.tsv"
COLLATERAL_NAME = "collateral_changes.tsv"
ASSAY_METRICS_NAME = "assay_metrics.tsv"
MANIFEST_NAME = "run_manifest.json"
CANDIDATE_PROFILE_NAME = "candidate_experimental_profile.csv"
SCENARIO_TABLE_NAME = "essentiality_per_gene.tsv"
SCENARIO_SUMMARY_NAME = "scenario_summary.json"


def _parse_bool(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Unrecognised Boolean value: {value!r}")


def _parse_growth_cutoffs(value: str) -> tuple[float, ...]:
    """Parse a deterministic, strictly increasing cutoff list."""
    try:
        cutoffs = tuple(float(item.strip()) for item in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "growth cutoffs must be comma-separated numbers"
        ) from exc
    if not cutoffs or any(not 0.0 < cutoff < 1.0 for cutoff in cutoffs):
        raise argparse.ArgumentTypeError(
            "growth cutoffs must be fractions strictly between 0 and 1"
        )
    if tuple(sorted(set(cutoffs))) != cutoffs:
        raise argparse.ArgumentTypeError(
            "growth cutoffs must be unique and strictly increasing"
        )
    return cutoffs


def _fraction_label(value: float) -> str:
    return f"{value:.3f}".replace(".", "p")


def expected_scenario_id(r4_fraction: float, r1846_fraction: float) -> str:
    """Return the canonical ID for one approved grid point."""
    return f"R4_{_fraction_label(r4_fraction)}__R1846_{_fraction_label(r1846_fraction)}"


def load_joint_capacity_grid(path: Path) -> pd.DataFrame:
    """Load and strictly validate the approved nine-point Cartesian grid."""
    grid = pd.read_csv(path, dtype=str).fillna("")
    missing = GRID_REQUIRED_COLUMNS - set(grid.columns)
    if missing:
        raise ValueError(f"Joint capacity grid is missing columns {sorted(missing)}")
    if grid.empty:
        raise ValueError(f"Joint capacity grid is empty: {path}")

    grid["scenario_id"] = grid["scenario_id"].str.strip()
    if grid["scenario_id"].eq("").any():
        raise ValueError("Joint capacity grid has an empty scenario_id")
    if grid["scenario_id"].duplicated().any():
        repeated = sorted(grid.loc[grid["scenario_id"].duplicated(), "scenario_id"])
        raise ValueError(f"Duplicate joint capacity scenario IDs: {repeated}")

    fraction_columns = (
        "r4_capacity_fraction_of_wt_flux",
        "r1846_capacity_fraction_of_wt_flux",
    )
    for column in fraction_columns:
        grid[column] = pd.to_numeric(grid[column], errors="raise")
        if (~grid[column].map(math.isfinite)).any() or (grid[column] <= 0).any():
            raise ValueError(f"{column} values must be finite and positive")

    grid["exploratory_only"] = grid["exploratory_only"].map(_parse_bool)
    if not grid["exploratory_only"].all():
        raise ValueError("Every joint capacity scenario must be exploratory_only=true")
    for column in ("basis", "rationale"):
        grid[column] = grid[column].str.strip()
        if grid[column].eq("").any():
            raise ValueError(f"Joint capacity grid has an empty {column}")

    expected_points = {
        (r4, r1846)
        for r4 in EXPECTED_R4_FRACTIONS
        for r1846 in EXPECTED_R1846_FRACTIONS
    }
    actual_points = {
        (
            float(row.r4_capacity_fraction_of_wt_flux),
            float(row.r1846_capacity_fraction_of_wt_flux),
        )
        for row in grid.itertuples(index=False)
    }
    if actual_points != expected_points or len(grid) != len(expected_points):
        raise ValueError(
            "Joint capacity grid must be exactly the approved 3x3 Cartesian "
            f"product; found={sorted(actual_points)}"
        )
    for row in grid.itertuples(index=False):
        expected = expected_scenario_id(
            float(row.r4_capacity_fraction_of_wt_flux),
            float(row.r1846_capacity_fraction_of_wt_flux),
        )
        if row.scenario_id != expected:
            raise ValueError(
                f"Non-canonical scenario_id {row.scenario_id!r}; expected {expected!r}"
            )

    return grid.sort_values("scenario_id").reset_index(drop=True)


def _load_assay_fitness(path: Path) -> pd.DataFrame:
    """Use the normalized assay loader while keeping this module importable alone."""
    try:
        from .validate_essential_genes import load_assay_fitness
    except ImportError as exc:  # pragma: no cover - protects partial installations
        raise RuntimeError(
            "The normalized Cas9/Cas12a assay loader is unavailable"
        ) from exc
    return load_assay_fitness(path)


def _normalized_assay_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def concordant_nonessential_proxies(assay_fitness: pd.DataFrame) -> set[str]:
    """Genes called nonessential by both Cas9 and Cas12a screens."""
    required = {"gene_id", "assay", "experimental_call"}
    missing = required - set(assay_fitness.columns)
    if missing:
        raise ValueError(f"Assay fitness table is missing columns {sorted(missing)}")

    table = assay_fitness.loc[:, sorted(required)].copy()
    table["gene_id"] = table["gene_id"].astype(str).str.strip()
    table["assay_key"] = table["assay"].map(_normalized_assay_name)
    table["call_key"] = table["experimental_call"].astype(str).str.strip().str.lower()
    required_assays = {"cas9", "cas12a"}
    present_assays = set(table["assay_key"])
    if not required_assays <= present_assays:
        raise ValueError(
            "Assay fitness input must contain both Cas9 and Cas12a; "
            f"found={sorted(present_assays)}"
        )

    proxies: set[str] = set()
    for gene_id, group in table.groupby("gene_id", sort=True):
        calls_by_assay = {
            assay: set(rows["call_key"])
            for assay, rows in group.groupby("assay_key", sort=False)
        }
        if all(
            calls_by_assay.get(assay) == {"nonessential"} for assay in required_assays
        ):
            proxies.add(str(gene_id))
    return proxies


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _atomic_write_bytes(path: Path, payload: bytes) -> None:
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
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def _atomic_write_json(path: Path, value: Any) -> None:
    _atomic_write_bytes(path, _json_text(value).encode("utf-8"))


def _dataframe_tsv_bytes(frame: pd.DataFrame) -> bytes:
    buffer = io.StringIO(newline="")
    frame.to_csv(buffer, sep="\t", index=False, na_rep="", lineterminator="\n")
    return buffer.getvalue().encode("utf-8")


def _atomic_write_tsv(path: Path, frame: pd.DataFrame) -> None:
    _atomic_write_bytes(path, _dataframe_tsv_bytes(frame))


def _profile_csv_bytes(
    rows: Iterable[dict[str, Any]],
    fieldnames: list[str],
) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        serialized = dict(row)
        for key, value in serialized.items():
            if isinstance(value, bool):
                serialized[key] = str(value).lower()
        writer.writerow({field: serialized.get(field, "") for field in fieldnames})
    return buffer.getvalue().encode("utf-8")


def _profile_fieldnames(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        fieldnames = list(csv.DictReader(handle).fieldnames or [])
    missing = PROFILE_REQUIRED_COLUMNS - set(fieldnames)
    if missing:
        raise ValueError(f"Capacity profile is missing columns {sorted(missing)}")
    return fieldnames


def _profile_rows_by_reaction(
    profile_rows: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    by_reaction = {str(row["source_reaction_id"]): row for row in profile_rows}
    if set(by_reaction) != set(TARGET_REACTIONS) or len(profile_rows) != 2:
        raise ValueError(
            "Joint calibration requires exactly the active R4 and R1846 profile rows; "
            f"found={sorted(by_reaction)}"
        )
    return by_reaction


def _genes(rule: str) -> set[str]:
    return set(GPR.from_string(rule).genes) if rule.strip() else set()


def _profile_gene_roles(
    profile_rows: list[dict[str, Any]],
) -> tuple[dict[str, str], set[str]]:
    target_by_reaction: dict[str, str] = {}
    backup_genes: set[str] = set()
    for row in profile_rows:
        primary = _genes(str(row["primary_gpr"]))
        backup = _genes(str(row["backup_gpr"]))
        if len(primary) != 1:
            raise ValueError(
                f"{row['capacity_id']} must have exactly one primary target gene"
            )
        target_by_reaction[str(row["source_reaction_id"])] = next(iter(primary))
        backup_genes.update(backup)
    return target_by_reaction, backup_genes


def _input_manifest(
    *,
    model_path: Path,
    profile_path: Path,
    grid_path: Path,
    experimental_path: Path,
    assay_fitness_path: Path,
    media_path: Path,
) -> dict[str, dict[str, str]]:
    inputs: dict[str, dict[str, str]] = {}
    for name, path in (
        ("model", model_path),
        ("profile", profile_path),
        ("grid", grid_path),
        ("experimental", experimental_path),
        ("assay_fitness", assay_fitness_path),
        ("media", media_path),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{name} input not found: {path}")
        inputs[name] = {
            "path": str(path.resolve()),
            "sha256": sha256_file(path),
        }
    return inputs


def _run_configuration(solver: str, cutoffs: tuple[float, ...]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "solver": solver,
        "workers": WORKERS,
        "cutoffs": list(cutoffs),
        "target_ratio_window": {
            "lower_inclusive": TARGET_LOWER_BOUND,
            "upper_exclusive": TARGET_UPPER_BOUND,
        },
        "wt_growth_delta_tolerance": WT_DELTA_TOLERANCE,
        "backup_gene_ko_ratio_strict_minimum": BACKUP_KO_RATIO_MINIMUM,
        "ratio_change_epsilon": RATIO_CHANGE_EPS,
        "blocked_control_reaction": BLOCKED_CONTROL_REACTION,
        "target_reactions": list(TARGET_REACTIONS),
        "global_capacity_split": True,
        "single_gene_knockout_scope": "all_model_genes",
    }


def _run_key(inputs: dict[str, Any], configuration: dict[str, Any]) -> str:
    return _sha256_text(
        canonical_json({"inputs": inputs, "configuration": configuration})
    )


def _scenario_configuration(row: Any) -> dict[str, Any]:
    return {
        "scenario_id": str(row.scenario_id),
        "r4_capacity_fraction_of_wt_flux": float(row.r4_capacity_fraction_of_wt_flux),
        "r1846_capacity_fraction_of_wt_flux": float(
            row.r1846_capacity_fraction_of_wt_flux
        ),
        "exploratory_only": bool(row.exploratory_only),
        "basis": str(row.basis),
        "rationale": str(row.rationale),
    }


def _scenario_key(run_key: str, scenario: dict[str, Any]) -> str:
    return _sha256_text(canonical_json({"run_key": run_key, "scenario": scenario}))


def _configure_solver(model, solver: str) -> None:
    model.solver = solver
    model.solver.configuration.verbosity = 0
    if solver.lower() == "gurobi":
        # COBRA's ``processes=1`` controls deletion workers.  This additionally
        # prevents Gurobi from introducing internal thread-level variability.
        model.solver.problem.Params.Threads = 1


def _load_model(path: Path, solver: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = read_sbml_model(str(path))
    _configure_solver(model, solver)
    return model


def _optimal_growth_and_fluxes(model) -> tuple[float, pd.Series]:
    solution = model.optimize()
    if solution.status != "optimal" or solution.objective_value is None:
        raise RuntimeError(f"Wild-type FBA is not optimal: {solution.status}")
    growth = float(solution.objective_value)
    if not (0.1 <= growth <= 2.0):
        raise RuntimeError(
            f"Wild-type growth {growth:.6g} h^-1 is outside the accepted "
            "0.1-2.0 h^-1 SD-Leu range"
        )
    return growth, solution.fluxes


def _scenario_profile_rows(
    profile_rows: list[dict[str, Any]],
    scenario: dict[str, Any],
    wt_fluxes: pd.Series,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    fractions = {
        "R4": float(scenario["r4_capacity_fraction_of_wt_flux"]),
        "R1846": float(scenario["r1846_capacity_fraction_of_wt_flux"]),
    }
    rows: list[dict[str, Any]] = []
    caps: dict[str, float] = {}
    for original in profile_rows:
        row = dict(original)
        reaction_id = str(row["source_reaction_id"])
        wt_flux = abs(float(wt_fluxes[reaction_id]))
        if wt_flux <= RATIO_CHANGE_EPS:
            raise ValueError(
                f"{reaction_id} has zero WT flux; a WT-flux-relative capacity is undefined"
            )
        cap = wt_flux * fractions[reaction_id]
        row["provisional_upper_bound"] = cap
        rows.append(row)
        caps[reaction_id] = cap
    return rows, caps


def _apply_scenario_profile(
    model,
    profile_rows: list[dict[str, Any]],
    fieldnames: list[str],
    reference_model_sha256: str,
    directory: Path,
) -> list[dict[str, Any]]:
    """Apply scenario rows through the production profile validator."""
    payload = _profile_csv_bytes(profile_rows, fieldnames)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=directory,
        prefix=".scenario-profile.",
        suffix=".csv",
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
        return apply_provisional_isozyme_capacities(
            model,
            temporary_path,
            reference_model_sha256=reference_model_sha256,
        )
    finally:
        temporary_path.unlink(missing_ok=True)


def _is_essential(ratio: float, cutoff: float) -> bool:
    return float(ratio) < float(cutoff)


def _build_scenario_per_gene(
    baseline: pd.DataFrame,
    scenario_predictions: pd.DataFrame,
    *,
    target_genes: set[str],
    backup_genes: set[str],
    proxy_genes: set[str],
    experimental_positive_genes: set[str],
    cutoffs: tuple[float, ...],
) -> pd.DataFrame:
    baseline_indexed = baseline.set_index("gene_id")
    scenario_indexed = scenario_predictions.set_index("gene_id")
    if set(baseline_indexed.index) != set(scenario_indexed.index):
        raise RuntimeError("Scenario and baseline model gene sets differ")

    rows: list[dict[str, Any]] = []
    for gene_id in sorted(baseline_indexed.index.astype(str)):
        base = baseline_indexed.loc[gene_id]
        current = scenario_indexed.loc[gene_id]
        baseline_ratio = float(base["ko_growth_ratio"])
        scenario_ratio = float(current["ko_growth_ratio"])
        stable_baseline_growth = round(float(base["ko_growth"]), RESULT_DECIMALS)
        stable_scenario_growth = round(float(current["ko_growth"]), RESULT_DECIMALS)
        stable_baseline_ratio = round(baseline_ratio, RESULT_DECIMALS)
        stable_scenario_ratio = round(scenario_ratio, RESULT_DECIMALS)
        stable_delta = round(
            stable_scenario_ratio - stable_baseline_ratio, RESULT_DECIMALS
        )
        record: dict[str, Any] = {
            "gene_id": gene_id,
            "baseline_ko_status": str(base["ko_status"]),
            "baseline_ko_growth": stable_baseline_growth,
            "baseline_ko_growth_ratio": stable_baseline_ratio,
            "scenario_ko_status": str(current["ko_status"]),
            "scenario_ko_growth": stable_scenario_growth,
            "scenario_ko_growth_ratio": stable_scenario_ratio,
            "ko_growth_ratio_delta": stable_delta,
            "abs_ko_growth_ratio_delta": abs(stable_delta),
            "is_target_gene": gene_id in target_genes,
            "is_backup_gene": gene_id in backup_genes,
            "is_concordant_nonessential_proxy": gene_id in proxy_genes,
            "experimental_positive": gene_id in experimental_positive_genes,
        }
        for cutoff in cutoffs:
            label = f"{cutoff * 100:g}pct".replace(".", "p")
            baseline_call = _is_essential(baseline_ratio, cutoff)
            scenario_call = _is_essential(scenario_ratio, cutoff)
            record[f"baseline_essential_at_{label}"] = baseline_call
            record[f"scenario_essential_at_{label}"] = scenario_call
            record[f"call_flip_at_{label}"] = baseline_call != scenario_call
            record[f"new_essential_at_{label}"] = scenario_call and not baseline_call
        rows.append(record)
    return pd.DataFrame(rows)


def _scenario_gates(
    per_gene: pd.DataFrame,
    *,
    target_by_reaction: dict[str, str],
    backup_genes: set[str],
    proxy_genes: set[str],
    baseline_wt_growth: float,
    scenario_wt_growth: float,
    r1843_bounds: tuple[float, float],
    family_totals: dict[str, tuple[float, float]],
    cutoffs: tuple[float, ...],
) -> dict[str, Any]:
    indexed = per_gene.set_index("gene_id")
    target_genes = set(target_by_reaction.values())
    target_ratios = {
        reaction_id: float(indexed.loc[gene_id, "scenario_ko_growth_ratio"])
        for reaction_id, gene_id in target_by_reaction.items()
    }
    target_window_pass = all(
        TARGET_LOWER_BOUND <= ratio < TARGET_UPPER_BOUND
        for ratio in target_ratios.values()
    )
    raw_wt_delta = abs(scenario_wt_growth - baseline_wt_growth)
    wt_delta = (
        0.0
        if raw_wt_delta <= RATIO_CHANGE_EPS
        else round(raw_wt_delta, RESULT_DECIMALS)
    )
    wt_pass = raw_wt_delta <= WT_DELTA_TOLERANCE
    r1843_pass = tuple(float(value) for value in r1843_bounds) == (0.0, 0.0)
    family_total_deltas = {}
    for reaction_id, (baseline, current) in family_totals.items():
        raw_delta = abs(float(current) - float(baseline))
        family_total_deltas[reaction_id] = (
            0.0
            if raw_delta <= RATIO_CHANGE_EPS
            else round(raw_delta, RESULT_DECIMALS)
        )
    family_total_pass = all(delta <= 1e-12 for delta in family_total_deltas.values())

    non_targets = per_gene.loc[~per_gene["gene_id"].isin(target_genes)]
    flip_columns = [
        f"call_flip_at_{cutoff * 100:g}pct".replace(".", "p") for cutoff in cutoffs
    ]
    flip_mask = non_targets[flip_columns].any(axis=1)
    non_target_flip_genes = sorted(non_targets.loc[flip_mask, "gene_id"].astype(str))
    no_non_target_flips = not non_target_flip_genes

    represented_backup = sorted(backup_genes & set(indexed.index.astype(str)))
    if represented_backup != sorted(backup_genes):
        missing = sorted(backup_genes - set(represented_backup))
        raise RuntimeError(f"Backup genes missing from deletion output: {missing}")
    backup_ratios = {
        gene_id: float(indexed.loc[gene_id, "scenario_ko_growth_ratio"])
        for gene_id in represented_backup
    }
    backup_pass = bool(backup_ratios) and all(
        ratio > BACKUP_KO_RATIO_MINIMUM for ratio in backup_ratios.values()
    )

    represented_proxies = proxy_genes & set(indexed.index.astype(str))
    proxy_new_essential_genes = sorted(
        gene_id
        for gene_id in represented_proxies
        if (
            float(indexed.loc[gene_id, "scenario_ko_growth_ratio"])
            < TARGET_UPPER_BOUND
            <= float(indexed.loc[gene_id, "baseline_ko_growth_ratio"])
        )
    )
    proxy_pass = not proxy_new_essential_genes

    boundary_margin = round(
        min(
            min(ratio - TARGET_LOWER_BOUND, TARGET_UPPER_BOUND - ratio)
            for ratio in target_ratios.values()
        ),
        RESULT_DECIMALS,
    )
    meaningful_changes = non_targets["abs_ko_growth_ratio_delta"].where(
        non_targets["abs_ko_growth_ratio_delta"] > RATIO_CHANGE_EPS,
        0.0,
    )
    non_target_change = round(float(meaningful_changes.max()), RESULT_DECIMALS)
    non_target_change_sum = round(float(meaningful_changes.sum()), RESULT_DECIMALS)
    feasible = all(
        (
            target_window_pass,
            wt_pass,
            r1843_pass,
            family_total_pass,
            no_non_target_flips,
            backup_pass,
            proxy_pass,
        )
    )
    return {
        "feasible": feasible,
        "target_window_pass": target_window_pass,
        "target_ko_ratios": target_ratios,
        "boundary_margin": boundary_margin,
        "wt_growth_delta": wt_delta,
        "wt_growth_gate_pass": wt_pass,
        "r1843_bounds": list(r1843_bounds),
        "r1843_bounds_gate_pass": r1843_pass,
        "reaction_family_total_upper_bound_deltas": family_total_deltas,
        "reaction_family_total_upper_bound_gate_pass": family_total_pass,
        "non_target_call_flip_count": len(non_target_flip_genes),
        "non_target_call_flip_genes": non_target_flip_genes,
        "no_non_target_call_flips_gate_pass": no_non_target_flips,
        "backup_gene_ko_ratios": backup_ratios,
        "backup_gene_ko_ratio_gate_pass": backup_pass,
        "concordant_nonessential_proxy_count": len(represented_proxies),
        "proxy_new_essential_count": len(proxy_new_essential_genes),
        "proxy_new_essential_genes": proxy_new_essential_genes,
        "concordant_nonessential_proxy_gate_pass": proxy_pass,
        "non_target_ratio_change_score": non_target_change,
        "non_target_ratio_change_sum": non_target_change_sum,
    }


def _assay_metrics(
    assay_fitness: pd.DataFrame,
    per_gene: pd.DataFrame,
    scenario_id: str,
) -> list[dict[str, Any]]:
    model_gene_ids = set(per_gene["gene_id"].astype(str))
    scenario_predictions = per_gene.rename(
        columns={
            "scenario_ko_status": "ko_status",
            "scenario_ko_growth": "ko_growth",
            "scenario_ko_growth_ratio": "ko_growth_ratio",
        }
    )[["gene_id", "ko_status", "ko_growth", "ko_growth_ratio"]]
    baseline_predictions = per_gene.rename(
        columns={
            "baseline_ko_status": "ko_status",
            "baseline_ko_growth": "ko_growth",
            "baseline_ko_growth_ratio": "ko_growth_ratio",
        }
    )[["gene_id", "ko_status", "ko_growth", "ko_growth_ratio"]]
    scenario_summary = make_assay_fitness_summary(
        build_assay_fitness_table(
            assay_fitness, scenario_predictions, TARGET_UPPER_BOUND
        ),
        model_gene_ids,
        TARGET_UPPER_BOUND,
        tuple(DEFAULT_CUTOFFS),
    )
    baseline_summary = make_assay_fitness_summary(
        build_assay_fitness_table(
            assay_fitness, baseline_predictions, TARGET_UPPER_BOUND
        ),
        model_gene_ids,
        TARGET_UPPER_BOUND,
        tuple(DEFAULT_CUTOFFS),
    )
    rows: list[dict[str, Any]] = []
    for assay in sorted(scenario_summary["per_assay"]):
        scenario_metrics = scenario_summary["per_assay"][assay]
        baseline_metrics = baseline_summary["per_assay"][assay]
        record: dict[str, Any] = {
            "scenario_id": scenario_id,
            "assay": assay,
        }
        for key, value in scenario_metrics.items():
            record[f"scenario_{key}"] = value
        for key, value in baseline_metrics.items():
            record[f"baseline_{key}"] = value
        proxy = scenario_summary["concordant_nonessential_proxy_safety"]
        baseline_proxy = baseline_summary["concordant_nonessential_proxy_safety"]
        record.update(
            {
                "concordant_nonessential_source_proxy_genes": proxy[
                    "source_proxy_genes"
                ],
                "concordant_nonessential_model_proxy_genes": proxy["model_proxy_genes"],
                "scenario_proxy_predicted_essential_at_10pct": proxy["primary"][
                    "predicted_essential_count"
                ],
                "baseline_proxy_predicted_essential_at_10pct": baseline_proxy[
                    "primary"
                ]["predicted_essential_count"],
                "scenario_proxy_safety_rate_at_10pct": proxy["primary"]["safety_rate"],
                "baseline_proxy_safety_rate_at_10pct": baseline_proxy["primary"][
                    "safety_rate"
                ],
            }
        )
        rows.append(record)
    return rows


def _evaluate_scenario(
    *,
    baseline_model,
    baseline_predictions: pd.DataFrame,
    baseline_wt_growth: float,
    baseline_family_upper_bounds: dict[str, float],
    wt_fluxes: pd.Series,
    profile_rows: list[dict[str, Any]],
    profile_fieldnames: list[str],
    scenario: dict[str, Any],
    target_by_reaction: dict[str, str],
    backup_genes: set[str],
    proxy_genes: set[str],
    experimental_positive_genes: set[str],
    assay_fitness: pd.DataFrame,
    cutoffs: tuple[float, ...],
    solver: str,
    model_sha256: str,
    scenario_directory: Path,
) -> tuple[pd.DataFrame, dict[str, Any], list[dict[str, Any]]]:
    scenario_directory.mkdir(parents=True, exist_ok=True)
    model = baseline_model.copy()
    _configure_solver(model, solver)
    scenario_rows, caps = _scenario_profile_rows(profile_rows, scenario, wt_fluxes)
    audit = _apply_scenario_profile(
        model,
        scenario_rows,
        profile_fieldnames,
        model_sha256,
        scenario_directory,
    )
    scenario_wt_growth, _ = _optimal_growth_and_fluxes(model)
    family_totals: dict[str, tuple[float, float]] = {}
    for row in scenario_rows:
        source_id = str(row["source_reaction_id"])
        backup_id = str(row["backup_reaction_id"])
        current_total = float(model.reactions.get_by_id(source_id).upper_bound) + float(
            model.reactions.get_by_id(backup_id).upper_bound
        )
        family_totals[source_id] = (
            baseline_family_upper_bounds[source_id],
            current_total,
        )

    scenario_predictions, deletion_wt_growth = run_single_gene_deletions(model, solver)
    if abs(deletion_wt_growth - scenario_wt_growth) > WT_DELTA_TOLERANCE:
        raise RuntimeError("Scenario WT changed between optimization and deletion run")
    per_gene = _build_scenario_per_gene(
        baseline_predictions,
        scenario_predictions,
        target_genes=set(target_by_reaction.values()),
        backup_genes=backup_genes,
        proxy_genes=proxy_genes,
        experimental_positive_genes=experimental_positive_genes,
        cutoffs=cutoffs,
    )
    gates = _scenario_gates(
        per_gene,
        target_by_reaction=target_by_reaction,
        backup_genes=backup_genes,
        proxy_genes=proxy_genes,
        baseline_wt_growth=baseline_wt_growth,
        scenario_wt_growth=scenario_wt_growth,
        r1843_bounds=tuple(model.reactions.get_by_id(BLOCKED_CONTROL_REACTION).bounds),
        family_totals=family_totals,
        cutoffs=cutoffs,
    )
    summary: dict[str, Any] = {
        "status": "complete",
        "scenario": scenario,
        "baseline_wt_growth": baseline_wt_growth,
        "scenario_wt_growth": scenario_wt_growth,
        "applied_backup_upper_bounds": caps,
        "capacity_audit": audit,
        "gates": gates,
        "workers": WORKERS,
        "solver": solver,
        "completed_at": utc_now(),
    }
    assay_rows = _assay_metrics(assay_fitness, per_gene, str(scenario["scenario_id"]))
    return per_gene, summary, assay_rows


def _load_checkpoint(
    scenario_directory: Path,
    *,
    run_key: str,
    scenario_key: str,
) -> tuple[pd.DataFrame, dict[str, Any]] | None:
    table_path = scenario_directory / SCENARIO_TABLE_NAME
    summary_path = scenario_directory / SCENARIO_SUMMARY_NAME
    if not table_path.exists() and not summary_path.exists():
        return None
    if not table_path.is_file() or not summary_path.is_file():
        return None
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("run_key") != run_key or summary.get("scenario_key") != scenario_key:
        raise RuntimeError(
            f"Checkpoint provenance mismatch in {scenario_directory}; refusing resume"
        )
    if summary.get("status") != "complete":
        return None
    if summary.get("per_gene_sha256") != sha256_file(table_path):
        raise RuntimeError(f"Checkpoint table hash mismatch: {table_path}")
    return pd.read_csv(table_path, sep="\t"), summary


def _write_checkpoint(
    scenario_directory: Path,
    per_gene: pd.DataFrame,
    summary: dict[str, Any],
    *,
    run_key: str,
    scenario_key: str,
) -> None:
    scenario_directory.mkdir(parents=True, exist_ok=True)
    table_path = scenario_directory / SCENARIO_TABLE_NAME
    table_payload = _dataframe_tsv_bytes(per_gene)
    _atomic_write_bytes(table_path, table_payload)
    checkpoint = dict(summary)
    checkpoint["run_key"] = run_key
    checkpoint["scenario_key"] = scenario_key
    checkpoint["per_gene_sha256"] = hashlib.sha256(table_payload).hexdigest()
    # Written last: the summary is the commit marker for an immediately
    # resumable, complete scenario.
    _atomic_write_json(scenario_directory / SCENARIO_SUMMARY_NAME, checkpoint)


def _summary_row(summary: dict[str, Any]) -> dict[str, Any]:
    scenario = summary["scenario"]
    gates = summary["gates"]
    return {
        "scenario_id": scenario["scenario_id"],
        "r4_capacity_fraction_of_wt_flux": scenario["r4_capacity_fraction_of_wt_flux"],
        "r1846_capacity_fraction_of_wt_flux": scenario[
            "r1846_capacity_fraction_of_wt_flux"
        ],
        "r4_applied_backup_upper_bound": summary["applied_backup_upper_bounds"]["R4"],
        "r1846_applied_backup_upper_bound": summary["applied_backup_upper_bounds"][
            "R1846"
        ],
        "r4_target_ko_ratio": gates["target_ko_ratios"]["R4"],
        "r1846_target_ko_ratio": gates["target_ko_ratios"]["R1846"],
        "feasible": gates["feasible"],
        "target_window_pass": gates["target_window_pass"],
        "boundary_margin": gates["boundary_margin"],
        "wt_growth_delta": gates["wt_growth_delta"],
        "wt_growth_gate_pass": gates["wt_growth_gate_pass"],
        "r1843_bounds_gate_pass": gates["r1843_bounds_gate_pass"],
        "reaction_family_total_upper_bound_gate_pass": gates[
            "reaction_family_total_upper_bound_gate_pass"
        ],
        "non_target_call_flip_count": gates["non_target_call_flip_count"],
        "no_non_target_call_flips_gate_pass": gates[
            "no_non_target_call_flips_gate_pass"
        ],
        "backup_gene_ko_ratio_gate_pass": gates["backup_gene_ko_ratio_gate_pass"],
        "proxy_new_essential_count": gates["proxy_new_essential_count"],
        "concordant_nonessential_proxy_gate_pass": gates[
            "concordant_nonessential_proxy_gate_pass"
        ],
        "non_target_ratio_change_score": gates["non_target_ratio_change_score"],
        "non_target_ratio_change_sum": gates["non_target_ratio_change_sum"],
    }


def rank_scenarios(summary_rows: Iterable[dict[str, Any]]) -> pd.DataFrame:
    """Apply the approved deterministic lexicographic ranking."""
    frame = pd.DataFrame(list(summary_rows))
    if frame.empty:
        raise ValueError("Cannot rank an empty scenario set")
    frame = frame.sort_values(
        [
            "feasible",
            "boundary_margin",
            "non_target_ratio_change_score",
            "scenario_id",
        ],
        ascending=[False, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    frame.insert(0, "rank", range(1, len(frame) + 1))
    return frame


def _collateral_rows(
    scenario_id: str,
    per_gene: pd.DataFrame,
    cutoffs: tuple[float, ...],
) -> list[dict[str, Any]]:
    target_mask = per_gene["is_target_gene"].map(_parse_bool)
    delta_mask = per_gene["abs_ko_growth_ratio_delta"] > RATIO_CHANGE_EPS
    flip_columns = [
        f"call_flip_at_{cutoff * 100:g}pct".replace(".", "p") for cutoff in cutoffs
    ]
    flip_mask = (
        per_gene[flip_columns].apply(lambda column: column.map(_parse_bool)).any(axis=1)
    )
    selected = per_gene.loc[~target_mask & (delta_mask | flip_mask)].copy()
    rows: list[dict[str, Any]] = []
    for row in selected.itertuples(index=False):
        record = {
            "scenario_id": scenario_id,
            "gene_id": row.gene_id,
            "baseline_ko_growth_ratio": row.baseline_ko_growth_ratio,
            "scenario_ko_growth_ratio": row.scenario_ko_growth_ratio,
            "ko_growth_ratio_delta": row.ko_growth_ratio_delta,
            "abs_ko_growth_ratio_delta": row.abs_ko_growth_ratio_delta,
            "is_backup_gene": row.is_backup_gene,
            "is_concordant_nonessential_proxy": row.is_concordant_nonessential_proxy,
        }
        for column in flip_columns:
            record[column] = getattr(row, column)
        rows.append(record)
    return rows


def _candidate_profile_rows(
    profile_rows: list[dict[str, Any]],
    selected: pd.Series,
    grid_sha256: str,
) -> list[dict[str, Any]]:
    bounds = {
        "R4": float(selected["r4_applied_backup_upper_bound"]),
        "R1846": float(selected["r1846_applied_backup_upper_bound"]),
    }
    rows: list[dict[str, Any]] = []
    for original in profile_rows:
        row = dict(original)
        reaction_id = str(row["source_reaction_id"])
        row["provisional_upper_bound"] = bounds[reaction_id]
        row["parameter_basis"] = (
            "joint_SD-Leu_grid_selected_not_measured;"
            f"scenario={selected['scenario_id']};grid_sha256={grid_sha256}"
        )
        row["rationale"] = (
            f"Selected by the provenance-locked R4/R1846 joint grid as "
            f"{selected['scenario_id']}; still exploratory and requires "
            "condition-matched abundance and kcat replacement."
        )
        rows.append(row)
    return rows


def calibrate_isozyme_capacities(
    *,
    model_path: Path = DEFAULT_MODEL,
    profile_path: Path = DEFAULT_PROFILE,
    grid_path: Path = DEFAULT_GRID,
    experimental_path: Path = DEFAULT_EXPERIMENTAL,
    assay_fitness_path: Path = DEFAULT_ASSAY_FITNESS,
    media_path: Path = DEFAULT_MEDIA,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    solver: str = "gurobi",
    growth_cutoffs: Iterable[float] = DEFAULT_CUTOFFS,
    workers: int = WORKERS,
    resume: bool = False,
) -> dict[str, Any]:
    """Run or exactly resume the nine-scenario joint capacity calibration."""
    if workers != WORKERS:
        raise ValueError("Joint capacity calibration is locked to workers=1")
    paths = [
        Path(model_path),
        Path(profile_path),
        Path(grid_path),
        Path(experimental_path),
        Path(assay_fitness_path),
        Path(media_path),
    ]
    (
        model_path,
        profile_path,
        grid_path,
        experimental_path,
        assay_fitness_path,
        media_path,
    ) = paths
    output_dir = Path(output_dir)
    candidate_path = output_dir / CANDIDATE_PROFILE_NAME
    if candidate_path.resolve() == profile_path.resolve():
        raise ValueError(
            "Candidate profile cannot overwrite the input capacity profile"
        )
    if output_dir.resolve() == model_path.resolve():
        raise ValueError("Output directory cannot be the canonical model path")

    inputs = _input_manifest(
        model_path=model_path,
        profile_path=profile_path,
        grid_path=grid_path,
        experimental_path=experimental_path,
        assay_fitness_path=assay_fitness_path,
        media_path=media_path,
    )
    cutoffs = tuple(float(value) for value in growth_cutoffs)
    if not cutoffs or any(not 0.0 < cutoff < 1.0 for cutoff in cutoffs):
        raise ValueError("growth cutoffs must be fractions strictly between 0 and 1")
    if tuple(sorted(set(cutoffs))) != cutoffs:
        raise ValueError("growth cutoffs must be unique and strictly increasing")
    configuration = _run_configuration(solver, cutoffs)
    run_key = _run_key(inputs, configuration)
    manifest_path = output_dir / MANIFEST_NAME

    existing_manifest: dict[str, Any] | None = None
    if manifest_path.exists():
        existing_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not resume:
            raise FileExistsError(
                f"Output manifest already exists; pass --resume for {output_dir}"
            )
        if existing_manifest.get("run_key") != run_key:
            raise RuntimeError(
                "Resume refused: configuration or input SHA differs from the manifest"
            )
    elif output_dir.exists() and any(output_dir.iterdir()):
        if resume:
            raise RuntimeError(
                "Resume refused: non-empty output has no provenance manifest"
            )
        raise FileExistsError(
            f"Output directory is not empty and has no manifest: {output_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "analysis_type": "exploratory_global_R4_R1846_joint_capacity_grid",
        "status": "running",
        "exploratory_only": True,
        "not_a_canonical_model_patch": True,
        "model_xml_written": False,
        "input_profile_overwritten": False,
        "inputs": inputs,
        "configuration": configuration,
        "run_key": run_key,
        "started_at": (
            existing_manifest.get("started_at")
            if existing_manifest is not None
            else utc_now()
        ),
        "updated_at": utc_now(),
    }
    _atomic_write_json(manifest_path, manifest)

    grid = load_joint_capacity_grid(grid_path)
    profile_rows = load_provisional_capacity_table(profile_path)
    profile_by_reaction = _profile_rows_by_reaction(profile_rows)
    profile_fieldnames = _profile_fieldnames(profile_path)
    target_by_reaction, backup_genes = _profile_gene_roles(profile_rows)
    assay_fitness = _load_assay_fitness(assay_fitness_path)
    proxy_genes = concordant_nonessential_proxies(assay_fitness)
    experimental = load_experimental(experimental_path, positive_only=True)
    experimental_positive_genes = set(
        experimental.loc[experimental["essential"], "gene_id"].astype(str)
    )

    baseline_model = _load_model(model_path, solver)
    apply_media(baseline_model, load_media(media_path))
    if BLOCKED_CONTROL_REACTION not in baseline_model.reactions:
        raise ValueError(
            f"Blocked control reaction not found: {BLOCKED_CONTROL_REACTION}"
        )
    if tuple(baseline_model.reactions.get_by_id(BLOCKED_CONTROL_REACTION).bounds) != (
        0.0,
        0.0,
    ):
        raise ValueError(f"{BLOCKED_CONTROL_REACTION} must start with bounds (0, 0)")
    for reaction_id, row in profile_by_reaction.items():
        if row["validated_model_sha256"] != inputs["model"]["sha256"]:
            raise ValueError(
                f"Profile {row['capacity_id']} model SHA is stale for {reaction_id}"
            )
    baseline_wt_growth, wt_fluxes = _optimal_growth_and_fluxes(baseline_model)
    baseline_predictions, deletion_wt_growth = run_single_gene_deletions(
        baseline_model, solver
    )
    if abs(deletion_wt_growth - baseline_wt_growth) > WT_DELTA_TOLERANCE:
        raise RuntimeError("Baseline WT changed between optimization and deletion run")
    baseline_family_upper_bounds = {
        reaction_id: float(baseline_model.reactions.get_by_id(reaction_id).upper_bound)
        for reaction_id in TARGET_REACTIONS
    }

    summaries: list[dict[str, Any]] = []
    per_gene_by_scenario: dict[str, pd.DataFrame] = {}
    assay_rows: list[dict[str, Any]] = []
    scenarios_root = output_dir / "scenarios"
    for grid_row in grid.itertuples(index=False):
        scenario = _scenario_configuration(grid_row)
        scenario_id = str(scenario["scenario_id"])
        scenario_directory = scenarios_root / scenario_id
        scenario_key = _scenario_key(run_key, scenario)
        checkpoint = (
            _load_checkpoint(
                scenario_directory,
                run_key=run_key,
                scenario_key=scenario_key,
            )
            if resume
            else None
        )
        if checkpoint is None:
            per_gene, summary, scenario_assay_rows = _evaluate_scenario(
                baseline_model=baseline_model,
                baseline_predictions=baseline_predictions,
                baseline_wt_growth=baseline_wt_growth,
                baseline_family_upper_bounds=baseline_family_upper_bounds,
                wt_fluxes=wt_fluxes,
                profile_rows=profile_rows,
                profile_fieldnames=profile_fieldnames,
                scenario=scenario,
                target_by_reaction=target_by_reaction,
                backup_genes=backup_genes,
                proxy_genes=proxy_genes,
                experimental_positive_genes=experimental_positive_genes,
                assay_fitness=assay_fitness,
                cutoffs=cutoffs,
                solver=solver,
                model_sha256=inputs["model"]["sha256"],
                scenario_directory=scenario_directory,
            )
            _write_checkpoint(
                scenario_directory,
                per_gene,
                summary,
                run_key=run_key,
                scenario_key=scenario_key,
            )
            summary["run_key"] = run_key
            summary["scenario_key"] = scenario_key
            summary["per_gene_sha256"] = sha256_file(
                scenario_directory / SCENARIO_TABLE_NAME
            )
        else:
            per_gene, summary = checkpoint
            scenario_assay_rows = _assay_metrics(assay_fitness, per_gene, scenario_id)
        summaries.append(summary)
        per_gene_by_scenario[scenario_id] = per_gene
        assay_rows.extend(scenario_assay_rows)

    if len(summaries) != 9:
        raise RuntimeError(
            f"Joint grid aggregate requires 9 scenarios, got {len(summaries)}"
        )
    ranked = rank_scenarios(_summary_row(summary) for summary in summaries)
    collateral_rows = [
        record
        for scenario_id, per_gene in sorted(per_gene_by_scenario.items())
        for record in _collateral_rows(scenario_id, per_gene, cutoffs)
    ]
    collateral_columns = [
        "scenario_id",
        "gene_id",
        "baseline_ko_growth_ratio",
        "scenario_ko_growth_ratio",
        "ko_growth_ratio_delta",
        "abs_ko_growth_ratio_delta",
        "is_backup_gene",
        "is_concordant_nonessential_proxy",
        *[f"call_flip_at_{cutoff * 100:g}pct".replace(".", "p") for cutoff in cutoffs],
    ]
    collateral = pd.DataFrame(collateral_rows, columns=collateral_columns)
    assay_metrics = (
        pd.DataFrame(assay_rows)
        .sort_values(["scenario_id", "assay"])
        .reset_index(drop=True)
    )

    grid_summary_path = output_dir / GRID_SUMMARY_NAME
    collateral_path = output_dir / COLLATERAL_NAME
    assay_metrics_path = output_dir / ASSAY_METRICS_NAME
    _atomic_write_tsv(grid_summary_path, ranked)
    _atomic_write_tsv(collateral_path, collateral)
    _atomic_write_tsv(assay_metrics_path, assay_metrics)

    feasible = ranked.loc[ranked["feasible"].map(_parse_bool)]
    candidate_created = False
    selected_scenario_id: str | None = None
    if not feasible.empty:
        selected = feasible.iloc[0]
        selected_scenario_id = str(selected["scenario_id"])
        candidate_rows = _candidate_profile_rows(
            profile_rows,
            selected,
            inputs["grid"]["sha256"],
        )
        _atomic_write_bytes(
            candidate_path,
            _profile_csv_bytes(candidate_rows, profile_fieldnames),
        )
        # The candidate must pass the production profile's static schema before
        # the run may be marked complete.
        load_provisional_capacity_table(candidate_path)
        candidate_created = True

    outputs: dict[str, Any] = {
        "joint_capacity_grid": {
            "path": str(grid_summary_path.resolve()),
            "sha256": sha256_file(grid_summary_path),
        },
        "collateral_changes": {
            "path": str(collateral_path.resolve()),
            "sha256": sha256_file(collateral_path),
        },
        "assay_metrics": {
            "path": str(assay_metrics_path.resolve()),
            "sha256": sha256_file(assay_metrics_path),
        },
        "scenario_checkpoints": {
            "count": len(summaries),
            "root": str(scenarios_root.resolve()),
        },
    }
    if candidate_created:
        outputs["candidate_experimental_profile"] = {
            "path": str(candidate_path.resolve()),
            "sha256": sha256_file(candidate_path),
        }
    manifest.update(
        {
            "status": "complete",
            "updated_at": utc_now(),
            "completed_at": utc_now(),
            "scenario_count": len(summaries),
            "feasible_scenario_count": len(feasible),
            "selected_scenario_id": selected_scenario_id,
            "candidate_profile_created": candidate_created,
            "outputs": outputs,
        }
    )
    _atomic_write_json(manifest_path, manifest)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the provenance-locked R4/R1846 joint capacity grid"
    )
    parser.add_argument("command", nargs="?", choices=("run",), default="run")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--profile", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--grid", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--experimental", type=Path, default=DEFAULT_EXPERIMENTAL)
    parser.add_argument("--assay-fitness", type=Path, default=DEFAULT_ASSAY_FITNESS)
    parser.add_argument("--media", type=Path, default=DEFAULT_MEDIA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--solver", default="gurobi")
    parser.add_argument(
        "--growth-cutoffs",
        type=_parse_growth_cutoffs,
        default=tuple(DEFAULT_CUTOFFS),
        help="Comma-separated fractions of WT (default: 0.01,0.05,0.10,0.15)",
    )
    parser.add_argument("--workers", type=int, default=WORKERS)
    parser.add_argument("--resume", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = calibrate_isozyme_capacities(
        model_path=args.model,
        profile_path=args.profile,
        grid_path=args.grid,
        experimental_path=args.experimental,
        assay_fitness_path=args.assay_fitness,
        media_path=args.media,
        output_dir=args.output_dir,
        solver=args.solver,
        growth_cutoffs=args.growth_cutoffs,
        workers=args.workers,
        resume=args.resume,
    )
    selected = manifest.get("selected_scenario_id") or "none"
    print(
        f"Completed {manifest['scenario_count']} scenarios; "
        f"feasible={manifest['feasible_scenario_count']}; selected={selected}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
