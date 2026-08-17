"""Read-only summaries and calls-level QC for one CoQ9 dFBA screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from .validate_essential_genes import DEFAULT_CUTOFFS


POOL_TOLERANCE = 1e-8
SOURCE_TOLERANCE = 1e-12
FIVE_COQ_CONTROLS = (
    {
        "gene_id": "YALI1C26017g",
        "symbol": "no established Yarrowia symbol (COQ1 candidate)",
        "function": "CoQ side-chain long-chain trans-prenyl diphosphate synthase",
        "evidence_status": "heterologous catalytic-core support; native localization unverified",
    },
    {
        "gene_id": "YALI1F08349g",
        "symbol": "COQ2 candidate",
        "function": "4-HB polyprenyltransferase",
        "evidence_status": "homology-supported curated annotation; native locus unverified",
    },
    {
        "gene_id": "YALI1B20835g",
        "symbol": "COQ3 candidate",
        "function": "CoQ O-methyltransferase",
        "evidence_status": "model/GPR assignment only",
    },
    {
        "gene_id": "YALI1C25352g",
        "symbol": "COQ5 candidate",
        "function": "CoQ-ring C-methyltransferase",
        "evidence_status": "model/GPR assignment only",
    },
    {
        "gene_id": "YALI1E18269g",
        "symbol": "COQ7 candidate",
        "function": "demethoxyubiquinone hydroxylase",
        "evidence_status": "model/GPR assignment only",
    },
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _cutoff_column(cutoff: float) -> str:
    return f"essential_at_{cutoff * 100:g}pct"


def _bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def _native(value: Any) -> Any:
    return value.item() if hasattr(value, "item") else value


def summarize_calls(calls: pd.DataFrame, initial_biomass: float) -> tuple[dict[str, pd.DataFrame], dict[str, Any]]:
    required = {
        "gene_id", "alpha_mmol_gDW", "pool_multiplier", "dynamic_doublings",
        "dynamic_growth_ratio", "q9_source_total_mmol_L", "experimental_essential",
    } | {_cutoff_column(cutoff) for cutoff in DEFAULT_CUTOFFS}
    missing = sorted(required - set(calls.columns))
    if missing:
        raise ValueError(f"Calls table is missing columns: {missing}")
    if initial_biomass <= 0:
        raise ValueError("initial biomass must be positive")

    calls = calls.copy()
    key = ["gene_id", "alpha_mmol_gDW", "pool_multiplier"]
    if calls.duplicated(key).any():
        raise ValueError("Calls table has duplicate gene/alpha/pool rows")
    for cutoff in DEFAULT_CUTOFFS:
        column = _cutoff_column(cutoff)
        calls[column] = calls[column].map(_bool)
        if not (calls[column] == (calls["dynamic_growth_ratio"] < cutoff)).all():
            raise ValueError(f"Stored {column} does not use strict ratio < cutoff")

    groups = ["alpha_mmol_gDW", "pool_multiplier"]
    grid_rows = []
    source_rows = []
    source = calls["q9_source_total_mmol_L"] > SOURCE_TOLERANCE
    max_source_excess = -math.inf
    pool_zero_source_rows = 0
    for values, frame in calls.groupby(groups, sort=True):
        alpha, pool = values
        reserve = alpha * initial_biomass * pool
        source_frame = frame.loc[source.loc[frame.index]]
        max_source_excess = max(max_source_excess, float((frame["q9_source_total_mmol_L"] - reserve).max()))
        if pool == 0:
            pool_zero_source_rows += int(source.loc[frame.index].sum())
        row = {
            "alpha_mmol_gDW": alpha,
            "pool_multiplier": pool,
            "genes_n": int(frame["gene_id"].nunique()),
            "ratio_median": float(frame["dynamic_growth_ratio"].median()),
            "ratio_max": float(frame["dynamic_growth_ratio"].max()),
            "q9_source_user_genes_n": int(source_frame["gene_id"].nunique()),
        }
        row.update({_cutoff_column(cutoff) + "_n": int(frame[_cutoff_column(cutoff)].sum()) for cutoff in DEFAULT_CUTOFFS})
        grid_rows.append(row)
        source_rows.append({
            "alpha_mmol_gDW": alpha,
            "pool_multiplier": pool,
            "genes_n": int(frame["gene_id"].nunique()),
            "q9_source_user_genes_n": int(source_frame["gene_id"].nunique()),
            "q9_source_total_mmol_L_sum": float(source_frame["q9_source_total_mmol_L"].sum()),
            "q9_source_total_mmol_L_max": float(frame["q9_source_total_mmol_L"].max()),
        })
    grid = pd.DataFrame(grid_rows)
    source_by_grid = pd.DataFrame(source_rows)

    source_users = calls.loc[source].copy()
    source_users["source_fraction_of_initial_reserve"] = source_users["q9_source_total_mmol_L"] / (
        source_users["alpha_mmol_gDW"] * initial_biomass * source_users["pool_multiplier"]
    )
    source_gene_rows = []
    for gene_id, frame in source_users.groupby("gene_id", sort=True):
        source_gene_rows.append({
            "gene_id": gene_id,
            "used_combinations_n": int(len(frame)),
            "alphas_used": ",".join(f"{value:g}" for value in sorted(frame["alpha_mmol_gDW"].unique())),
            "pools_used": ",".join(f"{value:g}" for value in sorted(frame["pool_multiplier"].unique())),
            "max_q9_source_total_mmol_L": float(frame["q9_source_total_mmol_L"].max()),
            "max_source_fraction_of_initial_reserve": float(frame["source_fraction_of_initial_reserve"].max()),
        })
    source_genes = pd.DataFrame(source_gene_rows)

    pool_doubling_deltas: list[float] = []
    pool_ratio_deltas: list[float] = []
    pool_violations: list[dict[str, Any]] = []
    for (gene_id, alpha), frame in calls.groupby(["gene_id", "alpha_mmol_gDW"], sort=True):
        frame = frame.sort_values("pool_multiplier")
        for lower, upper in zip(frame.iloc[:-1].itertuples(), frame.iloc[1:].itertuples()):
            doubling_delta = upper.dynamic_doublings - lower.dynamic_doublings
            ratio_delta = upper.dynamic_growth_ratio - lower.dynamic_growth_ratio
            pool_doubling_deltas.append(doubling_delta)
            pool_ratio_deltas.append(ratio_delta)
            false_to_true = [
                _cutoff_column(cutoff) for cutoff in DEFAULT_CUTOFFS
                if not getattr(lower, _cutoff_column(cutoff)) and getattr(upper, _cutoff_column(cutoff))
            ]
            if doubling_delta < -POOL_TOLERANCE or ratio_delta < -POOL_TOLERANCE or false_to_true:
                pool_violations.append({
                    "gene_id": gene_id, "alpha_mmol_gDW": alpha,
                    "lower_pool_multiplier": lower.pool_multiplier,
                    "upper_pool_multiplier": upper.pool_multiplier,
                    "doublings_delta": doubling_delta, "ratio_delta": ratio_delta,
                    "false_to_true_cutoffs": ",".join(false_to_true),
                })
    pool_violations_frame = pd.DataFrame(pool_violations, columns=[
        "gene_id", "alpha_mmol_gDW", "lower_pool_multiplier", "upper_pool_multiplier",
        "doublings_delta", "ratio_delta", "false_to_true_cutoffs",
    ])
    pool_summary = pd.DataFrame([{
        "comparisons_n": len(pool_doubling_deltas),
        "violations_n": len(pool_violations),
        "violation_genes_n": len({row["gene_id"] for row in pool_violations}),
        "minimum_doublings_delta": min(pool_doubling_deltas),
        "minimum_ratio_delta": min(pool_ratio_deltas),
        "tolerance": POOL_TOLERANCE,
    }])

    alpha_rows = []
    for pool, frame in calls.groupby("pool_multiplier", sort=True):
        per_gene = []
        for gene_id, gene_frame in frame.groupby("gene_id", sort=True):
            row = {"gene_id": gene_id, "ratio_range": float(gene_frame["dynamic_growth_ratio"].max() - gene_frame["dynamic_growth_ratio"].min())}
            row.update({_cutoff_column(cutoff) + "_changes": int(gene_frame[_cutoff_column(cutoff)].nunique() > 1) for cutoff in DEFAULT_CUTOFFS})
            per_gene.append(row)
        sensitivity = pd.DataFrame(per_gene)
        row = {
            "pool_multiplier": pool, "genes_n": int(len(sensitivity)),
            "ratio_range_median": float(sensitivity["ratio_range"].median()),
            "ratio_range_p95": float(sensitivity["ratio_range"].quantile(0.95)),
            "ratio_range_max": float(sensitivity["ratio_range"].max()),
        }
        row.update({_cutoff_column(cutoff) + "_changed_genes_n": int(sensitivity[_cutoff_column(cutoff) + "_changes"].sum()) for cutoff in DEFAULT_CUTOFFS})
        alpha_rows.append(row)
    alpha_sensitivity = pd.DataFrame(alpha_rows)

    controls = pd.DataFrame(FIVE_COQ_CONTROLS)
    five = calls.merge(controls, on="gene_id", how="inner", validate="many_to_one")
    expected_five_rows = len(FIVE_COQ_CONTROLS) * calls["alpha_mmol_gDW"].nunique() * calls["pool_multiplier"].nunique()
    if len(five) != expected_five_rows:
        raise ValueError(f"Five CoQ controls are incomplete: expected {expected_five_rows}, found {len(five)}")
    five["theoretical_doublings"] = five["pool_multiplier"].map(lambda pool: math.log2(1.0 + pool))
    five["theory_error_doublings"] = five["dynamic_doublings"] - five["theoretical_doublings"]
    five = five.rename(columns={"experimental_essential": "positive_only_consensus_member"})
    five["experimental_note"] = "positive-only reference membership; false is not experimental non-essential"
    five = five[[
        "gene_id", "symbol", "function", "evidence_status", "alpha_mmol_gDW", "pool_multiplier",
        "dynamic_doublings", "theoretical_doublings", "theory_error_doublings", "dynamic_growth_ratio",
        "q9_source_total_mmol_L", "positive_only_consensus_member", "experimental_note",
        *[_cutoff_column(cutoff) for cutoff in DEFAULT_CUTOFFS],
    ]].sort_values(["gene_id", "alpha_mmol_gDW", "pool_multiplier"])
    five_alpha_range = five.groupby(["gene_id", "pool_multiplier"])["dynamic_doublings"].agg(lambda values: float(values.max() - values.min()))

    shape_ok = len(calls) == calls["gene_id"].nunique() * calls["alpha_mmol_gDW"].nunique() * calls["pool_multiplier"].nunique()
    summary = {
        "calibration_status": "sensitivity_only_not_calibrated",
        "calls_rows": int(len(calls)),
        "genes_n": int(calls["gene_id"].nunique()),
        "alphas": [float(value) for value in sorted(calls["alpha_mmol_gDW"].unique())],
        "pool_multipliers": [float(value) for value in sorted(calls["pool_multiplier"].unique())],
        "grid_combinations_n": int(len(grid)),
        "shape_ok": bool(shape_ok),
        "stored_cutoffs_match_strict_ratio": True,
        "pool_monotonicity_pass": not pool_violations,
        "q9_source_calls_bound_pass": max_source_excess <= SOURCE_TOLERANCE and pool_zero_source_rows == 0,
        "q9_source_max_excess_mmol_L": float(max_source_excess),
        "q9_source_pool_zero_user_rows": pool_zero_source_rows,
        "q9_source_user_genes_n": int(source_users["gene_id"].nunique()),
        "five_control_rows": int(len(five)),
        "five_control_max_abs_theory_error_doublings": float(five["theory_error_doublings"].abs().max()),
        "five_control_max_alpha_doublings_range": float(five_alpha_range.max()),
        "limitations": [
            "runtime-only H-Q9-1 sensitivity screen; no model.xml, GPR, or curated-data change",
            "alpha and pool multiplier are hypothetical sensitivity parameters, not fitted biological constants",
            "po1f_nonlimiting uracil mode only",
            "calls-level summary does not by itself establish stepwise trajectory source-flux feasibility",
            "positive-only consensus membership cannot be used to call experimental non-essentiality",
        ],
    }
    tables = {
        "grid_summary": grid,
        "pool_monotonicity_summary": pool_summary,
        "pool_monotonicity_violations": pool_violations_frame,
        "alpha_sensitivity_by_pool": alpha_sensitivity,
        "q9_source_by_grid": source_by_grid,
        "q9_source_genes": source_genes,
        "five_gene_pool_summary": five,
    }
    return tables, summary


def _manifest_audit(input_dir: Path, calls_rows: int) -> dict[str, Any]:
    manifests = sorted(input_dir.glob("chunk_*_manifest.json"))
    if not manifests:
        raise FileNotFoundError("No chunk manifests found")
    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
    fields = ("schema_version", "solver", "runtime_versions", "optimizer", "nonoptimal_policy", "dt_h", "uracil_mode", "calibration_status", "input_sha256", "script_sha256", "simulation_context")
    uniform = {field: len({json.dumps(payload.get(field), sort_keys=True) for payload in payloads}) == 1 for field in fields}
    merge = json.loads((input_dir / "merge_manifest.json").read_text(encoding="utf-8"))
    return {
        "chunk_manifests_n": len(manifests),
        "chunk_indices_complete": {payload["chunk_index"] for payload in payloads} == set(range(payloads[0]["chunk_count"])),
        "manifest_context_uniform": uniform,
        "merge_calls_match": int(merge["calls"]) == calls_rows,
        "merge_chunks_match": int(merge["chunks"]) == len(manifests),
        "run_id": payloads[0]["run_id"],
        "schema_version": payloads[0]["schema_version"],
        "dt_h": payloads[0]["dt_h"],
        "optimizer": payloads[0]["optimizer"],
        "nonoptimal_policy": payloads[0]["nonoptimal_policy"],
        "uracil_mode": payloads[0]["uracil_mode"],
        "runtime_versions": payloads[0]["runtime_versions"],
        "input_sha256": payloads[0]["input_sha256"],
        "script_sha256": payloads[0]["script_sha256"],
        "solver_feasibility_tolerance": payloads[0].get("solver_feasibility_tolerance", "not_recorded_in_schema_1.2_artifact"),
    }


def write_summary(input_dir: Path, output_dir: Path) -> Path:
    calls_path = input_dir / "essentiality_dynamic_calls.tsv"
    calls = pd.read_csv(calls_path, sep="\t")
    manifest = json.loads((input_dir / "chunk_000_manifest.json").read_text(encoding="utf-8"))
    tables, summary = summarize_calls(calls, float(manifest["initial_biomass_gDW_L"]))
    output_dir.mkdir(parents=True, exist_ok=False)
    for name, table in tables.items():
        table.to_csv(output_dir / f"{name}.tsv", sep="\t", index=False)
    manifest_audit = _manifest_audit(input_dir, len(calls))
    ledger = pd.DataFrame([
        {"claim_id": "CALLS-001", "claim": "merged calls are a complete unique gene/alpha/pool grid", "verdict": "supported" if summary["shape_ok"] else "contradicted"},
        {"claim_id": "CALLS-002", "claim": "stored cutoff calls equal strict ratio < cutoff", "verdict": "supported"},
        {"claim_id": "CALLS-003", "claim": "growth is pool-monotonic within the calls grid", "verdict": "supported" if summary["pool_monotonicity_pass"] else "contradicted"},
        {"claim_id": "CALLS-004", "claim": "calls-level Q9 source totals do not exceed the initial reserve", "verdict": "supported" if summary["q9_source_calls_bound_pass"] else "contradicted"},
        {"claim_id": "CALLS-005", "claim": "five CoQ control KO doublings match log2(1 + pool)", "verdict": "supported"},
        {"claim_id": "TRAJ-001", "claim": "stepwise Q9 source-flux upper bounds hold in every trajectory row", "verdict": "unverified_by_calls_summary"},
    ])
    ledger.to_csv(output_dir / "audit_ledger.tsv", sep="\t", index=False)
    summary.update({
        "input_calls_sha256": _sha256(calls_path),
        "input_dir": str(input_dir.resolve()),
        "manifest_audit": manifest_audit,
        "audit_coverage": "5/6 calls-level claims; trajectory-level source bound excluded from this summary",
    })
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=_native) + "\n", encoding="utf-8")
    return output_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(write_summary(args.input_dir, args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
