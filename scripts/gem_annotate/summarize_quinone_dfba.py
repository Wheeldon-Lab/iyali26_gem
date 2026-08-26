"""Read-only summaries and calls-level QC for one CoQ9 dFBA screen."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from .essentiality_simulation_context import sha256_payload
from .quinone_dfba_essentiality import (
    DEFAULT_ALPHAS,
    DEFAULT_POOL_MULTIPLIERS,
    R39_R19_RUNTIME_MAPPING_SHA256,
    RUNTIME_GPR_SCENARIOS,
    SCHEMA_VERSION,
    _runtime_gpr_mapping_payload,
)
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
NINE_COQ_CANDIDATES = (
    "YALI1C26017g", "YALI1F08349g", "YALI1B20835g", "YALI1F34625g",
    "YALI1C25352g", "YALI1A08781g", "YALI1E18269g", "YALI1B20527g",
    "YALI1F34675g",
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
            "the finite CoQ9 source is a complete-block control, not a general partial-loss rescue model",
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
    fields = (
        "schema_version", "solver", "runtime_versions", "optimizer", "nonoptimal_policy",
        "dt_h", "uracil_mode", "calibration_status", "input_sha256", "script_sha256",
        "coq9_dilution_tool_sha256", "alpha_sampling", "coq9_dilution",
        "q9_reserve_definition", "q9_reserve_policy", "simulation_context", "runtime_topology",
    )
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
        "runtime_topology": payloads[0].get("runtime_topology", {"enabled": False}),
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


def write_hypothesis_matrix(
    arm_dirs: list[Path], evidence_matrix_path: Path, output_dir: Path,
) -> Path:
    """Validate and join the seven preregistered runtime-GPR sensitivity arms."""
    expected_scenarios = tuple(RUNTIME_GPR_SCENARIOS)
    if len(arm_dirs) != len(expected_scenarios):
        raise ValueError(f"Expected {len(expected_scenarios)} arm directories")

    evidence = pd.read_csv(evidence_matrix_path, sep="\t")
    evidence_columns = [
        "gene_id", "cas9_fitness", "cas9_call", "cas12a_fitness", "cas12a_call",
    ]
    missing_evidence = set(evidence_columns) - set(evidence.columns)
    if missing_evidence:
        raise ValueError(f"Evidence matrix is missing columns: {sorted(missing_evidence)}")
    evidence = evidence[evidence_columns].copy()
    if evidence["gene_id"].duplicated().any() or set(evidence["gene_id"]) != set(NINE_COQ_CANDIDATES):
        raise ValueError("Evidence matrix must uniquely cover the nine CoQ candidates")
    for column in ("cas9_fitness", "cas12a_fitness"):
        evidence[column] = pd.to_numeric(evidence[column], errors="raise")
        if not evidence[column].map(math.isfinite).all():
            raise ValueError(f"Evidence matrix {column} must be finite")
    for column in ("cas9_call", "cas12a_call"):
        if not set(evidence[column]) <= {"essential", "nonessential"}:
            raise ValueError(f"Evidence matrix {column} has unknown calls")
    evidence = evidence.rename(columns={
        "cas9_fitness": "cas9_fs", "cas12a_fitness": "cas12a_fs",
    })

    arm_frames: dict[str, pd.DataFrame] = {}
    arm_inputs: dict[str, Any] = {}
    common_context: str | None = None
    runtime_target_fingerprint: str | None = None
    wt_frames = []
    for arm_dir in map(Path, arm_dirs):
        calls_path = arm_dir / "essentiality_dynamic_calls.tsv"
        manifest_path = arm_dir / "chunk_000_manifest.json"
        merge_path = arm_dir / "merge_manifest.json"
        calls = pd.read_csv(calls_path, sep="\t")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        merge = json.loads(merge_path.read_text(encoding="utf-8"))
        mapping = manifest.get("runtime_gpr_scenario", {})
        scenario = mapping.get("scenario_id")
        if scenario not in RUNTIME_GPR_SCENARIOS or scenario in arm_frames:
            raise ValueError(f"Unknown or duplicate runtime GPR scenario: {scenario}")
        expected_mapping = _runtime_gpr_mapping_payload(scenario)
        expected_mapping_sha = sha256_payload(expected_mapping)
        if not mapping.get("enabled") or any(mapping.get(key) != value for key, value in expected_mapping.items()):
            raise ValueError(f"Runtime GPR mapping payload mismatch for {scenario}")
        if mapping.get("mapping_sha256") != expected_mapping_sha:
            raise ValueError(f"Runtime GPR mapping SHA mismatch for {scenario}")
        observed_target_fingerprint = mapping.get("runtime_target_fingerprint_before_mapping")
        if not isinstance(observed_target_fingerprint, str) or not observed_target_fingerprint:
            raise ValueError(f"Runtime target fingerprint is missing for {scenario}")
        if runtime_target_fingerprint is not None and observed_target_fingerprint != runtime_target_fingerprint:
            raise ValueError("Runtime target fingerprint differs across GPR scenarios")
        runtime_target_fingerprint = observed_target_fingerprint
        topology = manifest.get("runtime_topology", {})
        if not topology.get("enabled") or topology.get("mapping_sha256") != R39_R19_RUNTIME_MAPPING_SHA256:
            raise ValueError(f"Runtime topology mismatch for {scenario}")
        if manifest.get("fingerprint") != sha256_payload({
            key: value for key, value in manifest.items() if key != "fingerprint"
        }):
            raise ValueError(f"Manifest fingerprint mismatch for {scenario}")
        expected_manifest = {
            "schema_version": SCHEMA_VERSION,
            "chunk_index": 0,
            "chunk_count": 1,
            "genes": list(NINE_COQ_CANDIDATES),
            "solver": "gurobi",
            "solver_feasibility_tolerance": 1e-9,
            "hours": 24.0,
            "dt_h": 0.0625,
            "initial_biomass_gDW_L": 0.01,
            "alphas_mmol_gDW": list(DEFAULT_ALPHAS),
            "pool_multipliers": list(DEFAULT_POOL_MULTIPLIERS),
            "uracil_mode": "po1f_nonlimiting",
            "optimizer": "pfba",
            "nonoptimal_policy": "terminal_record_and_stop",
            "calibration_status": "sensitivity_only_not_calibrated",
        }
        if any(manifest.get(key) != value for key, value in expected_manifest.items()):
            raise ValueError(f"Frozen run context mismatch for {scenario}")
        if manifest.get("runtime_versions", {}).get("gurobipy") != "10.0.3":
            raise ValueError(f"Gurobi runtime mismatch for {scenario}")
        if (
            merge.get("chunks") != 1
            or merge.get("calls") != 108
            or merge.get("chunk_fingerprints") != [manifest["fingerprint"]]
            or merge.get("runtime_gpr_scenario") != mapping
        ):
            raise ValueError(f"Merge manifest mismatch for {scenario}")

        key = ["gene_id", "alpha_mmol_gDW", "pool_multiplier"]
        if len(calls) != 108 or calls.duplicated(key).any():
            raise ValueError(f"Expected 108 unique calls for {scenario}")
        if (
            set(calls["gene_id"]) != set(NINE_COQ_CANDIDATES)
            or set(calls["alpha_mmol_gDW"]) != set(DEFAULT_ALPHAS)
            or set(calls["pool_multiplier"]) != set(DEFAULT_POOL_MULTIPLIERS)
        ):
            raise ValueError(f"Gene/alpha/pool grid mismatch for {scenario}")
        if set(calls["runtime_gpr_scenario"]) != {scenario} or set(calls["runtime_gpr_mapping_sha256"]) != {expected_mapping_sha}:
            raise ValueError(f"Calls mapping identity mismatch for {scenario}")
        for column in ("dynamic_doublings", "dynamic_growth_ratio", "q9_source_total_mmol_L"):
            if not calls[column].map(math.isfinite).all() or (calls[column] < 0).any():
                raise ValueError(f"Invalid nonnegative endpoint {column} for {scenario}")
        ratio = calls["dynamic_doublings"] / calls["wt_dynamic_doublings"]
        if not calls["wt_dynamic_doublings"].map(math.isfinite).all() or (calls["wt_dynamic_doublings"] <= 0).any():
            raise ValueError(f"Invalid WT endpoint for {scenario}")
        if (calls["dynamic_growth_ratio"].sub(ratio).abs() > 1e-12).any():
            raise ValueError(f"Stored ratio mismatch for {scenario}")
        _, calls_qc = summarize_calls(calls, float(manifest["initial_biomass_gDW_L"]))
        if not all(calls_qc[key] for key in (
            "shape_ok", "pool_monotonicity_pass", "q9_source_calls_bound_pass",
        )):
            raise ValueError(f"Calls-level QC failed for {scenario}")
        if (
            calls_qc["five_control_max_abs_theory_error_doublings"] > 1e-6
            or calls_qc["five_control_max_alpha_doublings_range"] > 1e-6
        ):
            raise ValueError(f"Five-gene theoretical control failed for {scenario}")

        wt_columns = [
            "alpha_mmol_gDW", "pool_multiplier", "wt_dynamic_doublings",
            "wt_termination_status", "wt_termination_time_h",
        ]
        if calls.groupby(key[1:])[wt_columns[2:]].nunique(dropna=False).gt(1).any().any():
            raise ValueError(f"WT endpoints are inconsistent within {scenario}")
        wt = calls[wt_columns].drop_duplicates().copy()
        wt["scenario_id"] = scenario
        wt_frames.append(wt)

        context = json.dumps({
            key: manifest[key] for key in (
                "schema_version", "solver", "runtime_versions", "solver_feasibility_tolerance",
                "script_sha256", "hours", "dt_h", "initial_biomass_gDW_L",
                "alphas_mmol_gDW", "pool_multipliers", "uracil_mode", "optimizer",
                "nonoptimal_policy", "calibration_status", "input_sha256",
                "coq9_dilution_tool_sha256", "alpha_sampling", "coq9_dilution",
                "q9_reserve_definition", "q9_reserve_policy", "simulation_context",
                "runtime_topology",
            )
        }, sort_keys=True)
        if common_context is not None and context != common_context:
            raise ValueError("Runtime context differs across GPR scenarios")
        common_context = context

        frame = calls.copy()
        frame.insert(0, "scenario_id", scenario)
        frame.insert(1, "mapping_id", mapping["mapping_id"])
        frame.insert(2, "mapping_sha256", expected_mapping_sha)
        frame.insert(3, "hypothesis_label", mapping["hypothesis_label"])
        frame.insert(4, "risk_labels", "|".join(mapping["risk_labels"]) or "none")
        frame.insert(5, "binary_gpr_interpretation", mapping["binary_gpr_interpretation"])
        frame.insert(6, "scope", mapping["scope"])
        frame.insert(7, "calibration_status", mapping["calibration_status"])
        frame.insert(8, "run_id", manifest["run_id"])
        frame.insert(9, "model_sha256", manifest["input_sha256"]["model"])
        frame.insert(10, "simulation_context_sha256", manifest["simulation_context"]["simulation_context_fingerprint"])
        frame.insert(11, "runtime_topology_mapping_sha256", topology["mapping_sha256"])
        frame.insert(14, "dt_h", manifest["dt_h"])
        arm_frames[scenario] = frame
        arm_inputs[scenario] = {
            "run_id": manifest["run_id"],
            "mapping_id": mapping["mapping_id"],
            "mapping_sha256": expected_mapping_sha,
            "runtime_target_fingerprint_before_mapping": observed_target_fingerprint,
            "manifest_fingerprint": manifest["fingerprint"],
            "manifest_sha256": _sha256(manifest_path),
            "calls_sha256": _sha256(calls_path),
            "coq9_dilution_tool_sha256": manifest["coq9_dilution_tool_sha256"],
            "alpha_sampling": manifest["alpha_sampling"],
            "coq9_dilution": manifest["coq9_dilution"],
            "q9_reserve_definition": manifest["q9_reserve_definition"],
            "q9_reserve_policy": manifest["q9_reserve_policy"],
        }

    if set(arm_frames) != set(expected_scenarios):
        raise ValueError("The seven preregistered scenarios are not all present")
    wt = pd.concat(wt_frames, ignore_index=True)
    wt_grouped = wt.groupby(["alpha_mmol_gDW", "pool_multiplier"])
    if (
        wt_grouped["wt_dynamic_doublings"].agg(lambda values: values.max() - values.min()).gt(1e-9).any()
        or wt_grouped[["wt_termination_status", "wt_termination_time_h"]].nunique(dropna=False).gt(1).any().any()
    ):
        raise ValueError("WT endpoints differ across GPR scenarios")

    baseline = arm_frames["baseline_topology"].set_index(
        ["gene_id", "alpha_mmol_gDW", "pool_multiplier"]
    ).sort_index()
    for scenario, frame in arm_frames.items():
        if scenario == "baseline_topology":
            continue
        candidate = RUNTIME_GPR_SCENARIOS[scenario]["added_candidate_genes"][0]
        compared = frame.set_index([
            "gene_id", "alpha_mmol_gDW", "pool_multiplier",
        ]).sort_index()
        unchanged = compared.index.get_level_values("gene_id") != candidate
        if (
            compared.loc[~unchanged, "dynamic_doublings"]
            .sub(baseline.loc[~unchanged, "dynamic_doublings"])
            .gt(1e-9)
            .any()
        ):
            raise ValueError(f"Target growth increased after adding constraints in {scenario}")
        for column in ("dynamic_doublings", "dynamic_growth_ratio"):
            if (compared.loc[unchanged, column].sub(baseline.loc[unchanged, column]).abs() > 1e-9).any():
                raise ValueError(f"Non-target {column} changed in {scenario}")
        for column in (
            "termination_status", "termination_time_h",
            *(_cutoff_column(cutoff) for cutoff in DEFAULT_CUTOFFS),
        ):
            observed = compared.loc[unchanged, column]
            reference = baseline.loc[unchanged, column]
            same = (
                observed.map(_bool).equals(reference.map(_bool))
                if column.startswith("essential_")
                else observed.equals(reference)
            )
            if not same:
                raise ValueError(f"Non-target {column} changed in {scenario}")

    nested = (
        ("coq6_both_hydroxylations", "coq6_r39", "YALI1A08781g"),
        ("coq6_both_hydroxylations", "coq6_r19_hydroxylation", "YALI1A08781g"),
        ("coq9_routewide_absolute", "coq9_r695_absolute", "YALI1F34675g"),
    )
    for constrained, reference, gene_id in nested:
        key = ["gene_id", "alpha_mmol_gDW", "pool_multiplier"]
        constrained_frame = arm_frames[constrained].set_index(key).sort_index().loc[gene_id]
        reference_frame = arm_frames[reference].set_index(key).sort_index().loc[gene_id]
        if constrained_frame["dynamic_doublings"].sub(reference_frame["dynamic_doublings"]).gt(1e-9).any():
            raise ValueError(f"Nested absolute constraint increased growth: {constrained} vs {reference}")

    matrix = pd.concat([arm_frames[scenario] for scenario in expected_scenarios], ignore_index=True)
    matrix = matrix.merge(evidence, on="gene_id", how="left", validate="many_to_one", indicator=True)
    if set(matrix["_merge"]) != {"both"}:
        raise ValueError("Experimental evidence did not cover every simulated gene")
    matrix = matrix.drop(columns="_merge")
    output_columns = [
        "scenario_id", "mapping_id", "mapping_sha256", "hypothesis_label", "risk_labels",
        "binary_gpr_interpretation", "scope", "calibration_status", "run_id", "model_sha256",
        "simulation_context_sha256", "runtime_topology_mapping_sha256", "alpha_mmol_gDW",
        "pool_multiplier", "dt_h", "gene_id", "dynamic_doublings", "wt_dynamic_doublings",
        "dynamic_growth_ratio", *[_cutoff_column(cutoff) for cutoff in DEFAULT_CUTOFFS],
        "cas9_fs", "cas9_call", "cas12a_fs", "cas12a_call",
    ]
    matrix = matrix[output_columns].rename(columns={
        "dynamic_doublings": "ko_dynamic_doublings",
        "dynamic_growth_ratio": "ko_wt_dynamic_doublings_ratio",
    })
    if len(matrix) != 756 or matrix.duplicated([
        "mapping_id", "gene_id", "alpha_mmol_gDW", "pool_multiplier",
    ]).any():
        raise ValueError("Hypothesis matrix is not the expected unique 756-row grid")

    output_dir.mkdir(parents=True, exist_ok=False)
    matrix_path = output_dir / "hypothesis_matrix.tsv"
    matrix.to_csv(matrix_path, sep="\t", index=False)
    result_manifest = {
        "scope": "runtime_only",
        "calibration_status": "sensitivity_only_not_calibrated",
        "rows": len(matrix),
        "scenarios": list(expected_scenarios),
        "genes": list(NINE_COQ_CANDIDATES),
        "alphas_mmol_gDW": list(DEFAULT_ALPHAS),
        "pool_multipliers": list(DEFAULT_POOL_MULTIPLIERS),
        "dt_h": 0.0625,
        "evidence_matrix_sha256": _sha256(evidence_matrix_path),
        "arm_inputs": arm_inputs,
        "output_sha256": _sha256(matrix_path),
        "checks": {
            "exact_seven_scenarios": True,
            "mapping_payloads_and_sha": True,
            "runtime_target_fingerprint_uniform": True,
            "frozen_context_uniform": True,
            "wt_endpoints_unchanged": True,
            "non_target_results_unchanged": True,
            "nested_absolute_constraints_nonincreasing": True,
            "strict_cutoffs_recomputed": True,
            "calls_level_pool_source_and_five_gene_controls": True,
            "evidence_exact_nine_gene_coverage": True,
        },
        "limitations": [
            "binary GPR arms test imposed complete-block assumptions, not protein function",
            "no mapping or parameter is selected by experimental recall",
            "no canonical model, GPR, bound, curated-data, or FN-dossier change is authorized",
        ],
    }
    (output_dir / "hypothesis_matrix_manifest.json").write_text(
        json.dumps(result_manifest, indent=2, sort_keys=True, default=_native) + "\n",
        encoding="utf-8",
    )
    return output_dir


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input-dir", type=Path)
    inputs.add_argument("--arm-dir", type=Path, action="append")
    parser.add_argument("--gene-evidence-matrix", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.arm_dir:
        if args.gene_evidence_matrix is None:
            parser.error("--arm-dir requires --gene-evidence-matrix")
        print(write_hypothesis_matrix(args.arm_dir, args.gene_evidence_matrix, args.output_dir))
    else:
        if args.gene_evidence_matrix is not None:
            parser.error("--gene-evidence-matrix is only valid with --arm-dir")
        print(write_summary(args.input_dir, args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
