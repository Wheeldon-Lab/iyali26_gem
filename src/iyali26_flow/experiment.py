"""Checkpointed orchestration for the R4/R1846 phase-one experiment."""

from __future__ import annotations

import csv
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from . import __version__
from .analysis import (
    AffineDiagnostics,
    affine_holdout_diagnostics,
    grid_posterior_recovery,
    make_phase1_decision,
    posterior_summary,
    reference_lattice,
    sobol_parameter_points,
)
from .config import ExperimentConfig
from .core import Fidelity, ParameterPoint, Phase1Decision, SimulationResult
from .provenance import (
    atomic_write_json,
    atomic_write_tsv,
    output_inventory,
    sha256_file,
    sha256_json,
    utc_now,
)
from .simulator import R4R1846CapacitySimulator


MANIFEST_NAME = "run_manifest.json"
DECISION_NAME = "phase1_decision.json"
BASELINE_METRICS_NAME = "baseline_metrics.json"
POSTERIOR_NAME = "posterior_recovery.tsv"
AUDIT_COMPARISON_NAME = "full_audit_comparison.tsv"


class ExperimentTimeoutError(RuntimeError):
    """Base exception for a bounded run that preserved resumable outputs."""


class SoftTimeBudgetExceeded(ExperimentTimeoutError):
    pass


class HardTimeoutExceeded(ExperimentTimeoutError):
    pass


@dataclass(frozen=True, slots=True)
class ExperimentOutcome:
    status: str
    output_dir: Path
    manifest_path: Path
    decision: Phase1Decision | None

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "output_dir": str(self.output_dir),
            "manifest": str(self.manifest_path),
            "decision": None if self.decision is None else self.decision.to_dict(),
        }


def _optional_float(value: object) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def _optional_int(value: object) -> int | None:
    if value is None or str(value).strip() == "":
        return None
    return int(float(value))


def _optional_bool(value: object) -> bool | None:
    if value is None or str(value).strip() == "":
        return None
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError(f"Invalid Boolean checkpoint value: {value!r}")


def _result_from_record(
    record: dict[str, object],
    backup_genes: Sequence[str],
) -> SimulationResult:
    return SimulationResult(
        sample_id=str(record["sample_id"]),
        fidelity=str(record["fidelity"]),  # type: ignore[arg-type]
        point=ParameterPoint(
            float(record["r4_capacity_fraction_of_wt_flux"]),
            float(record["r1846_capacity_fraction_of_wt_flux"]),
        ),
        r4_applied_upper_bound=float(record["r4_applied_backup_upper_bound"]),
        r1846_applied_upper_bound=float(
            record["r1846_applied_backup_upper_bound"]
        ),
        r4_target_ko_ratio=float(record["r4_target_ko_ratio"]),
        r1846_target_ko_ratio=float(record["r1846_target_ko_ratio"]),
        backup_ko_ratios=tuple(
            (gene_id, float(record[f"backup_ko_ratio__{gene_id}"]))
            for gene_id in backup_genes
        ),
        baseline_wt_growth=float(record["baseline_wt_growth"]),
        scenario_wt_growth=float(record["scenario_wt_growth"]),
        wt_growth_delta=float(record["wt_growth_delta"]),
        target_window_pass=bool(_optional_bool(record["target_window_pass"])),
        backup_gate_pass=bool(
            _optional_bool(record["backup_gene_ko_ratio_gate_pass"])
        ),
        r1843_gate_pass=bool(_optional_bool(record["r1843_bounds_gate_pass"])),
        family_total_gate_pass=bool(
            _optional_bool(
                record["reaction_family_total_upper_bound_gate_pass"]
            )
        ),
        targeted_feasible=bool(_optional_bool(record["targeted_feasible"])),
        full_feasible=_optional_bool(record["full_feasible"]),
        non_target_call_flip_count=_optional_int(
            record["non_target_call_flip_count"]
        ),
        proxy_new_essential_count=_optional_int(
            record["proxy_new_essential_count"]
        ),
        no_non_target_call_flips_gate_pass=_optional_bool(
            record["no_non_target_call_flips_gate_pass"]
        ),
        proxy_gate_pass=_optional_bool(
            record["concordant_nonessential_proxy_gate_pass"]
        ),
        gene_count=int(float(record["gene_count"])),
        elapsed_seconds=float(record["elapsed_seconds"]),
    )


def _read_tsv(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _protected_paths(config: ExperimentConfig) -> tuple[Path, ...]:
    candidates = (
        config.input_map["model"].path,
        config.input_map["profile"].path,
        config.repo_root / "data" / "essentiality" / "curated_model_patches.csv",
        config.repo_root / "data" / "essentiality" / "curation_cases.csv",
    )
    return tuple(dict.fromkeys(path.resolve() for path in candidates))


def _hash_paths(paths: Iterable[Path]) -> dict[str, dict[str, str | None]]:
    return {
        str(path): {
            "sha256": sha256_file(path) if path.is_file() else None,
        }
        for path in paths
    }


def _code_provenance(repo_root: Path) -> dict[str, str]:
    paths = [
        *sorted((repo_root / "src" / "iyali26_flow").glob("*.py")),
        repo_root / "scripts" / "gem_annotate" / "calibrate_isozyme_capacities.py",
        repo_root / "scripts" / "gem_annotate" / "provisional_capacity.py",
        repo_root / "scripts" / "gem_annotate" / "validate_essential_genes.py",
        repo_root / "pyproject.toml",
        repo_root / "uv.lock",
    ]
    return {
        str(path.resolve()): sha256_file(path)
        for path in paths
        if path.is_file()
    }


class ExperimentRunner:
    """Run, checkpoint, resume and analyze the fixed phase-one protocol."""

    def __init__(
        self,
        config: ExperimentConfig,
        output_dir: Path,
        *,
        resume: bool = True,
    ) -> None:
        self.config = config
        self.output_dir = Path(output_dir).resolve()
        self.resume = resume
        self.manifest_path = self.output_dir / MANIFEST_NAME
        self.checkpoint_dir = self.output_dir / "checkpoints"
        self._invocation_started = time.monotonic()
        self._soft_deadline = (
            self._invocation_started + config.soft_time_budget_seconds
        )
        self._verified_inputs = config.verify_inputs()
        self._code_sources = _code_provenance(config.repo_root)
        self._run_key = sha256_json(
            {
                "experiment_id": config.experiment_id,
                "config_sha256": config.config_sha256,
                "inputs": self._verified_inputs,
                "package_version": __version__,
                "code_sources": self._code_sources,
            }
        )
        self._protected_before = _hash_paths(_protected_paths(config))
        self._manifest = self._prepare_manifest()
        self.simulator = R4R1846CapacitySimulator(config)
        self._manifest["target_fingerprints"] = dict(
            self.simulator.target_fingerprints
        )
        self._manifest["solver"] = self.simulator.solver_metadata
        self._manifest["baseline_invariants"] = self.simulator.baseline_invariants
        self._write_manifest()

    def _prepare_manifest(self) -> dict[str, object]:
        existing: dict[str, object] | None = None
        if self.manifest_path.is_file():
            import json

            existing = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            if existing.get("run_key") != self._run_key:
                raise RuntimeError(
                    "Resume refused: configuration, package version, or input SHA "
                    "differs from the existing manifest"
                )
            if not self.resume:
                raise FileExistsError(
                    f"Output directory already contains a run: {self.output_dir}"
                )
        elif self.output_dir.exists() and any(self.output_dir.iterdir()):
            raise RuntimeError(
                "Resume refused: non-empty output directory has no provenance manifest"
            )
        self.output_dir.mkdir(parents=True, exist_ok=True)
        previous_elapsed = (
            float(existing.get("elapsed_seconds", 0.0)) if existing else 0.0
        )
        return {
            "schema_version": 1,
            "analysis_type": "iyali26_R4_R1846_two_parameter_phase1",
            "experiment_id": self.config.experiment_id,
            "status": "running",
            "run_key": self._run_key,
            "package_version": __version__,
            "seed": self.config.seed,
            "config": {
                "path": str(self.config.source_path),
                "sha256": self.config.config_sha256,
            },
            "inputs": self._verified_inputs,
            "code_sources": self._code_sources,
            "material_passport": {
                "origin_skill": "experiment-agent",
                "origin_mode": "plan",
                "origin_date": "2026-07-20",
                "verification_status": "RUNNING",
                "version_label": "code_plan_v1",
            },
            "exploratory_only": True,
            "synthetic_truth_only": True,
            "model_xml_written": False,
            "media_modified": False,
            "gpr_modified": False,
            "candidate_profile_created": False,
            "accepted_case_created": False,
            "curated_patch_modified": False,
            "protected_files_before": self._protected_before,
            "started_at": existing.get("started_at", utc_now()) if existing else utc_now(),
            "updated_at": utc_now(),
            "elapsed_seconds": previous_elapsed,
            "invocation_count": int(existing.get("invocation_count", 0)) + 1
            if existing
            else 1,
        }

    def _write_manifest(self) -> None:
        self._manifest["updated_at"] = utc_now()
        atomic_write_json(self.manifest_path, self._manifest)

    def _checkpoint_paths(self, stage: str) -> tuple[Path, Path]:
        return (
            self.output_dir / f"{stage}.tsv",
            self.checkpoint_dir / f"{stage}.json",
        )

    def _stage_key(
        self,
        stage: str,
        points: Sequence[tuple[str, ParameterPoint]],
        fidelity: Fidelity,
    ) -> str:
        return sha256_json(
            {
                "run_key": self._run_key,
                "stage": stage,
                "fidelity": fidelity,
                "points": [
                    {"sample_id": sample_id, "theta": point.as_tuple()}
                    for sample_id, point in points
                ],
            }
        )

    def _load_stage(
        self,
        stage: str,
        points: Sequence[tuple[str, ParameterPoint]],
        fidelity: Fidelity,
    ) -> tuple[list[SimulationResult], bool]:
        import json

        table_path, checkpoint_path = self._checkpoint_paths(stage)
        if not table_path.exists() and not checkpoint_path.exists():
            return [], False
        if not table_path.is_file() or not checkpoint_path.is_file():
            raise RuntimeError(f"Incomplete checkpoint metadata for stage {stage}")
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        expected_stage_key = self._stage_key(stage, points, fidelity)
        if checkpoint.get("run_key") != self._run_key:
            raise RuntimeError(f"Run-key mismatch in {stage} checkpoint")
        if checkpoint.get("stage_key") != expected_stage_key:
            raise RuntimeError(f"Point-design mismatch in {stage} checkpoint")
        if checkpoint.get("table_sha256") != sha256_file(table_path):
            raise RuntimeError(f"Table SHA mismatch in {stage} checkpoint")
        records = _read_tsv(table_path)
        results = [
            _result_from_record(record, self.simulator.backup_genes)
            for record in records
        ]
        expected_ids = [sample_id for sample_id, _ in points]
        actual_ids = [result.sample_id for result in results]
        if actual_ids != expected_ids[: len(actual_ids)]:
            raise RuntimeError(f"Non-prefix sample sequence in {stage} checkpoint")
        complete = checkpoint.get("status") == "complete"
        if complete and len(results) != len(points):
            raise RuntimeError(f"Completed {stage} checkpoint has the wrong row count")
        return results, complete

    def _write_stage(
        self,
        stage: str,
        points: Sequence[tuple[str, ParameterPoint]],
        fidelity: Fidelity,
        results: Sequence[SimulationResult],
        *,
        complete: bool,
    ) -> None:
        table_path, checkpoint_path = self._checkpoint_paths(stage)
        atomic_write_tsv(table_path, (result.to_record() for result in results))
        atomic_write_json(
            checkpoint_path,
            {
                "schema_version": 1,
                "run_key": self._run_key,
                "stage": stage,
                "stage_key": self._stage_key(stage, points, fidelity),
                "fidelity": fidelity,
                "status": "complete" if complete else "running",
                "completed_count": len(results),
                "expected_count": len(points),
                "table_path": str(table_path),
                "table_sha256": sha256_file(table_path),
                "updated_at": utc_now(),
            },
        )

    def _check_soft_budget(self, stage: str, completed: int, total: int) -> None:
        if time.monotonic() >= self._soft_deadline:
            raise SoftTimeBudgetExceeded(
                f"Soft time budget reached in {stage} ({completed}/{total}); "
                "checkpoint is resumable"
            )

    def _run_stage(
        self,
        stage: str,
        points: Sequence[tuple[str, ParameterPoint]],
        *,
        fidelity: Fidelity,
    ) -> list[SimulationResult]:
        results, complete = self._load_stage(stage, points, fidelity)
        if complete:
            print(
                f"[iyali26-flow] {stage}: resumed {len(results)}/{len(points)} complete",
                flush=True,
            )
            return results
        print(
            f"[iyali26-flow] {stage}: starting at {len(results)}/{len(points)}",
            flush=True,
        )
        try:
            for sample_id, point in points[len(results) :]:
                self._check_soft_budget(stage, len(results), len(points))
                results.append(
                    self.simulator.simulate(
                        point,
                        sample_id=sample_id,
                        fidelity=fidelity,
                    )
                )
                if (
                    len(results) % self.config.design.checkpoint_every == 0
                    or len(results) == len(points)
                ):
                    self._write_stage(
                        stage,
                        points,
                        fidelity,
                        results,
                        complete=len(results) == len(points),
                    )
                    print(
                        f"[iyali26-flow] {stage}: checkpoint "
                        f"{len(results)}/{len(points)}",
                        flush=True,
                    )
        except BaseException:
            if results:
                self._write_stage(
                    stage,
                    points,
                    fidelity,
                    results,
                    complete=False,
                )
            raise
        return results

    def _reference_anchors(
        self,
    ) -> tuple[pd.DataFrame, list[tuple[str, ParameterPoint]]]:
        path = self.config.input_map["anchor_reference"].path
        reference = pd.read_csv(path, sep="\t")
        required = {
            "scenario_id",
            "r4_capacity_fraction_of_wt_flux",
            "r1846_capacity_fraction_of_wt_flux",
            "r4_applied_backup_upper_bound",
            "r1846_applied_backup_upper_bound",
            "r4_target_ko_ratio",
            "r1846_target_ko_ratio",
            "feasible",
            "target_window_pass",
            "wt_growth_delta",
            "wt_growth_gate_pass",
            "r1843_bounds_gate_pass",
            "reaction_family_total_upper_bound_gate_pass",
            "non_target_call_flip_count",
            "no_non_target_call_flips_gate_pass",
            "backup_gene_ko_ratio_gate_pass",
            "proxy_new_essential_count",
            "concordant_nonessential_proxy_gate_pass",
        }
        missing = required - set(reference.columns)
        if missing:
            raise ValueError(f"Reference anchor grid is missing {sorted(missing)}")
        if len(reference) != 9 or reference["scenario_id"].duplicated().any():
            raise ValueError("Reference anchor grid must contain nine unique scenarios")
        reference = reference.sort_values("scenario_id").reset_index(drop=True)
        anchor_grid = pd.read_csv(self.config.input_map["anchor_grid"].path)
        expected_points = {
            (
                str(row.scenario_id),
                float(row.r4_capacity_fraction_of_wt_flux),
                float(row.r1846_capacity_fraction_of_wt_flux),
            )
            for row in anchor_grid.itertuples(index=False)
        }
        reference_points = {
            (
                str(row.scenario_id),
                float(row.r4_capacity_fraction_of_wt_flux),
                float(row.r1846_capacity_fraction_of_wt_flux),
            )
            for row in reference.itertuples(index=False)
        }
        if expected_points != reference_points or len(expected_points) != 9:
            raise ValueError("Anchor design and frozen anchor reference do not match")
        points = [
            (
                f"anchor__{row.scenario_id}",
                ParameterPoint(
                    float(row.r4_capacity_fraction_of_wt_flux),
                    float(row.r1846_capacity_fraction_of_wt_flux),
                ),
            )
            for row in reference.itertuples(index=False)
        ]
        return reference, points

    def _compare_anchors(
        self,
        reference: pd.DataFrame,
        results: Sequence[SimulationResult],
    ) -> tuple[list[dict[str, object]], list[str]]:
        tolerance = self.config.thresholds.anchor_tolerance
        result_by_id = {
            result.sample_id.removeprefix("anchor__"): result for result in results
        }
        comparisons: list[dict[str, object]] = []
        invalid: list[str] = []
        scalar_pairs = {
            "r4_applied_backup_upper_bound": "r4_applied_upper_bound",
            "r1846_applied_backup_upper_bound": "r1846_applied_upper_bound",
            "r4_target_ko_ratio": "r4_target_ko_ratio",
            "r1846_target_ko_ratio": "r1846_target_ko_ratio",
            "wt_growth_delta": "wt_growth_delta",
        }
        bool_pairs = {
            "feasible": "full_feasible",
            "target_window_pass": "target_window_pass",
            "r1843_bounds_gate_pass": "r1843_gate_pass",
            "reaction_family_total_upper_bound_gate_pass": "family_total_gate_pass",
            "no_non_target_call_flips_gate_pass": (
                "no_non_target_call_flips_gate_pass"
            ),
            "backup_gene_ko_ratio_gate_pass": "backup_gate_pass",
            "concordant_nonessential_proxy_gate_pass": "proxy_gate_pass",
        }
        for row in reference.itertuples(index=False):
            scenario_id = str(row.scenario_id)
            result = result_by_id[scenario_id]
            differences = {
                reference_name: abs(
                    float(getattr(row, reference_name))
                    - float(getattr(result, result_name))
                )
                for reference_name, result_name in scalar_pairs.items()
            }
            boolean_mismatches = [
                reference_name
                for reference_name, result_name in bool_pairs.items()
                if bool(getattr(row, reference_name))
                != bool(getattr(result, result_name))
            ]
            count_mismatches = []
            for reference_name, result_name in (
                ("non_target_call_flip_count", "non_target_call_flip_count"),
                ("proxy_new_essential_count", "proxy_new_essential_count"),
            ):
                if int(getattr(row, reference_name)) != int(
                    getattr(result, result_name)
                ):
                    count_mismatches.append(reference_name)
            wt_gate_mismatch = bool(row.wt_growth_gate_pass) != (
                result.wt_growth_delta <= 1e-8
            )
            max_difference = max(differences.values())
            passed = (
                max_difference <= tolerance
                and not boolean_mismatches
                and not count_mismatches
                and not wt_gate_mismatch
            )
            comparisons.append(
                {
                    "scenario_id": scenario_id,
                    "maximum_absolute_numeric_difference": max_difference,
                    "boolean_mismatches": ";".join(boolean_mismatches),
                    "count_mismatches": ";".join(count_mismatches),
                    "wt_growth_gate_mismatch": wt_gate_mismatch,
                    "anchor_pass": passed,
                }
            )
            if not passed:
                invalid.append(
                    f"anchor {scenario_id} did not reproduce within {tolerance:g}"
                )
        atomic_write_tsv(self.output_dir / "anchor_comparison.tsv", comparisons)
        return comparisons, invalid

    def _determinism_check(
        self,
        anchor_reference: pd.DataFrame,
    ) -> tuple[list[dict[str, object]], list[str]]:
        feasible = anchor_reference.loc[anchor_reference["feasible"].astype(bool)]
        selected = feasible.iloc[0] if not feasible.empty else anchor_reference.iloc[0]
        point = ParameterPoint(
            float(selected["r4_capacity_fraction_of_wt_flux"]),
            float(selected["r1846_capacity_fraction_of_wt_flux"]),
        )
        points = [(f"determinism__{index:02d}", point) for index in range(3)]
        results = self._run_stage("determinism", points, fidelity="targeted")
        baseline = results[0]
        tolerance = self.config.thresholds.anchor_tolerance
        records: list[dict[str, object]] = []
        invalid: list[str] = []
        baseline_backups = dict(baseline.backup_ko_ratios)
        for result in results:
            differences = [
                abs(result.r4_target_ko_ratio - baseline.r4_target_ko_ratio),
                abs(result.r1846_target_ko_ratio - baseline.r1846_target_ko_ratio),
                abs(result.scenario_wt_growth - baseline.scenario_wt_growth),
            ]
            differences.extend(
                abs(dict(result.backup_ko_ratios)[gene] - baseline_backups[gene])
                for gene in self.simulator.backup_genes
            )
            gates_match = (
                result.target_window_pass == baseline.target_window_pass
                and result.backup_gate_pass == baseline.backup_gate_pass
                and result.r1843_gate_pass == baseline.r1843_gate_pass
                and result.family_total_gate_pass == baseline.family_total_gate_pass
                and result.targeted_feasible == baseline.targeted_feasible
            )
            maximum = max(differences)
            passed = maximum <= tolerance and gates_match
            records.append(
                {
                    "sample_id": result.sample_id,
                    "maximum_absolute_difference_from_first": maximum,
                    "gates_match_first": gates_match,
                    "deterministic_within_tolerance": passed,
                }
            )
            if not passed:
                invalid.append(
                    f"targeted simulator repeat {result.sample_id} is non-deterministic"
                )
        atomic_write_tsv(self.output_dir / "determinism_check.tsv", records)
        return records, invalid

    def _audit_selection(
        self,
        quick_results: Sequence[SimulationResult],
    ) -> tuple[list[tuple[str, ParameterPoint]], dict[str, SimulationResult]]:
        feasible = [result for result in quick_results if result.targeted_feasible]
        nonfeasible = [
            result for result in quick_results if not result.targeted_feasible
        ]
        cap = self.config.design.full_feasible_cap
        if len(feasible) > cap:
            rng = np.random.default_rng(self.config.seed + 401)
            indices = sorted(rng.choice(len(feasible), size=cap, replace=False))
            selected_feasible = [feasible[int(index)] for index in indices]
        else:
            selected_feasible = feasible

        def distance_to_window(result: SimulationResult) -> tuple[float, str]:
            distance_squared = 0.0
            for ratio in (
                result.r4_target_ko_ratio,
                result.r1846_target_ko_ratio,
            ):
                if ratio < 0.05:
                    distance_squared += (0.05 - ratio) ** 2
                elif ratio >= 0.10:
                    distance_squared += (ratio - 0.10) ** 2
            return math.sqrt(distance_squared), result.sample_id

        selected_nonfeasible = sorted(nonfeasible, key=distance_to_window)[
            : self.config.design.full_nonfeasible_count
        ]
        selected = [*selected_feasible, *selected_nonfeasible]
        source_by_audit = {
            f"audit__{result.sample_id}": result for result in selected
        }
        points = [
            (audit_id, result.point)
            for audit_id, result in source_by_audit.items()
        ]
        return points, source_by_audit

    def _compare_audits(
        self,
        audits: Sequence[SimulationResult],
        source_by_audit: dict[str, SimulationResult],
    ) -> tuple[list[dict[str, object]], list[str]]:
        tolerance = self.config.thresholds.anchor_tolerance
        records: list[dict[str, object]] = []
        invalid: list[str] = []
        for audit in audits:
            source = source_by_audit[audit.sample_id]
            differences = [
                abs(audit.r4_target_ko_ratio - source.r4_target_ko_ratio),
                abs(audit.r1846_target_ko_ratio - source.r1846_target_ko_ratio),
                abs(audit.scenario_wt_growth - source.scenario_wt_growth),
            ]
            audit_backups = dict(audit.backup_ko_ratios)
            differences.extend(
                abs(audit_backups[gene] - ratio)
                for gene, ratio in source.backup_ko_ratios
            )
            global_safety_pass = bool(
                audit.no_non_target_call_flips_gate_pass
                and audit.proxy_gate_pass
            )
            fidelity_match = audit.full_feasible == source.targeted_feasible
            maximum = max(differences)
            passed = (
                maximum <= tolerance
                and global_safety_pass
                and fidelity_match
                and audit.r1843_gate_pass
                and audit.family_total_gate_pass
            )
            records.append(
                {
                    "audit_sample_id": audit.sample_id,
                    "source_sample_id": source.sample_id,
                    "source_targeted_feasible": source.targeted_feasible,
                    "full_feasible": audit.full_feasible,
                    "maximum_targeted_full_absolute_difference": maximum,
                    "no_non_target_call_flips_gate_pass": (
                        audit.no_non_target_call_flips_gate_pass
                    ),
                    "concordant_nonessential_proxy_gate_pass": audit.proxy_gate_pass,
                    "fidelity_feasibility_match": fidelity_match,
                    "full_audit_pass": passed,
                }
            )
            if not passed:
                invalid.append(f"full audit failed for {source.sample_id}")
        atomic_write_tsv(self.output_dir / AUDIT_COMPARISON_NAME, records)
        return records, invalid

    def _write_analysis(
        self,
        *,
        design: Sequence[SimulationResult],
        holdout: Sequence[SimulationResult],
        lattice: Sequence[SimulationResult],
        invalid_reasons: Sequence[str],
    ) -> tuple[Phase1Decision, AffineDiagnostics, dict[str, object]]:
        affine = affine_holdout_diagnostics(design, holdout)
        truths = holdout[: self.config.design.synthetic_truths]
        posterior_records = grid_posterior_recovery(
            lattice,
            truths,
            grid_size=self.config.design.reference_grid_size,
            sigmas=self.config.observation_noise_sigma,
            seed=self.config.seed + 301,
            component_minimum_mass=(
                self.config.thresholds.component_minimum_mass
            ),
        )
        atomic_write_tsv(self.output_dir / POSTERIOR_NAME, posterior_records)
        posterior = posterior_summary(posterior_records)
        decision = make_phase1_decision(
            self.config,
            invalid_baseline_reasons=invalid_reasons,
            affine=affine,
            posterior=posterior,
        )
        metrics = {
            "schema_version": 1,
            "observation_vector": [
                "R4 target KO growth ratio",
                "R1846 target KO growth ratio",
            ],
            "affine_holdout": affine.to_dict(),
            "grid_posterior": posterior,
            "decision_thresholds": {
                "affine_max_normalized_error": (
                    self.config.thresholds.affine_max_normalized_error
                ),
                "multimodal_fraction": self.config.thresholds.multimodal_fraction,
                "component_minimum_mass": (
                    self.config.thresholds.component_minimum_mass
                ),
                "coverage_interval": [
                    self.config.thresholds.coverage_lower,
                    self.config.thresholds.coverage_upper,
                ],
            },
        }
        atomic_write_json(self.output_dir / BASELINE_METRICS_NAME, metrics)
        atomic_write_json(self.output_dir / DECISION_NAME, decision.to_dict())
        return decision, affine, posterior

    def _finalize(
        self,
        *,
        status: str,
        decision: Phase1Decision | None,
        error: str | None = None,
    ) -> ExperimentOutcome:
        protected_after = _hash_paths(_protected_paths(self.config))
        changed = [
            path
            for path, before in self._protected_before.items()
            if protected_after.get(path) != before
        ]
        elapsed = time.monotonic() - self._invocation_started
        self._manifest["elapsed_seconds"] = float(
            self._manifest.get("elapsed_seconds", 0.0)
        ) + elapsed
        self._manifest["last_invocation_elapsed_seconds"] = elapsed
        self._manifest["protected_files_after"] = protected_after
        self._manifest["protected_files_changed"] = changed
        self._manifest["model_sha256_before"] = self.simulator.model_sha256
        self._manifest["model_sha256_after"] = sha256_file(
            self.config.input_map["model"].path
        )
        if changed:
            status = "failed_protected_input_changed"
            error = f"Protected files changed during experiment: {changed}"
        self._manifest["status"] = status
        if decision is not None:
            self._manifest["phase1_decision"] = decision.to_dict()
        if error is not None:
            self._manifest["error"] = error
        self._manifest["outputs"] = output_inventory(
            path
            for path in self.output_dir.rglob("*")
            if path.is_file() and path != self.manifest_path
        )
        passport = dict(self._manifest["material_passport"])
        passport["verification_status"] = (
            "VERIFIED" if status == "complete" and not changed else "PARTIAL"
        )
        self._manifest["material_passport"] = passport
        if status == "complete":
            self._manifest["completed_at"] = utc_now()
        self._write_manifest()
        if changed:
            raise RuntimeError(error)
        return ExperimentOutcome(status, self.output_dir, self.manifest_path, decision)

    def run_phase1(self) -> ExperimentOutcome:
        if self._manifest.get("status") == "complete" and (
            self.output_dir / DECISION_NAME
        ).is_file():
            return self.analyze_only()
        try:
            reference, anchor_points = self._reference_anchors()
            anchors = self._run_stage("anchors", anchor_points, fidelity="full")
            _, anchor_invalid = self._compare_anchors(reference, anchors)
            _, determinism_invalid = self._determinism_check(reference)

            design_points = sobol_parameter_points(
                self.config,
                count=self.config.design.sobol_samples,
                seed=self.config.seed,
            )
            design_spec = [
                (f"design__{index:04d}", point)
                for index, point in enumerate(design_points)
            ]
            design = self._run_stage("design", design_spec, fidelity="targeted")

            holdout_points = sobol_parameter_points(
                self.config,
                count=self.config.design.holdout_samples,
                seed=self.config.seed + 1,
            )
            holdout_spec = [
                (f"holdout__{index:04d}", point)
                for index, point in enumerate(holdout_points)
            ]
            holdout = self._run_stage("holdout", holdout_spec, fidelity="targeted")

            lattice_points = reference_lattice(self.config)
            lattice_spec = [
                (f"lattice__{index:05d}", point)
                for index, point in enumerate(lattice_points)
            ]
            lattice = self._run_stage(
                "reference_lattice",
                lattice_spec,
                fidelity="targeted",
            )

            audit_points, source_by_audit = self._audit_selection(
                [*design, *holdout]
            )
            audits = self._run_stage("full_audit", audit_points, fidelity="full")
            _, audit_invalid = self._compare_audits(audits, source_by_audit)

            self.simulator.assert_input_model_unchanged()
            protected_now = _hash_paths(_protected_paths(self.config))
            protected_invalid = [
                f"protected file changed during phase one: {path}"
                for path, digest in self._protected_before.items()
                if protected_now.get(path) != digest
            ]
            decision, _, _ = self._write_analysis(
                design=design,
                holdout=holdout,
                lattice=lattice,
                invalid_reasons=[
                    *anchor_invalid,
                    *determinism_invalid,
                    *audit_invalid,
                    *protected_invalid,
                ],
            )
            return self._finalize(status="complete", decision=decision)
        except ExperimentTimeoutError as exc:
            return self._finalize(
                status="partial_timeout",
                decision=None,
                error=str(exc),
            )
        except BaseException as exc:
            self._finalize(status="failed", decision=None, error=repr(exc))
            raise

    def analyze_only(self) -> ExperimentOutcome:
        """Recompute phase-one diagnostics from completed, SHA-checked stages."""
        try:
            reference, anchor_points = self._reference_anchors()
            anchors, anchor_complete = self._load_stage(
                "anchors", anchor_points, "full"
            )
            if not anchor_complete:
                raise RuntimeError("Cannot analyze before anchors are complete")
            _, anchor_invalid = self._compare_anchors(reference, anchors)

            determinism_point_row = reference.loc[
                reference["feasible"].astype(bool)
            ]
            selected = (
                determinism_point_row.iloc[0]
                if not determinism_point_row.empty
                else reference.iloc[0]
            )
            determinism_point = ParameterPoint(
                float(selected["r4_capacity_fraction_of_wt_flux"]),
                float(selected["r1846_capacity_fraction_of_wt_flux"]),
            )
            determinism_spec = [
                (f"determinism__{index:02d}", determinism_point)
                for index in range(3)
            ]
            determinism, determinism_complete = self._load_stage(
                "determinism", determinism_spec, "targeted"
            )
            if not determinism_complete:
                raise RuntimeError("Cannot analyze before determinism stage is complete")
            baseline = determinism[0]
            determinism_invalid = []
            for result in determinism:
                if max(
                    abs(result.r4_target_ko_ratio - baseline.r4_target_ko_ratio),
                    abs(
                        result.r1846_target_ko_ratio
                        - baseline.r1846_target_ko_ratio
                    ),
                ) > self.config.thresholds.anchor_tolerance:
                    determinism_invalid.append(
                        f"targeted simulator repeat {result.sample_id} is non-deterministic"
                    )

            design_points = sobol_parameter_points(
                self.config,
                count=self.config.design.sobol_samples,
                seed=self.config.seed,
            )
            design_spec = [
                (f"design__{index:04d}", point)
                for index, point in enumerate(design_points)
            ]
            design, design_complete = self._load_stage(
                "design", design_spec, "targeted"
            )
            holdout_points = sobol_parameter_points(
                self.config,
                count=self.config.design.holdout_samples,
                seed=self.config.seed + 1,
            )
            holdout_spec = [
                (f"holdout__{index:04d}", point)
                for index, point in enumerate(holdout_points)
            ]
            holdout, holdout_complete = self._load_stage(
                "holdout", holdout_spec, "targeted"
            )
            lattice_points = reference_lattice(self.config)
            lattice_spec = [
                (f"lattice__{index:05d}", point)
                for index, point in enumerate(lattice_points)
            ]
            lattice, lattice_complete = self._load_stage(
                "reference_lattice", lattice_spec, "targeted"
            )
            if not all((design_complete, holdout_complete, lattice_complete)):
                raise RuntimeError("Cannot analyze before targeted stages are complete")
            audit_points, source_by_audit = self._audit_selection(
                [*design, *holdout]
            )
            audits, audit_complete = self._load_stage(
                "full_audit", audit_points, "full"
            )
            if not audit_complete:
                raise RuntimeError("Cannot analyze before full audits are complete")
            _, audit_invalid = self._compare_audits(audits, source_by_audit)
            decision, _, _ = self._write_analysis(
                design=design,
                holdout=holdout,
                lattice=lattice,
                invalid_reasons=[
                    *anchor_invalid,
                    *determinism_invalid,
                    *audit_invalid,
                ],
            )
            return self._finalize(status="complete", decision=decision)
        except BaseException as exc:
            self._finalize(status="failed", decision=None, error=repr(exc))
            raise
