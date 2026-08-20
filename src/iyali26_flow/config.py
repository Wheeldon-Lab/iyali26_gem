"""Strict JSON configuration for the phase-one experiment."""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from scripts.gem_annotate.config import load_project_paths

from .provenance import sha256_file


SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_INPUTS = {
    "model",
    "profile",
    "media",
    "experimental",
    "assay_fitness",
    "anchor_grid",
    "anchor_reference",
}


def _require_exact_keys(
    value: Mapping[str, Any], expected: set[str], label: str
) -> None:
    actual = set(value)
    missing = expected - actual
    unknown = actual - expected
    if missing or unknown:
        raise ValueError(
            f"{label} keys mismatch: missing={sorted(missing)}, "
            f"unknown={sorted(unknown)}"
        )


def _as_positive_int(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a positive integer")
    parsed = int(value)
    if parsed != value or parsed <= 0:
        raise ValueError(f"{label} must be a positive integer")
    return parsed


def _as_nonnegative_int(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a non-negative integer")
    parsed = int(value)
    if parsed != value or parsed < 0:
        raise ValueError(f"{label} must be a non-negative integer")
    return parsed


def _as_finite_float(value: object, label: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"{label} must be finite")
    return parsed


@dataclass(frozen=True, slots=True)
class InputSpec:
    path: Path
    sha256: str


@dataclass(frozen=True, slots=True)
class Bounds:
    lower: float
    upper: float

    def __post_init__(self) -> None:
        if not (math.isfinite(self.lower) and math.isfinite(self.upper)):
            raise ValueError("parameter bounds must be finite")
        if self.lower <= 0 or self.lower >= self.upper:
            raise ValueError("parameter bounds must satisfy 0 < lower < upper")


@dataclass(frozen=True, slots=True)
class DesignConfig:
    sobol_samples: int
    holdout_samples: int
    reference_grid_size: int
    synthetic_truths: int
    full_feasible_cap: int
    full_nonfeasible_count: int
    checkpoint_every: int


@dataclass(frozen=True, slots=True)
class ThresholdConfig:
    anchor_tolerance: float
    affine_max_normalized_error: float
    multimodal_fraction: float
    component_minimum_mass: float
    coverage_lower: float
    coverage_upper: float


@dataclass(frozen=True, slots=True)
class ExperimentConfig:
    source_path: Path
    repo_root: Path
    schema_version: int
    experiment_id: str
    seed: int
    solver: str
    inputs: tuple[tuple[str, InputSpec], ...]
    r4_bounds: Bounds
    r1846_bounds: Bounds
    design: DesignConfig
    observation_noise_sigma: tuple[float, ...]
    growth_cutoffs: tuple[float, ...]
    thresholds: ThresholdConfig
    soft_time_budget_seconds: int
    hard_timeout_seconds: int
    fmpe: tuple[tuple[str, object], ...]

    @property
    def input_map(self) -> dict[str, InputSpec]:
        return dict(self.inputs)

    @property
    def fmpe_map(self) -> dict[str, object]:
        return dict(self.fmpe)

    @property
    def config_sha256(self) -> str:
        return sha256_file(self.source_path)

    def verify_inputs(self) -> dict[str, dict[str, str]]:
        verified: dict[str, dict[str, str]] = {}
        for name, spec in self.inputs:
            if not spec.path.is_file():
                raise FileNotFoundError(f"configured {name} input not found: {spec.path}")
            actual = sha256_file(spec.path)
            if actual != spec.sha256:
                raise ValueError(
                    f"configured {name} SHA is stale: expected {spec.sha256}, got {actual}"
                )
            verified[name] = {"path": str(spec.path), "sha256": actual}
        return verified


TOP_LEVEL_KEYS = {
    "schema_version",
    "experiment_id",
    "seed",
    "solver",
    "inputs",
    "parameter_space",
    "design",
    "observation_noise_sigma",
    "growth_cutoffs",
    "thresholds",
    "soft_time_budget_seconds",
    "hard_timeout_seconds",
    "fmpe",
}


def load_experiment_config(
    path: Path,
    *,
    repo_root: Path | None = None,
    research_root: Path | None = None,
) -> ExperimentConfig:
    source_path = Path(path).resolve()
    raw = json.loads(source_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("experiment config must be a JSON object")
    _require_exact_keys(raw, TOP_LEVEL_KEYS, "experiment config")

    root = Path.cwd().resolve() if repo_root is None else Path(repo_root).resolve()
    project_paths = load_project_paths(research_root, required=True)
    inputs_raw = raw["inputs"]
    if not isinstance(inputs_raw, dict):
        raise ValueError("inputs must be an object")
    _require_exact_keys(inputs_raw, REQUIRED_INPUTS, "inputs")
    inputs: list[tuple[str, InputSpec]] = []
    for name in sorted(REQUIRED_INPUTS):
        item = inputs_raw[name]
        if not isinstance(item, dict):
            raise ValueError(f"inputs.{name} must be an object")
        _require_exact_keys(item, {"path", "sha256"}, f"inputs.{name}")
        digest = str(item["sha256"]).strip().lower()
        if not SHA256_PATTERN.fullmatch(digest):
            raise ValueError(f"inputs.{name}.sha256 must be 64 lowercase hex digits")
        input_path = Path(str(item["path"]))
        if not input_path.is_absolute():
            input_path = project_paths.resolve_legacy_path(input_path)
        inputs.append((name, InputSpec(input_path.resolve(), digest)))

    parameter_space = raw["parameter_space"]
    if not isinstance(parameter_space, dict):
        raise ValueError("parameter_space must be an object")
    _require_exact_keys(
        parameter_space,
        {"r4_capacity_fraction", "r1846_capacity_fraction"},
        "parameter_space",
    )

    def parse_bounds(name: str) -> Bounds:
        values = parameter_space[name]
        if not isinstance(values, list) or len(values) != 2:
            raise ValueError(f"parameter_space.{name} must contain [lower, upper]")
        return Bounds(float(values[0]), float(values[1]))

    design_raw = raw["design"]
    if not isinstance(design_raw, dict):
        raise ValueError("design must be an object")
    design_keys = {
        "sobol_samples",
        "holdout_samples",
        "reference_grid_size",
        "synthetic_truths",
        "full_feasible_cap",
        "full_nonfeasible_count",
        "checkpoint_every",
    }
    _require_exact_keys(design_raw, design_keys, "design")
    design = DesignConfig(
        **{
            key: _as_positive_int(design_raw[key], f"design.{key}")
            for key in design_keys
        }
    )
    for label, count in (
        ("sobol_samples", design.sobol_samples),
        ("holdout_samples", design.holdout_samples),
    ):
        if count & (count - 1):
            raise ValueError(f"design.{label} must be a power of two")
    if design.reference_grid_size < 3 or design.reference_grid_size % 2 == 0:
        raise ValueError("design.reference_grid_size must be odd and at least 3")
    if design.synthetic_truths > design.holdout_samples:
        raise ValueError("synthetic_truths cannot exceed holdout_samples")
    if design.full_feasible_cap + design.full_nonfeasible_count > 128:
        raise ValueError("configured full-gene audit count cannot exceed 128")

    if not isinstance(raw["observation_noise_sigma"], list):
        raise ValueError("observation_noise_sigma must be a list")
    sigmas = tuple(
        _as_finite_float(value, "observation_noise_sigma")
        for value in raw["observation_noise_sigma"]
    )
    if not sigmas or any(value <= 0 for value in sigmas):
        raise ValueError("observation_noise_sigma values must be positive")
    if tuple(sorted(set(sigmas))) != sigmas:
        raise ValueError("observation_noise_sigma must be unique and increasing")

    cutoffs = tuple(
        _as_finite_float(value, "growth_cutoffs") for value in raw["growth_cutoffs"]
    )
    if not cutoffs or any(not 0 < value < 1 for value in cutoffs):
        raise ValueError("growth_cutoffs must lie strictly between zero and one")
    if tuple(sorted(set(cutoffs))) != cutoffs:
        raise ValueError("growth_cutoffs must be unique and increasing")

    thresholds_raw = raw["thresholds"]
    if not isinstance(thresholds_raw, dict):
        raise ValueError("thresholds must be an object")
    threshold_keys = {
        "anchor_tolerance",
        "affine_max_normalized_error",
        "multimodal_fraction",
        "component_minimum_mass",
        "coverage_lower",
        "coverage_upper",
    }
    _require_exact_keys(thresholds_raw, threshold_keys, "thresholds")
    thresholds = ThresholdConfig(
        **{
            key: _as_finite_float(thresholds_raw[key], f"thresholds.{key}")
            for key in threshold_keys
        }
    )
    if thresholds.anchor_tolerance <= 0:
        raise ValueError("anchor_tolerance must be positive")
    if not 0 < thresholds.affine_max_normalized_error < 1:
        raise ValueError("affine_max_normalized_error must lie between zero and one")
    for label in (
        "multimodal_fraction",
        "component_minimum_mass",
        "coverage_lower",
        "coverage_upper",
    ):
        if not 0 <= getattr(thresholds, label) <= 1:
            raise ValueError(f"thresholds.{label} must lie between zero and one")
    if thresholds.coverage_lower >= thresholds.coverage_upper:
        raise ValueError("coverage_lower must be less than coverage_upper")

    fmpe_raw = raw["fmpe"]
    if not isinstance(fmpe_raw, dict):
        raise ValueError("fmpe must be an object")
    fmpe_keys = {
        "simulations",
        "hidden_features",
        "epochs",
        "patience",
        "batch_size",
        "learning_rate",
        "gradient_clip",
        "seeds",
    }
    _require_exact_keys(fmpe_raw, fmpe_keys, "fmpe")
    fmpe_simulations = _as_positive_int(fmpe_raw["simulations"], "fmpe.simulations")
    if fmpe_simulations & (fmpe_simulations - 1):
        raise ValueError("fmpe.simulations must be a power of two")
    hidden_features = fmpe_raw["hidden_features"]
    if not isinstance(hidden_features, list) or not hidden_features:
        raise ValueError("fmpe.hidden_features must be a non-empty list")
    normalized_hidden = [
        _as_positive_int(value, "fmpe.hidden_features") for value in hidden_features
    ]
    epochs = _as_positive_int(fmpe_raw["epochs"], "fmpe.epochs")
    patience = _as_positive_int(fmpe_raw["patience"], "fmpe.patience")
    batch_size = _as_positive_int(fmpe_raw["batch_size"], "fmpe.batch_size")
    if patience > epochs:
        raise ValueError("fmpe.patience cannot exceed fmpe.epochs")
    if batch_size > fmpe_simulations:
        raise ValueError("fmpe.batch_size cannot exceed fmpe.simulations")
    learning_rate = _as_finite_float(fmpe_raw["learning_rate"], "fmpe.learning_rate")
    gradient_clip = _as_finite_float(fmpe_raw["gradient_clip"], "fmpe.gradient_clip")
    if learning_rate <= 0 or gradient_clip <= 0:
        raise ValueError("FMPE learning rate and gradient clip must be positive")
    seeds_raw = fmpe_raw["seeds"]
    if not isinstance(seeds_raw, list) or len(seeds_raw) != 3:
        raise ValueError("fmpe.seeds must contain exactly three seeds")
    seeds = [_as_nonnegative_int(value, "fmpe.seeds") for value in seeds_raw]
    if len(set(seeds)) != 3:
        raise ValueError("fmpe.seeds must be unique")
    normalized_fmpe: dict[str, object] = {
        "simulations": fmpe_simulations,
        "hidden_features": normalized_hidden,
        "epochs": epochs,
        "patience": patience,
        "batch_size": batch_size,
        "learning_rate": learning_rate,
        "gradient_clip": gradient_clip,
        "seeds": seeds,
    }

    soft_budget = _as_positive_int(
        raw["soft_time_budget_seconds"], "soft_time_budget_seconds"
    )
    hard_timeout = _as_positive_int(raw["hard_timeout_seconds"], "hard_timeout_seconds")
    if hard_timeout <= soft_budget:
        raise ValueError("hard_timeout_seconds must exceed soft_time_budget_seconds")

    schema_version = _as_positive_int(raw["schema_version"], "schema_version")
    if schema_version != 1:
        raise ValueError(f"unsupported experiment schema_version: {schema_version}")
    experiment_id = str(raw["experiment_id"]).strip()
    if not experiment_id:
        raise ValueError("experiment_id cannot be empty")
    solver = str(raw["solver"]).strip()
    if not solver:
        raise ValueError("solver cannot be empty")

    return ExperimentConfig(
        source_path=source_path,
        repo_root=root,
        schema_version=schema_version,
        experiment_id=experiment_id,
        seed=_as_nonnegative_int(raw["seed"], "seed"),
        solver=solver,
        inputs=tuple(inputs),
        r4_bounds=parse_bounds("r4_capacity_fraction"),
        r1846_bounds=parse_bounds("r1846_capacity_fraction"),
        design=design,
        observation_noise_sigma=sigmas,
        growth_cutoffs=cutoffs,
        thresholds=thresholds,
        soft_time_budget_seconds=soft_budget,
        hard_timeout_seconds=hard_timeout,
        fmpe=tuple(sorted(normalized_fmpe.items())),
    )
