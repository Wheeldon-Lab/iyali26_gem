"""Deterministic baseline and grid-posterior analysis for phase one."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from scipy import ndimage
from scipy.stats import qmc

from .config import ExperimentConfig
from .core import ParameterPoint, Phase1Decision, SimulationResult


@dataclass(frozen=True, slots=True)
class AffineDiagnostics:
    coefficients: tuple[tuple[float, ...], ...]
    maximum_normalized_error: float
    median_normalized_error: float
    r_squared: tuple[float, float]

    def to_dict(self) -> dict[str, object]:
        return {
            "coefficients": [list(row) for row in self.coefficients],
            "maximum_normalized_error": self.maximum_normalized_error,
            "median_normalized_error": self.median_normalized_error,
            "r_squared": list(self.r_squared),
        }


def sobol_parameter_points(
    config: ExperimentConfig,
    *,
    count: int,
    seed: int,
) -> tuple[ParameterPoint, ...]:
    """Generate a scrambled, reproducible power-of-two Sobol design."""
    if count <= 0 or count & (count - 1):
        raise ValueError("Sobol sample count must be a positive power of two")
    exponent = int(math.log2(count))
    unit = qmc.Sobol(d=2, scramble=True, seed=seed).random_base2(exponent)
    lower = np.array(
        [config.r4_bounds.lower, config.r1846_bounds.lower],
        dtype=float,
    )
    upper = np.array(
        [config.r4_bounds.upper, config.r1846_bounds.upper],
        dtype=float,
    )
    scaled = qmc.scale(unit, lower, upper)
    return tuple(ParameterPoint(float(row[0]), float(row[1])) for row in scaled)


def reference_lattice(config: ExperimentConfig) -> tuple[ParameterPoint, ...]:
    """Return the configured Cartesian lattice in row-major R4/R1846 order."""
    size = config.design.reference_grid_size
    r4_values = np.linspace(config.r4_bounds.lower, config.r4_bounds.upper, size)
    r1846_values = np.linspace(
        config.r1846_bounds.lower,
        config.r1846_bounds.upper,
        size,
    )
    return tuple(
        ParameterPoint(float(r4), float(r1846))
        for r4 in r4_values
        for r1846 in r1846_values
    )


def _parameter_matrix(results: Sequence[SimulationResult]) -> np.ndarray:
    return np.asarray([result.point.as_tuple() for result in results], dtype=float)


def _observation_matrix(results: Sequence[SimulationResult]) -> np.ndarray:
    return np.asarray(
        [
            (result.r4_target_ko_ratio, result.r1846_target_ko_ratio)
            for result in results
        ],
        dtype=float,
    )


def affine_holdout_diagnostics(
    design: Sequence[SimulationResult],
    holdout: Sequence[SimulationResult],
) -> AffineDiagnostics:
    """Fit an affine simulator and report errors normalized by output ranges."""
    if len(design) < 3 or not holdout:
        raise ValueError("Affine diagnostics require design and holdout simulations")
    train_theta = _parameter_matrix(design)
    train_y = _observation_matrix(design)
    test_theta = _parameter_matrix(holdout)
    test_y = _observation_matrix(holdout)
    design_matrix = np.column_stack((np.ones(len(train_theta)), train_theta))
    coefficients, *_ = np.linalg.lstsq(design_matrix, train_y, rcond=None)
    predicted = np.column_stack((np.ones(len(test_theta)), test_theta)) @ coefficients
    ranges = np.ptp(train_y, axis=0)
    if np.any(ranges <= np.finfo(float).eps):
        raise ValueError("Cannot normalize affine error for a constant observation")
    normalized = np.abs(predicted - test_y) / ranges
    residual = test_y - predicted
    centered = test_y - np.mean(test_y, axis=0)
    denominator = np.sum(centered**2, axis=0)
    r_squared = np.where(
        denominator > 0,
        1.0 - np.sum(residual**2, axis=0) / denominator,
        1.0,
    )
    return AffineDiagnostics(
        coefficients=tuple(tuple(float(value) for value in row) for row in coefficients),
        maximum_normalized_error=float(np.max(normalized)),
        median_normalized_error=float(np.median(normalized)),
        r_squared=(float(r_squared[0]), float(r_squared[1])),
    )


def _hpd_mask(weights: np.ndarray, probability: float = 0.90) -> np.ndarray:
    flat = np.asarray(weights, dtype=float).ravel()
    order = np.argsort(-flat, kind="mergesort")
    cumulative = np.cumsum(flat[order])
    count = int(np.searchsorted(cumulative, probability, side="left")) + 1
    selected = np.zeros(flat.size, dtype=bool)
    selected[order[:count]] = True
    return selected.reshape(weights.shape)


def grid_posterior_recovery(
    lattice: Sequence[SimulationResult],
    truths: Sequence[SimulationResult],
    *,
    grid_size: int,
    sigmas: Iterable[float],
    seed: int,
    component_minimum_mass: float,
) -> list[dict[str, object]]:
    """Recover noisy synthetic truths with a uniform-prior lattice posterior."""
    if len(lattice) != grid_size * grid_size:
        raise ValueError("Lattice result count does not match grid_size squared")
    theta = _parameter_matrix(lattice)
    predictions = _observation_matrix(lattice)
    theta_ranges = np.ptp(theta, axis=0)
    if np.any(theta_ranges <= 0):
        raise ValueError("Reference lattice does not span both parameters")
    rng = np.random.default_rng(seed)
    records: list[dict[str, object]] = []
    structure = np.ones((3, 3), dtype=int)

    for sigma in sigmas:
        for truth in truths:
            truth_theta = np.asarray(truth.point.as_tuple(), dtype=float)
            truth_y = np.asarray(
                [truth.r4_target_ko_ratio, truth.r1846_target_ko_ratio],
                dtype=float,
            )
            observation = truth_y + rng.normal(0.0, float(sigma), size=2)
            log_weight = -0.5 * np.sum(
                ((predictions - observation) / float(sigma)) ** 2,
                axis=1,
            )
            log_weight -= np.max(log_weight)
            weights = np.exp(log_weight)
            weights /= np.sum(weights)
            weight_grid = weights.reshape((grid_size, grid_size))
            hpd = _hpd_mask(weight_grid, probability=0.90)
            labels, component_count = ndimage.label(hpd, structure=structure)
            component_masses = sorted(
                (
                    float(np.sum(weight_grid[labels == label]))
                    for label in range(1, component_count + 1)
                ),
                reverse=True,
            )
            material_components = sum(
                mass >= component_minimum_mass for mass in component_masses
            )
            nearest = int(np.argmin(np.sum((theta - truth_theta) ** 2, axis=1)))
            covered = bool(hpd.ravel()[nearest])
            posterior_mean = weights @ theta
            normalized_mean_error = np.abs(posterior_mean - truth_theta) / theta_ranges
            records.append(
                {
                    "truth_id": truth.sample_id,
                    "sigma": float(sigma),
                    "truth_r4_fraction": float(truth_theta[0]),
                    "truth_r1846_fraction": float(truth_theta[1]),
                    "truth_r4_ratio": float(truth_y[0]),
                    "truth_r1846_ratio": float(truth_y[1]),
                    "observed_r4_ratio": float(observation[0]),
                    "observed_r1846_ratio": float(observation[1]),
                    "posterior_mean_r4_fraction": float(posterior_mean[0]),
                    "posterior_mean_r1846_fraction": float(posterior_mean[1]),
                    "posterior_mean_normalized_error": float(
                        np.max(normalized_mean_error)
                    ),
                    "covered_by_90pct_hpd": covered,
                    "hpd_component_count": int(component_count),
                    "material_hpd_component_count": int(material_components),
                    "material_multimodal": material_components > 1,
                    "hpd_component_masses": ";".join(
                        f"{mass:.12g}" for mass in component_masses
                    ),
                    "nearest_lattice_index": nearest,
                }
            )
    return records


def posterior_summary(records: Sequence[dict[str, object]]) -> dict[str, object]:
    if not records:
        raise ValueError("Posterior summary requires at least one recovery record")
    sigmas = sorted({float(record["sigma"]) for record in records})
    coverage_by_sigma = {
        sigma: float(
            np.mean(
                [
                    bool(record["covered_by_90pct_hpd"])
                    for record in records
                    if float(record["sigma"]) == sigma
                ]
            )
        )
        for sigma in sigmas
    }
    multimodal_fraction = float(
        np.mean([bool(record["material_multimodal"]) for record in records])
    )
    return {
        "coverage_by_sigma": coverage_by_sigma,
        "multimodal_fraction": multimodal_fraction,
        "median_normalized_posterior_mean_error": float(
            np.median(
                [float(record["posterior_mean_normalized_error"]) for record in records]
            )
        ),
        "recovery_count": len(records),
    }


def make_phase1_decision(
    config: ExperimentConfig,
    *,
    invalid_baseline_reasons: Iterable[str],
    affine: AffineDiagnostics,
    posterior: dict[str, object],
) -> Phase1Decision:
    """Apply the plan's mandatory baseline, coverage and complexity gates."""
    reasons = list(invalid_baseline_reasons)
    coverage_raw = posterior["coverage_by_sigma"]
    if not isinstance(coverage_raw, dict):
        raise TypeError("coverage_by_sigma must be a mapping")
    coverage = tuple(
        sorted((float(sigma), float(value)) for sigma, value in coverage_raw.items())
    )
    for sigma, value in coverage:
        if not config.thresholds.coverage_lower <= value <= config.thresholds.coverage_upper:
            reasons.append(
                f"90% reference posterior coverage at sigma={sigma:g} is {value:.6f}, "
                f"outside [{config.thresholds.coverage_lower:.2f}, "
                f"{config.thresholds.coverage_upper:.2f}]"
            )
    multimodal_fraction = float(posterior["multimodal_fraction"])
    if reasons:
        status = "blocked_invalid_baseline"
    elif (
        multimodal_fraction >= config.thresholds.multimodal_fraction
        or affine.maximum_normalized_error
        > config.thresholds.affine_max_normalized_error
    ):
        status = "proceed_fmpe"
        if multimodal_fraction >= config.thresholds.multimodal_fraction:
            reasons.append(
                f"material multimodal fraction {multimodal_fraction:.6f} meets "
                f"threshold {config.thresholds.multimodal_fraction:.6f}"
            )
        if (
            affine.maximum_normalized_error
            > config.thresholds.affine_max_normalized_error
        ):
            reasons.append(
                f"affine normalized max error "
                f"{affine.maximum_normalized_error:.6f} exceeds "
                f"{config.thresholds.affine_max_normalized_error:.6f}"
            )
    else:
        status = "stop_simple_baseline_sufficient"
        reasons.append(
            "reference posteriors are effectively unimodal and the affine "
            "holdout error is within the 1% threshold"
        )
    return Phase1Decision(
        status=status,
        reasons=tuple(reasons),
        affine_max_normalized_error=affine.maximum_normalized_error,
        multimodal_fraction=multimodal_fraction,
        coverage_by_sigma=coverage,
    )
