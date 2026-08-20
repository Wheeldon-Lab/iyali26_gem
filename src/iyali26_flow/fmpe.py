"""Conditional LAMPE FMPE training, gated behind the phase-one decision."""

from __future__ import annotations

import copy
import csv
import importlib
import json
import math
import time
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from scipy import ndimage

from .analysis import _hpd_mask, sobol_parameter_points
from .config import ExperimentConfig
from .core import SimulationResult
from .experiment import _read_tsv, _result_from_record
from .provenance import atomic_write_json, atomic_write_tsv, sha256_file, sha256_json, utc_now
from .simulator import R4R1846CapacitySimulator


FMPE_EXTRA_ERROR = (
    "FMPE optional dependencies are not installed. Install them only after a "
    "phase-one proceed_fmpe decision with `uv sync --extra fmpe`."
)


def require_fmpe_dependencies(
    importer: Callable[[str], object] = importlib.import_module,
) -> tuple[object, object, object]:
    """Import the optional backend or raise one actionable, stable error."""
    try:
        torch = importer("torch")
        inference = importer("lampe.inference")
        nn = importer("torch.nn")
    except (ImportError, ModuleNotFoundError) as exc:
        raise RuntimeError(FMPE_EXTRA_ERROR) from exc
    return torch, inference, nn


def bounded_logit(
    values: np.ndarray | Sequence[float],
    lower: np.ndarray | Sequence[float],
    upper: np.ndarray | Sequence[float],
) -> np.ndarray:
    values_array = np.asarray(values, dtype=float)
    lower_array = np.asarray(lower, dtype=float)
    upper_array = np.asarray(upper, dtype=float)
    if np.any(upper_array <= lower_array):
        raise ValueError("Every upper bound must exceed its lower bound")
    unit = (values_array - lower_array) / (upper_array - lower_array)
    if np.any((unit <= 0.0) | (unit >= 1.0)):
        raise ValueError("bounded_logit values must lie strictly inside the bounds")
    return np.log(unit) - np.log1p(-unit)


def bounded_sigmoid(
    logits: np.ndarray | Sequence[float],
    lower: np.ndarray | Sequence[float],
    upper: np.ndarray | Sequence[float],
) -> np.ndarray:
    logits_array = np.asarray(logits, dtype=float)
    lower_array = np.asarray(lower, dtype=float)
    upper_array = np.asarray(upper, dtype=float)
    if np.any(upper_array <= lower_array):
        raise ValueError("Every upper bound must exceed its lower bound")
    positive = logits_array >= 0
    unit = np.empty_like(logits_array, dtype=float)
    unit[positive] = 1.0 / (1.0 + np.exp(-logits_array[positive]))
    exponent = np.exp(logits_array[~positive])
    unit[~positive] = exponent / (1.0 + exponent)
    return lower_array + unit * (upper_array - lower_array)


def _decision_allows_training(output_dir: Path) -> dict[str, object]:
    path = output_dir / "phase1_decision.json"
    if not path.is_file():
        raise RuntimeError("Phase-one decision is missing; run `iyali26-flow phase1` first")
    decision = json.loads(path.read_text(encoding="utf-8"))
    if decision.get("status") != "proceed_fmpe":
        raise RuntimeError(
            "FMPE training is gated off because phase one did not return proceed_fmpe"
        )
    return decision


def _prior_checkpoint(
    config: ExperimentConfig,
    output_dir: Path,
    simulator: R4R1846CapacitySimulator,
) -> list[SimulationResult]:
    fmpe = config.fmpe_map
    count = int(fmpe["simulations"])
    points = sobol_parameter_points(config, count=count, seed=config.seed + 701)
    table_path = output_dir / "fmpe" / "prior_simulations.tsv"
    checkpoint_path = output_dir / "fmpe" / "prior_checkpoint.json"
    stage_key = sha256_json(
        {
            "config_sha256": config.config_sha256,
            "input_model_sha256": simulator.model_sha256,
            "count": count,
            "seed": config.seed + 701,
        }
    )
    results: list[SimulationResult] = []
    if table_path.is_file() or checkpoint_path.is_file():
        if not table_path.is_file() or not checkpoint_path.is_file():
            raise RuntimeError("Incomplete FMPE prior checkpoint")
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if checkpoint.get("stage_key") != stage_key:
            raise RuntimeError("FMPE prior checkpoint provenance mismatch")
        if checkpoint.get("table_sha256") != sha256_file(table_path):
            raise RuntimeError("FMPE prior checkpoint table SHA mismatch")
        results = [
            _result_from_record(record, simulator.backup_genes)
            for record in _read_tsv(table_path)
        ]
        if checkpoint.get("status") == "complete" and len(results) == count:
            return results
    started = time.monotonic()
    for index in range(len(results), count):
        if time.monotonic() - started >= config.soft_time_budget_seconds:
            raise RuntimeError("FMPE prior simulation reached the soft time budget")
        result = simulator.simulate(
            points[index],
            sample_id=f"fmpe_prior__{index:05d}",
            fidelity="targeted",
        )
        results.append(result)
        if len(results) % config.design.checkpoint_every == 0 or len(results) == count:
            atomic_write_tsv(table_path, (item.to_record() for item in results))
            atomic_write_json(
                checkpoint_path,
                {
                    "schema_version": 1,
                    "stage_key": stage_key,
                    "status": "complete" if len(results) == count else "running",
                    "completed_count": len(results),
                    "expected_count": count,
                    "table_sha256": sha256_file(table_path),
                    "updated_at": utc_now(),
                },
            )
            print(
                f"[iyali26-flow] fmpe_prior: checkpoint {len(results)}/{count}",
                flush=True,
            )
    return results


def _make_dataset(
    config: ExperimentConfig,
    results: Sequence[SimulationResult],
) -> tuple[np.ndarray, np.ndarray, dict[str, list[float]]]:
    theta = np.asarray([result.point.as_tuple() for result in results], dtype=float)
    noiseless = np.asarray(
        [
            (result.r4_target_ko_ratio, result.r1846_target_ko_ratio)
            for result in results
        ],
        dtype=float,
    )
    lower = np.asarray(
        [config.r4_bounds.lower, config.r1846_bounds.lower], dtype=float
    )
    upper = np.asarray(
        [config.r4_bounds.upper, config.r1846_bounds.upper], dtype=float
    )
    transformed = bounded_logit(theta, lower, upper)
    rng = np.random.default_rng(config.seed + 702)
    sigmas = rng.choice(np.asarray(config.observation_noise_sigma), size=len(results))
    noisy = noiseless + rng.normal(0.0, sigmas[:, None], size=noiseless.shape)
    condition = np.column_stack((noisy, np.log(sigmas)))
    theta_mean = np.mean(transformed, axis=0)
    theta_std = np.std(transformed, axis=0)
    condition_mean = np.mean(condition, axis=0)
    condition_std = np.std(condition, axis=0)
    standardized_theta = (transformed - theta_mean) / theta_std
    standardized_condition = (condition - condition_mean) / condition_std
    normalization = {
        "theta_mean": theta_mean.tolist(),
        "theta_std": theta_std.tolist(),
        "condition_mean": condition_mean.tolist(),
        "condition_std": condition_std.tolist(),
        "lower": lower.tolist(),
        "upper": upper.tolist(),
    }
    return standardized_theta, standardized_condition, normalization


def _network_builder(nn: object, hidden_features: Sequence[int]):
    def build(input_features: int, output_features: int):
        layers = []
        previous = input_features
        for features in hidden_features:
            layers.extend((nn.Linear(previous, int(features)), nn.ELU()))
            previous = int(features)
        layers.append(nn.Linear(previous, output_features))
        return nn.Sequential(*layers)

    return build


def _train_seed(
    *,
    config: ExperimentConfig,
    theta: np.ndarray,
    condition: np.ndarray,
    seed: int,
    torch: object,
    inference: object,
    nn: object,
    output_dir: Path,
    normalization: dict[str, list[float]],
):
    fmpe_config = config.fmpe_map
    torch.manual_seed(seed)
    estimator = inference.FMPE(
        theta_dim=2,
        x_dim=3,
        build=_network_builder(nn, fmpe_config["hidden_features"]),
    ).to("cpu")
    loss_module = inference.FMPELoss(estimator)
    optimizer = torch.optim.Adam(
        estimator.parameters(),
        lr=float(fmpe_config["learning_rate"]),
    )
    theta_tensor = torch.as_tensor(theta, dtype=torch.float32, device="cpu")
    condition_tensor = torch.as_tensor(
        condition,
        dtype=torch.float32,
        device="cpu",
    )
    generator = torch.Generator(device="cpu").manual_seed(seed)
    permutation = torch.randperm(len(theta_tensor), generator=generator)
    validation_count = max(1, len(theta_tensor) // 10)
    validation_indices = permutation[:validation_count]
    training_indices = permutation[validation_count:]
    best_loss = math.inf
    best_state = None
    patience_left = int(fmpe_config["patience"])
    history: list[dict[str, object]] = []
    batch_size = int(fmpe_config["batch_size"])

    for epoch in range(int(fmpe_config["epochs"])):
        estimator.train()
        order = training_indices[
            torch.randperm(len(training_indices), generator=generator)
        ]
        training_loss = 0.0
        training_batches = 0
        for start in range(0, len(order), batch_size):
            indices = order[start : start + batch_size]
            optimizer.zero_grad(set_to_none=True)
            loss = loss_module(theta_tensor[indices], condition_tensor[indices])
            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                estimator.parameters(),
                float(fmpe_config["gradient_clip"]),
            )
            optimizer.step()
            training_loss += float(loss.detach())
            training_batches += 1
        estimator.eval()
        with torch.no_grad():
            validation_loss = float(
                loss_module(
                    theta_tensor[validation_indices],
                    condition_tensor[validation_indices],
                )
            )
        history.append(
            {
                "epoch": epoch + 1,
                "training_loss": training_loss / training_batches,
                "validation_loss": validation_loss,
            }
        )
        if validation_loss < best_loss - 1e-8:
            best_loss = validation_loss
            best_state = copy.deepcopy(estimator.state_dict())
            patience_left = int(fmpe_config["patience"])
        else:
            patience_left -= 1
            if patience_left == 0:
                break
    if best_state is None:
        raise RuntimeError(f"FMPE seed {seed} produced no finite validation state")
    estimator.load_state_dict(best_state)
    model_path = output_dir / "fmpe" / f"fmpe_seed_{seed}.pt"
    model_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dict": best_state,
            "seed": seed,
            "normalization": normalization,
            "configuration": fmpe_config,
        },
        model_path,
    )
    atomic_write_tsv(
        output_dir / "fmpe" / f"training_history_seed_{seed}.tsv",
        history,
    )
    return estimator, {
        "seed": seed,
        "best_validation_loss": best_loss,
        "epochs_completed": len(history),
        "model_path": str(model_path),
        "model_sha256": sha256_file(model_path),
    }


def _evaluate_seed(
    *,
    config: ExperimentConfig,
    estimator: object,
    normalization: dict[str, list[float]],
    posterior_records: Sequence[dict[str, str]],
    torch: object,
    seed: int,
) -> dict[str, object]:
    rng_seed = seed + 10_000
    torch.manual_seed(rng_seed)
    lower = np.asarray(normalization["lower"], dtype=float)
    upper = np.asarray(normalization["upper"], dtype=float)
    theta_mean = np.asarray(normalization["theta_mean"], dtype=float)
    theta_std = np.asarray(normalization["theta_std"], dtype=float)
    condition_mean = np.asarray(normalization["condition_mean"], dtype=float)
    condition_std = np.asarray(normalization["condition_std"], dtype=float)
    ranges = upper - lower
    size = config.design.reference_grid_size
    component_matches: list[bool] = []
    covered: list[bool] = []
    grid_covered: list[bool] = []
    mean_errors: list[float] = []
    out_of_bounds = 0
    estimator.eval()
    for record in posterior_records:
        condition = np.asarray(
            [
                float(record["observed_r4_ratio"]),
                float(record["observed_r1846_ratio"]),
                math.log(float(record["sigma"])),
            ],
            dtype=float,
        )
        standardized = (condition - condition_mean) / condition_std
        with torch.no_grad():
            draws = (
                estimator.flow(
                    torch.as_tensor(standardized, dtype=torch.float32)
                )
                .sample((2048,))
                .cpu()
                .numpy()
            )
        logits = draws * theta_std + theta_mean
        samples = bounded_sigmoid(logits, lower, upper)
        out_of_bounds += int(
            np.sum((samples < lower[None, :]) | (samples > upper[None, :]))
        )
        posterior_mean = np.mean(samples, axis=0)
        grid_mean = np.asarray(
            [
                float(record["posterior_mean_r4_fraction"]),
                float(record["posterior_mean_r1846_fraction"]),
            ]
        )
        mean_errors.append(float(np.max(np.abs(posterior_mean - grid_mean) / ranges)))
        histogram, r4_edges, r1846_edges = np.histogram2d(
            samples[:, 0],
            samples[:, 1],
            bins=size,
            range=((lower[0], upper[0]), (lower[1], upper[1])),
        )
        weights = histogram / np.sum(histogram)
        hpd = _hpd_mask(weights, probability=0.90)
        labels, component_count = ndimage.label(hpd, structure=np.ones((3, 3)))
        material_count = sum(
            float(np.sum(weights[labels == label]))
            >= config.thresholds.component_minimum_mass
            for label in range(1, component_count + 1)
        )
        component_matches.append(
            material_count == int(record["material_hpd_component_count"])
        )
        truth = np.asarray(
            [float(record["truth_r4_fraction"]), float(record["truth_r1846_fraction"])]
        )
        r4_index = int(np.clip(np.searchsorted(r4_edges, truth[0], side="right") - 1, 0, size - 1))
        r1846_index = int(
            np.clip(
                np.searchsorted(r1846_edges, truth[1], side="right") - 1,
                0,
                size - 1,
            )
        )
        covered.append(bool(hpd[r4_index, r1846_index]))
        grid_covered.append(record["covered_by_90pct_hpd"].lower() == "true")
    coverage = float(np.mean(covered))
    grid_coverage = float(np.mean(grid_covered))
    return {
        "seed": seed,
        "coverage": coverage,
        "grid_coverage": grid_coverage,
        "coverage_deviation": abs(coverage - grid_coverage),
        "median_normalized_posterior_mean_error": float(np.median(mean_errors)),
        "hpd_component_match_fraction": float(np.mean(component_matches)),
        "out_of_bounds_sample_count": out_of_bounds,
    }


def train_fmpe(config: ExperimentConfig, output_dir: Path) -> dict[str, object]:
    """Train three CPU FMPE models after, and only after, a proceed decision."""
    output_dir = Path(output_dir).resolve()
    decision = _decision_allows_training(output_dir)
    torch, inference, nn = require_fmpe_dependencies()
    simulator = R4R1846CapacitySimulator(config)
    prior_results = _prior_checkpoint(config, output_dir, simulator)
    theta, condition, normalization = _make_dataset(config, prior_results)
    with (output_dir / "posterior_recovery.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        posterior_records = list(csv.DictReader(handle, delimiter="\t"))
    training_runs: list[dict[str, object]] = []
    evaluations: list[dict[str, object]] = []
    for seed_value in config.fmpe_map["seeds"]:
        seed = int(seed_value)
        estimator, training = _train_seed(
            config=config,
            theta=theta,
            condition=condition,
            seed=seed,
            torch=torch,
            inference=inference,
            nn=nn,
            output_dir=output_dir,
            normalization=normalization,
        )
        training_runs.append(training)
        evaluations.append(
            _evaluate_seed(
                config=config,
                estimator=estimator,
                normalization=normalization,
                posterior_records=posterior_records,
                torch=torch,
                seed=seed,
            )
        )
    metric_names = (
        "coverage_deviation",
        "median_normalized_posterior_mean_error",
        "hpd_component_match_fraction",
    )
    between_seed_spreads = {
        metric: max(float(item[metric]) for item in evaluations)
        - min(float(item[metric]) for item in evaluations)
        for metric in metric_names
    }
    accepted = all(
        (
            all(float(item["coverage_deviation"]) <= 0.05 for item in evaluations),
            all(
                float(item["median_normalized_posterior_mean_error"]) <= 0.05
                for item in evaluations
            ),
            all(
                float(item["hpd_component_match_fraction"]) >= 0.95
                for item in evaluations
            ),
            all(int(item["out_of_bounds_sample_count"]) == 0 for item in evaluations),
            all(spread <= 0.05 for spread in between_seed_spreads.values()),
        )
    )
    report = {
        "schema_version": 1,
        "status": "passed" if accepted else "failed",
        "phase1_decision": decision,
        "device": "cpu",
        "training_runs": training_runs,
        "evaluations": evaluations,
        "between_seed_absolute_spreads": between_seed_spreads,
        "acceptance_pass": accepted,
        "completed_at": utc_now(),
    }
    atomic_write_json(output_dir / "fmpe" / "fmpe_evaluation.json", report)
    simulator.assert_input_model_unchanged()
    return report
