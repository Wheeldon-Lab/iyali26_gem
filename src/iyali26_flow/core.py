"""Public immutable types for iYali26 parameter-inference experiments."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Protocol, runtime_checkable


DecisionStatus = Literal[
    "blocked_invalid_baseline",
    "proceed_fmpe",
    "stop_simple_baseline_sufficient",
]
Fidelity = Literal["targeted", "full"]


@dataclass(frozen=True, slots=True)
class ParameterPoint:
    """The two exploratory backup-capacity fractions used by the first adapter."""

    r4_fraction: float
    r1846_fraction: float

    def __post_init__(self) -> None:
        for name, value in (
            ("r4_fraction", self.r4_fraction),
            ("r1846_fraction", self.r1846_fraction),
        ):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")

    def as_tuple(self) -> tuple[float, float]:
        return (float(self.r4_fraction), float(self.r1846_fraction))


@dataclass(frozen=True, slots=True)
class SimulationResult:
    """One provenance-ready targeted or full-fidelity simulation result."""

    sample_id: str
    fidelity: Fidelity
    point: ParameterPoint
    r4_applied_upper_bound: float
    r1846_applied_upper_bound: float
    r4_target_ko_ratio: float
    r1846_target_ko_ratio: float
    backup_ko_ratios: tuple[tuple[str, float], ...]
    baseline_wt_growth: float
    scenario_wt_growth: float
    wt_growth_delta: float
    target_window_pass: bool
    backup_gate_pass: bool
    r1843_gate_pass: bool
    family_total_gate_pass: bool
    targeted_feasible: bool
    full_feasible: bool | None
    non_target_call_flip_count: int | None
    proxy_new_essential_count: int | None
    no_non_target_call_flips_gate_pass: bool | None
    proxy_gate_pass: bool | None
    gene_count: int
    elapsed_seconds: float

    def to_record(self) -> dict[str, object]:
        record: dict[str, object] = {
            "sample_id": self.sample_id,
            "fidelity": self.fidelity,
            "r4_capacity_fraction_of_wt_flux": self.point.r4_fraction,
            "r1846_capacity_fraction_of_wt_flux": self.point.r1846_fraction,
            "r4_applied_backup_upper_bound": self.r4_applied_upper_bound,
            "r1846_applied_backup_upper_bound": self.r1846_applied_upper_bound,
            "r4_target_ko_ratio": self.r4_target_ko_ratio,
            "r1846_target_ko_ratio": self.r1846_target_ko_ratio,
            "baseline_wt_growth": self.baseline_wt_growth,
            "scenario_wt_growth": self.scenario_wt_growth,
            "wt_growth_delta": self.wt_growth_delta,
            "target_window_pass": self.target_window_pass,
            "backup_gene_ko_ratio_gate_pass": self.backup_gate_pass,
            "r1843_bounds_gate_pass": self.r1843_gate_pass,
            "reaction_family_total_upper_bound_gate_pass": (
                self.family_total_gate_pass
            ),
            "targeted_feasible": self.targeted_feasible,
            "full_feasible": self.full_feasible,
            "non_target_call_flip_count": self.non_target_call_flip_count,
            "proxy_new_essential_count": self.proxy_new_essential_count,
            "no_non_target_call_flips_gate_pass": (
                self.no_non_target_call_flips_gate_pass
            ),
            "concordant_nonessential_proxy_gate_pass": self.proxy_gate_pass,
            "gene_count": self.gene_count,
            "elapsed_seconds": self.elapsed_seconds,
        }
        for gene_id, ratio in self.backup_ko_ratios:
            record[f"backup_ko_ratio__{gene_id}"] = ratio
        return record


@dataclass(frozen=True, slots=True)
class Phase1Decision:
    """The mandatory gate between the baseline experiment and FMPE."""

    status: DecisionStatus
    reasons: tuple[str, ...]
    affine_max_normalized_error: float | None
    multimodal_fraction: float | None
    coverage_by_sigma: tuple[tuple[float, float], ...]

    @property
    def proceed_to_fmpe(self) -> bool:
        return self.status == "proceed_fmpe"

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "proceed_to_fmpe": self.proceed_to_fmpe,
            "reasons": list(self.reasons),
            "affine_max_normalized_error": self.affine_max_normalized_error,
            "multimodal_fraction": self.multimodal_fraction,
            "coverage_by_sigma": {
                str(sigma): coverage for sigma, coverage in self.coverage_by_sigma
            },
        }


@runtime_checkable
class Simulator(Protocol):
    """Minimal interface shared by present and future simulator adapters."""

    @property
    def parameter_bounds(self) -> tuple[tuple[float, float], tuple[float, float]]: ...

    @property
    def target_genes(self) -> tuple[str, ...]: ...

    @property
    def backup_genes(self) -> tuple[str, ...]: ...

    def simulate(
        self,
        point: ParameterPoint,
        *,
        sample_id: str,
        fidelity: Fidelity = "targeted",
    ) -> SimulationResult: ...
