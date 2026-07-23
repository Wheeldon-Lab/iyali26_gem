"""iYali26 adapter for the exploratory R4/R1846 capacity experiment.

The adapter deliberately reuses the validated SD-Leu and provisional-capacity
implementation already present in :mod:`scripts.gem_annotate`.  Every
parameter point is evaluated on an in-memory model copy; neither the SBML nor
the curated capacity table is written by this module.
"""

from __future__ import annotations

import platform
import time
from typing import Any

import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from scripts.gem_annotate.calibrate_isozyme_capacities import (
    BACKUP_KO_RATIO_MINIMUM,
    BLOCKED_CONTROL_REACTION,
    RATIO_CHANGE_EPS,
    TARGET_LOWER_BOUND,
    TARGET_UPPER_BOUND,
    WT_DELTA_TOLERANCE,
    _build_scenario_per_gene,
    _configure_solver,
    _load_assay_fitness,
    _load_model,
    _optimal_growth_and_fluxes,
    _profile_gene_roles,
    _scenario_gates,
    concordant_nonessential_proxies,
)
from scripts.gem_annotate.provisional_capacity import (
    apply_provisional_isozyme_capacities,
    load_provisional_capacity_table,
)
from scripts.gem_annotate.validate_essential_genes import (
    apply_media,
    load_experimental,
    load_media,
    run_single_gene_deletions,
)

from .config import ExperimentConfig
from .core import Fidelity, ParameterPoint, SimulationResult
from .provenance import sha256_file


TARGET_REACTIONS = ("R4", "R1846")
FAMILY_TOTAL_TOLERANCE = 1e-12


def _deletion_rows(model, gene_ids: tuple[str, ...]) -> tuple[pd.DataFrame, float]:
    """Run a deterministic deletion subset and normalize COBRA's row format."""
    solution = model.optimize()
    if solution.status != "optimal" or solution.objective_value is None:
        raise RuntimeError(f"Wild-type FBA is not optimal: {solution.status}")
    wt_growth = float(solution.objective_value)
    deletion = single_gene_deletion(
        model,
        gene_list=list(gene_ids),
        processes=1,
    )
    rows: list[dict[str, object]] = []
    for index, row in deletion.iterrows():
        ids = row.get("ids")
        if isinstance(ids, (set, frozenset)) and ids:
            gene_id = str(next(iter(ids)))
        elif isinstance(index, (set, frozenset)) and index:
            gene_id = str(next(iter(index)))
        else:
            gene_id = str(ids if ids is not None else index)
        status = str(row.get("status", "optimal"))
        value = row.get("growth")
        growth = (
            max(0.0, float(value))
            if status == "optimal" and value is not None and not pd.isna(value)
            else 0.0
        )
        rows.append(
            {
                "gene_id": gene_id,
                "ko_status": status,
                "ko_growth": growth,
                "ko_growth_ratio": growth / wt_growth,
            }
        )
    frame = pd.DataFrame(rows).sort_values("gene_id").reset_index(drop=True)
    missing = sorted(set(gene_ids) - set(frame["gene_id"].astype(str)))
    if missing:
        raise RuntimeError(f"Targeted deletion output omitted genes: {missing}")
    return frame, wt_growth


class R4R1846CapacitySimulator:
    """Deterministic two-parameter simulator backed by the current iYali26 GEM."""

    def __init__(self, config: ExperimentConfig) -> None:
        self.config = config
        self.verified_inputs = config.verify_inputs()
        inputs = config.input_map
        self.model_path = inputs["model"].path
        self.model_sha256 = sha256_file(self.model_path)

        profile_rows = load_provisional_capacity_table(inputs["profile"].path)
        by_reaction = {
            str(row["source_reaction_id"]): row for row in profile_rows
        }
        if set(by_reaction) != set(TARGET_REACTIONS) or len(profile_rows) != 2:
            raise ValueError(
                "R4/R1846 simulator requires exactly the two active capacity rows"
            )
        for reaction_id, row in by_reaction.items():
            if row["validated_model_sha256"] != self.model_sha256:
                raise ValueError(
                    f"Profile model SHA is stale for {reaction_id}: "
                    f"{row['validated_model_sha256']} != {self.model_sha256}"
                )
        target_by_reaction, backup_genes = _profile_gene_roles(profile_rows)
        self._target_by_reaction = dict(target_by_reaction)
        self._backup_genes = tuple(sorted(backup_genes))
        self._target_genes = tuple(
            self._target_by_reaction[reaction_id]
            for reaction_id in TARGET_REACTIONS
        )
        if len(self._target_genes) != 2 or len(self._backup_genes) != 3:
            raise ValueError(
                "Expected two target genes and three backup genes for R4/R1846"
            )
        self.target_fingerprints = {
            reaction_id: str(row["target_fingerprint"])
            for reaction_id, row in by_reaction.items()
        }

        baseline_model = _load_model(self.model_path, config.solver)
        apply_media(baseline_model, load_media(inputs["media"].path))
        if BLOCKED_CONTROL_REACTION not in baseline_model.reactions:
            raise ValueError(f"Blocked control reaction missing: {BLOCKED_CONTROL_REACTION}")
        r1843_bounds = tuple(
            float(value)
            for value in baseline_model.reactions.get_by_id(
                BLOCKED_CONTROL_REACTION
            ).bounds
        )
        if r1843_bounds != (0.0, 0.0):
            raise ValueError(
                f"{BLOCKED_CONTROL_REACTION} must be locked at (0, 0), "
                f"found {r1843_bounds}"
            )
        self._baseline_model = baseline_model
        self._baseline_wt_growth, self._wt_fluxes = _optimal_growth_and_fluxes(
            baseline_model
        )
        self._baseline_family_upper_bounds = {
            reaction_id: float(
                baseline_model.reactions.get_by_id(reaction_id).upper_bound
            )
            for reaction_id in TARGET_REACTIONS
        }
        for reaction_id in TARGET_REACTIONS:
            if abs(float(self._wt_fluxes[reaction_id])) <= RATIO_CHANGE_EPS:
                raise ValueError(
                    f"{reaction_id} has zero WT flux; relative capacity is undefined"
                )

        # Apply the existing, SHA-validated profile once.  Per-point copies only
        # change the source/backup upper-bound split while preserving its sum.
        template = baseline_model.copy()
        _configure_solver(template, config.solver)
        apply_provisional_isozyme_capacities(
            template,
            inputs["profile"].path,
            reference_model_sha256=self.model_sha256,
        )
        self._scenario_template = template
        self._profile_rows = by_reaction

        assay = _load_assay_fitness(inputs["assay_fitness"].path)
        self._proxy_genes = concordant_nonessential_proxies(assay)
        experimental = load_experimental(
            inputs["experimental"].path,
            positive_only=True,
        )
        self._experimental_positive_genes = set(
            experimental.loc[experimental["essential"], "gene_id"].astype(str)
        )
        self._baseline_predictions: pd.DataFrame | None = None

    @property
    def parameter_bounds(self) -> tuple[tuple[float, float], tuple[float, float]]:
        return (
            (self.config.r4_bounds.lower, self.config.r4_bounds.upper),
            (self.config.r1846_bounds.lower, self.config.r1846_bounds.upper),
        )

    @property
    def target_genes(self) -> tuple[str, ...]:
        return self._target_genes

    @property
    def backup_genes(self) -> tuple[str, ...]:
        return self._backup_genes

    @property
    def fast_gene_ids(self) -> tuple[str, ...]:
        return tuple(sorted((*self._target_genes, *self._backup_genes)))

    @property
    def baseline_wt_growth(self) -> float:
        return self._baseline_wt_growth

    @property
    def baseline_invariants(self) -> dict[str, object]:
        return {
            "model_sha256": self.model_sha256,
            "wild_type_growth": self._baseline_wt_growth,
            "family_total_upper_bounds": dict(self._baseline_family_upper_bounds),
            "r1843_bounds": list(
                self._baseline_model.reactions.get_by_id(
                    BLOCKED_CONTROL_REACTION
                ).bounds
            ),
            "target_genes": list(self.target_genes),
            "backup_genes": list(self.backup_genes),
            "target_fingerprints": dict(self.target_fingerprints),
        }

    @property
    def solver_metadata(self) -> dict[str, object]:
        metadata: dict[str, object] = {
            "requested": self.config.solver,
            "interface": self._baseline_model.solver.interface.__name__,
            "python": platform.python_version(),
        }
        try:
            import gurobipy

            metadata["gurobi"] = ".".join(
                str(value) for value in gurobipy.gurobi.version()
            )
        except (ImportError, AttributeError):
            pass
        return metadata

    def _validate_point(self, point: ParameterPoint) -> None:
        for name, value, bounds in (
            ("r4_fraction", point.r4_fraction, self.parameter_bounds[0]),
            ("r1846_fraction", point.r1846_fraction, self.parameter_bounds[1]),
        ):
            if not bounds[0] <= value <= bounds[1]:
                raise ValueError(
                    f"{name}={value} is outside configured bounds {bounds}"
                )

    def _scenario_model(
        self, point: ParameterPoint
    ) -> tuple[Any, dict[str, float], dict[str, tuple[float, float]]]:
        self._validate_point(point)
        # The runner is deliberately single-threaded.  Reusing this private
        # working model avoids a multi-second solver/model copy per lattice
        # point.  Both bounds are overwritten for every call, and COBRA's
        # deletion routines restore temporary gene knockouts before returning.
        model = self._scenario_template
        caps: dict[str, float] = {}
        family_totals: dict[str, tuple[float, float]] = {}
        fractions = {"R4": point.r4_fraction, "R1846": point.r1846_fraction}
        for reaction_id in TARGET_REACTIONS:
            row = self._profile_rows[reaction_id]
            backup_id = str(row["backup_reaction_id"])
            total = self._baseline_family_upper_bounds[reaction_id]
            cap = abs(float(self._wt_fluxes[reaction_id])) * fractions[reaction_id]
            if cap > total + FAMILY_TOTAL_TOLERANCE:
                raise ValueError(
                    f"{reaction_id} backup capacity {cap} exceeds family total {total}"
                )
            source = model.reactions.get_by_id(reaction_id)
            backup = model.reactions.get_by_id(backup_id)
            source.upper_bound = total - cap
            backup.upper_bound = cap
            caps[reaction_id] = cap
            family_totals[reaction_id] = (
                total,
                float(source.upper_bound) + float(backup.upper_bound),
            )
        return model, caps, family_totals

    def _ensure_baseline_predictions(self) -> pd.DataFrame:
        if self._baseline_predictions is None:
            predictions, deletion_wt = run_single_gene_deletions(
                self._baseline_model,
                self.config.solver,
            )
            if abs(deletion_wt - self._baseline_wt_growth) > WT_DELTA_TOLERANCE:
                raise RuntimeError(
                    "Baseline WT changed between optimization and deletion run"
                )
            self._baseline_predictions = predictions
        return self._baseline_predictions

    def simulate(
        self,
        point: ParameterPoint,
        *,
        sample_id: str,
        fidelity: Fidelity = "targeted",
    ) -> SimulationResult:
        if fidelity not in {"targeted", "full"}:
            raise ValueError(f"Unsupported fidelity: {fidelity}")
        started = time.monotonic()
        model, caps, family_totals = self._scenario_model(point)
        scenario_wt_growth, _ = _optimal_growth_and_fluxes(model)
        wt_delta = abs(scenario_wt_growth - self._baseline_wt_growth)
        r1843_pass = tuple(
            float(value)
            for value in model.reactions.get_by_id(BLOCKED_CONTROL_REACTION).bounds
        ) == (0.0, 0.0)
        family_total_pass = all(
            abs(current - baseline) <= FAMILY_TOTAL_TOLERANCE
            for baseline, current in family_totals.values()
        )

        full_feasible: bool | None = None
        non_target_flip_count: int | None = None
        proxy_new_essential_count: int | None = None
        no_non_target_flips: bool | None = None
        proxy_pass: bool | None = None

        if fidelity == "targeted":
            predictions, deletion_wt = _deletion_rows(model, self.fast_gene_ids)
            indexed = predictions.set_index("gene_id")
            target_ratios = {
                reaction_id: float(
                    indexed.loc[gene_id, "ko_growth_ratio"]
                )
                for reaction_id, gene_id in self._target_by_reaction.items()
            }
            backup_ratios = {
                gene_id: float(indexed.loc[gene_id, "ko_growth_ratio"])
                for gene_id in self.backup_genes
            }
            target_window_pass = all(
                TARGET_LOWER_BOUND <= ratio < TARGET_UPPER_BOUND
                for ratio in target_ratios.values()
            )
            backup_pass = bool(backup_ratios) and all(
                ratio > BACKUP_KO_RATIO_MINIMUM
                for ratio in backup_ratios.values()
            )
            gene_count = len(predictions)
        else:
            predictions, deletion_wt = run_single_gene_deletions(
                model,
                self.config.solver,
            )
            per_gene = _build_scenario_per_gene(
                self._ensure_baseline_predictions(),
                predictions,
                target_genes=set(self.target_genes),
                backup_genes=set(self.backup_genes),
                proxy_genes=self._proxy_genes,
                experimental_positive_genes=self._experimental_positive_genes,
                cutoffs=self.config.growth_cutoffs,
            )
            gates = _scenario_gates(
                per_gene,
                target_by_reaction=self._target_by_reaction,
                backup_genes=set(self.backup_genes),
                proxy_genes=self._proxy_genes,
                baseline_wt_growth=self._baseline_wt_growth,
                scenario_wt_growth=scenario_wt_growth,
                r1843_bounds=tuple(
                    model.reactions.get_by_id(BLOCKED_CONTROL_REACTION).bounds
                ),
                family_totals=family_totals,
                cutoffs=self.config.growth_cutoffs,
            )
            target_ratios = {
                key: float(value)
                for key, value in gates["target_ko_ratios"].items()
            }
            backup_ratios = {
                key: float(value)
                for key, value in gates["backup_gene_ko_ratios"].items()
            }
            target_window_pass = bool(gates["target_window_pass"])
            backup_pass = bool(gates["backup_gene_ko_ratio_gate_pass"])
            r1843_pass = bool(gates["r1843_bounds_gate_pass"])
            family_total_pass = bool(
                gates["reaction_family_total_upper_bound_gate_pass"]
            )
            full_feasible = bool(gates["feasible"])
            non_target_flip_count = int(gates["non_target_call_flip_count"])
            proxy_new_essential_count = int(gates["proxy_new_essential_count"])
            no_non_target_flips = bool(
                gates["no_non_target_call_flips_gate_pass"]
            )
            proxy_pass = bool(
                gates["concordant_nonessential_proxy_gate_pass"]
            )
            gene_count = len(predictions)

        if abs(deletion_wt - scenario_wt_growth) > WT_DELTA_TOLERANCE:
            raise RuntimeError("Scenario WT changed between optimization and deletion run")
        wt_pass = wt_delta <= WT_DELTA_TOLERANCE
        targeted_feasible = all(
            (
                target_window_pass,
                backup_pass,
                r1843_pass,
                family_total_pass,
                wt_pass,
            )
        )
        return SimulationResult(
            sample_id=sample_id,
            fidelity=fidelity,
            point=point,
            r4_applied_upper_bound=caps["R4"],
            r1846_applied_upper_bound=caps["R1846"],
            r4_target_ko_ratio=target_ratios["R4"],
            r1846_target_ko_ratio=target_ratios["R1846"],
            backup_ko_ratios=tuple(sorted(backup_ratios.items())),
            baseline_wt_growth=self._baseline_wt_growth,
            scenario_wt_growth=scenario_wt_growth,
            wt_growth_delta=0.0 if wt_delta <= RATIO_CHANGE_EPS else wt_delta,
            target_window_pass=target_window_pass,
            backup_gate_pass=backup_pass,
            r1843_gate_pass=r1843_pass,
            family_total_gate_pass=family_total_pass,
            targeted_feasible=targeted_feasible,
            full_feasible=full_feasible,
            non_target_call_flip_count=non_target_flip_count,
            proxy_new_essential_count=proxy_new_essential_count,
            no_non_target_call_flips_gate_pass=no_non_target_flips,
            proxy_gate_pass=proxy_pass,
            gene_count=gene_count,
            elapsed_seconds=time.monotonic() - started,
        )

    def assert_input_model_unchanged(self) -> None:
        actual = sha256_file(self.model_path)
        if actual != self.model_sha256:
            raise RuntimeError(
                f"Input model SHA changed during experiment: {self.model_sha256} -> {actual}"
            )
