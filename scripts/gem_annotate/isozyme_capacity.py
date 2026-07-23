"""Exploratory KO-specific capacity scans for putative isozymes.

This module deliberately does *not* patch the GEM.  A scan temporarily caps a
reaction only after a named target gene has been knocked out, representing a
hypothesis about the combined residual capacity of the remaining OR partners.
The uncapped wild-type solution remains the denominator, so a shared global
reaction bound cannot create an artificial gene-specific phenotype.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
from cobra.manipulation.delete import knock_out_model_genes

from .essentiality_evidence import target_fingerprint


REQUIRED_COLUMNS = {
    "scenario_id",
    "case_id",
    "target_gene",
    "reaction_id",
    "capacity_fraction_of_wt_flux",
    "model_sha256",
    "media_sha256",
    "target_fingerprint",
    "exploratory_only",
    "evidence_level",
    "basis",
    "rationale",
}
FLUX_EPS = 1e-12


def _parse_bool(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Unrecognised Boolean value: {value!r}")


def load_isozyme_capacity_scan(path: Path) -> pd.DataFrame:
    """Load a model-SHA-pinned, exploratory-only capacity scan definition."""
    scan = pd.read_csv(path, dtype=str).fillna("")
    missing = REQUIRED_COLUMNS - set(scan.columns)
    if missing:
        raise ValueError(
            f"Isozyme capacity scan is missing columns {sorted(missing)}: {path}"
        )
    if scan.empty:
        raise ValueError(f"Isozyme capacity scan is empty: {path}")

    for column in ("scenario_id", "case_id", "target_gene", "reaction_id"):
        scan[column] = scan[column].str.strip()
        if scan[column].eq("").any():
            raise ValueError(f"Isozyme capacity scan has an empty {column}")

    scan["capacity_fraction_of_wt_flux"] = pd.to_numeric(
        scan["capacity_fraction_of_wt_flux"], errors="raise"
    )
    if (~scan["capacity_fraction_of_wt_flux"].map(math.isfinite)).any() or (
        scan["capacity_fraction_of_wt_flux"] < 0
    ).any():
        raise ValueError("Capacity fractions must be finite and non-negative")

    scan["exploratory_only"] = scan["exploratory_only"].map(_parse_bool)
    if not scan["exploratory_only"].all():
        raise ValueError(
            "Every isozyme capacity row must be exploratory_only=true; "
            "measured capacities require the evidence-gated model pipeline"
        )

    duplicated = scan.duplicated(
        subset=["scenario_id", "target_gene", "reaction_id"], keep=False
    )
    if duplicated.any():
        rows = scan.loc[
            duplicated, ["scenario_id", "target_gene", "reaction_id"]
        ].to_dict("records")
        raise ValueError(f"Duplicate isozyme capacity scan targets: {rows[:5]}")

    targets_per_scenario = scan.groupby("scenario_id")["target_gene"].nunique()
    if (targets_per_scenario != 1).any():
        invalid = targets_per_scenario[targets_per_scenario != 1].index.tolist()
        raise ValueError(
            "Each scenario_id must describe exactly one target gene; invalid: "
            f"{invalid}"
        )
    return scan.sort_values(
        ["target_gene", "reaction_id", "capacity_fraction_of_wt_flux", "scenario_id"]
    ).reset_index(drop=True)


def validate_isozyme_capacity_scan(
    model,
    scan: pd.DataFrame,
    model_sha256: str,
    media_sha256: str,
) -> None:
    """Reject stale, non-isozyme or directionally ambiguous scan definitions."""
    declared_shas = {str(value).strip() for value in scan["model_sha256"]}
    if declared_shas != {model_sha256}:
        raise ValueError(
            "Isozyme capacity scan model SHA mismatch: "
            f"declared={sorted(declared_shas)}, current={model_sha256}"
        )
    declared_media_shas = {str(value).strip() for value in scan["media_sha256"]}
    if declared_media_shas != {media_sha256}:
        raise ValueError(
            "Isozyme capacity scan medium SHA mismatch: "
            f"declared={sorted(declared_media_shas)}, current={media_sha256}"
        )

    model_reactions = {reaction.id for reaction in model.reactions}
    model_genes = {gene.id for gene in model.genes}
    for row in scan.itertuples(index=False):
        if row.target_gene not in model_genes:
            raise ValueError(f"Capacity target gene not found: {row.target_gene}")
        if row.reaction_id not in model_reactions:
            raise ValueError(f"Capacity target reaction not found: {row.reaction_id}")
        reaction = model.reactions.get_by_id(row.reaction_id)
        reaction_genes = {gene.id for gene in reaction.genes}
        if row.target_gene not in reaction_genes:
            raise ValueError(
                f"{row.target_gene} is not in the GPR for {row.reaction_id}: "
                f"{reaction.gene_reaction_rule}"
            )
        if " or " not in reaction.gene_reaction_rule.lower():
            raise ValueError(
                f"{row.reaction_id} is not an OR-isozyme reaction: "
                f"{reaction.gene_reaction_rule}"
            )
        if reaction.lower_bound < 0:
            raise ValueError(
                f"{row.reaction_id} is reversible; a one-sided exploratory "
                "capacity requires explicit forward and reverse bounds"
            )
        if reaction.upper_bound <= 0:
            raise ValueError(
                f"{row.reaction_id} is blocked in the current model and must not "
                "be opened by an exploratory capacity scan"
            )

    for scenario_id, group in scan.groupby("scenario_id", sort=True):
        contexts = []
        for reaction_id in sorted(set(group["reaction_id"].astype(str))):
            reaction = model.reactions.get_by_id(reaction_id)
            contexts.append(
                {
                    "reaction_id": reaction.id,
                    "stoichiometry": {
                        metabolite.id: float(coefficient)
                        for metabolite, coefficient in reaction.metabolites.items()
                    },
                    "lower_bound": float(reaction.lower_bound),
                    "upper_bound": float(reaction.upper_bound),
                    "gpr": reaction.gene_reaction_rule,
                }
            )
        current_fingerprint = target_fingerprint(contexts)
        declared_fingerprints = {
            str(value).strip() for value in group["target_fingerprint"]
        }
        if declared_fingerprints != {current_fingerprint}:
            raise ValueError(
                f"Isozyme capacity scan target fingerprint mismatch for {scenario_id}: "
                f"declared={sorted(declared_fingerprints)}, "
                f"current={current_fingerprint}"
            )


def _cutoff_label(cutoff: float) -> str:
    percent = cutoff * 100
    return f"essential_at_{percent:g}pct".replace(".", "p")


def _strictly_below_cutoff(value: float, cutoff: float) -> bool:
    """Apply a strict cutoff without classifying floating-point ties as below."""
    tolerance = max(1e-12, abs(cutoff) * 1e-9)
    return value < cutoff - tolerance


def run_isozyme_capacity_scan(
    model,
    scan: pd.DataFrame,
    baseline_predictions: pd.DataFrame,
    experimental: pd.DataFrame,
    wt_growth: float,
    cutoffs: tuple[float, ...],
    solver: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Evaluate KO-specific residual-capacity hypotheses without mutating *model*."""
    model.solver = solver
    wt_solution = model.optimize()
    if wt_solution.status != "optimal" or wt_solution.objective_value is None:
        raise RuntimeError(f"Wild-type FBA is not optimal: {wt_solution.status}")
    if abs(float(wt_solution.objective_value) - wt_growth) > 1e-8:
        raise RuntimeError("Capacity scan WT does not match the baseline WT solution")

    baseline = baseline_predictions.set_index("gene_id")
    experimental_positive = set(
        experimental.loc[experimental["essential"], "gene_id"].astype(str)
    )
    model_gene_ids = set(baseline.index.astype(str))
    evaluated_positive_ids = experimental_positive & model_gene_ids
    rows: list[dict[str, Any]] = []

    for scenario_id, group in scan.groupby("scenario_id", sort=True):
        target_gene = str(group["target_gene"].iloc[0])
        if target_gene not in baseline.index:
            raise ValueError(f"Capacity target has no baseline prediction: {target_gene}")

        reaction_ids: list[str] = []
        wt_fluxes: list[float] = []
        fractions: list[float] = []
        applied_caps: list[float] = []
        with model:
            knock_out_model_genes(model, [target_gene])
            for item in group.itertuples(index=False):
                reaction = model.reactions.get_by_id(item.reaction_id)
                wt_flux = abs(float(wt_solution.fluxes[item.reaction_id]))
                if wt_flux <= FLUX_EPS:
                    raise ValueError(
                        f"{item.reaction_id} has zero WT flux; a WT-flux-relative "
                        "capacity is undefined"
                    )
                fraction = float(item.capacity_fraction_of_wt_flux)
                cap = wt_flux * fraction
                reaction.upper_bound = min(float(reaction.upper_bound), cap)
                reaction_ids.append(str(item.reaction_id))
                wt_fluxes.append(wt_flux)
                fractions.append(fraction)
                applied_caps.append(cap)

            solution = model.optimize()
            status = str(solution.status)
            if status != "optimal" or solution.objective_value is None:
                scenario_growth = 0.0
            else:
                scenario_growth = max(0.0, float(solution.objective_value))

        baseline_growth = float(baseline.loc[target_gene, "ko_growth"])
        baseline_ratio = float(baseline.loc[target_gene, "ko_growth_ratio"])
        scenario_ratio = scenario_growth / wt_growth
        result: dict[str, Any] = {
            "scenario_id": str(scenario_id),
            "case_id": ";".join(sorted(set(group["case_id"].astype(str)))),
            "target_gene": target_gene,
            "reaction_ids": ";".join(reaction_ids),
            "remaining_isozymes": ";".join(
                sorted(
                    {
                        gene.id
                        for reaction_id in reaction_ids
                        for gene in model.reactions.get_by_id(reaction_id).genes
                        if gene.id != target_gene
                    }
                )
            ),
            "wt_growth": wt_growth,
            "baseline_ko_growth": baseline_growth,
            "baseline_ko_growth_ratio": baseline_ratio,
            "scenario_ko_status": status,
            "scenario_ko_growth": scenario_growth,
            "scenario_ko_growth_ratio": scenario_ratio,
            "wt_reaction_fluxes": ";".join(f"{value:.12g}" for value in wt_fluxes),
            "capacity_fractions_of_wt_flux": ";".join(
                f"{value:.12g}" for value in fractions
            ),
            "capacity_fraction_of_wt_flux": (
                fractions[0] if len(fractions) == 1 else float("nan")
            ),
            "applied_upper_bounds": ";".join(
                f"{value:.12g}" for value in applied_caps
            ),
            "applied_upper_bound": (
                applied_caps[0] if len(applied_caps) == 1 else float("nan")
            ),
            "evidence_level": ";".join(
                sorted(set(group["evidence_level"].astype(str)))
            ),
            "basis": ";".join(sorted(set(group["basis"].astype(str)))),
            "rationale": ";".join(sorted(set(group["rationale"].astype(str)))),
            "experimental_positive": target_gene in experimental_positive,
            "exploratory_only": True,
        }

        for cutoff in cutoffs:
            label = _cutoff_label(cutoff)
            baseline_call = _strictly_below_cutoff(baseline_ratio, cutoff)
            scenario_call = _strictly_below_cutoff(scenario_ratio, cutoff)
            base_tp = sum(
                _strictly_below_cutoff(
                    float(baseline.loc[gene_id, "ko_growth_ratio"]), cutoff
                )
                for gene_id in evaluated_positive_ids
            )
            scenario_tp = base_tp
            if target_gene in evaluated_positive_ids:
                scenario_tp += int(scenario_call) - int(baseline_call)
            result[f"baseline_{label}"] = baseline_call
            result[f"scenario_{label}"] = scenario_call
            result[f"gained_{label}"] = bool(scenario_call and not baseline_call)
            result[f"lost_{label}"] = bool(baseline_call and not scenario_call)
            result[f"scenario_TP_at_{cutoff * 100:g}pct".replace(".", "p")] = scenario_tp
            result[f"scenario_FN_at_{cutoff * 100:g}pct".replace(".", "p")] = (
                len(evaluated_positive_ids) - scenario_tp
            )
            result[f"scenario_recall_at_{cutoff * 100:g}pct".replace(".", "p")] = (
                scenario_tp / len(evaluated_positive_ids)
                if evaluated_positive_ids
                else 0.0
            )
        rows.append(result)

    results = pd.DataFrame(rows).sort_values(
        ["target_gene", "capacity_fraction_of_wt_flux", "scenario_id"]
    ).reset_index(drop=True)
    summary = {
        "analysis_type": "exploratory_ko_specific_isozyme_capacity_sensitivity",
        "exploratory_only": True,
        "not_a_model_patch": True,
        "wild_type_is_uncapped": True,
        "wt_growth": wt_growth,
        "num_scenarios": len(results),
        "targets": sorted(results["target_gene"].unique()),
        "cutoffs": list(cutoffs),
        "warning": (
            "Bounds are sensitivity hypotheses, not measured enzyme capacities. "
            "They must not be promoted to model patches or agent evidence."
        ),
        "scenarios": json.loads(results.to_json(orient="records")),
    }
    return results, summary
