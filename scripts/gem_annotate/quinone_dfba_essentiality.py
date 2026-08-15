"""Runtime-only CoQ9 reserve dFBA screen for the PO1f SD-Leu model.

This is a sensitivity experiment, not a calibrated biomass change.  It adds a
temporary CoQ9 demand proportional to biomass and a finite, one-way CoQ9
reserve.  A knockout can grow from that reserve only until it is exhausted.
Neither the canonical SBML nor the essentiality dossier is changed.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from cobra import Reaction
from cobra.manipulation.delete import knock_out_model_genes

from .config import load_project_paths
from .essentiality_simulation_context import (
    load_effective_simulation_context,
    sha256_file,
    sha256_payload,
)
from .validate_essential_genes import DEFAULT_CUTOFFS, load_experimental


WORKFLOW = "quinone_dfba_essentiality"
SCHEMA_VERSION = "1.1"
Q9_METABOLITE_ID = "m468[C_mi]"
BIOMASS_ID = "biomass_C"
POOL_SOURCE_ID = "DFBA_Q9_POOL_SOURCE"
DEMAND_ID = "DFBA_Q9_GROWTH_DEMAND"
URACIL_MODES = ("finite_batch", "po1f_nonlimiting")

# Values are the concentration statements already recorded in the SD-Leu
# medium comments.  They are finite batch inventories, not uptake kinetics.
INITIAL_POOLS_MMOL_L = {
    "R1070": 111.0,
    "R1003": 0.054,
    "R1202": 0.287,
    "R1204": 0.601,
    "R1215": 0.095,
    "R1217": 0.381,
    "R1220": 0.274,
    "R1222": 0.134,
    "R1223": 0.303,
    "R1231": 0.840,
    "R1232": 0.245,
    "R1233": 0.276,
    "R1234": 1.195,
    "R1354": 0.178428,
}
DEFAULT_ALPHAS = (1e-6, 1e-4, 1e-3)
DEFAULT_POOL_MULTIPLIERS = (0.0, 0.5, 1.0, 2.0)


def _add_runtime_q9(model, alpha: float) -> tuple[Reaction, Reaction]:
    """Add a growth-coupled CoQ9 drain and a finite source to one model copy."""
    if alpha < 0:
        raise ValueError("alpha must be non-negative")
    q9 = model.metabolites.get_by_id(Q9_METABOLITE_ID)
    biomass = model.reactions.get_by_id(BIOMASS_ID)
    source = Reaction(POOL_SOURCE_ID, lower_bound=0.0, upper_bound=0.0)
    source.add_metabolites({q9: 1.0})
    demand = Reaction(DEMAND_ID, lower_bound=0.0, upper_bound=1000.0)
    demand.add_metabolites({q9: -1.0})
    model.add_reactions([source, demand])
    constraint = model.problem.Constraint(
        demand.flux_expression - alpha * biomass.flux_expression,
        lb=0.0,
        ub=0.0,
        name="dfba_q9_growth_coupling",
    )
    model.add_cons_vars(constraint)
    return source, demand


def _finite_medium(base_medium: dict[str, float], pools: dict[str, float], biomass: float, dt: float) -> dict[str, float]:
    medium = dict(base_medium)
    for reaction_id, amount in pools.items():
        if reaction_id not in medium:
            continue
        medium[reaction_id] = min(float(medium[reaction_id]), max(0.0, amount / (biomass * dt)))
    return medium


def _optimize_minimal_pool(model, source: Reaction) -> tuple[float, Any]:
    """Maximize growth first, then use the reserve only if synthesis cannot serve it."""
    model.objective = model.reactions.get_by_id(BIOMASS_ID)
    source_limit = source.upper_bound
    source.upper_bound = 0.0
    primary = model.optimize()
    source.upper_bound = source_limit
    if primary.status != "optimal" or primary.objective_value is None:
        return 0.0, primary
    growth = max(0.0, float(primary.objective_value))
    if growth > 1e-9:
        return growth, primary
    primary = model.optimize()
    if primary.status != "optimal" or primary.objective_value is None:
        return 0.0, primary
    growth = max(0.0, float(primary.objective_value))
    biomass = model.reactions.get_by_id(BIOMASS_ID)
    with model:
        growth_lock = model.problem.Constraint(
            biomass.flux_expression, lb=max(0.0, growth - 1e-9), name="dfba_growth_lock"
        )
        model.add_cons_vars(growth_lock)
        model.objective = source
        model.objective_direction = "min"
        secondary = model.optimize()
        if secondary.status != "optimal":
            raise RuntimeError(f"Reserve minimization failed: {secondary.status}")
        return growth, secondary


def _software_versions(solver: str) -> dict[str, str]:
    versions = {"python": sys.version.split()[0], "cobra": __import__("cobra").__version__}
    if solver.lower() == "gurobi":
        import gurobipy

        versions["gurobipy"] = ".".join(map(str, gurobipy.gurobi.version()))
    return versions


def simulate_gene(
    base_model,
    *,
    gene_id: str | None,
    alpha: float,
    pool_multiplier: float,
    hours: float,
    dt: float,
    initial_biomass: float,
    uracil_mode: str = "finite_batch",
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Run a disposable Euler batch simulation for WT or one single-gene KO."""
    if hours <= 0 or dt <= 0 or initial_biomass <= 0:
        raise ValueError("hours, dt, and initial_biomass must be positive")
    if uracil_mode not in URACIL_MODES:
        raise ValueError(f"uracil_mode must be one of {URACIL_MODES}")
    if not math.isclose(hours / dt, round(hours / dt), rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("hours must be an exact multiple of dt")
    with base_model as model:
        if gene_id is not None:
            knock_out_model_genes(model, [gene_id])
        source, _ = _add_runtime_q9(model, alpha)
        base_medium = dict(model.medium)
        pools = dict(INITIAL_POOLS_MMOL_L)
        if uracil_mode == "po1f_nonlimiting":
            pools.pop("R1354")
        initial_q9_pool = alpha * initial_biomass * pool_multiplier
        q9_pool = initial_q9_pool
        q9_tolerance = max(1e-12, initial_q9_pool * 1e-9)
        biomass = initial_biomass
        depleted_at: float | None = 0.0 if q9_pool == 0 else None
        trajectory: list[dict[str, Any]] = []
        for step in range(int(round(hours / dt))):
            time = step * dt
            model.medium = _finite_medium(base_medium, pools, biomass, dt)
            source.upper_bound = q9_pool / (biomass * dt) if q9_pool > 0 else 0.0
            growth, solution = _optimize_minimal_pool(model, source)
            source_flux = max(0.0, float(solution.fluxes[source.id]))
            consumed_q9 = min(q9_pool, source_flux * biomass * dt)
            q9_pool = max(0.0, q9_pool - consumed_q9)
            if q9_pool <= q9_tolerance:
                q9_pool = 0.0
            if q9_pool == 0.0 and depleted_at is None:
                depleted_at = time + dt
            for reaction_id in pools:
                uptake = max(0.0, -float(solution.fluxes[reaction_id]))
                pools[reaction_id] = max(0.0, pools[reaction_id] - uptake * biomass * dt)
            trajectory.append(
                {
                    "gene_id": gene_id or "WT",
                    "uracil_mode": uracil_mode,
                    "time_h": time,
                    "biomass_gDW_L": biomass,
                    "growth_h-1": growth,
                    "q9_pool_mmol_L": q9_pool,
                    "q9_source_flux_mmol_gDW_h": source_flux,
                    "uracil_mmol_L": pools.get("R1354", math.nan),
                    "glucose_mmol_L": pools["R1070"],
                    "status": solution.status,
                }
            )
            biomass *= 1.0 + growth * dt
        return (
            {
                "gene_id": gene_id or "WT",
                "uracil_mode": uracil_mode,
                "final_biomass_gDW_L": biomass,
                "dynamic_doublings": math.log2(biomass / initial_biomass),
                "initial_growth_h-1": trajectory[0]["growth_h-1"],
                "q9_pool_depleted_h": depleted_at,
                "q9_source_total_mmol_L": initial_q9_pool - q9_pool,
                "final_glucose_mmol_L": pools["R1070"],
                "final_uracil_mmol_L": pools.get("R1354", math.nan),
            },
            trajectory,
        )


def _gene_ids(model, requested: str | None, excluded: tuple[str, ...], chunk_index: int, chunk_count: int) -> list[str]:
    if chunk_count <= 0 or not 0 <= chunk_index < chunk_count:
        raise ValueError("chunk-index must be in [0, chunk-count)")
    available = sorted(gene.id for gene in model.genes if gene.id not in set(excluded))
    selected = [item.strip() for item in requested.split(",")] if requested else available
    missing = sorted(set(selected) - set(available))
    if missing:
        raise ValueError(f"Requested genes are not screenable: {missing}")
    return selected[chunk_index::chunk_count]


def _run_chunk(args: argparse.Namespace) -> Path:
    paths = load_project_paths(args.research_root, required=True)
    media_path = paths.media / "sd_leu.csv"
    profile_path = paths.strain_profiles / "po1f_sd_leu.json"
    experimental_path = paths.essentiality / "consensus_essential_genes.csv"
    paths.require(media_path, profile_path, experimental_path, paths.output_model)
    context = load_effective_simulation_context(
        model_path=paths.output_model, media_path=media_path, strain_profile_path=profile_path
    )
    context.model.solver = args.solver
    base_r1354_bound = float(context.model.medium["R1354"])
    genes = _gene_ids(context.model, args.genes, context.excluded_runtime_genes, args.chunk_index, args.chunk_count)
    out = Path(args.output_dir or paths.results / WORKFLOW / args.run_id).resolve()
    out.mkdir(parents=True, exist_ok=True)
    experimental = load_experimental(experimental_path, positive_only=True)
    rows: list[dict[str, Any]] = []
    trajectories: list[dict[str, Any]] = []
    for alpha in args.alphas:
        for multiplier in args.pool_multipliers:
            wt, wt_trace = simulate_gene(context.model, gene_id=None, alpha=alpha, pool_multiplier=multiplier, hours=args.hours, dt=args.dt, initial_biomass=args.initial_biomass, uracil_mode=args.uracil_mode)
            trajectories.extend([{**item, "alpha_mmol_gDW": alpha, "pool_multiplier": multiplier} for item in wt_trace])
            for gene_id in genes:
                result, trace = simulate_gene(context.model, gene_id=gene_id, alpha=alpha, pool_multiplier=multiplier, hours=args.hours, dt=args.dt, initial_biomass=args.initial_biomass, uracil_mode=args.uracil_mode)
                ratio = result["dynamic_doublings"] / wt["dynamic_doublings"] if wt["dynamic_doublings"] > 0 else math.nan
                rows.append({
                    **result, "alpha_mmol_gDW": alpha, "pool_multiplier": multiplier,
                    "wt_dynamic_doublings": wt["dynamic_doublings"], "dynamic_growth_ratio": ratio,
                    "experimental_essential": gene_id in set(experimental["gene_id"]),
                    **{f"essential_at_{cutoff * 100:g}pct": bool(ratio < cutoff) for cutoff in DEFAULT_CUTOFFS},
                })
                trajectories.extend([{**item, "alpha_mmol_gDW": alpha, "pool_multiplier": multiplier} for item in trace])
    pd.DataFrame(rows).to_csv(out / f"chunk_{args.chunk_index:03d}_calls.tsv", sep="\t", index=False)
    pd.DataFrame(trajectories).to_csv(out / f"chunk_{args.chunk_index:03d}_trajectory.tsv", sep="\t", index=False)
    manifest = {
        "workflow": WORKFLOW, "schema_version": SCHEMA_VERSION, "run_id": args.run_id,
        "chunk_index": args.chunk_index, "chunk_count": args.chunk_count, "genes": genes,
        "solver": args.solver, "runtime_versions": _software_versions(args.solver),
        "script_sha256": sha256_file(Path(__file__)), "hours": args.hours, "dt_h": args.dt,
        "initial_biomass_gDW_L": args.initial_biomass, "alphas_mmol_gDW": args.alphas,
        "pool_multipliers": args.pool_multipliers,
        "uracil_mode": args.uracil_mode,
        "initial_pools_mmol_L": {**INITIAL_POOLS_MMOL_L, "R1354": None} if args.uracil_mode == "po1f_nonlimiting" else INITIAL_POOLS_MMOL_L,
        "base_r1354_bound_mmol_gDW_h": base_r1354_bound,
        "q9_reserve_definition": "alpha * initial_biomass * pool_multiplier mmol/L",
        "calibration_status": "sensitivity_only_not_calibrated",
        "simulation_context": context.provenance(),
        "input_sha256": {"model": context.canonical_model_sha256, "medium": sha256_file(media_path), "profile": sha256_file(profile_path), "experimental": sha256_file(experimental_path)},
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest["fingerprint"] = sha256_payload(manifest)
    (out / f"chunk_{args.chunk_index:03d}_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def _merge(args: argparse.Namespace) -> Path:
    out = Path(args.output_dir).resolve()
    manifests = sorted(out.glob("chunk_*_manifest.json"))
    if not manifests:
        raise FileNotFoundError(f"No chunk manifests in {out}")
    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
    fingerprints = {payload["fingerprint"] for payload in payloads}
    comparable = [{key: value for key, value in payload.items() if key not in {"chunk_index", "genes", "created_at_utc", "fingerprint"}} for payload in payloads]
    if len({sha256_payload(item) for item in comparable}) != 1:
        raise ValueError("Chunk manifests do not share one simulation context")
    expected = set(range(payloads[0]["chunk_count"]))
    observed = {payload["chunk_index"] for payload in payloads}
    if observed != expected:
        raise ValueError(f"Incomplete chunks: missing {sorted(expected - observed)}")
    calls = pd.concat([pd.read_csv(path, sep="\t") for path in sorted(out.glob("chunk_*_calls.tsv"))], ignore_index=True)
    calls.to_csv(out / "essentiality_dynamic_calls.tsv", sep="\t", index=False)
    summary = {"workflow": WORKFLOW, "chunks": len(manifests), "chunk_fingerprints": sorted(fingerprints), "calls": len(calls), "calibration_status": "sensitivity_only_not_calibrated"}
    (out / "merge_manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-root", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--run-id", default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    parser.add_argument("--solver", default="glpk")
    parser.add_argument("--genes", help="comma-separated pilot genes; defaults to every model gene")
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--chunk-count", type=int, default=1)
    parser.add_argument("--hours", type=float, default=24.0)
    parser.add_argument("--dt", type=float, default=0.25)
    parser.add_argument("--initial-biomass", type=float, default=0.01)
    parser.add_argument("--alphas", type=float, nargs="+", default=list(DEFAULT_ALPHAS))
    parser.add_argument("--pool-multipliers", type=float, nargs="+", default=list(DEFAULT_POOL_MULTIPLIERS))
    parser.add_argument("--uracil-mode", choices=URACIL_MODES, default="finite_batch")
    parser.add_argument("--merge", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.merge:
        if not args.output_dir:
            raise ValueError("--merge requires --output-dir")
        print(_merge(args))
    else:
        print(_run_chunk(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
