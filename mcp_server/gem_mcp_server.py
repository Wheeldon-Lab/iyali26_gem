"""GEM MCP server — diagnose tools for genome-scale metabolic models via COBRApy.

This is the *diagnose* half of an eventual diagnose -> repair -> validate loop for
autonomous model curation. It exposes read-only COBRApy diagnostics as MCP tools so
an agent (e.g. Claude Desktop) can inspect a model, find mass/charge imbalances,
drill into individual reactions and metabolites, and run FBA.

The loaded model is cached in module-level server state so tools operate on the live
cobra.Model object rather than re-parsing the SBML on every call.

Tools:
    load_model(path)                 -- load an SBML/JSON model and cache it; return a summary
    list_mass_charge_imbalances(...) -- every reaction whose mass or charge does not balance
    get_reaction(reaction_id)        -- full detail on one reaction
    get_metabolite(metabolite_id)    -- formula/charge/compartment + participating reactions
    run_fba()                        -- flux balance analysis: objective value + nonzero fluxes
"""

from __future__ import annotations

import os
from typing import Any, Optional

from fastmcp import FastMCP

# COBRApy is imported lazily inside load_model so the server can still start (and
# report a clean error) if the cobra install is broken, rather than failing at import.

mcp = FastMCP("gem-curation")

# ---------------------------------------------------------------------------
# Server state
# ---------------------------------------------------------------------------

# Directory of the working model, so a bare filename like "model.xml" resolves
# against the project root regardless of the server's working directory.
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MODEL = "model.xml"

# Flux magnitudes below this are treated as zero when reporting FBA / imbalances.
ZERO_TOL = 1e-9

_STATE: dict[str, Any] = {
    "model": None,   # the live cobra.Model, or None if nothing loaded yet
    "path": None,    # absolute path it was loaded from
}


def _require_model() -> Optional[dict]:
    """Return an error dict if no model is loaded, else None.

    Every tool calls this first so the "no model loaded yet" case comes back as
    clean structured data instead of an exception.
    """
    if _STATE["model"] is None:
        return {
            "error": "no_model_loaded",
            "message": "No model is loaded yet. Call load_model(path) first.",
            "hint": f"Default working model is '{DEFAULT_MODEL}' in {PROJECT_DIR}.",
        }
    return None


def _resolve_path(path: str) -> str:
    """Resolve a possibly-relative model path against the project directory."""
    if os.path.isabs(path):
        return path
    return os.path.join(PROJECT_DIR, path)


def _mass_balance(reaction) -> dict[str, float]:
    """check_mass_balance() with exceptions turned into structured data.

    Returns {} when the reaction balances, otherwise a dict mixing element keys
    ('C', 'H', 'O', ...) and 'charge'. On error returns {'_error': msg}.
    """
    try:
        return reaction.check_mass_balance()
    except Exception as exc:  # e.g. a metabolite with no formula
        return {"_error": str(exc)}


def _split_imbalance(balance: dict[str, float]) -> dict[str, Any]:
    """Split a check_mass_balance() dict into charge vs. per-element mass parts."""
    charge = balance.get("charge")
    mass = {k: v for k, v in balance.items() if k not in ("charge", "_error")}
    out: dict[str, Any] = {
        "charge_imbalance": charge,  # None if charge balances
        "mass_imbalance": mass,      # {} if all elements balance
    }
    if "_error" in balance:
        out["balance_error"] = balance["_error"]
    return out


def _annotation_coverage(model) -> dict[str, int]:
    """Count how many metabolites carry a formula / charge.

    Mass balance can only be computed when every metabolite in a reaction has a
    formula; charge balance needs every charge. This exposes whether the model is
    even *diagnosable* so callers don't mistake "can't check" for "balanced".
    """
    with_formula = sum(1 for m in model.metabolites if m.formula)
    with_charge = sum(1 for m in model.metabolites if m.charge is not None)
    total = len(model.metabolites)
    return {
        "num_metabolites": total,
        "metabolites_with_formula": with_formula,
        "metabolites_without_formula": total - with_formula,
        "metabolites_with_charge": with_charge,
        "metabolites_without_charge": total - with_charge,
    }


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------


@mcp.tool()
def load_model(path: str = DEFAULT_MODEL) -> dict:
    """Load a metabolic model, cache it in server state, and return a summary.

    Args:
        path: Path to the model file. Relative paths resolve against the project
            root. SBML (.xml) and JSON (.json) are supported. Defaults to the
            canonical working model 'model.xml'.

    Returns:
        A summary dict (id, path, counts, objective, compartments) or an error dict.
    """
    abs_path = _resolve_path(path)
    if not os.path.exists(abs_path):
        return {
            "error": "file_not_found",
            "message": f"No such file: {abs_path}",
            "project_dir": PROJECT_DIR,
        }

    try:
        # Import here so a broken cobra install surfaces as a clean tool error.
        from cobra.io import load_json_model, read_sbml_model

        if abs_path.endswith(".json"):
            model = load_json_model(abs_path)
        else:
            model = read_sbml_model(abs_path)
    except Exception as exc:
        return {"error": "load_failed", "message": str(exc), "path": abs_path}

    _STATE["model"] = model
    _STATE["path"] = abs_path

    try:
        objective = str(model.objective.expression)
    except Exception:
        objective = None

    return {
        "id": model.id,
        "path": abs_path,
        "num_reactions": len(model.reactions),
        "num_metabolites": len(model.metabolites),
        "num_genes": len(model.genes),
        "objective": objective,
        "compartments": dict(model.compartments),
        "annotation_coverage": _annotation_coverage(model),
    }


@mcp.tool()
def list_mass_charge_imbalances(
    include_boundary: bool = False,
    limit: Optional[int] = None,
) -> dict:
    """List every reaction whose mass or charge does not balance.

    This is the core diagnostic. Boundary/exchange reactions are excluded by
    default because they are intentionally open (they add or remove mass at the
    system boundary and are expected to be "imbalanced").

    Args:
        include_boundary: If True, also check boundary (exchange/sink/demand)
            reactions. Usually leave False.
        limit: If set, return at most this many imbalanced reactions (they are
            still fully counted in the summary). Useful to cap payload size.

    Returns:
        {count, num_charge_imbalanced, num_mass_imbalanced, num_boundary_skipped,
         truncated, diagnosability, reactions: [...]}, where each reaction has its
         charge and per-element discrepancies. `diagnosability` reports how many
         reactions can actually be checked -- a reaction with a formula-less
         metabolite silently reports "balanced" even though nothing was verified,
         so callers must not read count==0 as "clean" without checking it.
    """
    err = _require_model()
    if err:
        return err
    model = _STATE["model"]

    imbalanced: list[dict] = []
    num_charge = 0
    num_mass = 0
    boundary_skipped = 0
    num_checked = 0
    num_mass_checkable = 0
    num_charge_checkable = 0

    for rxn in model.reactions:
        if rxn.boundary and not include_boundary:
            boundary_skipped += 1
            continue

        # Diagnosability: mass balance needs every metabolite to have a formula;
        # charge balance needs every charge. Track these so a model missing
        # annotations can't masquerade as "balanced".
        num_checked += 1
        mets = list(rxn.metabolites)
        if all(m.formula for m in mets):
            num_mass_checkable += 1
        if all(m.charge is not None for m in mets):
            num_charge_checkable += 1

        balance = _mass_balance(rxn)
        if not balance:
            continue  # balances (or an empty dict) -> skip

        parts = _split_imbalance(balance)
        has_charge = parts["charge_imbalance"] is not None
        has_mass = bool(parts["mass_imbalance"])
        if not (has_charge or has_mass or "balance_error" in parts):
            continue

        if has_charge:
            num_charge += 1
        if has_mass:
            num_mass += 1

        imbalanced.append(
            {
                "reaction_id": rxn.id,
                "name": rxn.name,
                "reaction_string": rxn.reaction,
                "is_boundary": rxn.boundary,
                **parts,
            }
        )

    count = len(imbalanced)
    truncated = False
    if limit is not None and count > limit:
        imbalanced = imbalanced[:limit]
        truncated = True

    num_mass_uncheckable = num_checked - num_mass_checkable
    num_charge_uncheckable = num_checked - num_charge_checkable
    warning = None
    if num_checked and num_mass_uncheckable == num_checked:
        warning = (
            f"Mass-balance diagnosis is meaningless for this model: all "
            f"{num_checked} checked reactions contain a metabolite with no "
            f"formula, so check_mass_balance() returns empty. A count of "
            f"{num_mass} mass imbalances does NOT mean the model is balanced -- "
            f"the formulas are missing. Annotate metabolite formulas first."
        )
    elif num_mass_uncheckable:
        warning = (
            f"{num_mass_uncheckable} of {num_checked} checked reactions contain "
            f"a formula-less metabolite and cannot be mass-checked; they are "
            f"absent from the imbalance count but are not verified as balanced."
        )

    return {
        "count": count,
        "num_charge_imbalanced": num_charge,
        "num_mass_imbalanced": num_mass,
        "num_boundary_skipped": boundary_skipped,
        "truncated": truncated,
        "diagnosability": {
            "num_reactions_checked": num_checked,
            "num_mass_checkable": num_mass_checkable,
            "num_mass_uncheckable": num_mass_uncheckable,
            "num_charge_checkable": num_charge_checkable,
            "num_charge_uncheckable": num_charge_uncheckable,
            "warning": warning,
        },
        "reactions": imbalanced,
    }


@mcp.tool()
def get_reaction(reaction_id: str) -> dict:
    """Return full detail on one reaction: stoichiometry, metabolites, bounds, GPR.

    Args:
        reaction_id: The reaction identifier (e.g. 'R5').

    Returns:
        A detail dict, or an error dict if the reaction is not found.
    """
    err = _require_model()
    if err:
        return err
    model = _STATE["model"]

    try:
        rxn = model.reactions.get_by_id(reaction_id)
    except KeyError:
        return {
            "error": "reaction_not_found",
            "message": f"No reaction with id '{reaction_id}'.",
        }

    metabolites = []
    for met, coeff in rxn.metabolites.items():
        metabolites.append(
            {
                "id": met.id,
                "name": met.name,
                "coefficient": coeff,
                "role": "reactant" if coeff < 0 else "product",
                "formula": met.formula,
                "charge": met.charge,
                "compartment": met.compartment,
            }
        )

    parts = _split_imbalance(_mass_balance(rxn))

    return {
        "id": rxn.id,
        "name": rxn.name,
        "reaction_string": rxn.reaction,
        "lower_bound": rxn.lower_bound,
        "upper_bound": rxn.upper_bound,
        "reversibility": rxn.reversibility,
        "gene_reaction_rule": rxn.gene_reaction_rule,
        "genes": sorted(g.id for g in rxn.genes),
        "subsystem": rxn.subsystem,
        "metabolites": metabolites,
        "balances": not (
            parts["charge_imbalance"] is not None or parts["mass_imbalance"]
        ),
        **parts,
    }


@mcp.tool()
def get_metabolite(metabolite_id: str) -> dict:
    """Return a metabolite's formula, charge, compartment, and participating reactions.

    Args:
        metabolite_id: The metabolite identifier (e.g. 'm1[C_cy]').

    Returns:
        A detail dict, or an error dict if the metabolite is not found.
    """
    err = _require_model()
    if err:
        return err
    model = _STATE["model"]

    try:
        met = model.metabolites.get_by_id(metabolite_id)
    except KeyError:
        return {
            "error": "metabolite_not_found",
            "message": f"No metabolite with id '{metabolite_id}'.",
        }

    reactions = []
    for rxn in met.reactions:
        coeff = rxn.get_coefficient(met.id)
        reactions.append(
            {
                "id": rxn.id,
                "name": rxn.name,
                "coefficient": coeff,
                "role": "reactant" if coeff < 0 else "product",
                "reaction_string": rxn.reaction,
            }
        )

    return {
        "id": met.id,
        "name": met.name,
        "formula": met.formula,
        "charge": met.charge,
        "compartment": met.compartment,
        "num_reactions": len(reactions),
        "reactions": reactions,
    }


@mcp.tool()
def run_fba() -> dict:
    """Run flux balance analysis and return the objective value and nonzero fluxes.

    Returns:
        {status, objective_value, objective_expression, num_nonzero_fluxes,
         fluxes: {reaction_id: flux}} sorted by descending flux magnitude, or an
        error dict. On an infeasible/unbounded solve, status reflects that and
        fluxes is empty.
    """
    err = _require_model()
    if err:
        return err
    model = _STATE["model"]

    try:
        solution = model.optimize()
    except Exception as exc:
        return {"error": "fba_failed", "message": str(exc)}

    try:
        objective_expression = str(model.objective.expression)
    except Exception:
        objective_expression = None

    nonzero: dict[str, float] = {}
    if solution.status == "optimal":
        # solution.fluxes is a pandas Series keyed by reaction id.
        for rxn_id, flux in solution.fluxes.items():
            if abs(flux) > ZERO_TOL:
                nonzero[rxn_id] = float(flux)

    # Sort by descending magnitude so the agent sees the dominant fluxes first.
    nonzero = dict(
        sorted(nonzero.items(), key=lambda kv: abs(kv[1]), reverse=True)
    )

    return {
        "status": solution.status,
        "objective_value": (
            float(solution.objective_value)
            if solution.objective_value is not None
            else None
        ),
        "objective_expression": objective_expression,
        "num_nonzero_fluxes": len(nonzero),
        "fluxes": nonzero,
    }


if __name__ == "__main__":
    # stdio transport — what Claude Desktop launches and speaks to.
    mcp.run()
