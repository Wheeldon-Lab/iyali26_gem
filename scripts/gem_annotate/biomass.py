"""
biomass.py — biomass reaction diagnosis and fixing.
"""

import logging
import re

from .metabolites import _bigg_id

logger = logging.getLogger(__name__)


def _formula_mw(formula: str) -> float:
    """
    Compute molecular weight (g/mol) from a molecular formula string.
    Handles simple formulas like C6H12O6N2P (no parentheses, no charge).
    Returns 0.0 if formula is empty or unparseable.
    """
    atomic_weights = {
        "C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007,
        "P": 30.974, "S": 32.06,  "Fe": 55.845, "Mg": 24.305,
        "Zn": 65.38, "Cu": 63.546, "Mn": 54.938, "Co": 58.933,
        "Mo": 95.96, "Se": 78.96, "Ca": 40.078, "Na": 22.990,
        "K": 39.098, "Cl": 35.45, "R": 0.0,  # R = unknown residue → 0
    }
    if not formula:
        return 0.0
    mw = 0.0
    for element, count_str in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        count = int(count_str) if count_str else 1
        mw += atomic_weights.get(element, 0.0) * count
    return mw


def fix_biomass_reaction(model) -> None:
    """
    Diagnose and partially fix the biomass reaction R1372.

    Step 1 — Open bounds (was [0,0]).
    Step 2 — Mass normalisation: if the sum of reactant masses deviates more
              than 1 % from 1.0 g/gDW, rescale all stoichiometric coefficients
              proportionally so the reaction represents exactly 1 g/gDW biomass.
              This is a destructive operation; all changes are logged.
    Step 3 — Blocked-precursor diagnosis: temporarily open all exchange
              reactions, then check whether each reactant can carry any flux.
              Blocked precursors are reported (not auto-fixed).
    Step 4 — GAM check: verify that the ATP-hydrolysis maintenance term
              (atp + h2o → adp + pi + h) is present; report if missing.
              A typical value is ~30 mmol/gDW/h but requires biological
              validation before adding, so this step only reports.
    """
    if "R1372" not in [r.id for r in model.reactions]:
        logger.warning("R1372 not found in model")
        return

    rxn = model.reactions.get_by_id("R1372")
    logger.info(f"R1372 current bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")

    # ── Step 1: open bounds ──────────────────────────────────────────────────
    if rxn.lower_bound == 0 and rxn.upper_bound == 0:
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
        logger.info("R1372: bounds opened to [0, 1000]")

    reactants = [(m, abs(v)) for m, v in rxn.metabolites.items() if v < 0]
    products  = [(m, abs(v)) for m, v in rxn.metabolites.items() if v > 0]
    logger.info(f"R1372 composition: {len(reactants)} reactants, {len(products)} products")

    # ── Step 2: mass normalisation ───────────────────────────────────────────
    total_mass = sum(
        coeff * _formula_mw(met.formula)
        for met, coeff in reactants
        if met.formula
    )
    logger.info(f"R1372 reactant mass sum: {total_mass:.4f} g/gDW (target ≈ 1.0)")

    if total_mass > 0 and abs(total_mass - 1.0) > 0.01:
        scale = 1.0 / total_mass
        logger.warning(
            f"R1372 mass deviates {total_mass:.4f} g/gDW from 1.0; "
            f"rescaling all stoichiometric coefficients by {scale:.6f}"
        )
        new_stoich: dict = {}
        for met, v in rxn.metabolites.items():
            new_stoich[met] = v * scale
            logger.info(
                f"  {met.id}: {v:.6f} → {v * scale:.6f}"
            )
        rxn.add_metabolites(new_stoich, combine=False)
        logger.info("R1372: stoichiometry rescaled")
    else:
        logger.info("R1372: mass sum within 1 % of 1.0 — no rescaling needed")

    # Check mass balance after potential rescaling
    try:
        balance = rxn.check_mass_balance()
        if balance:
            logger.warning(f"R1372 residual mass imbalance: {balance}")
        else:
            logger.info("R1372 mass balance: OK")
    except Exception as e:
        logger.warning(f"R1372 mass balance check failed: {e}")

    # ── Step 3: blocked-precursor diagnosis ─────────────────────────────────
    logger.info("R1372: checking which precursors can carry flux …")
    blocked_precursors = []

    with model:
        model.solver = "glpk"
        # Temporarily open all exchange reactions so precursor availability is
        # not limited by missing transport — we test the network topology only.
        for ex in model.exchanges:
            ex.lower_bound = -1000
            ex.upper_bound = 1000

        for met, coeff in reactants:
            # Create a temporary demand reaction for this metabolite
            with model:
                demand_id = f"_tmp_demand_{met.id}"
                demand = model.add_boundary(met, type="demand", reaction_id=demand_id)
                model.objective = demand
                sol = model.optimize()
                flux = sol.objective_value if sol.status == "optimal" else 0.0
                if flux < 1e-6:
                    blocked_precursors.append((met.id, met.name))

    if blocked_precursors:
        logger.warning(
            f"R1372: {len(blocked_precursors)} precursors cannot carry flux "
            f"(network gaps — need manual curation):"
        )
        for mid, mname in blocked_precursors:
            logger.warning(f"  BLOCKED: {mid}  ({mname})")
    else:
        logger.info("R1372: all precursors can carry flux")

    # ── Step 4: GAM (growth-associated maintenance) check ───────────────────
    # GAM term: atp + h2o → adp + pi + h
    gam_reactants = {"atp", "h2o"}
    gam_products  = {"adp", "pi", "h"}

    rxn_reactant_bigg = {
        _bigg_id(m) for m, v in rxn.metabolites.items() if v < 0
    }
    rxn_product_bigg = {
        _bigg_id(m) for m, v in rxn.metabolites.items() if v > 0
    }

    missing_r = gam_reactants - rxn_reactant_bigg
    missing_p = gam_products  - rxn_product_bigg

    if missing_r or missing_p:
        logger.warning(
            "R1372: GAM term (atp + h2o → adp + pi + h) appears INCOMPLETE. "
            "Missing reactants: %s  Missing products: %s. "
            "Typical GAM coefficient is ~30 mmol/gDW/h — add manually after "
            "biological validation.", missing_r or "none", missing_p or "none"
        )
    else:
        logger.info("R1372: GAM term (atp + h2o → adp + pi + h) is present")

    logger.info("R1372: diagnosis complete — review warnings above for manual curation")
