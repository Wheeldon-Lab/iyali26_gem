"""
exchange.py — exchange bound calibration.
"""

import logging
import re

from .config import MINIMAL_MEDIUM_BIGG, MINIMAL_MEDIUM_NAMES
from .metabolites import _parse_name_formula

logger = logging.getLogger(__name__)


def set_exchange_bounds(model, medium_bigg: dict[str, float] | None = None,
                        medium_names: dict[str, float] | None = None) -> None:
    """
    Calibrate exchange-reaction bounds to a defined minimal medium.

    Problem: by default many exchange reactions are wide open (lb=-1000 or
    ub=1000), creating unbounded flux capacity that makes FBA/FVA physiologically
    meaningless — memote reports ~20% of reactions as having unbounded flux.

    Matching strategy (two-tier, first match wins):
      Tier 1 — bigg.metabolite annotation (exact, reliable):
        Strip the compartment suffix from the BiGG ID (e.g. "glc__D_e" → "glc__D")
        and look up in medium_bigg.
      Tier 2 — chemical name fallback (for metabolites without BiGG annotation):
        Parse the name from the iYli21 "name_FORMULA" convention and look up
        in medium_names (case-insensitive).

    If either tier matches  → set lb = medium value (negative = uptake allowed).
    Otherwise              → set lb = 0 (uptake blocked, secretion still free).

    Secretion is left open (ub = 1000) everywhere because closing it would
    introduce artificial growth blocks.
    """
    if medium_bigg is None:
        medium_bigg = MINIMAL_MEDIUM_BIGG
    if medium_names is None:
        medium_names = MINIMAL_MEDIUM_NAMES

    opened_bigg  = 0
    opened_name  = 0
    closed       = 0
    unchanged    = 0

    for ex in model.exchanges:
        if len(ex.metabolites) != 1:
            continue    # defensive: exchange reactions should have exactly 1 met
        met = next(iter(ex.metabolites))

        # ── Tier 1: bigg.metabolite annotation ────────────────────────────
        new_lb = None
        match_tier = None
        bigg_raw = met.annotation.get("bigg.metabolite")
        if bigg_raw:
            bigg_id = bigg_raw[0] if isinstance(bigg_raw, list) else str(bigg_raw)
            # Strip compartment suffix: "glc__D_e" → "glc__D"
            bigg_base = re.sub(r"_[a-z]$", "", bigg_id)
            if bigg_base in medium_bigg:
                new_lb = medium_bigg[bigg_base]
                match_tier = "bigg"

        # ── Tier 2: chemical name fallback ────────────────────────────────
        if new_lb is None:
            chem_name, _ = _parse_name_formula(met.name)
            key = chem_name.lower().strip()
            if key in medium_names:
                new_lb = medium_names[key]
                match_tier = "name"

        # ── Apply bound ────────────────────────────────────────────────────
        if new_lb is not None:
            if ex.lower_bound != new_lb:
                ex.lower_bound = new_lb
                if match_tier == "bigg":
                    opened_bigg += 1
                else:
                    opened_name += 1
            else:
                unchanged += 1
        else:
            if ex.lower_bound < 0:
                ex.lower_bound = 0.0
                closed += 1
            else:
                unchanged += 1

        # Ensure secretion is allowed (don't over-constrain)
        if ex.upper_bound < 1000:
            ex.upper_bound = 1000.0

    logger.info(
        f"Exchange bounds: {opened_bigg} via BiGG annotation, "
        f"{opened_name} via name fallback, "
        f"{closed} uptake closed, {unchanged} unchanged"
    )


# Reaction IDs for mineral salts + vitamins components not covered by the
# BiGG/name-based tier-1/2 matching in set_exchange_bounds.
# Y. lipolytica is a biotin auxotroph — biotin must be supplied exogenously.
_MINERAL_SALTS_VITAMINS: dict[str, float] = {
    "R2061": -1000.0,   # Magnesium
    "R1298": -1000.0,   # Potassium
    "R1323": -1000.0,   # Sodium
    "R1029": -1000.0,   # biotin (auxotroph — must supply)
    "R1340": -1000.0,   # thiamine
    "R1305": -1000.0,   # pyridoxine
}


def configure_medium(model,
                     extra_exchanges: dict[str, float] | None = None) -> None:
    """
    Extend the current medium to a mineral salts + vitamins formulation.

    Opens the six exchange reactions listed in _MINERAL_SALTS_VITAMINS
    (Mg, K, Na, biotin, thiamine, pyridoxine) that are not matched by the
    BiGG/name tiers in set_exchange_bounds.  Carbon source and other
    nutrients set by set_exchange_bounds are left unchanged.

    After updating bounds, syncs model.medium so COBRApy's medium dict
    reflects the actual open uptakes.

    Parameters
    ----------
    model         : cobra.Model (modified in-place)
    extra_exchanges : optional additional {rxn_id: lb} overrides
    """
    targets = dict(_MINERAL_SALTS_VITAMINS)
    if extra_exchanges:
        targets.update(extra_exchanges)

    opened = []
    missing = []
    for rxn_id, lb in targets.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
        except KeyError:
            missing.append(rxn_id)
            continue
        if rxn.lower_bound != lb:
            rxn.lower_bound = lb
            opened.append(rxn_id)

    if missing:
        logger.warning(f"configure_medium: reaction IDs not found: {missing}")

    model.medium = {
        ex.id: abs(ex.lower_bound)
        for ex in model.exchanges
        if ex.lower_bound < 0
    }

    logger.info(
        f"configure_medium: opened {len(opened)} mineral/vitamin exchanges "
        f"({', '.join(opened)}); medium now has {len(model.medium)} components"
    )
