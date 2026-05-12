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


# BiGG metabolite base IDs for mineral salts + vitamins, with uptake bounds.
# Y. lipolytica is a biotin auxotroph — biotin must be supplied exogenously.
_MINERAL_SALTS_METS: dict[str, float] = {
    "mg2":   -1000.0,   # Magnesium
    "k":     -1000.0,   # Potassium (also in MINIMAL_MEDIUM_BIGG; harmless duplicate)
    "na1":   -1000.0,   # Sodium
    "btn":   -1000.0,   # biotin (auxotroph — must supply)
    "thm":   -1000.0,   # thiamine
    "pydx":  -1000.0,   # pyridoxine
}

# Reaction-ID fallback for metabolites that lack bigg.metabolite annotation.
# Only used when the BiGG-based lookup above fails.
_MINERAL_SALTS_FALLBACK_IDS: dict[str, float] = {
    "R2061": -1000.0,   # Magnesium
    "R1298": -1000.0,   # Potassium
    "R1323": -1000.0,   # Sodium
    "R1029": -1000.0,   # biotin
    "R1340": -1000.0,   # thiamine
    "R1305": -1000.0,   # pyridoxine
}


def configure_medium(model,
                     extra_exchanges: dict[str, float] | None = None) -> None:
    """
    Extend the current medium to a mineral salts + vitamins formulation.

    Tries to match exchange reactions via their metabolite's bigg.metabolite
    annotation against _MINERAL_SALTS_METS.  Falls back to reaction-ID lookup
    (_MINERAL_SALTS_FALLBACK_IDS) for exchanges that lack BiGG annotation.

    After updating bounds, syncs model.medium so COBRApy's medium dict
    reflects the actual open uptakes.

    Parameters
    ----------
    model           : cobra.Model (modified in-place)
    extra_exchanges : optional additional {rxn_id: lb} overrides (by reaction ID)
    """
    # Build bigg_base → exchange_rxn map for metabolite-based matching
    bigg_to_ex: dict[str, object] = {}
    for ex in model.exchanges:
        if len(ex.metabolites) != 1:
            continue
        met = next(iter(ex.metabolites))
        raw = met.annotation.get("bigg.metabolite") if isinstance(met.annotation, dict) else None
        if raw:
            bigg_id = raw[0] if isinstance(raw, list) else str(raw)
            bigg_base = re.sub(r"_[a-z]$", "", bigg_id)
            if bigg_base not in bigg_to_ex:
                bigg_to_ex[bigg_base] = ex

    opened = []
    matched_by_bigg: set[str] = set()   # reaction IDs matched via BiGG metabolite

    # ── Primary: metabolite BiGG ID matching ─────────────────────────────
    for bigg_base, lb in _MINERAL_SALTS_METS.items():
        ex = bigg_to_ex.get(bigg_base)
        if ex is None:
            continue
        if ex.lower_bound != lb:
            ex.lower_bound = lb
        matched_by_bigg.add(ex.id)
        opened.append(ex.id)
        logger.debug(f"  configure_medium: {ex.id} opened via BiGG metabolite '{bigg_base}'")

    # ── Fallback: reaction ID for anything not matched above ──────────────
    for rxn_id, lb in _MINERAL_SALTS_FALLBACK_IDS.items():
        if rxn_id in matched_by_bigg:
            continue   # already handled via metabolite lookup
        try:
            rxn = model.reactions.get_by_id(rxn_id)
        except KeyError:
            continue
        if rxn.lower_bound != lb:
            rxn.lower_bound = lb
        opened.append(rxn_id)
        logger.debug(f"  configure_medium: {rxn_id} opened via reaction-ID fallback")

    # ── Extra overrides (by reaction ID) ─────────────────────────────────
    if extra_exchanges:
        for rxn_id, lb in extra_exchanges.items():
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                if rxn.lower_bound != lb:
                    rxn.lower_bound = lb
                opened.append(rxn_id)
            except KeyError:
                logger.warning(f"configure_medium: extra exchange '{rxn_id}' not found")

    model.medium = {
        ex.id: abs(ex.lower_bound)
        for ex in model.exchanges
        if ex.lower_bound < 0
    }

    logger.info(
        f"configure_medium: opened {len(opened)} mineral/vitamin exchanges "
        f"({', '.join(opened)}); medium now has {len(model.medium)} components"
    )
