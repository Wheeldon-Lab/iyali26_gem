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
