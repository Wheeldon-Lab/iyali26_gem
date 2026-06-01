"""
patches.py — known data-bug patches for the iYli21 model.

Each patch fixes a discrete, documented error in the source SBML.  Patches
are applied early in the pipeline, immediately after metabolite annotation,
so that downstream mass-balance, FBA, and gap analyses operate on a
chemically correct model.

Currently applied patches:

  1. NADP+ formula fix
     iYli21 stores NADP+ as C21H28N7O17P3 (one H short of the KEGG C00006
     reference C21H29N7O17P3).  This affects 6 compartmental copies and
     ~126 reactions that use NADP+/NADPH, which had cascading effects on
     ceramide synthesis (see patch 2).

  2. Ceramide formula corrections
     iYli21 stores ceramide-1-(C24/C26) with formula C63H125NO6, which is
     inconsistent with the LIPID MAPS / Yeast-GEM reference value
     C42H85NO3 (C24) / C44H89NO3 (C26).  ceramide-2/2'/3/4 are simply
     missing formulas.  These were left wrong/missing because the upstream
     NADP+ error made the relevant reactions unbalanceable — fixing patch 1
     now makes patch 2 work.

     Naming map (verified against LIPID MAPS):
       ceramide-1  = Cer(d18:0/24:0)         = C42H85NO3  / C44H89NO3
       ceramide-2  = Cer(t18:0/24:0)         = C42H85NO4  / C44H89NO4
       ceramide-2' = Cer(d18:0/24:0(2OH))    = C42H85NO4  / C44H89NO4
       ceramide-3  = Cer(t18:0/24:0(2OH))    = C42H85NO5  / C44H89NO5
       ceramide-4  = (extrapolated +1 O from ceramide-3)
                                              C42H85NO6  / C44H89NO6
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


# ── Patch 1: NADP+ formula fix ────────────────────────────────────────────

# iYli21 NADP+ copies (all compartments). Identified by metabolite name
# prefix "NADP(+)_" which is the iYli21 naming convention.
# Formula correction: C21H28N7O17P3 → C21H29N7O17P3 (KEGG C00006).
_NADP_PLUS_OLD_FORMULA = "C21H28N7O17P3"
_NADP_PLUS_NEW_FORMULA = "C21H29N7O17P3"


def _is_nadp_plus(met) -> bool:
    """Identify NADP+ by name (handles iYli21 'NADP(+)_C21H28N7O17P3' naming)."""
    name = (met.name or "")
    return name.startswith("NADP(+)") or name.lower().startswith("nadp(+)")


def fix_nadp_plus_formula(model) -> int:
    """
    Fix the NADP+ formula bug across all compartments.

    Returns the number of metabolites patched.
    """
    fixed = 0
    for met in model.metabolites:
        if not _is_nadp_plus(met):
            continue
        if met.formula == _NADP_PLUS_OLD_FORMULA:
            met.formula = _NADP_PLUS_NEW_FORMULA
            fixed += 1
            logger.debug(f"  NADP+ patch: {met.id} formula → {_NADP_PLUS_NEW_FORMULA}")
        elif met.formula == _NADP_PLUS_NEW_FORMULA:
            pass   # already correct
        else:
            logger.warning(
                f"  NADP+ patch: {met.id} has unexpected formula {met.formula!r} "
                f"— skipped (expected {_NADP_PLUS_OLD_FORMULA})"
            )
    return fixed


# ── Patch 2: Ceramide formula corrections ─────────────────────────────────

# Map: lowercased base name (compartment-independent) → (formula_C24, formula_C26)
# Base name is what comes BEFORE the trailing "_" or "_FORMULA" in iYli21 names.
# E.g. "ceramide-1-(C24)_" matches base "ceramide-1-(c24)".
_CERAMIDE_FORMULAS_C24 = {
    "ceramide-1-(c24)":  "C42H85NO3",   # Cer(d18:0/24:0)        LMSP02020012
    "ceramide-2-(c24)":  "C42H85NO4",   # Cer(t18:0/24:0)        LMSP02030004
    "ceramide-2'-(c24)": "C42H85NO4",   # Cer(d18:0/24:0(2OH))   (α-OH dihydroceramide)
    "ceramide-3-(c24)":  "C42H85NO5",   # Cer(t18:0/24:0(2OH))   LMSP02030002
    "ceramide-4-(c24)":  "C42H85NO6",   # (extrapolated, R208 product)
}
_CERAMIDE_FORMULAS_C26 = {
    "ceramide-1-(c26)":  "C44H89NO3",
    "ceramide-2-(c26)":  "C44H89NO4",
    "ceramide-2'-(c26)": "C44H89NO4",
    "ceramide-3-(c26)":  "C44H89NO5",
    "ceramide-4-(c26)":  "C44H89NO6",
}
_CERAMIDE_FORMULAS = {**_CERAMIDE_FORMULAS_C24, **_CERAMIDE_FORMULAS_C26}

# Old (wrong) formula on ceramide-1 — explicitly remove this value when we
# overwrite, so the log warns if it's anything else (defensive).
_CERAMIDE_1_KNOWN_BAD_FORMULA = "C63H125NO6"


def _base_name(met) -> str:
    """Return the iYli21 base name (lowercased, strip trailing _ and _FORMULA)."""
    n = (met.name or "").lower().rstrip("_").strip()
    # Strip trailing "_<FORMULA>" if present (e.g. "ceramide-1-(c24)_c63h125no6")
    if "_" in n:
        parts = n.rsplit("_", 1)
        last = parts[-1]
        # heuristic: looks like a chemical formula → strip
        if last and last[0] in "cChH" and any(c.isdigit() for c in last):
            n = parts[0].rstrip("_")
    return n


def fix_ceramide_formulas(model) -> int:
    """
    Fix ceramide-1 wrong formula + fill ceramide-2/2'/3/4 missing formulas.

    Authoritative source: LIPID MAPS LMSD (verified manually) cross-checked
    against Yeast-GEM (https://github.com/SysBioChalmers/yeast-GEM).

    Applies to all compartmental copies (C_er, C_mi, C_go, etc).
    Returns number of metabolites patched.
    """
    fixed = 0
    seen_bases = set()
    for met in model.metabolites:
        base = _base_name(met)
        if base not in _CERAMIDE_FORMULAS:
            continue
        new_formula = _CERAMIDE_FORMULAS[base]
        old_formula = met.formula or ""
        if old_formula == new_formula:
            continue   # already correct
        if old_formula and old_formula != _CERAMIDE_1_KNOWN_BAD_FORMULA:
            logger.warning(
                f"  Ceramide patch: {met.id} ({base!r}) has unexpected existing "
                f"formula {old_formula!r}, overwriting with {new_formula!r}"
            )
        met.formula = new_formula
        fixed += 1
        seen_bases.add(base)
        logger.debug(
            f"  Ceramide patch: {met.id} [{met.compartment}] {old_formula or '(none)'} "
            f"→ {new_formula}  ({base})"
        )
    # Report which base names did not appear in the model
    missing = set(_CERAMIDE_FORMULAS) - seen_bases
    if missing:
        logger.info(
            f"  Ceramide patch: {len(missing)} mapped base name(s) not found in model "
            f"(possibly C24 or C26 variants absent): {sorted(missing)[:5]}…"
        )
    return fixed


# ── Patch 3: free CoA charge fix ──────────────────────────────────────────

# iYli21 stores free coenzyme A (KEGG C00010) with charge 0, but the correct
# physiological charge is -4.  Only the 9 free-CoA copies are affected; the
# acyl-CoA species (tetracosanoyl-CoA, acetyl-CoA, …) already carry charge -4.
# Identified by: name starts with "coenzyme A", formula == C21H36N7O16P3S, charge == 0.
_COA_FREE_FORMULA   = "C21H36N7O16P3S"
_COA_CORRECT_CHARGE = -4


def fix_coa_charge(model) -> int:
    """
    Fix the free-CoA charge bug (KEGG C00010: charge should be -4, not 0).

    Patches only free coenzyme A; acyl-CoA species are left untouched because
    they already carry the correct charge.  Returns the number patched.
    """
    fixed = 0
    for met in model.metabolites:
        name = met.name or ""
        if not name.startswith("coenzyme A"):
            continue
        if met.formula == _COA_FREE_FORMULA and met.charge == 0:
            met.charge = _COA_CORRECT_CHARGE
            fixed += 1
            logger.debug(f"  CoA patch: {met.id} charge → {_COA_CORRECT_CHARGE}")
    return fixed


# ── Top-level driver ──────────────────────────────────────────────────────

def apply_all_patches(model) -> dict:
    """
    Apply all known model patches. Returns a dict with per-patch counts.
    Safe to call multiple times (idempotent — already-correct values are skipped).
    """
    logger.info("Applying iYli21 known-bug patches …")
    counts = {
        "nadp_plus_fixed":   fix_nadp_plus_formula(model),
        "ceramide_fixed":    fix_ceramide_formulas(model),
        "coa_charge_fixed":  fix_coa_charge(model),
    }
    logger.info(
        f"  patches applied: NADP+={counts['nadp_plus_fixed']} copies, "
        f"ceramide={counts['ceramide_fixed']} copies, "
        f"CoA-charge={counts['coa_charge_fixed']} copies"
    )
    return counts
