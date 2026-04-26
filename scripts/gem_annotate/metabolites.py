"""
metabolites.py — metabolite annotation and H+/H2O balancing.
"""

import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

_FORMULA_RE = re.compile(r"^[A-Z][A-Za-z0-9]*$")   # no *, +, -, spaces


def _parse_name_formula(raw_name: str) -> tuple[str, str]:
    """
    iYli21 name format: 'chemical name_FORMULA'
    Returns (chem_name, formula).  formula may be '' if not parseable.
    """
    if "_" not in raw_name:
        return raw_name.strip(), ""
    chem_name, maybe_formula = raw_name.rsplit("_", 1)
    maybe_formula = maybe_formula.strip()
    # Accept as formula if it looks like a molecular formula
    # or is a charge token like 'p+1', 'H2O', etc.
    if maybe_formula and _FORMULA_RE.match(maybe_formula):
        return chem_name.strip(), maybe_formula
    # charge tokens like 'p+1' — not a formula
    return raw_name.strip(), ""


def annotate_metabolites(model, chem_xref: dict, chem_prop_data: dict) -> None:
    """
    For every metabolite:
      - Extract chemical formula from name (2a)
      - Lookup MNXM via name → get cross-refs (1) + InChI/InChIKey (2a)
    Modifies model in-place.
    """
    by_source = chem_xref["by_source"]
    by_mnxid = chem_xref["by_mnxid"]
    prop = chem_prop_data["prop"]
    name_index = chem_prop_data["name_index"]

    formula_from_name = 0
    formula_from_prop = 0
    charge_from_prop = 0
    annotated = 0
    no_match = 0

    for met in model.metabolites:
        chem_name, formula_str = _parse_name_formula(met.name)

        # 2a: set formula from name if not already set
        if formula_str and not met.formula:
            met.formula = formula_str
            formula_from_name += 1

        # 1: look up MNXM
        # Strategy A: bigg.metabolite already in annotation (future-proof)
        mnxm_id = None
        bigg_raw = met.annotation.get("bigg.metabolite")
        if bigg_raw:
            bigg_id = bigg_raw[0] if isinstance(bigg_raw, list) else bigg_raw
            mnxm_id = by_source.get(f"bigg.metabolite:{bigg_id}")

        # Strategy B: match by chemical name (case-insensitive)
        if mnxm_id is None:
            mnxm_id = name_index.get(chem_name.lower())

        if mnxm_id is None:
            no_match += 1
            continue

        # Build new annotation dict
        new_ann: dict[str, list] = defaultdict(list)
        for db_prefix, db_id in by_mnxid.get(mnxm_id, []):
            new_ann[db_prefix].append(db_id)
        new_ann["metanetx.chemical"] = [mnxm_id]

        # Enrich with InChI, InChIKey, formula from chem_prop
        if mnxm_id in prop:
            p = prop[mnxm_id]
            if p["inchi"]:
                new_ann["inchi"] = p["inchi"]
            if p["inchikey"]:
                new_ann["inchikey"] = p["inchikey"]
            if p["formula"] and not met.formula and _FORMULA_RE.match(p["formula"]):
                met.formula = p["formula"]
                formula_from_prop += 1
            # Write MNX charge back to met.charge. This exposes the real charge
            # balance (memote previously showed a spurious 100% green because
            # no metabolite had a charge value).
            if p["charge"] not in ("", "NA"):
                try:
                    met.charge = int(float(p["charge"]))
                    charge_from_prop += 1
                except (ValueError, TypeError):
                    pass

        # Merge (keep existing, add new)
        merged = dict(met.annotation)
        merged.update(new_ann)
        met.annotation = merged
        annotated += 1

    logger.info(
        f"Metabolites: {annotated} annotated, {no_match} unmatched | "
        f"Formulas: {formula_from_name} from name, {formula_from_prop} from prop | "
        f"Charges: {charge_from_prop} from prop"
    )


def _bigg_id(met) -> str:
    """Return the bigg.metabolite annotation value as a plain string, or ''."""
    raw = met.annotation.get("bigg.metabolite")
    if not raw:
        return ""
    return raw[0] if isinstance(raw, list) else str(raw)


def _build_compartment_met_index(model) -> dict[str, dict[str, object]]:
    """
    Pre-build per-compartment lookup for H+ and H2O, keyed by bigg.metabolite.
    Returns {compartment: {"h": met_or_None, "h2o": met_or_None}}
    """
    index: dict[str, dict] = {}
    for met in model.metabolites:
        bid = _bigg_id(met)
        if bid not in ("h", "h2o"):
            continue
        comp = met.compartment
        if comp not in index:
            index[comp] = {"h": None, "h2o": None}
        if bid == "h" and index[comp]["h"] is None:
            index[comp]["h"] = met
        elif bid == "h2o" and index[comp]["h2o"] is None:
            index[comp]["h2o"] = met
    return index


def fix_proton_water_balance(model) -> None:
    """
    For each internal reaction unbalanced only in H and/or O, add H+ / H2O.

    Uses bigg.metabolite annotation ("h", "h2o") for exact metabolite lookup
    instead of name-string matching, avoiding wrong-compartment assignments.

    For transport reactions spanning multiple compartments the overall H/O
    imbalance is distributed evenly across compartments that have both H+ and
    H2O available; if no compartment qualifies the reaction is skipped.
    """
    comp_index = _build_compartment_met_index(model)
    fixed = 0
    skipped_no_formula = 0

    for rxn in model.reactions:
        # skip exchange/demand/sink reactions (≤1 metabolite)
        if len(rxn.metabolites) <= 1:
            continue
        # skip if any metabolite missing formula
        if any(not met.formula for met in rxn.metabolites):
            skipped_no_formula += 1
            continue

        try:
            balance = rxn.check_mass_balance()
        except Exception:
            continue

        if not balance:
            continue

        elements = set(balance.keys())
        # Only fix if imbalance is purely H and/or O (charge ignored)
        if not elements.issubset({"H", "O", "charge"}):
            continue

        h_imb = balance.get("H", 0)
        o_imb = balance.get("O", 0)
        if h_imb == 0 and o_imb == 0:
            continue

        # Collect compartments present in this reaction
        rxn_comps = {met.compartment for met in rxn.metabolites}

        # Choose target compartment: prefer one that already has both h and h2o
        # in comp_index; fall back to majority compartment.
        def has_needed(comp):
            ci = comp_index.get(comp, {})
            need_h2o = o_imb != 0
            need_h   = (h_imb != 0) or (o_imb != 0 and h_imb - 2 * o_imb != 0)
            if need_h2o and ci.get("h2o") is None:
                return False
            if need_h and ci.get("h") is None and not need_h2o:
                return False
            return True

        candidates = [c for c in rxn_comps if has_needed(c)]
        if candidates:
            # prefer majority compartment among candidates
            comps_list = [met.compartment for met in rxn.metabolites]
            target_comp = max(candidates, key=comps_list.count)
        else:
            # fall back to overall majority compartment
            comps_list = [met.compartment for met in rxn.metabolites]
            target_comp = max(rxn_comps, key=comps_list.count)

        ci = comp_index.get(target_comp, {})

        additions: dict = {}

        if o_imb != 0:
            water = ci.get("h2o")
            if water is None:
                continue
            additions[water] = rxn.metabolites.get(water, 0) - o_imb
            # After adding o_imb H2O molecules, remaining H imbalance:
            h_residual = h_imb - 2 * o_imb
        else:
            h_residual = h_imb

        if h_residual != 0:
            proton = ci.get("h")
            if proton is None and o_imb != 0:
                # H2O already staged — apply it even without proton fix
                pass
            elif proton is None:
                continue
            else:
                additions[proton] = rxn.metabolites.get(proton, 0) - h_residual

        if additions:
            rxn.add_metabolites(additions)
            fixed += 1

    logger.info(
        f"Proton/water balance: {fixed} reactions fixed, "
        f"{skipped_no_formula} skipped (missing formula)"
    )
