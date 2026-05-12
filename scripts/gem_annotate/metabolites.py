"""
metabolites.py — metabolite annotation and H+/H2O balancing.
"""

import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

# Strategy BD: direct MNXM ID lookup for metabolites whose name_index entry is
# absent or unreliable (e.g. WATER is a special non-MNXM ID).
_DIRECT_MNXM_TABLE: dict[str, str] = {
    "h2o": "WATER",
    "water": "WATER",
    "h+": "MNXM1",
    "h(+)": "MNXM1",
    "h+_p+1": "MNXM1",   # iYli21 charge-notation name for proton
    "proton": "MNXM1",
    "o2": "MNXM735438",
    "oxygen": "MNXM735438",
    "dioxygen": "MNXM735438",
    "co2": "MNXM13",
    "carbon dioxide": "MNXM13",
    "iron": "MNXM147471",
    "2-oxoadipic acid": "MNXM263",
    "2-oxoadipate": "MNXM263",
    "aminoacetone": "MNXM1106",
    "5-phospho-alpha-d-ribose 1-diphosphate": "MNXM1104453",
}

# Strategy B1: common synonym → MetaNetX canonical name (must exist in name_index)
_SYNONYM_TABLE: dict[str, str] = {
    "h2o": "water",
    "oxygen": "dioxygen",
    "o2": "dioxygen",
    "h+": "hydron",
    "proton": "hydron",
    "co2": "carbon dioxide",
    "coa": "coenzyme A",
    "prpp": "5-phospho-alpha-D-ribose 1-diphosphate",
    "adenylo-succinate": "adenylosuccinate",
    "s-adenosylmethioninamine": "S-adenosylmethioninamine",
    "aminoacetone": "aminoacetone",
    "homocitrate": "homocitrate",
    "2-oxoadipic acid": "2-oxoadipate",
    "adenosine 3',5'-bismonophosphate": "adenosine 3',5'-bisphosphate",
    "5'-adenylyl sulfate": "adenosine 5'-phosphosulfate",
    "2-isopropylmaleic acid": "2-isopropylmaleate",
}

_FORMULA_RE  = re.compile(r"^[A-Z][A-Za-z0-9\+\-\(\)\.]*$")
_ELEMENT_RE  = re.compile(r"[CHONSP]")
_STRIP_RE    = re.compile(r"[^a-z0-9]")
# Charge-notation suffix: p+1, p-2, n+1, charge+1, … — not a formula
_CHARGE_RE   = re.compile(r"^[pnq][\+\-]\d+$|^charge[\+\-]\d+$", re.IGNORECASE)


def _parse_name_formula(raw_name: str) -> tuple[str, str]:
    """
    iYli21 name format: 'chemical name_FORMULA'
    Returns (chem_name, formula).  formula may be '' if not parseable.

    Charge notation suffixes like 'p+1', 'p-2' are not formulas; strip them
    and return only the chemical name portion.
    """
    if "_" not in raw_name:
        return raw_name.strip(), ""
    chem_name, maybe_formula = raw_name.rsplit("_", 1)
    maybe_formula = maybe_formula.strip()
    # Reject charge-notation tokens like 'p+1' before formula check
    if _CHARGE_RE.match(maybe_formula):
        return chem_name.strip().strip("_"), ""
    # Require: starts with uppercase AND contains at least one element letter.
    if (maybe_formula
            and _FORMULA_RE.match(maybe_formula)
            and _ELEMENT_RE.search(maybe_formula)):
        return chem_name.strip().strip("_"), maybe_formula
    return raw_name.strip().strip("_"), ""


# HMDB00xxxxx — 7-digit padded form where the first two digits are "00"
# (safe to truncate to 5-digit form HMDB00171 = HMDB + last 5 digits)
_HMDB_PADDED_RE = re.compile(r"^HMDB00(\d{5})$")


def _normalize_annotation(ann: dict) -> dict:
    """
    Normalise annotation values in-place and return the same dict.

    chebi  — ensure every value is "CHEBI:<digits>"; deduplicate.
    hmdb   — canonicalise to Memote's ^HMDB\\d{5}$ format (5-digit suffix).
             HMDB0000171 (7-digit padded, leading "00") → HMDB00171.
             Genuine 7-digit IDs whose digits 5-6 are not "00" are kept as-is.
    others — deduplicate list values.
    """
    for key, val in ann.items():
        if not isinstance(val, list):
            continue

        if key == "chebi":
            seen: set[str] = set()
            for v in val:
                v = str(v)
                if v.upper().startswith("CHEBI:"):
                    seen.add(v.upper())
                elif v.isdigit():
                    seen.add(f"CHEBI:{v}")
            ann[key] = sorted(seen)

        elif key == "hmdb":
            seen: set[str] = set()
            for v in val:
                v = str(v)
                m = _HMDB_PADDED_RE.match(v.upper())
                if m:
                    v = f"HMDB{m.group(1)}"   # strip the leading "00" padding
                seen.add(v)
            ann[key] = sorted(seen)

        else:
            # basic dedup, preserve order
            seen_vals: list = []
            seen_set: set = set()
            for v in val:
                if v not in seen_set:
                    seen_set.add(v)
                    seen_vals.append(v)
            ann[key] = seen_vals

    return ann


def _apply_mnxm(met, mnxm_id: str, by_mnxid: dict, prop: dict) -> dict:
    """Build annotation dict for a matched MNXM ID and update met formula/charge."""
    new_ann: dict = defaultdict(list)
    for db_prefix, db_id in by_mnxid.get(mnxm_id, []):
        new_ann[db_prefix].append(db_id)
    new_ann["metanetx.chemical"] = [mnxm_id]

    if mnxm_id in prop:
        p = prop[mnxm_id]
        if p["inchi"]:
            new_ann["inchi"] = p["inchi"]
        if p["inchikey"]:
            new_ann["inchikey"] = p["inchikey"]
        if p["formula"] and not met.formula and _ELEMENT_RE.search(p["formula"]):
            met.formula = p["formula"]
        if p["charge"] not in ("", "NA"):
            try:
                met.charge = int(float(p["charge"]))
            except (ValueError, TypeError):
                pass
    return new_ann


def annotate_metabolites(model, chem_xref: dict, chem_prop_data: dict) -> None:
    """
    For every metabolite:
      - Extract chemical formula from name (2a)
      - Lookup MNXM via name → get cross-refs (1) + InChI/InChIKey (2a)
    Modifies model in-place.

    Match strategies (in priority order):
      A   — bigg.metabolite annotation → by_source lookup
      B   — exact chemical name (case-insensitive) → name_index
      B2a — normalized name (strip non-alnum) → normalized_name_index
      B2b — normalized prefix match (length ≥ 8, within ±3 chars)
      C   — formula + optional charge disambiguation → formula_index
    """
    by_source = chem_xref["by_source"]
    by_mnxid  = chem_xref["by_mnxid"]
    prop      = chem_prop_data["prop"]
    name_index = chem_prop_data["name_index"]

    # Validate _SYNONYM_TABLE values exist in name_index
    for syn_key, canonical in _SYNONYM_TABLE.items():
        if canonical.lower() not in name_index:
            logger.warning(
                f"_SYNONYM_TABLE value not found in name_index: {canonical!r} "
                f"(mapped from {syn_key!r})"
            )

    # Pre-build normalized name index for B2a / B2b
    normalized_name_index: dict[str, str] = {}
    for k, v in name_index.items():
        nk = _STRIP_RE.sub("", k.lower())
        if nk not in normalized_name_index:
            normalized_name_index[nk] = v

    # Pre-build formula index for Strategy C: formula → [mnxm_id, ...]
    formula_index: dict[str, list[str]] = defaultdict(list)
    # Also build mnxm_id → name for C name-disambiguation
    mnxm_name: dict[str, str] = {}
    for mnxm_id, p in prop.items():
        f = p.get("formula", "").strip()
        if f and _ELEMENT_RE.search(f):
            formula_index[f].append(mnxm_id)
        if p.get("name"):
            mnxm_name[mnxm_id] = p["name"].lower()

    formula_from_name = 0
    formula_from_prop = 0
    charge_from_prop  = 0
    hit_A = hit_B = hit_BD = hit_B0 = hit_B1 = hit_B2a = hit_B2b = hit_C = no_match = 0

    for met in model.metabolites:
        chem_name, formula_str = _parse_name_formula(met.name or "")

        # 2a: set formula from name if not already set
        if formula_str and not met.formula:
            met.formula = formula_str
            formula_from_name += 1

        # ── Strategy A ────────────────────────────────────────────────────
        mnxm_id = None
        strategy = None
        bigg_raw = met.annotation.get("bigg.metabolite")
        if bigg_raw:
            bigg_id = bigg_raw[0] if isinstance(bigg_raw, list) else bigg_raw
            mnxm_id = by_source.get(f"bigg.metabolite:{bigg_id}")
            if mnxm_id:
                strategy = "A"

        # ── Strategy BD: direct MNXM table (covers WATER + other specials) ─
        # Runs before fuzzy strategies so known-safe mappings are never overridden.
        if mnxm_id is None:
            mnxm_id = _DIRECT_MNXM_TABLE.get(chem_name.lower())
            if mnxm_id:
                strategy = "BD"

        # ── Strategy B: exact name ────────────────────────────────────────
        if mnxm_id is None:
            mnxm_id = name_index.get(chem_name.lower())
            if mnxm_id:
                strategy = "B"

        # ── Strategy B0: Excel corruption fix ────────────────────────────
        # Names like "foo_ActiveX VT_ERROR:" — strip the garbage suffix
        if mnxm_id is None and ("ActiveX" in chem_name or "VT_ERROR" in chem_name):
            clean = chem_name.split("_")[0].strip()
            mnxm_id = name_index.get(clean.lower())
            if mnxm_id:
                chem_name = clean
                strategy = "B0"

        # ── Strategy B1: synonym table ────────────────────────────────────
        if mnxm_id is None:
            canonical = _SYNONYM_TABLE.get(chem_name.lower())
            if canonical:
                mnxm_id = name_index.get(canonical.lower())
                if mnxm_id:
                    strategy = "B1"

        # ── Strategy B2a: normalized name ────────────────────────────────
        norm = _STRIP_RE.sub("", chem_name.lower())
        if mnxm_id is None:
            mnxm_id = normalized_name_index.get(norm)
            if mnxm_id:
                strategy = "B2a"

        # ── Strategy B2b: normalized prefix match ─────────────────────────
        if mnxm_id is None and len(norm) >= 8:
            for k, v in normalized_name_index.items():
                if len(k) < 8:
                    continue
                diff = len(k) - len(norm)
                if -3 <= diff <= 3:
                    if k.startswith(norm) or norm.startswith(k):
                        mnxm_id = v
                        strategy = "B2b"
                        break

        # ── Strategy C: formula match ─────────────────────────────────────
        if mnxm_id is None and met.formula:
            candidates = formula_index.get(met.formula.strip(), [])
            if len(candidates) == 1:
                mnxm_id = candidates[0]
                strategy = "C"
            elif len(candidates) > 1:
                # Try charge disambiguation with the model's charge value
                if met.charge is not None:
                    charge_str = str(int(met.charge))
                    disamb = [c for c in candidates if prop[c].get("charge") == charge_str]
                    if len(disamb) == 1:
                        mnxm_id = disamb[0]
                        strategy = "C"
                # If model charge is 0 (likely unset default) and all candidates
                # share the same charge, the formula is unambiguous in practice
                if mnxm_id is None and (met.charge == 0 or met.charge is None):
                    candidate_charges = {prop[c].get("charge", "") for c in candidates}
                    if len(candidate_charges) == 1:
                        mnxm_id = candidates[0]
                        strategy = "C"
                # Name token-overlap disambiguation
                if mnxm_id is None:
                    query_tokens = set(re.split(r"[^a-z0-9]", chem_name.lower())) - {""}
                    if query_tokens:
                        best_score, best_cand = -1, None
                        for cand in candidates:
                            cand_name = mnxm_name.get(cand, "")
                            cand_tokens = set(re.split(r"[^a-z0-9]", cand_name)) - {""}
                            score = len(query_tokens & cand_tokens)
                            if score > best_score:
                                best_score, best_cand = score, cand
                        if best_score > 0:
                            mnxm_id = best_cand
                            strategy = "C"

        if mnxm_id is None:
            no_match += 1
            continue

        new_ann = _normalize_annotation(_apply_mnxm(met, mnxm_id, by_mnxid, prop))

        # Track per-strategy formula/charge enrichment
        if met.formula and formula_str and met.formula != formula_str:
            formula_from_prop += 1
        p = prop.get(mnxm_id, {})
        if p.get("charge") not in ("", "NA", None):
            charge_from_prop += 1

        # Merge: keep existing keys, fill in missing ones only
        merged = dict(met.annotation)
        for key, val in new_ann.items():
            if key not in merged:
                merged[key] = val
        met.annotation = merged

        if strategy == "A":
            hit_A += 1
        elif strategy == "B":
            hit_B += 1
        elif strategy == "BD":
            hit_BD += 1
        elif strategy == "B0":
            hit_B0 += 1
        elif strategy == "B1":
            hit_B1 += 1
        elif strategy == "B2a":
            hit_B2a += 1
        elif strategy == "B2b":
            hit_B2b += 1
        else:
            hit_C += 1

    total_hit = hit_A + hit_B + hit_BD + hit_B0 + hit_B1 + hit_B2a + hit_B2b + hit_C
    logger.info(
        f"Metabolites annotated: {total_hit} | "
        f"A={hit_A}  B(exact)={hit_B}  BD(direct)={hit_BD}  B0(excel)={hit_B0}  "
        f"B1(synonym)={hit_B1}  B2a(norm)={hit_B2a}  B2b(prefix)={hit_B2b}  "
        f"C(formula)={hit_C}  unmatched={no_match}"
    )
    logger.info(
        f"Formulas: {formula_from_name} from name, {formula_from_prop} from prop | "
        f"Charges: {charge_from_prop} from prop"
    )


def normalize_all_annotations(model) -> None:
    """
    Walk every metabolite, reaction, and gene annotation and:
      - Merge keys that differ only in case into their canonical lowercase form.
      - Apply chebi / hmdb value normalisation via _normalize_annotation.
    Modifies model in-place.  Safe to call multiple times (idempotent).

    Also acts as a safety net for scalar→list normalization: any annotation
    value that slipped through as a plain string is caught by _normalize_annotation.
    """
    met_fixed = rxn_fixed = gene_fixed = 0

    def _clean(obj) -> bool:
        ann = obj.annotation
        if not isinstance(ann, dict) or not ann:
            return False
        changed = False

        # ── Step 1: merge case-variant keys into lowercase canonical ─────
        lower_map: dict[str, list[str]] = {}
        for k in list(ann.keys()):
            lower_map.setdefault(k.lower(), []).append(k)

        for canonical, variants in lower_map.items():
            if len(variants) == 1 and variants[0] == canonical:
                continue
            # collect all values under all variant keys
            combined: list = []
            for v in variants:
                val = ann.pop(v)
                if isinstance(val, list):
                    combined.extend(val)
                else:
                    combined.append(val)
            ann[canonical] = combined
            changed = True

        # ── Step 2: normalise chebi / hmdb values; dedup all lists ───────
        _normalize_annotation(ann)
        return changed

    for met in model.metabolites:
        if _clean(met):
            met_fixed += 1
    for rxn in model.reactions:
        if _clean(rxn):
            rxn_fixed += 1
    for gene in model.genes:
        if _clean(gene):
            gene_fixed += 1

    logger.info(
        f"normalize_all_annotations: fixed keys/values in "
        f"{met_fixed} metabolites, {rxn_fixed} reactions, {gene_fixed} genes"
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
            # add_metabolites uses combine=True, so pass only the delta to add
            additions[water] = -o_imb
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
                # add_metabolites uses combine=True, so pass only the delta to add
                additions[proton] = -h_residual

        if additions:
            rxn.add_metabolites(additions)
            fixed += 1

    logger.info(
        f"Proton/water balance: {fixed} reactions fixed, "
        f"{skipped_no_formula} skipped (missing formula)"
    )
