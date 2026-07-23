"""
update_model.py — full annotation pipeline for iYali26 GEM.

Priority order:
  1. Metabolite cross-database annotation via MetaNetX chem_xref + chem_prop
  2a. Fill missing chemical formulas from metabolite names / chem_prop
  2b. Balance H+/H2O in reactions
  3. Diagnose & fix biomass reaction R1372
  4a. Reaction annotation via MetaNetX reac_xref
  4b. Gene annotation via UniProt REST API
"""

import logging
import re
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
from cobra.io import read_sbml_model, write_sbml_model

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
STARTING_MODEL_PATH = REPO_ROOT / "data" / "iyali26.xml"
OUTPUT_MODEL_PATH = REPO_ROOT / "model.xml"
MNX_DIR = REPO_ROOT / "data" / "metanetx"

# ─────────────────────────────────────────────
# MetaNetX loaders  (polars-fast if available, else pandas)
# ─────────────────────────────────────────────

def _read_tsv(path: Path, names: list[str]) -> pd.DataFrame:
    """Read a MetaNetX TSV (skip # comment lines) into a pandas DataFrame."""
    return pd.read_csv(
        path, sep="\t", comment="#", header=None,
        names=names, dtype=str, low_memory=False,
    ).fillna("")


def load_chem_xref(path: Path) -> dict[str, list[tuple[str, str]]]:
    """
    Returns two indexes from chem_xref.tsv:
      by_source[source_key]  → list of mnx_ids
      by_mnxid[mnx_id]      → list of (db_prefix, db_id)
    Packed into a single dict with those two keys.
    """
    logger.info("Loading chem_xref.tsv …")
    df = _read_tsv(path, ["source", "mnx_id", "description"])
    # only keep rows that map to real MNXM IDs (not BIOMASS / EMPTY)
    df = df[df["mnx_id"].str.startswith("MNXM")]

    by_source: dict[str, str] = {}          # source_key → mnx_id (first hit)
    by_mnxid: dict[str, list] = defaultdict(list)  # mnx_id → [(prefix, id), …]

    for source, mnx_id, _ in df.itertuples(index=False):
        if source not in by_source:
            by_source[source] = mnx_id
        if ":" in source:
            prefix, sid = source.split(":", 1)
            by_mnxid[mnx_id].append((prefix, sid))

    logger.info(f"  {len(by_source):,} source entries, {len(by_mnxid):,} MNXM IDs")
    return {"by_source": by_source, "by_mnxid": dict(by_mnxid)}


def load_chem_prop(path: Path) -> dict[str, dict]:
    """
    Returns prop[mnx_id] = {name, formula, charge, inchi, inchikey, smiles}
    Also builds name_index: lower(name) → mnx_id  (first hit per name)
    """
    logger.info("Loading chem_prop.tsv …")
    df = _read_tsv(path, ["mnx_id", "name", "reference", "formula", "charge",
                           "mass", "inchi", "inchikey", "smiles"])
    df = df[df["mnx_id"].str.startswith("MNXM")]

    prop: dict[str, dict] = {}
    name_index: dict[str, str] = {}   # lower(name) → mnx_id

    for row in df.itertuples(index=False):
        prop[row.mnx_id] = {
            "name": row.name, "formula": row.formula, "charge": row.charge,
            "inchi": row.inchi, "inchikey": row.inchikey, "smiles": row.smiles,
        }
        key = row.name.lower().strip()
        if key and key not in name_index:
            name_index[key] = row.mnx_id

    logger.info(f"  {len(prop):,} compounds, {len(name_index):,} unique names")
    return {"prop": prop, "name_index": name_index}


def load_reac_xref(path: Path) -> dict[str, list[tuple[str, str]]]:
    """
    Returns:
      by_mnxr[mnx_id]         → list of (db_prefix, db_id)
      bigg_to_mnxr[bigg_id]   → mnx_id
      desc_index[lower(name)] → mnx_id   (name = part before '||' in description)
    """
    logger.info("Loading reac_xref.tsv …")
    df = _read_tsv(path, ["source", "mnx_id", "description"])
    df = df[df["mnx_id"].str.startswith("MNXR")]

    by_mnxr: dict[str, list] = defaultdict(list)
    bigg_to_mnxr: dict[str, str] = {}
    desc_index: dict[str, str] = {}    # lower(short name) → mnx_id

    for source, mnx_id, desc in df.itertuples(index=False):
        if ":" in source:
            prefix, sid = source.split(":", 1)
            by_mnxr[mnx_id].append((prefix, sid))
            if prefix == "bigg.reaction" and sid not in bigg_to_mnxr:
                bigg_to_mnxr[sid] = mnx_id

        # description format: "Short name||equation" — index the short name
        short = desc.split("||")[0].strip().lower() if "||" in desc else desc.strip().lower()
        if short and short not in desc_index:
            desc_index[short] = mnx_id

    logger.info(f"  {len(by_mnxr):,} MNXR IDs, {len(bigg_to_mnxr):,} BiGG reactions, {len(desc_index):,} names indexed")
    return {"by_mnxr": dict(by_mnxr), "bigg_to_mnxr": bigg_to_mnxr, "desc_index": desc_index}


_MNXM_IN_EQ = re.compile(r"(MNXM\d+)@")   # extracts MNXM IDs from equation strings


def load_reac_prop(path: Path) -> dict:
    """
    Parse reac_prop.tsv and build a stoichiometric fingerprint index.

    Equation column format:
        1 MNXM01@MNXD1 + 2 MNXM02@MNXD1 = 1 MNXM03@MNXD1

    Returns:
      fingerprint_index : frozenset(MNXM_ids) → list[mnxr_id]
        Maps the *set* of all MNXM participants in a reaction to the MNXR(s)
        that use exactly that set.  Most sets map to one MNXR; collisions
        (different stoichiometries or compartments for the same metabolite set)
        are kept as lists so callers can flag ambiguous matches.
    """
    logger.info("Loading reac_prop.tsv …")
    df = _read_tsv(path, ["mnx_id", "equation", "reference", "classifs", "is_balanced", "is_transport"])
    df = df[df["mnx_id"].str.startswith("MNXR")]

    fingerprint_index: dict[frozenset, list[str]] = defaultdict(list)

    for mnx_id, equation, *_ in df.itertuples(index=False):
        mnxm_ids = frozenset(_MNXM_IN_EQ.findall(equation))
        if len(mnxm_ids) >= 2:          # skip degenerate single-metabolite entries
            fingerprint_index[mnxm_ids].append(mnx_id)

    logger.info(f"  {len(fingerprint_index):,} unique metabolite-set fingerprints indexed")
    return {"fingerprint_index": dict(fingerprint_index)}


# ─────────────────────────────────────────────
# Helper: parse name_formula from iYali26 name strings
# ─────────────────────────────────────────────

_FORMULA_RE = re.compile(r"^[A-Z][A-Za-z0-9]*$")   # no *, +, -, spaces

def _parse_name_formula(raw_name: str) -> tuple[str, str]:
    """
    iYali26 name format: 'chemical name_FORMULA'
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


# ─────────────────────────────────────────────
# Priority 1 + 2a: metabolite annotation + formula fill
# ─────────────────────────────────────────────

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


# ─────────────────────────────────────────────
# Priority 2b: balance H+ / H2O
# ─────────────────────────────────────────────

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


# ─────────────────────────────────────────────
# Priority 2c: exchange-reaction bound calibration
# ─────────────────────────────────────────────

# Minimal aerobic medium for Y. lipolytica W29 on glucose.
# Primary key: BiGG metabolite ID (no compartment suffix).
# Fallback key: lowercased chemical name (before trailing "_FORMULA" in iYali26 names).
# Value: uptake lower bound (mmol / gDW / h).  -1000 = effectively unlimited.
MINIMAL_MEDIUM_BIGG: dict[str, float] = {
    # carbon
    "glc__D":   -10.0,    # D-glucose, sole carbon source
    # nitrogen
    "nh4":      -1000.0,  # ammonium
    # oxygen / aerobic
    "o2":       -1000.0,  # oxygen
    # phosphorus
    "pi":       -1000.0,  # inorganic phosphate
    "h2po4":    -1000.0,  # dihydrogen phosphate (alt BiGG form)
    # sulfur
    "so4":      -1000.0,  # sulphate / sulfate
    # proton and water (exchange reactions are always in extracellular compartment)
    "h":        -1000.0,  # proton
    "h2o":      -1000.0,  # water
    # trace metals
    "fe2":      -1000.0,  # Fe2+
    "fe3":      -1000.0,  # Fe3+
    # macroions
    "k":        -1000.0,  # potassium
    "na1":      -1000.0,  # sodium
    # CO2 (can be re-assimilated)
    "co2":      -1000.0,
}

# Fallback name-based keys for metabolites that lack bigg.metabolite annotation.
# Lowercased chemical name extracted from iYali26 "name_FORMULA" convention.
MINIMAL_MEDIUM_NAMES: dict[str, float] = {
    "d-glucose":      -10.0,
    "glucose":        -10.0,
    "ammonium":       -1000.0,
    "oxygen":         -1000.0,
    "phosphate":      -1000.0,
    "sulphate":       -1000.0,
    "sulfate":        -1000.0,
    "h+":             -1000.0,
    "h2o":            -1000.0,
    "water":          -1000.0,
    "iron":           -1000.0,
    "potassium":      -1000.0,
    "sodium":         -1000.0,
    "carbon dioxide": -1000.0,
    "co2":            -1000.0,
}


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
        Parse the name from the iYali26 "name_FORMULA" convention and look up
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


# ─────────────────────────────────────────────
# Priority 3: diagnose & fix biomass reaction R1372
# ─────────────────────────────────────────────

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


# ─────────────────────────────────────────────
# Priority 4a: reaction annotation via reac_xref
# ─────────────────────────────────────────────

def annotate_reactions(model, reac_xref: dict, reac_prop: dict | None = None) -> None:
    """
    Annotate reactions with MNXR IDs and cross-database identifiers.

    Three strategies, applied in order until one succeeds:

    A — bigg.reaction annotation → bigg_to_mnxr lookup  (fast, high-precision)
    B — reaction name → MetaNetX description short-name  (string match)
    C — metabolite-set fingerprint in reac_prop.tsv      (structure-based)
        Requires reac_prop dict from load_reac_prop().
        Uses the frozenset of MNXM IDs annotated on the reaction's metabolites.
        Unique match → high-confidence annotation.
        Multiple matches → logged as ambiguous, skipped.

    Modifies model in-place.
    """
    by_mnxr = reac_xref["by_mnxr"]
    desc_index = reac_xref["desc_index"]
    bigg_to_mnxr = reac_xref["bigg_to_mnxr"]
    fingerprint_index = (reac_prop or {}).get("fingerprint_index", {})

    annotated_A = 0
    annotated_B = 0
    annotated_C = 0
    ambiguous_C = 0
    no_match = 0

    def _apply_annotation(rxn, mnxr_id: str) -> None:
        new_ann: dict[str, list] = defaultdict(list)
        for db_prefix, db_id in by_mnxr.get(mnxr_id, []):
            new_ann[db_prefix].append(db_id)
        new_ann["metanetx.reaction"] = [mnxr_id]
        merged = dict(rxn.annotation)
        merged.update(new_ann)
        rxn.annotation = merged

    for rxn in model.reactions:
        mnxr_id = None
        strategy = None

        # Strategy A: existing bigg.reaction annotation
        bigg_raw = rxn.annotation.get("bigg.reaction")
        if bigg_raw:
            bid = bigg_raw[0] if isinstance(bigg_raw, list) else bigg_raw
            mnxr_id = bigg_to_mnxr.get(bid)
            if mnxr_id:
                strategy = "A"

        # Strategy B: reaction name → MetaNetX description short-name
        if mnxr_id is None:
            mnxr_id = desc_index.get(rxn.name.lower().strip())
            if mnxr_id:
                strategy = "B"

        # Strategy C: stoichiometric fingerprint via metabolite MNXM set
        if mnxr_id is None and fingerprint_index:
            mnxm_set = frozenset(
                ann[0] if isinstance(ann, list) else ann
                for met in rxn.metabolites
                for ann in [met.annotation.get("metanetx.chemical")]
                if ann
            )
            if len(mnxm_set) >= 2:
                candidates = fingerprint_index.get(mnxm_set, [])
                if len(candidates) == 1:
                    mnxr_id = candidates[0]
                    strategy = "C"
                elif len(candidates) > 1:
                    logger.debug(
                        f"R:{rxn.id} ({rxn.name!r}): ambiguous fingerprint match "
                        f"→ {candidates} — skipping"
                    )
                    ambiguous_C += 1

        if mnxr_id is None:
            no_match += 1
            continue

        _apply_annotation(rxn, mnxr_id)
        if strategy == "A":
            annotated_A += 1
        elif strategy == "B":
            annotated_B += 1
        else:
            annotated_C += 1

    total = annotated_A + annotated_B + annotated_C
    logger.info(
        f"Reactions annotated: {total} total  "
        f"(A={annotated_A}, B={annotated_B}, C={annotated_C}), "
        f"ambiguous-C={ambiguous_C}, unmatched={no_match}"
    )


# ─────────────────────────────────────────────
# Priority 4b: gene annotation via UniProt REST
# ─────────────────────────────────────────────

# ── Gene annotation constants ─────────────────────────────────────────────────
_UNIPROT_SEARCH_URL  = "https://rest.uniprot.org/uniprotkb/search"
_KEGG_CONV_URL       = "https://rest.kegg.jp/conv/uniprot/yli"
_TIER_B_LIMIT        = 1100  # per-gene UniProt search cap (~3 min at 0.15 s/call)
# W29/CLIB89 = UP000182444 (Other proteome); CLIB122 reference = UP000001300
_PROTEOME_IDS        = ("UP000182444", "UP000001300")
# NCBI E-utilities — Tier A-prime (bulk NCBI Gene, covers all YALI1 locus tags)
_NCBI_ESEARCH_URL    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_EFETCH_URL     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_NCBI_EFETCH_BATCH   = 400   # IDs per efetch request (NCBI recommends ≤500)


def _normalise_locus_tag(raw: str) -> set[str]:
    """
    Return all four candidate forms for a Y. lipolytica locus tag so that
    matching works regardless of strain (YALI0/YALI1) or underscore convention.

    Examples
    --------
    "YALI1C08548g"  →  {"YALI1C08548g", "YALI1_C08548g",
                         "YALI0C08548g", "YALI0_C08548g"}
    "YALI0_C08548g" →  same set
    """
    # strip optional underscore after prefix, normalise to bare form
    m = re.match(r"(YALI[01])_?([A-Za-z]\d+g)$", raw, re.IGNORECASE)
    if not m:
        return {raw}
    suffix = m.group(2)
    return {
        f"YALI0{suffix}",
        f"YALI0_{suffix}",
        f"YALI1{suffix}",
        f"YALI1_{suffix}",
    }


def _merge_gene_annotation(gene, new_data: dict) -> None:
    """
    Merge new_data into gene.annotation without overwriting existing values.
    List values are extended and de-duplicated; scalar values are kept as lists.
    """
    merged = dict(gene.annotation)
    for key, val in new_data.items():
        new_vals = val if isinstance(val, list) else [val]
        existing = merged.get(key)
        if existing is None:
            merged[key] = new_vals
        else:
            existing_list = existing if isinstance(existing, list) else [existing]
            merged[key] = list(dict.fromkeys(existing_list + new_vals))
    gene.annotation = merged


def _parse_uniprot_entry(entry: dict) -> tuple[set[str], dict]:
    """
    Parse one UniProt JSON entry into (locus_tag_candidates, annotation_dict).
    locus_tag_candidates: all normalised forms of every OLN in this entry.
    annotation_dict: ready to pass to _merge_gene_annotation().
    """
    ann: dict = {}

    acc = entry.get("primaryAccession", "")
    if acc:
        ann["uniprot"] = [acc]

    # KEGG cross-ref
    kegg_refs = entry.get("uniProtKBCrossReferences", [])
    for xref in kegg_refs:
        db = xref.get("database", "")
        xid = xref.get("id", "")
        if not xid:
            continue
        if db == "KEGG":
            ann.setdefault("kegg.genes", []).append(xid)
        elif db == "RefSeq":
            ann.setdefault("refseq", []).append(xid)
        elif db == "GeneID":
            ann.setdefault("ncbigene", []).append(xid)
        elif db == "EnsemblFungi":
            ann.setdefault("ensembl", []).append(xid)

    # Collect all Ordered Locus Names (OLNs) as candidate keys
    candidates: set[str] = set()
    gene_obj = entry.get("genes", [])
    for g in gene_obj:
        for oln in g.get("orderedLocusNames", []):
            raw = oln.get("value", "")
            if raw:
                candidates |= _normalise_locus_tag(raw)
        # Also index gene synonyms (catches YALI1 aliases stored as geneName)
        for syn in g.get("synonyms", []):
            raw = syn.get("value", "")
            if raw:
                candidates |= _normalise_locus_tag(raw)
        gn = g.get("geneName", {}).get("value", "")
        if gn:
            candidates |= _normalise_locus_tag(gn)

    return candidates, ann


def _fetch_proteome(proteome_id: str) -> list[dict]:
    """
    Download all entries for a UniProt proteome via paginated /search.
    Returns full JSON entries (including uniProtKBCrossReferences).
    """
    logger.info(f"  Fetching proteome {proteome_id} from UniProt …")
    results: list[dict] = []
    url = _UNIPROT_SEARCH_URL
    params: dict = {
        "query": f"proteome:{proteome_id}",
        "format": "json",
        "size": 500,
    }
    while url:
        try:
            resp = requests.get(url, params=params, timeout=120)
            resp.raise_for_status()
        except Exception as e:
            logger.warning(f"  Proteome fetch failed for {proteome_id}: {e}")
            break
        data = resp.json()
        page = data.get("results", [])
        results.extend(page)
        # Follow Link: rel="next" header for subsequent pages
        url = None
        params = {}
        for part in resp.headers.get("Link", "").split(","):
            part = part.strip()
            if 'rel="next"' in part:
                m = re.search(r"<([^>]+)>", part)
                if m:
                    url = m.group(1)
    logger.info(f"    {len(results)} entries retrieved")
    return results


def _tier_a(gene_ids: list[str]) -> dict[str, dict]:
    """
    Tier A: bulk proteome download.

    Downloads W29/CLIB89 proteome (UP000182444) first, then CLIB122 reference
    proteome (UP000001300).  Matches entries to model gene IDs via normalised OLN lookup
    (strip underscore, case-insensitive, YALI0↔YALI1 interchangeable).

    Returns {model_gene_id: annotation_dict}
    """
    # Build a lookup: normalised_candidate → model_gene_id
    candidate_to_model: dict[str, str] = {}
    for gid in gene_ids:
        for cand in _normalise_locus_tag(gid):
            candidate_to_model[cand.lower()] = gid

    mapping: dict[str, dict] = {}

    for proteome_id in _PROTEOME_IDS:
        entries = _fetch_proteome(proteome_id)
        for i, e in enumerate(entries[:3]):
            gene_info = e.get("genes", [])
            logger.debug(f"    [debug] entry {i} genes field: {gene_info}")
        for entry in entries:
            candidates, ann = _parse_uniprot_entry(entry)
            for cand in candidates:
                model_id = candidate_to_model.get(cand.lower())
                if model_id and model_id not in mapping:
                    mapping[model_id] = ann
        logger.info(f"    Tier A running total after {proteome_id}: {len(mapping)} mapped")

    return mapping


def _tier_b(unmapped: list[str]) -> dict[str, dict]:
    """
    Tier B: per-gene UniProt search for Tier A misses.

    Queries  gene_exact:<id> AND organism_id:4952  for each unmapped gene,
    trying all four normalised locus-tag forms.  Stops after _TIER_B_LIMIT
    API calls to avoid rate-limit issues.

    Returns {model_gene_id: annotation_dict}
    """
    mapping: dict[str, dict] = {}
    calls = 0

    for gid in unmapped:
        if calls >= _TIER_B_LIMIT:
            logger.warning(
                f"Tier B: reached {_TIER_B_LIMIT}-query limit; "
                f"{len(unmapped) - len(mapping)} genes still unmapped"
            )
            break

        candidates = _normalise_locus_tag(gid)
        hit_ann: dict | None = None

        for cand in sorted(candidates):   # deterministic order
            query = f'gene_exact:"{cand}" AND taxonomy_id:4952'
            try:
                resp = requests.get(
                    _UNIPROT_SEARCH_URL,
                    params={"query": query, "format": "json", "size": 1},
                    timeout=30,
                )
                calls += 1
                resp.raise_for_status()
                results = resp.json().get("results", [])
                if results:
                    _, hit_ann = _parse_uniprot_entry(results[0])
                    break
            except Exception as e:
                logger.debug(f"Tier B query failed for {cand}: {e}")
                calls += 1

        if hit_ann:
            mapping[gid] = hit_ann

    logger.info(f"  Tier B: {len(mapping)}/{len(unmapped)} mapped ({calls} API calls)")
    return mapping


def _tier_ncbi(unmapped: list[str]) -> dict[str, dict]:
    """
    Tier A-prime: bulk NCBI Gene lookup for genes missed by UniProt proteome download.

    NCBI Gene has 8 689 complete records for Y. lipolytica W29 including every
    YALI1 locus tag, so this tier should cover the vast majority of unmapped genes.

    Strategy
    --------
    1. esearch  — one query: "YALI1[All Fields] AND 4952[Taxonomy ID]"
                  returns the full list of NCBI Gene IDs for the organism.
    2. efetch   — fetch those IDs in batches of _NCBI_EFETCH_BATCH, XML format.
    3. Parse    — extract per-gene: Locus_tag, GeneID, RefSeq protein, UniProt.
    4. Match    — normalise locus tags and look up in candidate_to_model.

    XML fields harvested
    --------------------
    - Entrezgene_locus / Gene-ref_locus-tag  → locus tag for matching
    - Entrezgene_track-info / Gene-track_geneid → ncbigene
    - Entrezgene_locus[]/Gene-commentary_products[]/Gene-commentary_accession
      where type=protein → refseq
    - Entrezgene_comments[]/Gene-commentary_refs[]/Dbtag where db=UniProtKB → uniprot

    Returns {model_gene_id: annotation_dict}
    """
    import xml.etree.ElementTree as ET

    if not unmapped:
        return {}

    # ── Build candidate lookup ────────────────────────────────────────────
    candidate_to_model: dict[str, str] = {}
    for gid in unmapped:
        for cand in _normalise_locus_tag(gid):
            candidate_to_model[cand.lower()] = gid

    # ── Step 1: esearch — get all NCBI Gene IDs for Y. lipolytica ─────────
    # Search for all YALI locus-tagged genes under taxid 4952.
    # "txid4952[Organism]" covers W29 and all Y. lipolytica strains.
    gene_ids_ncbi: list[str] = []
    try:
        resp = requests.get(
            _NCBI_ESEARCH_URL,
            params={
                "db":      "gene",
                "term":    "txid4952[Organism] AND alive[property]",
                "retmax":  "15000",
                "retmode": "json",
            },
            timeout=60,
        )
        resp.raise_for_status()
        data = resp.json()
        gene_ids_ncbi = data.get("esearchresult", {}).get("idlist", [])
        logger.info(f"  Tier A′: NCBI esearch returned {len(gene_ids_ncbi)} Gene IDs")
    except Exception as e:
        logger.warning(f"  Tier A′: NCBI esearch failed: {e}")
        return {}

    if not gene_ids_ncbi:
        logger.warning("  Tier A′: esearch returned 0 IDs — skipping")
        return {}

    # ── Step 2+3: efetch in batches, parse XML ────────────────────────────
    mapping: dict[str, dict] = {}

    for batch_start in range(0, len(gene_ids_ncbi), _NCBI_EFETCH_BATCH):
        batch = gene_ids_ncbi[batch_start: batch_start + _NCBI_EFETCH_BATCH]
        try:
            resp = requests.get(
                _NCBI_EFETCH_URL,
                params={
                    "db":      "gene",
                    "id":      ",".join(batch),
                    "retmode": "xml",
                },
                timeout=120,
            )
            resp.raise_for_status()
        except Exception as e:
            logger.warning(f"  Tier A′: efetch batch {batch_start} failed: {e}")
            continue

        try:
            root = ET.fromstring(resp.content)
        except ET.ParseError as e:
            logger.warning(f"  Tier A′: XML parse error at batch {batch_start}: {e}")
            continue

        # Each <Entrezgene> element is one gene record
        for eg in root.iter("Entrezgene"):
            ann: dict = {}

            # ── ncbigene: Gene-track/Gene-track_geneid ───────────────────
            geneid_el = eg.find(".//Gene-track/Gene-track_geneid")
            if geneid_el is not None and geneid_el.text:
                ann["ncbigene"] = [geneid_el.text.strip()]

            # ── locus tag: Gene-ref/Gene-ref_locus-tag ───────────────────
            locus_tag_el = eg.find(".//Gene-ref/Gene-ref_locus-tag")
            locus_tag = locus_tag_el.text.strip() if (
                locus_tag_el is not None and locus_tag_el.text
            ) else ""

            # Also try Gene-ref_locus (gene symbol) as fallback
            if not locus_tag:
                sym_el = eg.find(".//Gene-ref/Gene-ref_locus")
                if sym_el is not None and sym_el.text:
                    locus_tag = sym_el.text.strip()

            if not locus_tag:
                continue

            # ── refseq: protein accessions from Gene-commentary ──────────
            # Type 8 = protein in NCBI Gene XML
            for gc in eg.iter("Gene-commentary"):
                type_el = gc.find("Gene-commentary_type")
                acc_el  = gc.find("Gene-commentary_accession")
                if (
                    type_el is not None
                    and type_el.get("value") == "protein"
                    and acc_el is not None
                    and acc_el.text
                ):
                    ann.setdefault("refseq", []).append(acc_el.text.strip())

            # ── uniprot: Dbtag where db="UniProtKB" ──────────────────────
            for dbtag in eg.iter("Dbtag"):
                db_el  = dbtag.find("Dbtag_db")
                tag_el = dbtag.find("Dbtag_tag/Object-id/Object-id_str")
                if (
                    db_el is not None
                    and db_el.text
                    and db_el.text.strip().lower() in ("uniprotkb", "uniprot")
                    and tag_el is not None
                    and tag_el.text
                ):
                    ann.setdefault("uniprot", []).append(tag_el.text.strip())

            # ── Match locus tag → model gene ─────────────────────────────
            for cand in _normalise_locus_tag(locus_tag):
                model_id = candidate_to_model.get(cand.lower())
                if model_id and model_id not in mapping:
                    mapping[model_id] = ann
                    break

        logger.info(
            f"  Tier A′: batch {batch_start}–{batch_start + len(batch) - 1} "
            f"processed → {len(mapping)} mapped so far"
        )

    logger.info(f"  Tier A′: {len(mapping)}/{len(unmapped)} mapped")
    return mapping


def _build_kegg_uniprot_index() -> dict[str, str]:
    """
    Download KEGG conv/uniprot/yli and return a reverse index:
        {uniprot_accession → kegg_gene_id}  (e.g. "Q6CFR0" → "yli:YALI_C00198g")

    KEGG's Y. lipolytica genome was re-annotated (DSM 3286 / GCF_014490615.1),
    so KEGG locus tags no longer match model YALI0/YALI1 IDs.  The only reliable
    bridge is via UniProt accession, which is shared across databases.
    """
    acc_to_kegg: dict[str, str] = {}
    try:
        resp = requests.get(_KEGG_CONV_URL, timeout=60)
        resp.raise_for_status()
        for line in resp.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                kegg_gene_id = parts[0].strip()          # e.g. "yli:YALI_C00198g"
                uniprot_acc  = parts[1].split(":")[-1].strip()  # strip "up:" prefix
                acc_to_kegg[uniprot_acc] = kegg_gene_id
        logger.info(f"  KEGG conv: {len(acc_to_kegg)} UniProt→KEGG entries loaded")
    except Exception as e:
        logger.warning(f"  KEGG conv fetch failed: {e}")
    return acc_to_kegg


def _enrich_kegg_genes(model, acc_to_kegg: dict[str, str]) -> int:
    """
    For every model gene that already has a 'uniprot' annotation, look up
    its accession(s) in acc_to_kegg and append the KEGG gene ID if found.

    Returns the number of genes enriched.
    """
    enriched = 0
    for gene in model.genes:
        ann = gene.annotation
        accessions = ann.get("uniprot", [])
        if isinstance(accessions, str):
            accessions = [accessions]
        kegg_hits = [acc_to_kegg[a] for a in accessions if a in acc_to_kegg]
        if kegg_hits:
            existing = ann.get("kegg.genes", [])
            if isinstance(existing, str):
                existing = [existing]
            merged = list(dict.fromkeys(existing + kegg_hits))
            gene.annotation = {**ann, "kegg.genes": merged}
            enriched += 1
    return enriched


def annotate_genes(model) -> None:
    """
    Map YALI1* locus IDs to UniProt accessions and cross-database identifiers.

    Tier A  — Bulk UniProt proteome download (W29/CLIB89: UP000182444,
              CLIB122: UP000001300).  Matches via OLN (YALI0↔YALI1, underscore
              optional).  Coverage limited by UniProt's incomplete W29 proteome.
    Tier A′ — Bulk NCBI Gene download (txid4952, ~8 700 records).  NCBI has
              complete YALI1 locus-tag coverage; this tier is the main driver
              toward 80 %+ annotation rate.
    Tier B  — Per-gene UniProt search (gene_exact + taxonomy_id:4952) for any
              genes still unmapped after A + A′, capped at _TIER_B_LIMIT queries.
    KEGG enrichment — After all tiers, appends kegg.genes via UniProt accession
              bridge (KEGG locus tags no longer match YALI0/YALI1).

    Writes uniprot, kegg.genes, refseq, ncbigene, ensembl into gene.annotation,
    extending (not overwriting) any values already present.
    """
    gene_ids = [g.id for g in model.genes]
    logger.info(f"Annotating {len(gene_ids)} genes …")

    # ── Tier A: UniProt proteome bulk download ────────────────────────────
    logger.info("=== Gene annotation Tier A: UniProt proteome download ===")
    tier_a_map = _tier_a(gene_ids)
    logger.info(f"Tier A: {len(tier_a_map)}/{len(gene_ids)} mapped")

    unmapped_a = [gid for gid in gene_ids if gid not in tier_a_map]

    # ── Tier A′: NCBI Gene bulk download ─────────────────────────────────
    tier_ap_map: dict[str, dict] = {}
    if unmapped_a:
        logger.info(f"=== Gene annotation Tier A′: NCBI Gene bulk ({len(unmapped_a)} remaining) ===")
        tier_ap_map = _tier_ncbi(unmapped_a)

    unmapped_ap = [gid for gid in unmapped_a if gid not in tier_ap_map]

    # ── Tier B: per-gene UniProt search for residual misses ───────────────
    tier_b_map: dict[str, dict] = {}
    if unmapped_ap:
        logger.info(f"=== Gene annotation Tier B: per-gene UniProt search ({len(unmapped_ap)} remaining) ===")
        tier_b_map = _tier_b(unmapped_ap)

    # ── Apply all tiers to model ──────────────────────────────────────────
    annotated = 0
    for gene in model.genes:
        ann = (
            tier_a_map.get(gene.id)
            or tier_ap_map.get(gene.id)
            or tier_b_map.get(gene.id)
        )
        if ann:
            _merge_gene_annotation(gene, ann)
            annotated += 1

    still_unmapped = len(gene_ids) - annotated
    logger.info(
        f"Genes annotated: {annotated}/{len(gene_ids)} "
        f"(A={len(tier_a_map)}, A′={len(tier_ap_map)}, B={len(tier_b_map)}), "
        f"{still_unmapped} unmapped"
    )

    # ── KEGG cross-ref enrichment ─────────────────────────────────────────
    logger.info("=== KEGG cross-ref enrichment: appending kegg.genes via UniProt bridge ===")
    acc_to_kegg = _build_kegg_uniprot_index()
    if acc_to_kegg:
        enriched = _enrich_kegg_genes(model, acc_to_kegg)
        logger.info(f"  KEGG enrichment: kegg.genes added to {enriched} genes")


# ─────────────────────────────────────────────
# Priority 5: gap analysis via FVA
# ─────────────────────────────────────────────

def find_gaps(model) -> dict:
    """
    Identify blocked reactions and classify their metabolites as orphan/dead-end.

    A reaction is "blocked" if FVA shows max ≈ 0 AND min ≈ 0
    (i.e. it can carry no flux in any feasible steady state).

    For each blocked reaction's metabolites:
      orphan    — metabolite has NO reaction that can produce it
                  (positive stoichiometric coefficient in any reaction)
      dead_end  — metabolite has NO reaction that can consume it
                  (negative stoichiometric coefficient in any reaction)

    Returns
    -------
    dict with keys:
      "blocked_reactions"    : list of reaction IDs
      "orphan_metabolites"   : list of metabolite IDs
      "dead_end_metabolites" : list of metabolite IDs
      "fva_result"           : cobra FVA DataFrame (full, for inspection)
    """
    from cobra.flux_analysis import flux_variability_analysis

    logger.info("Running FVA to identify blocked reactions …")
    logger.info(f"  Model: {len(model.reactions)} reactions, solver=glpk")

    model.solver = "glpk"

    # Run FVA on all reactions; loopless=False for speed
    # fraction_of_optimum=0 means we don't constrain the objective —
    # we want to know which reactions are STRUCTURALLY blocked (no feasible flux
    # regardless of the objective), so we use the full feasible polytope.
    fva = flux_variability_analysis(
        model,
        reaction_list=model.reactions,
        fraction_of_optimum=0.0,
        loopless=False,
    )

    TOL = 1e-6

    # Precompute per-reaction flux reachability from FVA:
    #   can_carry_forward[rxn_id]  = fva.maximum > TOL   (net positive flux possible)
    #   can_carry_backward[rxn_id] = fva.minimum < -TOL  (net negative flux possible)
    # A reaction is "blocked" when neither direction can carry flux.
    can_forward  = fva["maximum"] >  TOL    # Series indexed by rxn_id
    can_backward = fva["minimum"] < -TOL

    blocked_mask = ~can_forward & ~can_backward
    blocked_ids  = fva.index[blocked_mask].tolist()

    logger.info(f"  Blocked reactions: {len(blocked_ids)} / {len(model.reactions)}")

    # Classify metabolites using FVA-derived reachability, not static graph edges.
    #
    # A metabolite M is a "functional orphan" (no reachable producer) if:
    #   every reaction R that could produce M (stoich coeff > 0 in forward, or < 0
    #   in reverse) has zero reachable flux in that direction.
    #
    # Concretely, for each reaction R involving M:
    #   - Forward direction (coeff > 0) contributes production  if can_forward[R]
    #   - Reverse direction (coeff < 0) contributes production  if can_backward[R]
    #     (running backward turns a consumer into a producer)
    #
    # Symmetrically for dead-end (no reachable consumer).
    #
    # We only evaluate metabolites that appear in at least one blocked reaction
    # (otherwise they are clearly reachable).

    orphan_mets: set[str] = set()
    dead_end_mets: set[str] = set()

    candidate_mets: set = set()
    for rxn_id in blocked_ids:
        candidate_mets.update(model.reactions.get_by_id(rxn_id).metabolites)

    for met in candidate_mets:
        can_produce = False
        can_consume = False
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            fwd = bool(can_forward.get(rxn.id, False))
            rev = bool(can_backward.get(rxn.id, False))
            # Forward run produces M if coeff > 0; reverse run produces M if coeff < 0
            if (coeff > 0 and fwd) or (coeff < 0 and rev):
                can_produce = True
            # Forward run consumes M if coeff < 0; reverse run consumes M if coeff > 0
            if (coeff < 0 and fwd) or (coeff > 0 and rev):
                can_consume = True
            if can_produce and can_consume:
                break   # no need to check further

        if not can_produce:
            orphan_mets.add(met.id)
        if not can_consume:
            dead_end_mets.add(met.id)

    logger.info(f"  Orphan metabolites (no reachable producer):  {len(orphan_mets)}")
    logger.info(f"  Dead-end metabolites (no reachable consumer): {len(dead_end_mets)}")

    return {
        "blocked_reactions": blocked_ids,
        "orphan_metabolites": sorted(orphan_mets),
        "dead_end_metabolites": sorted(dead_end_mets),
        "fva_result": fva,
    }


def report_gaps(gaps: dict) -> None:
    """Log a human-readable gap analysis summary."""
    fva = gaps["fva_result"]
    blocked = gaps["blocked_reactions"]
    orphans = gaps["orphan_metabolites"]
    dead_ends = gaps["dead_end_metabolites"]

    logger.info("─── Gap Analysis Report ───────────────────────────────")
    logger.info(f"  Blocked reactions : {len(blocked)}")
    logger.info(f"  Orphan metabolites (no producer)  : {len(orphans)}")
    logger.info(f"  Dead-end metabolites (no consumer): {len(dead_ends)}")

    # Show worst blocked subsystems if reaction subsystem info available
    # (iYali26 may not have subsystems — fall back to first 20 IDs)
    logger.info(f"  First 20 blocked: {blocked[:20]}")
    if orphans:
        logger.info(f"  Orphans (first 10): {orphans[:10]}")
    if dead_ends:
        logger.info(f"  Dead-ends (first 10): {dead_ends[:10]}")
    logger.info("────────────────────────────────────────────────────────")


# ─────────────────────────────────────────────
# Stoichiometric consistency: merge duplicate metabolites
# ─────────────────────────────────────────────

def audit_mis_reactions(model, mis_metabolite_ids: list[str]) -> None:
    """
    Targeted stoichiometric audit for metabolites identified by MIS analysis.

    For each metabolite in mis_metabolite_ids, prints every reaction that
    involves it along with:
      - reaction ID, name
      - per-element imbalance (C / H / N / O / P / S / charge)
      - reaction equation string (for manual inspection)

    Use this after `merge_duplicate_metabolites()` to verify that the
    duplicates were the cause of inconsistency, and to find any residual
    stoichiometric errors.
    """
    def _element_balance(rxn) -> dict[str, float]:
        """Return per-element net balance across a reaction (should be 0)."""
        balance: dict[str, float] = {}
        for met, coeff in rxn.metabolites.items():
            formula = met.formula or ""
            for elem, cnt_str in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
                cnt = int(cnt_str) if cnt_str else 1
                balance[elem] = balance.get(elem, 0.0) + coeff * cnt
        return balance

    seen_rxns: set[str] = set()
    any_found = False

    for mid in mis_metabolite_ids:
        try:
            met = model.metabolites.get_by_id(mid)
        except KeyError:
            logger.warning(f"audit_mis_reactions: metabolite '{mid}' not in model")
            continue

        for rxn in met.reactions:
            if rxn.id in seen_rxns:
                continue
            seen_rxns.add(rxn.id)
            any_found = True

            balance = _element_balance(rxn)
            imbalanced_elems = {e: v for e, v in balance.items() if abs(v) > 1e-6}

            # Charge balance
            charge_bal = sum(
                coeff * (met_.charge or 0)
                for met_, coeff in rxn.metabolites.items()
            )

            status = "OK" if (not imbalanced_elems and abs(charge_bal) < 1e-6) else "IMBALANCED"
            imb_str = ", ".join(
                f"{e}:{v:+.3g}" for e, v in sorted(imbalanced_elems.items())
            )
            if abs(charge_bal) > 1e-6:
                imb_str += f", charge:{charge_bal:+.3g}"

            logger.info(
                f"  [{status}] {rxn.id}  ({rxn.name or 'no name'})  "
                f"imbalance={imb_str or 'none'}  "
                f"eq={rxn.reaction}"
            )

    if not any_found:
        logger.info("audit_mis_reactions: no reactions found for the given metabolite IDs")


def merge_duplicate_metabolites(
    model,
    pairs: list[tuple[str, str]],
) -> None:
    """
    Merge pairs of metabolites that represent the same chemical species
    but exist as separate nodes in the model (causing stoichiometric
    inconsistency).

    For each (keep_id, drop_id) pair:
      1. In every reaction that uses drop_id, replace drop_id with keep_id,
         summing coefficients if both appear in the same reaction.
      2. Remove drop_id from the model.
      3. Copy annotations from drop to keep (don't overwrite existing keys).

    pairs : list of (keep_id, drop_id) — keep_id survives, drop_id is removed.

    After merging, run `cobra.util.check_mass_balance()` manually to verify
    the reactions involving these metabolites are now balanced.
    """
    merged_count  = 0
    skipped_count = 0

    for keep_id, drop_id in pairs:
        try:
            keep_met = model.metabolites.get_by_id(keep_id)
        except KeyError:
            logger.warning(f"merge_duplicate_metabolites: keep ID '{keep_id}' not found — skipping pair ({keep_id}, {drop_id})")
            skipped_count += 1
            continue
        try:
            drop_met = model.metabolites.get_by_id(drop_id)
        except KeyError:
            logger.warning(f"merge_duplicate_metabolites: drop ID '{drop_id}' not found — skipping pair ({keep_id}, {drop_id})")
            skipped_count += 1
            continue

        # ── Replace drop_met with keep_met in all reactions ───────────────
        for rxn in list(drop_met.reactions):
            drop_coeff = rxn.metabolites[drop_met]
            # If keep_met already appears in this reaction, sum the coefficients
            existing_coeff = rxn.metabolites.get(keep_met, 0.0)
            new_coeff = existing_coeff + drop_coeff
            # Remove drop_met from reaction
            rxn.subtract_metabolites({drop_met: drop_coeff})
            if new_coeff != 0.0:
                rxn.add_metabolites({keep_met: new_coeff - existing_coeff
                                     if existing_coeff else new_coeff})

        # ── Copy annotations (don't overwrite existing keys in keep) ──────
        for key, val in (drop_met.annotation or {}).items():
            if key not in keep_met.annotation:
                keep_met.annotation[key] = val

        # ── Remove drop_met from model ────────────────────────────────────
        model.remove_metabolites([drop_met])
        logger.info(f"  Merged {drop_id} → {keep_id}")
        merged_count += 1

    logger.info(
        f"merge_duplicate_metabolites: {merged_count} pairs merged, "
        f"{skipped_count} skipped"
    )


# ─────────────────────────────────────────────
# main
# ─────────────────────────────────────────────

def _legacy_main():
    if not STARTING_MODEL_PATH.exists():
        logger.error(f"Could not find starting model at {STARTING_MODEL_PATH}")
        return

    logger.info(f"Loading raw model: {STARTING_MODEL_PATH.name}")
    model = read_sbml_model(str(STARTING_MODEL_PATH))

    mnx_ok = MNX_DIR.exists() and (MNX_DIR / "chem_xref.tsv").exists()

    if mnx_ok:
        # Load MetaNetX tables once, reuse across functions
        chem_xref = load_chem_xref(MNX_DIR / "chem_xref.tsv")
        chem_prop_data = load_chem_prop(MNX_DIR / "chem_prop.tsv")
        reac_xref = load_reac_xref(MNX_DIR / "reac_xref.tsv")

        reac_prop_path = MNX_DIR / "reac_prop.tsv"
        reac_prop = load_reac_prop(reac_prop_path) if reac_prop_path.exists() else None
        if reac_prop is None:
            logger.warning("reac_prop.tsv not found — Strategy C (fingerprint) disabled")

        # Priority 1 + 2a
        logger.info("=== Priority 1+2a: metabolite annotation + formulas ===")
        annotate_metabolites(model, chem_xref, chem_prop_data)

        # Priority 2b
        logger.info("=== Priority 2b: H+/H2O balance ===")
        fix_proton_water_balance(model)

        # Priority 4a  (Strategy C needs metabolite annotations from 1+2a above)
        logger.info("=== Priority 4a: reaction annotation ===")
        annotate_reactions(model, reac_xref, reac_prop)
    else:
        logger.warning(
            f"MetaNetX files not found in {MNX_DIR}. "
            "Download chem_xref.tsv, chem_prop.tsv, reac_xref.tsv from "
            "https://www.metanetx.org/mnxdoc/mnxref.html"
        )

    # Stoichiometric consistency: merge known duplicate metabolite pairs
    # These pairs were identified by MIS analysis (metabolites appearing in
    # stoichiometrically inconsistent reactions with identical formulas/charges
    # but separate model IDs).  Merge before any FBA/FVA to make the network solvable.
    logger.info("=== Stoichiometric consistency: merging duplicate metabolites ===")
    DUPLICATE_PAIRS: list[tuple[str, str]] = [
        # (keep_id, drop_id) — keep_id survives; drop_id is removed
        ("m200",  "m1855"),
        ("m772",  "m2043"),
        ("m878",  "m1963"),
    ]
    merge_duplicate_metabolites(model, DUPLICATE_PAIRS)

    # Priority 2c: exchange-bound calibration (must run before any FBA/FVA)
    logger.info("=== Priority 2c: exchange bounds (minimal medium) ===")
    set_exchange_bounds(model)   # uses MINIMAL_MEDIUM_BIGG (Tier 1) + MINIMAL_MEDIUM_NAMES (Tier 2)

    # Priority 3 (always run — independent of MetaNetX)
    logger.info("=== Priority 3: biomass reaction R1372 ===")
    fix_biomass_reaction(model)

    # Priority 4b — network required, skip if offline
    logger.info("=== Priority 4b: gene annotation via UniProt ===")
    annotate_genes(model)

    # Priority 5: gap analysis
    logger.info("=== Priority 5: gap analysis (FVA) ===")
    gaps = find_gaps(model)
    report_gaps(gaps)

    logger.info(f"Saving updated model to: {OUTPUT_MODEL_PATH.name}")
    write_sbml_model(model, str(OUTPUT_MODEL_PATH))
    logger.info("Model build complete.")


def main():
    """Delegate the retired monolith to the canonical, evidence-gated pipeline.

    The historical implementation above predates the explicit YALI0/YALI1
    crosswalk and must not be used for a formal rebuild.  Keeping this wrapper
    preserves the old command while ensuring it executes the maintained
    package implementation instead of the unsafe same-suffix locus matching.
    """
    logger.warning(
        "scripts/update_model.py is a compatibility entry point; delegating "
        "to scripts.gem_annotate.main"
    )
    from scripts.gem_annotate.main import main as canonical_main

    return canonical_main()


if __name__ == "__main__":
    main()
