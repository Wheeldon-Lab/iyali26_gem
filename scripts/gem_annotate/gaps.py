"""
gaps.py — gap analysis, MIS audit, metabolite merging, and gap-fill reaction insertion.
"""

import json
import logging
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# ── BiGG reaction-ID suffix → primary model compartment ──────────────────────
# Derived from the existing annotated reactions in iYli21.
# Used to choose a default compartment when an equation only has MNXD1.
_BIGG_SUFFIX_TO_COMPARTMENT: dict[str, str] = {
    "m":  "C_mi",   # mitochondria
    "mi": "C_mi",
    "mr": "C_mi",
    "pm": "C_mi",
    "x":  "C_pe",   # peroxisome
    "p":  "C_pe",
    "pp": "C_cy",   # cytoplasm (periplasm in bacteria → cytoplasm in yeast)
    "n":  "C_nu",   # nucleus
    "g":  "C_go",   # Golgi
    "e":  "C_en",   # endosome
    "er": "C_er",   # ER
    "r":  "C_er",
    "vm": "C_vm",   # vacuolar membrane
    "va": "C_va",   # vacuole
    "lp": "C_lp",   # lipid particle
    "c":  "C_cy",   # cytoplasm (default)
    "cy": "C_cy",
    "i":  "C_cy",
}
_DEFAULT_COMPARTMENT = "C_cy"

# Known duplicate metabolite pairs identified by MIS analysis.
# (keep_id, drop_id) — keep_id survives; drop_id is removed.
DUPLICATE_PAIRS: list[tuple[str, str]] = [
    ("m200",  "m1855"),
    ("m772",  "m2043"),
    ("m878",  "m1963"),
]


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
    #   every reaction R that could produce M (positive stoichiometric coefficient in any reaction)
    #   has zero reachable flux in that direction.
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
    # (iYli21 may not have subsystems — fall back to first 20 IDs)
    logger.info(f"  First 20 blocked: {blocked[:20]}")
    if orphans:
        logger.info(f"  Orphans (first 10): {orphans[:10]}")
    if dead_ends:
        logger.info(f"  Dead-ends (first 10): {dead_ends[:10]}")
    logger.info("────────────────────────────────────────────────────────")


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


# ── Gap-fill reaction insertion ───────────────────────────────────────────────

_EQ_TOKEN = re.compile(r"(\d+(?:\.\d+)?)\s+([\w]+)@(MNXD\d+)")


def _infer_compartment(bigg_id: str) -> str:
    """
    Guess the primary model compartment from a BiGG reaction ID suffix.
    Returns a C_* compartment string, defaulting to C_cy.
    """
    clean = bigg_id.upper()
    if clean.startswith("R_"):
        clean = clean[2:]
    # Walk suffixes from longest to shortest (up to 3 chars)
    for length in (3, 2, 1):
        suffix = clean[-length:].lower()
        if suffix in _BIGG_SUFFIX_TO_COMPARTMENT:
            return _BIGG_SUFFIX_TO_COMPARTMENT[suffix]
    return _DEFAULT_COMPARTMENT


def _build_mnxm_to_model_met(model) -> dict[tuple[str, str], object]:
    """
    Return a dict mapping (mnxm_id, compartment) → cobra Metabolite.
    Also includes a compartment-free fallback key (mnxm_id, None).
    For WATER and MNXM1 (H+) special tokens, uses formula-based lookup.
    """
    mnxm_comp: dict[tuple[str, str], object] = {}
    mnxm_any: dict[str, object] = {}     # first hit per MNXM regardless of compartment

    for met in model.metabolites:
        mnxm = met.annotation.get("metanetx.chemical")
        if not mnxm:
            continue
        if isinstance(mnxm, list):
            mnxm = mnxm[0]
        comp = met.compartment
        if (mnxm, comp) not in mnxm_comp:
            mnxm_comp[(mnxm, comp)] = met
        if mnxm not in mnxm_any:
            mnxm_any[mnxm] = met

    # Fallback: also index by (mnxm, None) pointing to the first hit
    for mnxm, met in mnxm_any.items():
        mnxm_comp[(mnxm, None)] = met

    return mnxm_comp


def _build_special_token_map(model) -> dict[str, dict[str, object]]:
    """
    Build compartment-keyed maps for the special MetaNetX tokens
    WATER (H2O, formula='H2O') and MNXM1/H+ (name contains 'H+').

    Returns {'WATER': {comp: met}, 'MNXM1': {comp: met}}.
    """
    water_by_comp: dict[str, object] = {}
    proton_by_comp: dict[str, object] = {}

    for met in model.metabolites:
        comp = met.compartment
        f = met.formula or ""
        name = met.name or ""
        if f == "H2O" and comp not in water_by_comp:
            water_by_comp[comp] = met
        if "H+" in name and comp not in proton_by_comp:
            proton_by_comp[comp] = met

    return {"WATER": water_by_comp, "MNXM1": proton_by_comp}


def _parse_equation(equation: str) -> list[tuple[float, str, str]]:
    """
    Parse a MetaNetX equation string.

    Returns list of (coefficient, mnxm_token, mnxd_compartment) tuples.
    Substrates have negative coefficients; products have positive.
    """
    left, _, right = equation.partition("=")
    result = []
    for side, sign in ((left, -1), (right, +1)):
        for m in _EQ_TOKEN.finditer(side):
            coeff = float(m.group(1)) * sign
            token = m.group(2)   # e.g. MNXM5 or WATER
            mnxd  = m.group(3)   # e.g. MNXD1 or MNXD2
            result.append((coeff, token, mnxd))
    return result


def _load_mnxm_cache(cache_path: Path, mnx_xref_path: Path | None) -> dict[str, str]:
    """
    Load or build the MNXM→bigg.metabolite cache from chem_xref.tsv.
    Returns mnxm_id → bigg_metabolite_id (bare, no compartment suffix).
    """
    if cache_path.exists():
        with open(cache_path) as fh:
            return json.load(fh)

    if mnx_xref_path is None or not mnx_xref_path.exists():
        logger.warning("chem_xref.tsv not available — MNXM→BiGG cache cannot be built")
        return {}

    logger.info("Building MNXM→BiGG metabolite cache from chem_xref.tsv …")
    df = pd.read_csv(mnx_xref_path, sep="\t", comment="#", header=None,
                     names=["source", "mnx_id", "desc"], dtype=str).fillna("")
    df = df[df["mnx_id"].str.startswith("MNXM")]

    cache: dict[str, str] = {}
    for source, mnx_id, _ in df.itertuples(index=False):
        if source.startswith("bigg.metabolite:") and mnx_id not in cache:
            bigg_id = source.split(":", 1)[1]
            cache[mnx_id] = bigg_id

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as fh:
        json.dump(cache, fh)
    logger.info(f"  Cached {len(cache):,} MNXM→BiGG entries to {cache_path.name}")
    return cache


def _load_gene_cache(cache_path: Path, model) -> dict[str, str]:
    """
    Build and cache a gene_id → model gene ID lookup so that GPR rules
    use the exact gene IDs present in the model.  The CSV uses YALI1* locus
    tags; the model may use YALI0* variants.

    Returns locus_tag (normalised) → model_gene_id.
    """
    if cache_path.exists():
        with open(cache_path) as fh:
            return json.load(fh)

    logger.info("Building gene locus-tag → model gene ID cache …")
    _TAG = re.compile(r"(YALI[01])_?([A-Za-z]\d+g)$", re.IGNORECASE)
    lookup: dict[str, str] = {}
    for gene in model.genes:
        m = _TAG.match(gene.id)
        if not m:
            lookup[gene.id.lower()] = gene.id
            continue
        suffix = m.group(2).lower()
        # Index all four variant forms
        for prefix in ("YALI0", "YALI0_", "YALI1", "YALI1_"):
            lookup[f"{prefix}{suffix}"] = gene.id

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as fh:
        json.dump(lookup, fh)
    logger.info(f"  Cached {len(lookup):,} gene tag entries to {cache_path.name}")
    return lookup


def add_gap_fill_reactions(model, csv_path: str | Path,
                           mnx_dir: Path | None = None,
                           cache_dir: Path | None = None) -> dict:
    """
    Read gap_fill_prioritized.csv and insert P0 reactions that are absent
    from the model.

    Parameters
    ----------
    model     : cobra.Model (modified in-place)
    csv_path  : path to gap_fill_prioritized.csv
    mnx_dir   : directory containing chem_xref.tsv (for MNXM→BiGG mapping)
    cache_dir : directory for persisted JSON caches (default: data/cache)

    Returns
    -------
    dict with keys: added (list of reaction IDs), skipped_duplicate,
                    skipped_missing_mets, imbalanced
    """
    from cobra import Gene, Reaction

    csv_path = Path(csv_path)
    if cache_dir is None:
        cache_dir = csv_path.parent / "cache"

    mnx_xref_path = (mnx_dir / "chem_xref.tsv") if mnx_dir else None

    # ── Load / build caches ───────────────────────────────────────────────────
    mnxm_bigg_cache = _load_mnxm_cache(
        cache_dir / "mnxm_to_bigg_metabolite.json", mnx_xref_path
    )
    gene_cache = _load_gene_cache(cache_dir / "gene_locus_tag_map.json", model)

    # ── Build MNXM → model metabolite lookup (per compartment) ───────────────
    mnxm_comp_to_met = _build_mnxm_to_model_met(model)
    special_tokens   = _build_special_token_map(model)

    # ── Load the CSV ─────────────────────────────────────────────────────────
    df = pd.read_csv(csv_path, dtype=str).fillna("")
    p0 = df[df["priority"] == "P0"].copy()
    logger.info(f"add_gap_fill_reactions: {len(p0)} P0 rows in {csv_path.name}")

    # Collect all genes per MNXR (for isozyme OR-GPR)
    mnxr_genes: dict[str, list[str]] = defaultdict(list)
    mnxr_rows: dict[str, dict] = {}
    for _, row in p0.iterrows():
        mnxr = row.get("mnxr_id", "")
        bigg = row.get("bigg_reaction", "").strip()
        if not mnxr or not bigg:
            continue
        gene_id = row.get("gene_id", "").strip()
        if gene_id and gene_id not in mnxr_genes[mnxr]:
            mnxr_genes[mnxr].append(gene_id)
        if mnxr not in mnxr_rows:
            mnxr_rows[mnxr] = row.to_dict()

    stats = {
        "added": [],
        "skipped_duplicate": [],
        "skipped_missing_mets": [],
        "imbalanced": [],
    }

    # ── Process each unique (bigg_reaction, MNXR) pair ───────────────────────
    seen_bigg: set[str] = set()   # deduplicate when same BiGG appears with multiple MNXR

    for mnxr, row in mnxr_rows.items():
        bigg_id = row.get("bigg_reaction", "").strip()
        equation = row.get("equation", "").strip()
        ec        = row.get("ec_number", "").strip()
        kegg_rxn  = row.get("kegg_reaction", "").strip()

        if not bigg_id or not equation:
            logger.debug(f"  Skipping {mnxr}: missing bigg_reaction or equation")
            continue

        # Skip if this BiGG ID was already handled (two MNXR share same BiGG)
        if bigg_id in seen_bigg:
            continue
        seen_bigg.add(bigg_id)

        # Skip reactions already in the model (check both raw ID and R_-prefixed)
        existing_ids = {r.id for r in model.reactions}
        bigg_ids_to_check = [bigg_id]
        if not bigg_id.startswith("R_"):
            bigg_ids_to_check.append("R_" + bigg_id)
        # Also check existing bigg.reaction annotations
        bigg_in_model = any(
            any(bid in (r.annotation.get("bigg.reaction") or []) for r in model.reactions)
            for bid in bigg_ids_to_check
        )
        id_in_model = any(bid in existing_ids for bid in bigg_ids_to_check)
        if id_in_model or bigg_in_model:
            logger.debug(f"  Skipping {bigg_id} ({mnxr}): already in model")
            stats["skipped_duplicate"].append(bigg_id)
            continue

        # ── Infer compartment from BiGG suffix ───────────────────────────────
        default_comp = _infer_compartment(bigg_id)

        # ── Resolve metabolites ───────────────────────────────────────────────
        stoich: dict[object, float] = {}
        missing_tokens: list[str] = []

        parsed = _parse_equation(equation)
        if not parsed:
            logger.warning(f"  {bigg_id}: could not parse equation '{equation}'")
            continue

        # Determine compartment: if equation uses MNXD2 (transport), we map
        # MNXD1 → default_comp and MNXD2 → extracellular
        has_mnxd2 = any(mnxd == "MNXD2" for _, _, mnxd in parsed)
        mnxd_map = {
            "MNXD1": default_comp,
            "MNXD2": "C_ex" if has_mnxd2 else default_comp,
        }

        for coeff, token, mnxd in parsed:
            comp = mnxd_map.get(mnxd, default_comp)

            # Special tokens: WATER and MNXM1 (H+) indexed by formula/name
            if token in special_tokens:
                met = special_tokens[token].get(comp) or special_tokens[token].get(
                    _DEFAULT_COMPARTMENT
                )
            else:
                # Try exact (mnxm, comp), then any-compartment fallback
                met = mnxm_comp_to_met.get((token, comp)) or mnxm_comp_to_met.get(
                    (token, None)
                )

            if met is None:
                missing_tokens.append(f"{token}@{comp}")
                continue

            stoich[met] = stoich.get(met, 0.0) + coeff

        if missing_tokens:
            logger.warning(
                f"  {bigg_id} ({mnxr}): {len(missing_tokens)} metabolites not found "
                f"in model — skipping. Missing: {missing_tokens}"
            )
            stats["skipped_missing_mets"].append(bigg_id)
            continue

        # Remove any zero-net stoichiometry entries
        stoich = {met: c for met, c in stoich.items() if abs(c) > 1e-9}
        if not stoich:
            logger.warning(f"  {bigg_id}: empty stoichiometry after parsing — skipping")
            continue

        # ── Build GPR (isozymes → OR) ─────────────────────────────────────────
        raw_genes = mnxr_genes.get(mnxr, [])
        resolved_genes: list[str] = []
        for g in raw_genes:
            canon = gene_cache.get(g) or gene_cache.get(g.lower())
            if canon:
                resolved_genes.append(canon)
            else:
                # Gene not yet in model — add it
                resolved_genes.append(g)
        gpr = " or ".join(resolved_genes) if resolved_genes else ""

        # ── Create reaction ───────────────────────────────────────────────────
        rxn = Reaction(bigg_id)
        rxn.name = bigg_id
        rxn.lower_bound = -1000.0
        rxn.upper_bound =  1000.0
        rxn.add_metabolites(stoich)
        if gpr:
            rxn.gene_reaction_rule = gpr

        # ── Annotations ──────────────────────────────────────────────────────
        ann: dict[str, object] = {
            "bigg.reaction":       [bigg_id],
            "metanetx.reaction":   [mnxr],
            "sbo":                 ["SBO:0000176"],   # biochemical reaction
        }
        if kegg_rxn:
            ann["kegg.reaction"] = [kegg_rxn]
        if ec:
            ann["ec-code"] = [ec]
        rxn.annotation = ann

        model.add_reactions([rxn])
        stats["added"].append(bigg_id)
        logger.info(f"  Added {bigg_id} ({mnxr})  GPR='{gpr}'  mets={len(stoich)}")

        # ── Mass-balance check ────────────────────────────────────────────────
        imbalance = rxn.check_mass_balance()
        if imbalance:
            logger.warning(
                f"  Mass imbalance in {bigg_id}: "
                + ", ".join(f"{e}:{v:+.3g}" for e, v in sorted(imbalance.items()))
            )
            stats["imbalanced"].append(bigg_id)

    logger.info(
        f"add_gap_fill_reactions: added={len(stats['added'])}, "
        f"skipped_duplicate={len(stats['skipped_duplicate'])}, "
        f"skipped_missing_mets={len(stats['skipped_missing_mets'])}, "
        f"imbalanced={len(stats['imbalanced'])}"
    )
    return stats
