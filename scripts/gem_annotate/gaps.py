"""
gaps.py — gap analysis, MIS audit, metabolite merging, and gap-fill reaction insertion.
"""

import json
import logging
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from .gap_fill_direction import load_gap_fill_direction_curation
from .locus_resolver import (
    LocusCrosswalk,
    canonical_locus_key,
    load_default_locus_crosswalk,
    model_gene_fingerprint,
)

logger = logging.getLogger(__name__)

# ── BiGG reaction-ID suffix → primary model compartment ──────────────────────
# Derived from the existing annotated reactions in iYali26.
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

# Some P0 candidates are already represented by a curated reaction under a
# different historical ID.  These mappings are intentionally explicit: a
# generic same-EC or same-GPR heuristic could suppress legitimate paralogous
# reactions.  The target reactions are reviewed iYali26 representations.
CURATED_EXISTING_GAP_FILL_REACTIONS: dict[tuple[str, str], str] = {
    ("SPHPL", "MNXR188844"): "R730",
    ("R_PSPHPL", "MNXR146152"): "R663",
}


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


def _resolve_met_id(model, bare_id: str) -> "list":
    """
    Resolve a bare metabolite ID (e.g. "m200") to all matching COBRApy
    Metabolite objects, covering three ID forms COBRApy may produce after an
    SBML round-trip:

      1. Exact match (bare_id as-is).
      2. bare_id + "[compartment]"  — bracket-compartment suffix.
      3. "M_" + bare_id  or  "M_" + bare_id + "[compartment]".

    Returns a list of all matches (may span multiple compartments).
    Returns an empty list if nothing matches.
    """
    hits: list = []

    # 1. Exact
    try:
        hits.append(model.metabolites.get_by_id(bare_id))
        return hits   # exact is unambiguous
    except KeyError:
        pass

    # 2 & 3. Prefix / M_-prefix scans
    bracket_prefixes = (bare_id + "[", "M_" + bare_id + "[")
    exact_prefixes   = ("M_" + bare_id,)

    for met in model.metabolites:
        mid = met.id
        if any(mid.startswith(p) for p in bracket_prefixes):
            hits.append(met)
        elif any(mid == p for p in exact_prefixes):
            hits.append(met)

    return hits


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

    pairs : list of (keep_id, drop_id) bare IDs — compartment suffixes are
    resolved automatically.  Merges only when both metabolites are in the same
    compartment (guards against cross-compartment collisions).

    After merging, run `cobra.util.check_mass_balance()` manually to verify
    the reactions involving these metabolites are now balanced.
    """
    # Build a quick bare-ID → [full_id, ...] map for debug logging
    bare_to_full: dict[str, list[str]] = {}
    for met in model.metabolites:
        mid = met.id
        # Strip M_ prefix
        core = mid[2:] if mid.startswith("M_") else mid
        # Strip [compartment] suffix
        bracket = core.find("[")
        bare = core[:bracket] if bracket != -1 else core
        bare_to_full.setdefault(bare, []).append(mid)

    sample = list(bare_to_full.items())[:10]
    logger.info(f"merge_duplicate_metabolites: bare-ID map sample (first 10): {sample}")

    merged_count  = 0
    skipped_count = 0

    for keep_id, drop_id in pairs:
        keep_hits = _resolve_met_id(model, keep_id)
        drop_hits = _resolve_met_id(model, drop_id)

        if not keep_hits:
            logger.warning(
                f"merge_duplicate_metabolites: keep ID '{keep_id}' not found "
                f"(known variants: {bare_to_full.get(keep_id, [])}) — "
                f"skipping pair ({keep_id}, {drop_id})"
            )
            skipped_count += 1
            continue
        if not drop_hits:
            logger.warning(
                f"merge_duplicate_metabolites: drop ID '{drop_id}' not found "
                f"(known variants: {bare_to_full.get(drop_id, [])}) — "
                f"skipping pair ({keep_id}, {drop_id})"
            )
            skipped_count += 1
            continue

        # Match keep/drop candidates that share the same compartment.
        # Build {compartment: met} maps for fast intersection.
        keep_by_comp = {m.compartment: m for m in keep_hits}
        drop_by_comp = {m.compartment: m for m in drop_hits}
        shared_comps = sorted(set(keep_by_comp) & set(drop_by_comp))

        if not shared_comps:
            logger.warning(
                f"merge_duplicate_metabolites: no shared compartment for "
                f"keep {[m.id for m in keep_hits]} / drop {[m.id for m in drop_hits]} — "
                f"skipping pair ({keep_id}, {drop_id})"
            )
            skipped_count += 1
            continue

        for comp in shared_comps:
            keep_met = keep_by_comp[comp]
            drop_met = drop_by_comp[comp]

            logger.info(
                f"  Merging {drop_met.id} → {keep_met.id} "
                f"(compartment {comp})"
            )

            # ── Replace drop_met with keep_met in all reactions ───────────
            for rxn in list(drop_met.reactions):
                drop_coeff = rxn.metabolites[drop_met]
                rxn.subtract_metabolites({drop_met: drop_coeff})
                if drop_coeff != 0.0:
                    rxn.add_metabolites({keep_met: drop_coeff})

            # ── Copy annotations (don't overwrite existing keys in keep) ──
            for key, val in (drop_met.annotation or {}).items():
                if key not in keep_met.annotation:
                    keep_met.annotation[key] = val

            # ── Remove drop_met from model ────────────────────────────────
            model.remove_metabolites([drop_met])
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


def _load_gene_cache(
    cache_path: Path, model, resolver: LocusCrosswalk | None = None
) -> dict[str, str]:
    """
    Build and cache a safe locus tag → exact model gene ID lookup.

    Keys are canonical lowercase spellings.  Same-assembly case/underscore
    variants collapse to one key; cross-assembly keys are added only through
    the explicit one-to-one crosswalk.  Legacy unversioned caches are ignored.
    """
    resolver = resolver or load_default_locus_crosswalk()
    gene_ids = [gene.id for gene in model.genes]
    expected_meta = {
        "schema": "safe-locus-cache-v1",
        "crosswalk_fingerprint": resolver.fingerprint,
        "model_gene_fingerprint": model_gene_fingerprint(gene_ids),
    }

    if cache_path.exists():
        try:
            with cache_path.open(encoding="utf-8") as handle:
                payload = json.load(handle)
            if payload.get("_meta") == expected_meta and isinstance(
                payload.get("lookup"), dict
            ):
                return {
                    canonical_locus_key(key): value
                    for key, value in payload["lookup"].items()
                }
            logger.info("Ignoring stale or legacy gene locus cache: %s", cache_path)
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            logger.warning("Ignoring unreadable gene locus cache: %s", cache_path)

    logger.info("Building gene locus-tag → model gene ID cache …")
    lookup = resolver.build_lookup(gene_ids)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", encoding="utf-8") as handle:
        json.dump(
            {"_meta": expected_meta, "lookup": lookup},
            handle,
            indent=2,
            sort_keys=True,
        )
    logger.info(f"  Cached {len(lookup):,} gene tag entries to {cache_path.name}")
    return lookup


def add_gap_fill_reactions(
    model,
    csv_path: str | Path,
    mnx_dir: Path | None = None,
    cache_dir: Path | None = None,
    direction_curation_path: str | Path | None = None,
) -> dict:
    """
    Read gap_fill_prioritized.csv and insert P0 reactions that are absent
    from the model.

    Parameters
    ----------
    model     : cobra.Model (modified in-place)
    csv_path  : path to gap_fill_prioritized.csv
    mnx_dir   : directory containing chem_xref.tsv (for MNXM→BiGG mapping)
    cache_dir : directory for persisted JSON caches (default: data/cache)
    direction_curation_path : optional durable table of evidence-backed equation
                              orientation and flux bounds

    Returns
    -------
    dict with keys: added (list of reaction IDs), skipped_duplicate,
                    skipped_curated_existing, skipped_missing_mets, skipped_unresolved_genes,
                    direction_curated, uncurated_direction, imbalanced
    """
    from cobra import Reaction

    csv_path = Path(csv_path)
    if cache_dir is None:
        cache_dir = csv_path.parent / "cache"

    mnx_xref_path = (mnx_dir / "chem_xref.tsv") if mnx_dir else None

    # ── Load / build caches ───────────────────────────────────────────────────
    _load_mnxm_cache(
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

    direction_rows = {}
    if direction_curation_path is not None:
        direction_rows = load_gap_fill_direction_curation(direction_curation_path)
        p0_pairs = {
            (row.get("bigg_reaction", "").strip(), row.get("mnxr_id", "").strip())
            for _, row in p0.iterrows()
        }
        stale_rows = sorted(
            row.bigg_reaction
            for row in direction_rows.values()
            if (row.bigg_reaction, row.mnxr_id) not in p0_pairs
        )
        if stale_rows:
            raise ValueError(
                "active gap-fill direction curation does not match a P0 "
                "reaction/MNXR pair: " + ", ".join(stale_rows)
            )

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
        "skipped_curated_existing": [],
        "skipped_missing_mets": [],
        "skipped_unresolved_genes": [],
        "direction_curated": [],
        "uncurated_direction": [],
        "imbalanced": [],
    }

    # ── Pre-build model reaction ID and BiGG annotation sets (updated incrementally) ─
    existing_ids: set[str] = {r.id for r in model.reactions}
    bigg_annotated: set[str] = set()
    for r in model.reactions:
        raw = (r.annotation or {}).get("bigg.reaction", [])
        if isinstance(raw, list):
            bigg_annotated.update(raw)
        elif raw:
            bigg_annotated.add(raw)

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

        # Some external candidates are known aliases of existing model
        # reactions even when a later annotation pass changes their BiGG/MNXR
        # labels.  Do not reinsert a duplicate merely because those labels
        # drifted; fail closed if the curated target is unexpectedly absent.
        curated_existing_id = CURATED_EXISTING_GAP_FILL_REACTIONS.get(
            (bigg_id, mnxr)
        )
        if curated_existing_id is not None:
            if curated_existing_id not in existing_ids:
                raise ValueError(
                    "curated existing-reaction target is absent: "
                    f"{bigg_id} ({mnxr}) -> {curated_existing_id}"
                )
            logger.info(
                "  Skipping %s (%s): represented by curated existing reaction %s",
                bigg_id,
                mnxr,
                curated_existing_id,
            )
            stats["skipped_curated_existing"].append(bigg_id)
            continue

        # Skip reactions already in the model (check both raw ID and R_-prefixed)
        bigg_ids_to_check = [bigg_id]
        if not bigg_id.startswith("R_"):
            bigg_ids_to_check.append("R_" + bigg_id)
        id_in_model   = any(bid in existing_ids    for bid in bigg_ids_to_check)
        bigg_in_model = any(bid in bigg_annotated  for bid in bigg_ids_to_check)
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
        unresolved_genes: list[str] = []
        for g in raw_genes:
            canon = gene_cache.get(canonical_locus_key(g))
            if canon:
                if canon not in resolved_genes:
                    resolved_genes.append(canon)
            else:
                unresolved_genes.append(g)
        if unresolved_genes:
            logger.warning(
                "  %s (%s): refusing gap-fill GPR because genes cannot be "
                "resolved exactly or through the explicit crosswalk: %s",
                bigg_id,
                mnxr,
                unresolved_genes,
            )
            stats["skipped_unresolved_genes"].append(bigg_id)
            continue
        gpr = " or ".join(resolved_genes) if resolved_genes else ""

        # ── Apply evidence-backed direction curation ─────────────────────────
        direction_row = direction_rows.get(bigg_id)
        direction_is_active = (
            direction_row is not None and direction_row.status == "active"
        )
        if direction_is_active:
            if direction_row.stoichiometry_action == "reverse":
                stoich = {met: -coefficient for met, coefficient in stoich.items()}
            lower_bound = direction_row.lower_bound
            upper_bound = direction_row.upper_bound
            stats["direction_curated"].append(bigg_id)
        else:
            # Backward-compatible only for legacy P0 reactions whose direction
            # has not yet been reviewed.  Keep these visible in the build audit
            # rather than silently treating reversibility as evidence.
            lower_bound = -1000.0
            upper_bound = 1000.0
            stats["uncurated_direction"].append(bigg_id)
            logger.warning(
                "  %s (%s): no active direction curation; retaining legacy "
                "reversible bounds",
                bigg_id,
                mnxr,
            )

        # ── Create reaction ───────────────────────────────────────────────────
        rxn = Reaction(bigg_id)
        rxn.name = bigg_id
        rxn.lower_bound = lower_bound
        rxn.upper_bound = upper_bound
        rxn.add_metabolites(stoich)
        if gpr:
            rxn.gene_reaction_rule = gpr
        if direction_is_active:
            rxn.notes = {
                "gap_fill_direction_status": "active",
                "gap_fill_stoichiometry_action": direction_row.stoichiometry_action,
                "gap_fill_direction_evidence": direction_row.evidence_url,
                "gap_fill_direction_rationale": direction_row.rationale,
            }
        else:
            rxn.notes = {"gap_fill_direction_status": "legacy_unreviewed"}

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
        existing_ids.add(bigg_id)
        bigg_annotated.add(bigg_id)
        stats["added"].append(bigg_id)
        logger.info(
            "  Added %s (%s)  GPR='%s'  mets=%d bounds=(%g, %g)",
            bigg_id,
            mnxr,
            gpr,
            len(stoich),
            lower_bound,
            upper_bound,
        )

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
        f"skipped_curated_existing={len(stats['skipped_curated_existing'])}, "
        f"skipped_missing_mets={len(stats['skipped_missing_mets'])}, "
        f"skipped_unresolved_genes={len(stats['skipped_unresolved_genes'])}, "
        f"direction_curated={len(stats['direction_curated'])}, "
        f"uncurated_direction={len(stats['uncurated_direction'])}, "
        f"imbalanced={len(stats['imbalanced'])}"
    )
    return stats
