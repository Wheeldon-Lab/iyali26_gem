"""
main.py — orchestration entry point for the iYli21 annotation pipeline.
"""

import logging

from cobra.io import read_sbml_model, write_sbml_model

from .biomass import fix_biomass_reaction
from .config import CACHE_DIR, MNX_DIR, OUTPUT_MODEL_PATH, REPO_ROOT, STARTING_MODEL_PATH
from .exchange import configure_medium, set_exchange_bounds
from .gaps import DUPLICATE_PAIRS, add_gap_fill_reactions, find_gaps, merge_duplicate_metabolites, report_gaps
from .annotate_reactions_extended import annotate_remaining_reactions
from .ec_annotation import enrich_genes_with_ec
from .genes import annotate_genes
from .idmapping import _enrich_via_idmapping
from .io import load_chem_prop, load_chem_xref, load_mnxm_depr, load_reac_prop, load_reac_xref
from .metabolites import annotate_metabolites, fix_proton_water_balance, normalize_all_annotations
from .reactions import annotate_reactions

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main():
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

        mnxm_depr = load_mnxm_depr(MNX_DIR / "chem_depr.tsv")
        if not mnxm_depr:
            logger.info(
                "  chem_depr.tsv not found in data/metanetx/ — "
                "download from https://www.metanetx.org/mnxdoc/mnxref.html "
                "to improve fingerprint match rates for deprecated MNXM IDs"
            )

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
    merge_duplicate_metabolites(model, DUPLICATE_PAIRS)

    # Priority 2c: exchange-bound calibration (must run before any FBA/FVA)
    logger.info("=== Priority 2c: exchange bounds (minimal medium) ===")
    set_exchange_bounds(model)   # uses MINIMAL_MEDIUM_BIGG (Tier 1) + MINIMAL_MEDIUM_NAMES (Tier 2)

    logger.info("=== Priority 2c+: mineral salts + vitamins medium extension ===")
    configure_medium(model)      # adds Mg, K, Na, biotin, thiamine, pyridoxine

    # Priority 3 (always run — independent of MetaNetX)
    logger.info("=== Priority 3: biomass reaction R1372 ===")
    fix_biomass_reaction(model)

    # Priority 4b — network required, skip if offline
    logger.info("=== Priority 4b: gene annotation via UniProt ===")
    annotate_genes(model)

    # Priority 4c — ncbigene → UniProt ID-mapping for genes still missing uniprot
    logger.info("=== Priority 4c: UniProt ID-mapping (ncbigene → UniProtKB) ===")
    _enrich_via_idmapping(model)

    # Priority 4d — enrich genes with EC numbers via UniProt stream API
    logger.info("=== Priority 4d: gene EC number enrichment via UniProt ===")
    enrich_genes_with_ec(model)

    # Priority 4e — extended reaction annotation (exchange / transport / EC→MNXR)
    # Runs after 4d so gene EC numbers are already populated.
    if mnx_ok:
        logger.info("=== Priority 4e: extended reaction annotation ===")
        annotate_remaining_reactions(model, reac_xref, reac_prop, mnxm_depr=mnxm_depr)

    # === EC backfill: copy gene EC numbers to reaction annotations ===
    logger.info("=== EC backfill: gene ec-code → reaction annotation ===")
    ec_backfill_count = 0
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        if "ec-code" in ann:
            continue
        ec_set = set()
        for gene in rxn.genes:
            g_ann = gene.annotation if isinstance(gene.annotation, dict) else {}
            ec_raw = g_ann.get("ec-code", [])
            if isinstance(ec_raw, str):
                ec_raw = [ec_raw]
            for ec in ec_raw:
                ec = ec.strip()
                if ec:
                    ec_set.add(ec)
        if ec_set:
            if not isinstance(rxn.annotation, dict):
                rxn.annotation = {}
            rxn.annotation["ec-code"] = sorted(ec_set)
            ec_backfill_count += 1
    logger.info(f"  EC backfill: ec-code added to {ec_backfill_count} reactions")

    # Priority 5: gap analysis — FVA before gap-fill
    logger.info("=== Priority 5: gap analysis (FVA, post-medium) ===")
    gaps = find_gaps(model)
    report_gaps(gaps)
    blocked_before_medium = len(gaps["blocked_reactions"])
    logger.info(f"  Blocked reactions after medium extension: {blocked_before_medium}")

    # Priority 6: gap-fill — insert P0 reactions from gap_fill_prioritized.csv
    gap_fill_csv = REPO_ROOT / "data" / "gap_fill_prioritized.csv"
    if gap_fill_csv.exists():
        logger.info("=== Priority 6: gap-fill reaction insertion (P0) ===")
        add_gap_fill_reactions(
            model,
            csv_path=gap_fill_csv,
            mnx_dir=MNX_DIR if mnx_ok else None,
            cache_dir=CACHE_DIR,
        )
        logger.info("=== Priority 6b: post-gap-fill FVA ===")
        gaps_after = find_gaps(model)
        before = len(gaps["blocked_reactions"])
        after  = len(gaps_after["blocked_reactions"])
        logger.info(
            f"  Blocked reactions: {before} → {after}  "
            f"(Δ {before - after:+d} unblocked)"
        )
        report_gaps(gaps_after)
    else:
        logger.warning(f"gap_fill_prioritized.csv not found at {gap_fill_csv} — skipping")

    # Migrate c_va (lowercase typo) → C_va before compartment naming
    cva_migrated = 0
    for met in model.metabolites:
        if met.compartment == "c_va":
            met.compartment = "C_va"
            cva_migrated += 1
    if cva_migrated:
        logger.info(f"Compartment fix: migrated {cva_migrated} metabolite(s) from c_va → C_va")

    # Fix compartment names so COBRApy can identify the external compartment
    _COMPARTMENT_NAMES = {
        "C_cy": "cytoplasm",
        "C_ex": "extracellular",
        "C_mi": "mitochondria",
        "C_nu": "nucleus",
        "C_er": "endoplasmic reticulum",
        "C_go": "Golgi apparatus",
        "C_va": "vacuole",
        "C_lp": "lipid particle",
        "C_pe": "peroxisome",
        "C_em": "endosomal membrane",
        "C_en": "endosome",
        "C_gm": "Golgi membrane",
        "C_mm": "mitochondrial membrane",
        "C_vm": "vacuolar membrane",
    }
    existing_comps = set(met.compartment for met in model.metabolites)
    model.compartments = {c: _COMPARTMENT_NAMES.get(c, c) for c in existing_comps}
    logger.info(f"Compartments set: {model.compartments}")

    # Verify COBRApy can resolve the external compartment via the "extracellular" name.
    # COBRApy checks compartment names for "extra" substring to identify the boundary.
    _ex_name = model.compartments.get("C_ex", "")
    if _ex_name != "extracellular":
        logger.warning(
            f"C_ex compartment name is {_ex_name!r}, expected 'extracellular' — "
            "COBRApy may not auto-detect the external compartment"
        )
    else:
        logger.info("  C_ex compartment verified as 'extracellular'")

    # === SBO annotation for pseudo-reactions ===
    logger.info("=== SBO annotation: pseudo-reactions and boundary reactions ===")
    _SBO_BIOMASS     = ["SBO:0000629"]
    _SBO_MAINTENANCE = ["SBO:0000630"]
    _SBO_ENCAPSULATE = ["SBO:0000395"]
    _SBO_EXCHANGE    = ["SBO:0000627"]
    _SBO_DEMAND      = ["SBO:0000628"]

    sbo_counts = {"biomass": 0, "maintenance": 0, "pool": 0, "exchange": 0, "demand": 0}

    _exchanges = set(model.exchanges)
    _demands   = set(model.demands)

    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        if "sbo" in ann:
            continue
        if not isinstance(rxn.annotation, dict):
            rxn.annotation = {}

        rid = rxn.id
        if rid.startswith(("xBIOMASS", "newBiom", "biomass_C")):
            rxn.annotation["sbo"] = _SBO_BIOMASS
            sbo_counts["biomass"] += 1
        elif rid.startswith("xMAINTENANCE"):
            rxn.annotation["sbo"] = _SBO_MAINTENANCE
            sbo_counts["maintenance"] += 1
        elif rid.startswith(("xLIPID", "xAMINOACID", "xPOOL_")):
            rxn.annotation["sbo"] = _SBO_ENCAPSULATE
            sbo_counts["pool"] += 1
        elif rxn in _exchanges:
            rxn.annotation["sbo"] = _SBO_EXCHANGE
            sbo_counts["exchange"] += 1
        elif rxn in _demands:
            rxn.annotation["sbo"] = _SBO_DEMAND
            sbo_counts["demand"] += 1

    logger.info(
        f"  SBO added: biomass={sbo_counts['biomass']}  "
        f"maintenance={sbo_counts['maintenance']}  "
        f"pool/encapsulating={sbo_counts['pool']}  "
        f"exchange={sbo_counts['exchange']}  demand={sbo_counts['demand']}"
    )

    ec_check = sum(1 for r in model.reactions if isinstance(r.annotation, dict) and 'ec-code' in r.annotation)
    logger.info(f"  DEBUG: reactions with ec-code BEFORE normalize: {ec_check}")

    logger.info("=== Final pass: normalise all annotation keys/values ===")
    normalize_all_annotations(model)

    ec_check2 = sum(1 for r in model.reactions if isinstance(r.annotation, dict) and 'ec-code' in r.annotation)
    logger.info(f"  DEBUG: reactions with ec-code AFTER normalize: {ec_check2}")

    logger.info(f"Saving updated model to: {OUTPUT_MODEL_PATH.name}")
    write_sbml_model(model, str(OUTPUT_MODEL_PATH))
    logger.info("Model build complete.")


if __name__ == "__main__":
    main()
