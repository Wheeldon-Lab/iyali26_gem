"""
config.py — all constants and path definitions for the gem_annotate pipeline.
"""

from pathlib import Path

REPO_ROOT            = Path(__file__).resolve().parent.parent.parent
STARTING_MODEL_PATH  = REPO_ROOT / "data" / "iyli21.xml"
OUTPUT_MODEL_PATH    = REPO_ROOT / "model.xml"
MNX_DIR              = REPO_ROOT / "data" / "metanetx"
CACHE_DIR            = REPO_ROOT / "data" / "cache"

# ── UniProt / NCBI / KEGG URL constants ──────────────────────────────────────
_UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
_KEGG_CONV_URL      = "https://rest.kegg.jp/conv/uniprot/yli"
_TIER_B_LIMIT       = 1100   # per-gene UniProt search cap (~3 min at 0.15 s/call)
# W29/CLIB89 = UP000182444 (Other proteome); CLIB122 reference = UP000001300
_PROTEOME_IDS       = ("UP000182444", "UP000001300")
# NCBI E-utilities — Tier A-prime (bulk NCBI Gene, covers all YALI1 locus tags)
_NCBI_ESEARCH_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_EFETCH_URL    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_NCBI_EFETCH_BATCH  = 400    # IDs per efetch request (NCBI recommends ≤500)

# ── Minimal aerobic medium for Y. lipolytica W29 on glucose ──────────────────
# Primary key: BiGG metabolite ID (no compartment suffix).
# Fallback key: lowercased chemical name (before trailing "_FORMULA" in iYli21 names).
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
# Lowercased chemical name extracted from iYli21 "name_FORMULA" convention.
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
