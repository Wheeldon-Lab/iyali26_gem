"""
gap_fill_prioritize.py — Layered filtering and prioritization of gap-fill candidates.

Steps:
  1. Filter promiscuous-enzyme noise (EC with >30 candidates: keep only BiGG-mapped ones)
  2. Fix MNXM → model metabolite mapping (WATER/MNXM2/MNXM1 and via chem_xref)
  3. Prioritize candidates into P0–P3 tiers
  4. Output gap_fill_prioritized.csv
"""

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
METANETX = DATA / "metanetx"

# ── 1. Load gap-fill candidates ──────────────────────────────────────────────

df = pd.read_csv(DATA / "gap_fill_candidates.csv")
print(f"Loaded {len(df)} candidates")

# ── 2. Load MetaNetX chem_xref: external_id → MNXM ──────────────────────────

print("Loading chem_xref.tsv …")
xref = pd.read_csv(
    METANETX / "chem_xref.tsv",
    sep="\t", comment="#", header=None,
    names=["xref_id", "mnx_id", "description"],
)

# Build MNXM → set of BiGG IDs (strip "bigg.metabolite:" prefix)
bigg_xref = xref[xref["xref_id"].str.startswith("bigg.metabolite:", na=False)].copy()
bigg_xref["bigg_id"] = bigg_xref["xref_id"].str.removeprefix("bigg.metabolite:")
mnxm_to_bigg: dict[str, set[str]] = (
    bigg_xref.groupby("mnx_id")["bigg_id"].apply(set).to_dict()
)

# Build MNXM → set of KEGG IDs (strip "kegg.compound:" prefix)
kegg_xref = xref[xref["xref_id"].str.startswith("kegg.compound:", na=False)].copy()
kegg_xref["kegg_id"] = kegg_xref["xref_id"].str.removeprefix("kegg.compound:")
mnxm_to_kegg: dict[str, set[str]] = (
    kegg_xref.groupby("mnx_id")["kegg_id"].apply(set).to_dict()
)

# Also build a direct MNXM alias map (some MNX IDs redirect to canonical IDs)
# xref rows where xref_id starts with "MNXM" and mnx_id differs are aliases
alias_rows = xref[
    xref["xref_id"].str.match(r"^MNXM\d+$", na=False)
    & (xref["xref_id"] != xref["mnx_id"])
].copy()
mnxm_alias: dict[str, str] = dict(zip(alias_rows["xref_id"], alias_rows["mnx_id"]))

# ── 3. Load chem_prop: MNXM → formula, charge ────────────────────────────────

print("Loading chem_prop.tsv …")
prop = pd.read_csv(
    METANETX / "chem_prop.tsv",
    sep="\t", comment="#", header=None,
    names=["mnx_id", "name", "reference", "formula", "charge", "mass", "inchi", "inchikey", "smiles"],
    dtype={"formula": str, "charge": str},
)
prop = prop.set_index("mnx_id")

def has_formula(mnxm_id: str) -> bool:
    """Return True if the MNXM entry has a non-empty formula in chem_prop."""
    mnxm_id = mnxm_alias.get(mnxm_id, mnxm_id)
    if mnxm_id not in prop.index:
        return False
    row = prop.loc[mnxm_id]
    formula = row["formula"] if isinstance(row, pd.Series) else row.iloc[0]["formula"]
    return bool(formula and str(formula).strip() not in ("", "nan"))

# ── 4. Load model metabolites ─────────────────────────────────────────────────

print("Loading model metabolites from iyli21.xml …")
tree = ET.parse(DATA / "iyli21.xml")
root = tree.getroot()
ns_sbml = "http://www.sbml.org/sbml/level3/version1/core"

model_mets = []
for model_elem in root.iter("{%s}model" % ns_sbml):
    for sp in model_elem.iter("{%s}species" % ns_sbml):
        model_mets.append({"id": sp.get("id"), "name": sp.get("name", "")})

df_model = pd.DataFrame(model_mets)

# Build a set of "base names" (strip compartment suffix and formula suffix) for fuzzy matching
# Model names have format: "compound name_FORMULA" (e.g. "H2O_H2O", "H+_p+1")
# Extract the part before the last underscore as the compound key
def base_name(full_name: str) -> str:
    """Strip trailing _FORMULA or trailing underscore from model name."""
    parts = full_name.rsplit("_", 1)
    return parts[0].strip().lower() if len(parts) > 1 else full_name.strip().lower()

df_model["base_name"] = df_model["name"].apply(base_name)
model_base_names = set(df_model["base_name"])

# Build formula → model metabolite set from chem_prop for cross-referencing
model_formulas: dict[str, set[str]] = {}
for _, row in df_model.iterrows():
    formula_part = row["name"].rsplit("_", 1)[-1].strip()
    if formula_part and formula_part not in ("nan", ""):
        model_formulas.setdefault(formula_part, set()).add(row["id"])

# ── 5. Build MNXM → "in model?" mapping using chem_xref ─────────────────────
# Strategy:
#   a. WATER → model has H2O metabolites (formula H2O)
#   b. MNXM2 → hydroxide (not water!) — keep as missing unless model has it
#   c. MNXM1 / MNXM01 → proton H+ — model does have these
#   d. For other MNXM IDs: check if any BiGG alias exists in model reactions
#      (we don't have direct BiGG ↔ model mapping, so we use formula matching)

# Build formula → canonical MNXM from prop
formula_to_mnxm: dict[str, str] = {}
for mnxm_id, row in prop.iterrows():
    formula = str(row.get("formula", "") or "").strip()
    if formula and formula != "nan" and mnxm_id not in formula_to_mnxm:
        formula_to_mnxm[formula] = mnxm_id

# Known special MNXM → model formula mappings
SPECIAL_MNXM_IN_MODEL = {
    "WATER": True,   # H2O — definitely in model
    "MNXM2": False,  # OH- (hydroxide, not water)
    "MNXM1": True,   # H+ (proton) — definitely in model
    "MNXM01": True,  # PMF/translocated proton — treat as in model (cosmetic variant)
}

def mnxm_in_model(mnxm_id: str) -> bool:
    """
    Check if a MetaNetX metabolite ID corresponds to a metabolite
    already in the iYli21 model.
    Strategy: formula matching via chem_prop, with special-case overrides.
    """
    if mnxm_id in SPECIAL_MNXM_IN_MODEL:
        return SPECIAL_MNXM_IN_MODEL[mnxm_id]

    # Resolve any alias (e.g. MNXM739538 → canonical form)
    canonical = mnxm_alias.get(mnxm_id, mnxm_id)

    if canonical not in prop.index:
        return False

    row = prop.loc[canonical]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    formula = str(row.get("formula", "") or "").strip()
    if not formula or formula == "nan":
        return False

    return formula in model_formulas

# ── 6. Parse missing_metabolites and recount after MNXM fix ──────────────────

def parse_missing(missing_str) -> list[str]:
    """Return list of missing MNXM IDs from the pipe-separated field."""
    if pd.isna(missing_str) or str(missing_str).strip() == "":
        return []
    return [m.strip() for m in str(missing_str).split("|") if m.strip()]

def recount_missing(missing_list: list[str]) -> tuple[int, list[str]]:
    """
    After checking each MNXM ID against the model, return
    (count_still_missing, list_still_missing).
    """
    still_missing = [m for m in missing_list if not mnxm_in_model(m)]
    return len(still_missing), still_missing

df["missing_list"] = df["missing_metabolites"].apply(parse_missing)
df[["missing_count_fixed", "missing_fixed_list"]] = df["missing_list"].apply(
    lambda ml: pd.Series(recount_missing(ml))
)
df["missing_fixed_str"] = df["missing_fixed_list"].apply(lambda lst: "|".join(lst))

# ── 7. Step 1 — Filter promiscuous EC numbers (>30 candidates) ───────────────

ec_counts = df["ec_number"].value_counts()
noisy_ecs = ec_counts[ec_counts > 30].index.tolist()

print(f"\nStep 1 — Promiscuous EC numbers (>30 candidates): {noisy_ecs}")
for ec in noisy_ecs:
    total = (df["ec_number"] == ec).sum()
    with_bigg = ((df["ec_number"] == ec) & df["bigg_reaction"].notna() & (df["bigg_reaction"] != "")).sum()
    print(f"  {ec}: {total} total → {with_bigg} with BiGG ID (kept)")

# For noisy ECs: keep only rows with a BiGG reaction ID
noisy_mask = df["ec_number"].isin(noisy_ecs)
bigg_mask = df["bigg_reaction"].notna() & (df["bigg_reaction"] != "")

keep_mask = ~noisy_mask | bigg_mask   # keep non-noisy OR (noisy but has BiGG)
df_filtered = df[keep_mask].copy()

print(f"  Rows after Step 1 filter: {len(df_filtered)} (removed {len(df) - len(df_filtered)})")

# Save the noise-filtered + MNXM-fixed mapping for reuse
mapping_out = df_filtered[
    ["gene_id", "ec_number", "mnxr_id", "bigg_reaction", "kegg_reaction",
     "metacyc_reaction", "in_model", "missing_metabolites",
     "missing_count_fixed", "missing_fixed_str", "equation"]
].copy()
mapping_out.to_csv(DATA / "gap_fill_filtered_mapping.csv", index=False)
print(f"  Saved filtered mapping → data/gap_fill_filtered_mapping.csv")

# ── 8. Step 3 — Assign priority ──────────────────────────────────────────────

def has_bigg(row) -> bool:
    return bool(row.get("bigg_reaction") and str(row["bigg_reaction"]).strip() not in ("", "nan"))

def has_kegg(row) -> bool:
    return bool(row.get("kegg_reaction") and str(row["kegg_reaction"]).strip() not in ("", "nan"))

def assign_priority(row) -> str:
    n_miss = row["missing_count_fixed"]
    missing_lst = row["missing_fixed_list"] if isinstance(row["missing_fixed_list"], list) else []

    # All missing metabolites have formula info in chem_prop?
    all_have_formula = all(has_formula(m) for m in missing_lst)

    if has_bigg(row) and n_miss == 0:
        return "P0"
    if has_bigg(row) and 1 <= n_miss <= 2 and all_have_formula:
        return "P1"
    if has_kegg(row) and n_miss <= 3:
        return "P2"
    return "P3"

df_filtered["priority"] = df_filtered.apply(assign_priority, axis=1)

priority_counts = df_filtered["priority"].value_counts().sort_index()
print("\nStep 3 — Priority distribution:")
for p, cnt in priority_counts.items():
    print(f"  {p}: {cnt}")

# ── 9. Step 4 — Output prioritized CSV ───────────────────────────────────────

out_cols = [
    "priority", "gene_id", "ec_number", "mnxr_id",
    "bigg_reaction", "kegg_reaction",
    "missing_count_fixed", "missing_fixed_str", "equation",
]

out = (
    df_filtered[out_cols]
    .rename(columns={
        "missing_count_fixed": "missing_metabolites_count",
        "missing_fixed_str": "missing_metabolites",
    })
    .sort_values(["priority", "gene_id", "ec_number"])
    .reset_index(drop=True)
)

out_path = DATA / "gap_fill_prioritized.csv"
out.to_csv(out_path, index=False)
print(f"\nStep 4 — Saved {len(out)} prioritized candidates → {out_path}")

# ── 10. Summary report ────────────────────────────────────────────────────────

print("\n=== Summary ===")
print(f"Original candidates : {len(df)}")
print(f"After Step 1 filter : {len(df_filtered)}")
print(f"  P0 (direct add)   : {(out['priority']=='P0').sum()}")
print(f"  P1 (high conf)    : {(out['priority']=='P1').sum()}")
print(f"  P2 (medium conf)  : {(out['priority']=='P2').sum()}")
print(f"  P3 (manual review): {(out['priority']=='P3').sum()}")

print("\nStep 2 — MNXM fix statistics:")
original_missing = df_filtered["missing_list"].apply(len)
fixed_missing = df_filtered["missing_count_fixed"]
improved = (original_missing > fixed_missing).sum()
print(f"  Rows where MNXM fix reduced missing count: {improved}")
print(f"  WATER/MNXM1 etc. resolved as in-model: see gap_fill_filtered_mapping.csv")

print("\nTop P0 candidates (first 10):")
p0 = out[out["priority"] == "P0"][["gene_id","ec_number","bigg_reaction","equation"]].head(10)
print(p0.to_string(index=False))

print("\nDone.")
