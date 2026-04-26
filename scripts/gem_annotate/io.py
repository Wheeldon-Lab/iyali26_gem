"""
io.py — MetaNetX TSV loaders.
"""

import logging
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

_MNXM_IN_EQ = re.compile(r"(MNXM\d+)@")   # extracts MNXM IDs from equation strings


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
