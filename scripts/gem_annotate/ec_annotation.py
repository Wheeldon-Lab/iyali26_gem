"""
ec_annotation.py — Enrich model genes with EC numbers via UniProt stream API.
"""

import logging
import re

import requests

from .config import CACHE_DIR
from .http_utils import _request_with_retry

logger = logging.getLogger(__name__)

_UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
_BATCH_SIZE = 50   # accessions per request — keeps query string well under URL limits


def _fetch_ec_for_accessions(accessions: list[str]) -> dict[str, list[str]]:
    """
    Query UniProt stream API for a batch of accessions.
    Returns {accession: [ec_number, ...]} for every accession that has EC data.
    """
    query = " OR ".join(f"accession:{a}" for a in accessions)
    try:
        resp = _request_with_retry(
            "GET", _UNIPROT_STREAM_URL,
            params={"query": f"({query})", "fields": "accession,ec", "format": "tsv"},
            timeout=60,
        )
        resp.raise_for_status()
    except Exception as e:
        logger.warning(f"UniProt stream request failed: {e}")
        return {}

    result: dict[str, list[str]] = {}
    for line in resp.text.strip().splitlines():
        if line.startswith("Entry"):   # header
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        acc = parts[0].strip()
        ec_raw = parts[1].strip()
        if not ec_raw:
            continue
        # EC numbers are semicolon-separated; strip trailing semicolons/spaces
        ec_list = [e.strip() for e in ec_raw.split(";") if e.strip()]
        if ec_list:
            result[acc] = ec_list

    return result


def enrich_genes_with_ec(model) -> None:
    """
    For every model gene that has a 'uniprot' annotation, batch-fetch EC numbers
    from UniProt and write them into gene.annotation["ec-code"].

    Results are cached to data/cache/gene_ec_cache.json so subsequent runs skip
    the network entirely for already-queried accessions.
    """
    import json

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / "gene_ec_cache.json"

    # {accession: [ec_number, ...] | null}  — null = queried, no EC found
    ec_cache: dict[str, list[str] | None] = {}
    if cache_file.exists():
        with cache_file.open() as f:
            ec_cache = json.load(f)

    # Collect all uniprot accessions that need EC lookup
    acc_to_genes: dict[str, list] = {}   # accession → list of gene objects
    for gene in model.genes:
        ann = gene.annotation if isinstance(gene.annotation, dict) else {}
        raw = ann.get("uniprot", [])
        accessions = raw if isinstance(raw, list) else [raw]
        for acc in accessions:
            if acc:
                acc_to_genes.setdefault(acc, []).append(gene)

    if not acc_to_genes:
        logger.info("enrich_genes_with_ec: no genes have uniprot annotations — skipping")
        return

    to_fetch = [a for a in acc_to_genes if a not in ec_cache]
    logger.info(
        f"enrich_genes_with_ec: {len(acc_to_genes)} unique accessions, "
        f"{len(to_fetch)} need fetching ({len(ec_cache)} cached)"
    )

    newly_fetched: dict[str, list[str]] = {}   # only accessions that have EC numbers
    total_hits = 0
    for i in range(0, len(to_fetch), _BATCH_SIZE):
        batch = to_fetch[i: i + _BATCH_SIZE]
        batch_result = _fetch_ec_for_accessions(batch)
        newly_fetched.update(batch_result)
        total_hits += len(batch_result)
        logger.debug(
            f"  Batch {i}–{i + len(batch) - 1}: "
            f"{len(batch_result)}/{len(batch)} accessions have EC numbers"
        )

    if newly_fetched:
        ec_cache.update(newly_fetched)
        with cache_file.open("w") as f:
            json.dump(ec_cache, f)
        logger.info(
            f"  EC cache updated: {total_hits} new EC hits cached, "
            f"{len(to_fetch) - total_hits} accessions had no EC (not cached — will retry next run)"
        )

    # Apply EC numbers to gene annotations
    enriched = 0
    for acc, genes in acc_to_genes.items():
        ec_list = ec_cache.get(acc)
        if not ec_list:
            continue
        for gene in genes:
            ann = dict(gene.annotation) if isinstance(gene.annotation, dict) else {}
            existing = ann.get("ec-code", [])
            if isinstance(existing, str):
                existing = [existing]
            merged = list(dict.fromkeys(existing + ec_list))
            ann["ec-code"] = merged
            gene.annotation = ann
            enriched += 1

    logger.info(
        f"enrich_genes_with_ec: ec-code written to {enriched} gene-accession pairs "
        f"({sum(1 for g in model.genes if g.annotation.get('ec-code'))} genes total with EC)"
    )
