"""
idmapping.py — UniProt search-based enrichment step.

Converts NCBI Gene IDs (already present in gene.annotation["ncbigene"]) to
UniProt accessions via the UniProt search API using xref:geneid- queries,
then merges the result back into model gene annotations without overwriting
existing values.

Uses synchronous GET requests (no job/poll), batched at 30 Gene IDs per
request to stay well under URL-length limits.
"""

import json
import logging
import re
import time

import requests

from .config import CACHE_DIR
from .http_utils import _request_with_retry

logger = logging.getLogger(__name__)

_UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
_SEARCH_BATCH_SIZE  = 30    # Gene IDs per request (~600 chars of query string)
_RATE_LIMIT_SLEEP   = 0.2   # seconds between requests (UniProt rate limit)


def _search_uniprot_by_geneids(ncbi_ids: list[str]) -> dict[str, dict]:
    """
    Query UniProt search API for a batch of NCBI Gene IDs using xref:geneid- syntax.

    Returns {ncbi_gene_id: annotation_dict} for every ID that resolves to a
    UniProt entry.  When multiple entries map to one Gene ID, the reviewed
    (Swiss-Prot) entry is preferred; unreviewed (TrEMBL) is the fallback.

    annotation_dict keys: uniprot (str), and optionally ncbigene, kegg.genes,
    refseq — all as single strings for consistency with _merge_str_annotation.
    """
    if not ncbi_ids:
        return {}

    query = " OR ".join(f"xref:geneid-{nid}" for nid in ncbi_ids)
    params = {
        "query":  f"({query})",
        "fields": "accession,reviewed,xref_geneid,xref_kegg,xref_refseq",
        "format": "tsv",
        "size":   500,
    }

    rows: list[str] = []
    url: str | None = _UNIPROT_SEARCH_URL

    while url:
        try:
            resp = _request_with_retry("GET", url, params=params, timeout=60)
            resp.raise_for_status()
        except Exception as e:
            logger.warning(f"  ID-mapping: search request failed: {e}")
            return {}

        lines = resp.text.strip().splitlines()
        if not rows:
            # First page: keep header to identify columns
            rows.extend(lines)
        else:
            # Subsequent pages: skip header line
            rows.extend(lines[1:])

        # Follow Link: rel="next" for pagination
        url = None
        params = {}
        for part in resp.headers.get("Link", "").split(","):
            part = part.strip()
            if 'rel="next"' in part:
                m = re.search(r"<([^>]+)>", part)
                if m:
                    url = m.group(1)

    if not rows:
        return {}

    # Parse TSV: header is first row
    header = [h.strip() for h in rows[0].split("\t")]
    try:
        col_acc      = header.index("Entry")
        col_reviewed = header.index("Reviewed")
        col_geneid   = header.index("Gene Names (ordered locus)")  # fallback name
    except ValueError:
        # Column names vary slightly; use positional fallback
        col_acc, col_reviewed, col_geneid = 0, 1, 2

    # Also find optional cross-ref columns by scanning header
    def _col(name: str) -> int | None:
        for i, h in enumerate(header):
            if name.lower() in h.lower():
                return i
        return None

    col_geneid_xref = _col("GeneID")   # "Gene Names (ordered locus)" vs "GeneID; ..."
    col_kegg        = _col("KEGG")
    col_refseq      = _col("RefSeq")

    # {ncbi_id: {"reviewed": bool, "acc": str, "kegg": str, "refseq": str}}
    best: dict[str, dict] = {}

    for line in rows[1:]:
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) <= col_acc:
            continue
        acc      = parts[col_acc]
        reviewed = parts[col_reviewed].lower() == "reviewed" if col_reviewed < len(parts) else False

        # Extract all GeneIDs from the xref_geneid column (semicolon-separated)
        gene_ids_in_row: list[str] = []
        if col_geneid_xref is not None and col_geneid_xref < len(parts):
            raw = parts[col_geneid_xref]
            # UniProt TSV format: "12345; 67890;" or just "12345"
            gene_ids_in_row = [g.strip().rstrip(";") for g in raw.split(";") if g.strip().rstrip(";")]

        def _first(raw: str) -> str:
            idx = raw.find(";")
            return raw[:idx].strip() if idx != -1 else raw.strip()

        kegg_val   = _first(parts[col_kegg])   if col_kegg   and col_kegg   < len(parts) and parts[col_kegg]   else ""
        refseq_val = _first(parts[col_refseq]) if col_refseq and col_refseq < len(parts) and parts[col_refseq] else ""

        for nid in gene_ids_in_row:
            if nid not in ncbi_ids:
                continue  # result row for a Gene ID we didn't ask for (cross-ref)
            existing = best.get(nid)
            # Prefer reviewed entry; replace unreviewed with reviewed
            if existing is None or (reviewed and not existing["reviewed"]):
                best[nid] = {
                    "reviewed": reviewed,
                    "acc":      acc,
                    "kegg":     kegg_val,
                    "refseq":   refseq_val,
                    "ncbigene": nid,
                }

    # Convert to annotation dicts; all values are list[str] per SBML/Memote spec
    result: dict[str, dict] = {}
    for nid, data in best.items():
        ann: dict = {"uniprot": [data["acc"]], "ncbigene": [data["ncbigene"]]}
        if data["kegg"]:
            ann["kegg.genes"] = [data["kegg"]]
        if data["refseq"]:
            ann["refseq"] = [data["refseq"]]
        result[nid] = ann

    return result


def _merge_str_annotation(gene, new_data: dict) -> None:
    """
    Merge new_data into gene.annotation without overwriting existing keys.
    Scalar string values are wrapped in lists to satisfy SBML/Memote spec.
    """
    merged = dict(gene.annotation)
    for key, val in new_data.items():
        if key not in merged or merged[key] is None:
            merged[key] = [val] if isinstance(val, str) else val
    gene.annotation = merged


def _enrich_via_idmapping(model) -> None:
    """
    Enrich model genes that have ncbigene but no uniprot annotation by querying
    the UniProt search API (xref:geneid-<id> → UniProtKB).

    Steps
    -----
    1. Collect genes with ncbigene and no uniprot.
    2. Check cache (data/cache/idmapping_results.json); query only uncached IDs.
    3. Batch-query UniProt search API (30 IDs/request, synchronous GET).
    4. Persist new hits to cache.
    5. Merge all results (cache + newly fetched) into model gene annotations.
    """
    # ── 1. Collect targets ────────────────────────────────────────────────
    targets: dict[str, str] = {}  # {ncbi_gene_id → model_gene_id}
    for gene in model.genes:
        ann = gene.annotation
        has_uniprot = bool(ann.get("uniprot"))
        ncbigene_val = ann.get("ncbigene")
        if has_uniprot or not ncbigene_val:
            continue
        # ncbigene may be str or list
        if isinstance(ncbigene_val, list):
            ncbi_id = ncbigene_val[0] if ncbigene_val else None
        else:
            ncbi_id = str(ncbigene_val)
        if ncbi_id:
            targets[ncbi_id] = gene.id

    if not targets:
        logger.info("  ID-mapping: no genes need ncbigene→uniprot enrichment")
        return

    # ── 2. Load cache; only query IDs not already cached ─────────────────
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / "idmapping_results.json"

    cache: dict[str, dict | None] = {}
    if cache_file.exists():
        with cache_file.open() as f:
            cache = json.load(f)

    to_query = [nid for nid in targets if nid not in cache]
    cached_hits = {nid: ann for nid, ann in cache.items() if ann is not None and nid in targets}
    logger.info(
        f"  ID-mapping: {len(cached_hits)} hits from cache, "
        f"{len(to_query)}/{len(targets)} IDs need querying"
    )

    # ── 3. Batch-query UniProt search API ─────────────────────────────────
    newly_cached: dict[str, dict] = {}

    if to_query:
        n_batches = (len(to_query) + _SEARCH_BATCH_SIZE - 1) // _SEARCH_BATCH_SIZE
        logger.info(
            f"=== UniProt search: querying {len(to_query)} NCBI Gene IDs "
            f"in {n_batches} batch(es) of {_SEARCH_BATCH_SIZE} ==="
        )

        for batch_start in range(0, len(to_query), _SEARCH_BATCH_SIZE):
            batch = to_query[batch_start: batch_start + _SEARCH_BATCH_SIZE]
            batch_num = batch_start // _SEARCH_BATCH_SIZE + 1

            batch_result = _search_uniprot_by_geneids(batch)
            newly_cached.update(batch_result)
            logger.info(
                f"  ID-mapping: batch {batch_num}/{n_batches} "
                f"({len(batch)} queried, {len(batch_result)} hits)"
            )

            if batch_num < n_batches:
                time.sleep(_RATE_LIMIT_SLEEP)

        # ── 4. Persist new hits ───────────────────────────────────────────
        if newly_cached:
            cache.update(newly_cached)
            with cache_file.open("w") as f:
                json.dump(cache, f)
            no_result_count = len(to_query) - len(newly_cached)
            logger.info(
                f"  ID-mapping: cache updated ({len(cache)} total entries). "
                f"{len(newly_cached)} new hits, "
                f"{no_result_count} IDs had no result (not cached — will retry next run)"
            )

    # ── 5. Apply all results to model ─────────────────────────────────────
    gene_lookup = {g.id: g for g in model.genes}
    merged_count = 0
    enriched_keys: dict[str, int] = {}

    all_results = {**cached_hits, **newly_cached}
    for ncbi_id, ann in all_results.items():
        if not ann:
            continue
        model_gene_id = targets.get(ncbi_id)
        if not model_gene_id:
            continue
        gene = gene_lookup.get(model_gene_id)
        if gene is None:
            continue
        _merge_str_annotation(gene, ann)
        merged_count += 1
        for k in ann:
            enriched_keys[k] = enriched_keys.get(k, 0) + 1

    keys_summary = ", ".join(f"{k}×{v}" for k, v in sorted(enriched_keys.items()))
    logger.info(
        f"  ID-mapping: {merged_count}/{len(targets)} genes enriched. "
        f"Keys added — {keys_summary or 'none'}"
    )
