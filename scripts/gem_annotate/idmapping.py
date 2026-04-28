"""
idmapping.py — UniProt ID-mapping enrichment step.

Converts NCBI Gene IDs (already present in gene.annotation["ncbigene"]) to
UniProt accessions via the UniProt ID-mapping REST API, then merges the result
back into model gene annotations without overwriting existing values.
"""

import logging
import time

import requests

logger = logging.getLogger(__name__)

_IDMAPPING_RUN_URL    = "https://rest.uniprot.org/idmapping/run"
_IDMAPPING_STATUS_URL = "https://rest.uniprot.org/idmapping/status/{jobId}"
_IDMAPPING_RESULTS_URL = "https://rest.uniprot.org/idmapping/uniprotkb/results/{jobId}"

_POLL_INTERVAL   = 5   # seconds between status polls
_POLL_TIMEOUT    = 300 # seconds before giving up


def _submit_idmapping(ncbi_ids: list[str]) -> str | None:
    """POST to UniProt ID-mapping and return the jobId, or None on failure."""
    try:
        resp = requests.post(
            _IDMAPPING_RUN_URL,
            data={
                "from": "GeneID",
                "to":   "UniProtKB",
                "ids":  ",".join(ncbi_ids),
            },
            timeout=60,
        )
        resp.raise_for_status()
        job_id = resp.json().get("jobId")
        if not job_id:
            logger.warning("  ID-mapping: no jobId in response")
        return job_id
    except Exception as e:
        logger.warning(f"  ID-mapping: submission failed: {e}")
        return None


def _poll_idmapping(job_id: str) -> bool:
    """Poll status endpoint until the job is FINISHED. Returns True on success."""
    url = _IDMAPPING_STATUS_URL.format(jobId=job_id)
    elapsed = 0
    while elapsed < _POLL_TIMEOUT:
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            status = data.get("jobStatus", "")
            if status == "FINISHED":
                return True
            if status == "FAILED":
                logger.warning(f"  ID-mapping job {job_id} FAILED")
                return False
            # RUNNING or NEW — keep polling
        except Exception as e:
            logger.debug(f"  ID-mapping poll error: {e}")
        time.sleep(_POLL_INTERVAL)
        elapsed += _POLL_INTERVAL
    logger.warning(f"  ID-mapping: job {job_id} timed out after {_POLL_TIMEOUT}s")
    return False


def _fetch_idmapping_results(job_id: str) -> list[dict]:
    """
    Fetch all paginated results from the uniprotkb results endpoint.
    Returns a list of {from: ncbi_gene_id, to: {UniProt entry}} dicts.
    """
    results: list[dict] = []
    url: str | None = _IDMAPPING_RESULTS_URL.format(jobId=job_id)
    params: dict = {"format": "json", "size": 500}
    import re

    while url:
        try:
            resp = requests.get(url, params=params, timeout=120)
            resp.raise_for_status()
        except Exception as e:
            logger.warning(f"  ID-mapping: results fetch failed: {e}")
            break
        data = resp.json()
        results.extend(data.get("results", []))
        # Follow Link: rel="next" for subsequent pages
        url = None
        params = {}
        for part in resp.headers.get("Link", "").split(","):
            part = part.strip()
            if 'rel="next"' in part:
                m = re.search(r"<([^>]+)>", part)
                if m:
                    url = m.group(1)

    return results


def _parse_idmapping_entry(result: dict) -> tuple[str, dict]:
    """
    Parse one ID-mapping result row.

    Returns (ncbi_gene_id, annotation_dict) where annotation_dict contains
    uniprot, kegg.genes, refseq, ncbigene with single string values
    (most-reliable pick per key, consistent with how write_sbml_model serialises).
    """
    ncbi_id = str(result.get("from", ""))
    entry   = result.get("to", {})
    ann: dict = {}

    acc = entry.get("primaryAccession", "")
    if acc:
        ann["uniprot"] = acc  # str, not list — see note in module docstring

    for xref in entry.get("uniProtKBCrossReferences", []):
        db  = xref.get("database", "")
        xid = xref.get("id", "")
        if not xid:
            continue
        if db == "KEGG" and "kegg.genes" not in ann:
            ann["kegg.genes"] = xid
        elif db == "RefSeq" and "refseq" not in ann:
            ann["refseq"] = xid
        elif db == "GeneID" and "ncbigene" not in ann:
            ann["ncbigene"] = xid

    return ncbi_id, ann


def _pick_best_entry(entries: list[dict]) -> dict | None:
    """
    From a list of UniProt entries mapped from one NCBI Gene ID, prefer a
    reviewed (Swiss-Prot) entry; fall back to the first unreviewed (TrEMBL).
    Returns the chosen entry dict, or None if the list is empty.
    """
    if not entries:
        return None
    reviewed = [e for e in entries if e.get("to", {}).get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"]
    return reviewed[0] if reviewed else entries[0]


def _merge_str_annotation(gene, new_data: dict) -> None:
    """
    Merge new_data (string values) into gene.annotation without overwriting.
    If the key already exists, skip — existing data takes priority.
    """
    merged = dict(gene.annotation)
    for key, val in new_data.items():
        if key not in merged or merged[key] is None:
            merged[key] = val
        # key already present — do not overwrite
    gene.annotation = merged


def _enrich_via_idmapping(model) -> None:
    """
    Enrich model genes that have ncbigene but no uniprot annotation by querying
    the UniProt ID-mapping API (GeneID → UniProtKB).

    Steps
    -----
    1. Collect genes with ncbigene and no uniprot.
    2. Submit all NCBI Gene IDs in one batch request.
    3. Poll until the job finishes.
    4. Fetch paginated results, group by NCBI Gene ID, pick best UniProt entry.
    5. Merge annotation back into model genes (no-overwrite policy).
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

    ncbi_ids = list(targets.keys())
    logger.info(
        f"=== UniProt ID-mapping: submitting {len(ncbi_ids)} NCBI Gene IDs "
        f"(GeneID → UniProtKB) ==="
    )

    # ── 2. Submit ─────────────────────────────────────────────────────────
    job_id = _submit_idmapping(ncbi_ids)
    if not job_id:
        logger.warning("  ID-mapping: skipped (submission failed)")
        return
    logger.info(f"  ID-mapping: job submitted, jobId={job_id}")

    # ── 3. Poll ───────────────────────────────────────────────────────────
    if not _poll_idmapping(job_id):
        logger.warning("  ID-mapping: skipped (job did not finish)")
        return
    logger.info("  ID-mapping: job finished, fetching results …")

    # ── 4. Fetch + parse ──────────────────────────────────────────────────
    raw_results = _fetch_idmapping_results(job_id)
    logger.info(f"  ID-mapping: {len(raw_results)} result rows received")

    # Group multiple UniProt entries per NCBI Gene ID, then pick best
    from collections import defaultdict
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in raw_results:
        ncbi_id = str(row.get("from", ""))
        if ncbi_id:
            grouped[ncbi_id].append(row)

    gene_lookup = {g.id: g for g in model.genes}
    merged_count    = 0
    enriched_keys: dict[str, int] = {}  # key → count of genes that gained it

    for ncbi_id, rows in grouped.items():
        model_gene_id = targets.get(ncbi_id)
        if not model_gene_id:
            continue
        gene = gene_lookup.get(model_gene_id)
        if gene is None:
            continue

        best = _pick_best_entry(rows)
        if best is None:
            continue

        _, ann = _parse_idmapping_entry(best)
        if not ann:
            continue

        _merge_str_annotation(gene, ann)
        merged_count += 1
        for k in ann:
            enriched_keys[k] = enriched_keys.get(k, 0) + 1

    # ── 5. Report ─────────────────────────────────────────────────────────
    keys_summary = ", ".join(f"{k}×{v}" for k, v in sorted(enriched_keys.items()))
    logger.info(
        f"  ID-mapping: {merged_count}/{len(targets)} genes enriched. "
        f"Keys added — {keys_summary or 'none'}"
    )
