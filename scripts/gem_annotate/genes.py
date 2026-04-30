"""
genes.py — gene annotation via UniProt REST API and NCBI E-utilities.
"""

import json
import logging
import re

import requests

from .config import (
    _KEGG_CONV_URL,
    _NCBI_EFETCH_BATCH,
    _NCBI_EFETCH_URL,
    _NCBI_ESEARCH_URL,
    _PROTEOME_IDS,
    _TIER_B_LIMIT,
    _UNIPROT_SEARCH_URL,
    CACHE_DIR,
)

logger = logging.getLogger(__name__)


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
    Results are cached to data/cache/uniprot_<proteome_id>.json so subsequent
    runs skip the network entirely.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / f"uniprot_{proteome_id}.json"

    if cache_file.exists():
        logger.info(f"  Loading proteome {proteome_id} from cache ({cache_file.name}) …")
        with cache_file.open() as f:
            results = json.load(f)
        logger.info(f"    {len(results)} entries loaded from cache")
        return results

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

    with cache_file.open("w") as f:
        json.dump(results, f)
    logger.info(f"    Cached to {cache_file}")

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

    Results are cached to data/cache/tier_b_results.json keyed by gene ID.
    On subsequent runs, already-queried genes (hit or miss) are skipped entirely.

    Returns {model_gene_id: annotation_dict}
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / "tier_b_results.json"

    # Load existing cache: {gene_id: ann_dict | null}
    # null means we queried before and found nothing (avoid re-querying misses)
    cache: dict[str, dict | None] = {}
    if cache_file.exists():
        with cache_file.open() as f:
            cache = json.load(f)

    to_query = [gid for gid in unmapped if gid not in cache]
    cached_hits = {gid: ann for gid, ann in cache.items() if ann is not None}
    logger.info(
        f"  Tier B: {len(cached_hits)} hits from cache, "
        f"{len(to_query)}/{len(unmapped)} genes need querying"
    )

    mapping: dict[str, dict] = dict(cached_hits)
    calls = 0
    newly_cached: dict[str, dict | None] = {}

    for gid in to_query:
        if calls >= _TIER_B_LIMIT:
            logger.warning(
                f"Tier B: reached {_TIER_B_LIMIT}-query limit; "
                f"{len(to_query) - len(newly_cached)} genes still unmapped"
            )
            break

        candidates = _normalise_locus_tag(gid)
        hit_ann: dict | None = None

        # Try gene_exact first (all candidates), then fall back to gene: (loose)
        query_templates = [
            lambda c: f'gene_exact:"{c}" AND taxonomy_id:4952',
            lambda c: f'gene:"{c}" AND taxonomy_id:4952',
        ]
        for make_query in query_templates:
            if hit_ann:
                break
            for cand in sorted(candidates):   # deterministic order
                query = make_query(cand)
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

        newly_cached[gid] = hit_ann  # store None for misses to skip next time
        if hit_ann:
            mapping[gid] = hit_ann

    # Persist updated cache
    if newly_cached:
        cache.update(newly_cached)
        with cache_file.open("w") as f:
            json.dump(cache, f)
        logger.info(f"  Tier B: cache updated ({len(cache)} total entries)")

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

    Results are cached to data/cache/ncbi_gene_mapping.json keyed by model gene ID.
    On subsequent runs the network is skipped entirely.

    Returns {model_gene_id: annotation_dict}
    """
    import xml.etree.ElementTree as ET

    if not unmapped:
        return {}

    # ── Cache: {model_gene_id: ann_dict} — covers the full organism bulk ──
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / "ncbi_gene_mapping.json"

    if cache_file.exists():
        with cache_file.open() as f:
            full_cache: dict[str, dict] = json.load(f)
        mapping = {gid: full_cache[gid] for gid in unmapped if gid in full_cache}
        logger.info(
            f"  Tier A′: {len(mapping)}/{len(unmapped)} loaded from cache "
            f"({cache_file.name})"
        )
        return mapping

    # ── Build candidate lookup ────────────────────────────────────────────
    candidate_to_model: dict[str, str] = {}
    for gid in unmapped:
        for cand in _normalise_locus_tag(gid):
            candidate_to_model[cand.lower()] = gid

    # ── Step 1: esearch — get all NCBI Gene IDs for Y. lipolytica ─────────
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

    with cache_file.open("w") as f:
        json.dump(mapping, f)
    logger.info(f"  Tier A′: cached to {cache_file.name}")

    return mapping


def _build_kegg_uniprot_index() -> dict[str, str]:
    """
    Download KEGG conv/uniprot/yli and return a reverse index:
        {uniprot_accession → kegg_gene_id}  (e.g. "Q6CFR0" → "yli:YALI_C00198g")

    KEGG's Y. lipolytica genome was re-annotated (DSM 3286 / GCF_014490615.1),
    so KEGG locus tags no longer match model YALI0/YALI1 IDs.  The only reliable
    bridge is via UniProt accession, which is shared across databases.

    Result is cached to data/cache/kegg_uniprot_index.json.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / "kegg_uniprot_index.json"

    if cache_file.exists():
        with cache_file.open() as f:
            acc_to_kegg = json.load(f)
        logger.info(
            f"  KEGG conv: {len(acc_to_kegg)} entries loaded from cache "
            f"({cache_file.name})"
        )
        return acc_to_kegg

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
        with cache_file.open("w") as f:
            json.dump(acc_to_kegg, f)
        logger.info(f"  KEGG conv: cached to {cache_file.name}")
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
