"""
annotate_reactions_extended.py — Supplementary reaction annotation strategies
applied after the primary annotate_reactions() pass.

Strategy A  — Exchange reaction auto-annotation via metabolite BiGG ID / MNXM
Strategy B  — Transport reaction name matching via MetaNetX desc_index / BiGG suffixes
Strategy B5 — Transport fingerprint: MNXM frozenset lookup in fingerprint_index / single_fingerprint_index,
              filtered to is_transport MNXR, then name-token Jaccard when multiple hits
Strategy C  — GPR → gene EC numbers → MetaNetX MNXR, disambiguated by metabolite fingerprint
Strategy D  — Fingerprint re-match for any remaining reaction with ≥2 MNXM annotations
              (catches reactions whose metabolites were annotated after the primary pass)
"""

import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

_STRIP_RE       = re.compile(r"[^a-z0-9]")
_COMP_SUFFIX_RE = re.compile(r"_[a-z]$")   # strips compartment suffix from BiGG IDs

# Longest/most-specific suffixes first; "r" (reverse/ER) and "p" (too ambiguous) excluded
_TRANSPORT_SUFFIXES = ("tex", "tpp", "abc", "t2", "t3", "pp", "t")

# Ubiquitous metabolites that dilute Jaccard similarity in Strategies C and D.
# MNXM1=H+, MNXM3=ATP, MNXM5=NAD+/NADH, MNXM9=ADP, WATER=H2O
_UBIQUITOUS_MNXM = frozenset({"MNXM1", "WATER", "MNXM3", "MNXM9", "MNXM5"})


def _apply_reaction_annotation(rxn, mnxr_id: str, by_mnxr: dict) -> None:
    """Merge MNXR cross-refs into rxn.annotation without overwriting existing keys."""
    new_ann: dict[str, list] = defaultdict(list)
    for db_prefix, db_id in by_mnxr.get(mnxr_id, []):
        new_ann[db_prefix].append(db_id)
    new_ann["metanetx.reaction"] = [mnxr_id]
    merged = dict(rxn.annotation)
    for key, val in new_ann.items():
        if key not in merged:
            merged[key] = val
    rxn.annotation = merged


def _get_met_bigg(met) -> str:
    """Return the bigg.metabolite annotation as a bare string (no compartment suffix)."""
    ann = met.annotation if isinstance(met.annotation, dict) else {}
    raw = ann.get("bigg.metabolite", "")
    if isinstance(raw, list):
        val = raw[0] if raw else ""
    else:
        val = str(raw) if raw else ""
    return _COMP_SUFFIX_RE.sub("", val) if val else ""


def _get_met_mnxm(met) -> str:
    """Return the metanetx.chemical annotation as a plain string."""
    ann = met.annotation if isinstance(met.annotation, dict) else {}
    raw = ann.get("metanetx.chemical", "")
    if isinstance(raw, list):
        return raw[0] if raw else ""
    return str(raw) if raw else ""


def _is_unannotated(rxn) -> bool:
    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    return not ann


def annotate_remaining_reactions(
    model,
    reac_xref: dict,
    reac_prop: dict | None = None,
    mnxm_depr: dict[str, str] | None = None,
) -> None:
    """
    Apply supplementary annotation strategies to reactions still unannotated after
    the primary annotate_reactions() pass.

    Parameters
    ----------
    model      : COBRApy model (modified in-place)
    reac_xref  : dict returned by load_reac_xref() — must include ec_to_mnxr key
    reac_prop  : dict returned by load_reac_prop() — used for fingerprint disambiguation
    mnxm_depr  : optional dict returned by load_mnxm_depr() — maps deprecated MNXM IDs
                 to their current canonical IDs; improves fingerprint match rates
    """
    _depr = mnxm_depr or {}

    def _normalise_mnxm(mnxm_id: str) -> str:
        """Resolve a potentially deprecated MNXM ID to its current canonical form."""
        seen: set[str] = set()
        while mnxm_id in _depr and mnxm_id not in seen:
            seen.add(mnxm_id)
            mnxm_id = _depr[mnxm_id]
        return mnxm_id

    def _get_met_mnxm_norm(met) -> str:
        """Return the metanetx.chemical annotation, normalised through the deprecation map."""
        raw = _get_met_mnxm(met)
        return _normalise_mnxm(raw) if raw else ""

    by_mnxr      = reac_xref["by_mnxr"]
    bigg_to_mnxr = reac_xref["bigg_to_mnxr"]
    desc_index   = reac_xref["desc_index"]

    fingerprint_index        = (reac_prop or {}).get("fingerprint_index", {})
    single_fingerprint_index = (reac_prop or {}).get("single_fingerprint_index", {})
    transport_mnxr           = (reac_prop or {}).get("transport_mnxr", set())

    # ── One-time reverse indexes ──────────────────────────────────────────────

    # mnxr_to_fingerprint: each MNXR → its metabolite frozenset (from fingerprint_index).
    # Used in Layer 2 Jaccard to look up a candidate's fingerprint in O(1) rather
    # than scanning all fingerprint_index entries.
    mnxr_to_fingerprint: dict[str, frozenset] = {}
    for fp, mnxr_list in fingerprint_index.items():
        for mnxr in mnxr_list:
            mnxr_to_fingerprint[mnxr] = fp

    # single_met_mnxm_to_mnxr: MNXM → unique MNXR for exchange-like reactions
    # that appear in exactly one MetaNetX single-metabolite equation.
    single_met_mnxm_to_mnxr: dict[str, str] = {}
    for mnxm, mnxr_list in single_fingerprint_index.items():
        if len(mnxr_list) == 1:
            single_met_mnxm_to_mnxr[mnxm] = mnxr_list[0]

    # ── EC→MNXR merge with source priority ───────────────────────────────────
    # reac_xref hits come from explicit ec-code: cross-refs and EC: description
    # headers — they are more specific than reac_prop classifs column entries.
    # Per EC: if reac_xref has ≥1 hit → use only those; else fall back to reac_prop.
    ec_from_xref: dict[str, list[str]] = reac_xref.get("ec_to_mnxr", {})
    ec_from_prop: dict[str, list[str]] = (reac_prop or {}).get("ec_to_mnxr", {})

    ec_to_mnxr: dict[str, list[str]] = {}
    all_ecs = set(ec_from_xref) | set(ec_from_prop)
    for ec in all_ecs:
        xref_hits = ec_from_xref.get(ec, [])
        if xref_hits:
            ec_to_mnxr[ec] = list(dict.fromkeys(xref_hits))
        else:
            ec_to_mnxr[ec] = list(dict.fromkeys(ec_from_prop.get(ec, [])))

    logger.info(f"  Strategy C: {len(ec_to_mnxr):,} EC numbers available for EC→MNXR lookup")

    # Pre-build a normalised desc_index for fuzzy matching
    norm_desc_index: dict[str, str] = {}
    for k, v in desc_index.items():
        nk = _STRIP_RE.sub("", k)
        if nk not in norm_desc_index:
            norm_desc_index[nk] = v

    # Build a reverse MNXR → short description map for Strategy C name disambiguation
    mnxr_to_desc: dict[str, str] = {}
    for short_name, mnxr in desc_index.items():
        if mnxr not in mnxr_to_desc:
            mnxr_to_desc[mnxr] = short_name

    # ── Strategy D pre-computation ────────────────────────────────────────────
    # Inverted index: MNXM → list of frozensets that contain it.
    # Lets Strategy D skip entries that share zero metabolites with the query,
    # reducing the scan from 77k fingerprints to only overlapping ones.
    mnxm_to_fingerprints: dict[str, list[frozenset]] = defaultdict(list)
    for fp in fingerprint_index:
        for mnxm in fp:
            mnxm_to_fingerprints[mnxm].append(fp)

    # MNXR → set of EC numbers (from reac_prop classifs column).
    # Used in Strategy D step 4 to break name-score ties by EC agreement.
    mnxr_to_ecs: dict[str, set[str]] = defaultdict(set)
    for ec, mnxr_list in ec_from_prop.items():
        for mnxr in mnxr_list:
            mnxr_to_ecs[mnxr].add(ec)

    hit_A = hit_B = hit_B5 = hit_C = hit_D = ambig_C = ambig_D = no_match = 0
    d_new_at_03 = 0   # Strategy D: candidates added by lowering threshold from 0.4 → 0.3
    unmatched_exchanges: list[tuple[str, str, str]] = []   # (rxn.id, rxn.name, met info)

    exchanges = set(model.exchanges)
    demands   = set(model.demands)

    for rxn in model.reactions:
        if not _is_unannotated(rxn):
            continue

        mnxr_id = None
        strategy = None

        # ── Strategy A: exchange reactions ────────────────────────────────
        if rxn in exchanges or (len(rxn.metabolites) == 1):
            if len(rxn.metabolites) == 1:
                met = next(iter(rxn.metabolites))
                bigg_met = _get_met_bigg(met)
                if bigg_met:
                    for cid in (f"EX_{bigg_met}_e", f"EX_{bigg_met}", f"EX_{bigg_met}(e)"):
                        mnxr_id = bigg_to_mnxr.get(cid)
                        if mnxr_id:
                            break
                # Fallback A1: single-MNXM fingerprint lookup
                if mnxr_id is None:
                    mnxm = _get_met_mnxm_norm(met)
                    if mnxm:
                        mnxr_id = single_met_mnxm_to_mnxr.get(mnxm)
                # Fallback A2: scan bigg_to_mnxr for EX_/DM_/SK_ BiGG IDs whose
                # single-metabolite equation matches this metabolite's MNXM.
                if mnxr_id is None:
                    mnxm = _get_met_mnxm_norm(met)
                    if mnxm and single_fingerprint_index:
                        mnxr_candidates = single_fingerprint_index.get(mnxm, [])
                        ex_candidates = [
                            m for m in mnxr_candidates
                            if any(
                                sid.startswith(("EX_", "DM_", "SK_"))
                                for pf, sid in by_mnxr.get(m, [])
                                if pf == "bigg.reaction"
                            )
                        ]
                        if len(ex_candidates) == 1:
                            mnxr_id = ex_candidates[0]
                if mnxr_id:
                    strategy = "A"
                else:
                    # Record for diagnostic logging
                    met_info = (
                        f"bigg={_get_met_bigg(met) or '—'} "
                        f"mnxm={_get_met_mnxm(met) or '—'}"
                    )
                    unmatched_exchanges.append((rxn.id, rxn.name or "", met_info))

        # ── Strategy B: transport reactions ───────────────────────────────
        if mnxr_id is None:
            comps = {m.compartment for m in rxn.metabolites}
            is_transport = len(comps) >= 2

            if is_transport:
                name_lower = (rxn.name or "").lower().strip()

                # B1: exact desc_index match on reaction name
                if name_lower:
                    mnxr_id = desc_index.get(name_lower)
                    if mnxr_id:
                        strategy = "B"

                # B2: normalised name match
                if mnxr_id is None and name_lower:
                    norm_name = _STRIP_RE.sub("", name_lower)
                    mnxr_id = norm_desc_index.get(norm_name)
                    if mnxr_id:
                        strategy = "B"

                # B3: extract metabolite name before "transport" and try desc_index
                if mnxr_id is None and name_lower and "transport" in name_lower:
                    before_transport = name_lower.split("transport")[0].strip().rstrip("-").strip()
                    if before_transport:
                        mnxr_id = desc_index.get(before_transport)
                        if mnxr_id:
                            strategy = "B"

                # B4: met BiGG ID + transport suffixes in bigg_to_mnxr
                if mnxr_id is None:
                    for met in rxn.metabolites:
                        bigg_met = _get_met_bigg(met)
                        if not bigg_met:
                            continue
                        for suf in _TRANSPORT_SUFFIXES:
                            cid = f"{bigg_met}{suf}"
                            mnxr_id = bigg_to_mnxr.get(cid)
                            if mnxr_id:
                                strategy = "B"
                                break
                        if mnxr_id:
                            break

                # B5: MNXM frozenset fingerprint lookup for transport reactions.
                # Normalise MNXM IDs through the deprecation map before lookup.
                # For single-metabolite transports (same compound, two compartments)
                # the deduped frozenset has size 1 → use single_fingerprint_index.
                if mnxr_id is None and (fingerprint_index or single_fingerprint_index):
                    mnxm_set = frozenset(
                        _get_met_mnxm_norm(m) for m in rxn.metabolites if _get_met_mnxm_norm(m)
                    )

                    def _b5_resolve(candidates: list[str], source_label: str) -> str | None:
                        """
                        Filter by transport flag, then name-token Jaccard (gap ≥ 0.1).
                        Logs ambiguous cases; returns the unique winner or None.
                        """
                        if not candidates:
                            return None
                        if len(candidates) == 1:
                            return candidates[0]
                        # Stage 1: filter to transport reactions
                        if transport_mnxr:
                            t_hits = [m for m in candidates if m in transport_mnxr]
                            if len(t_hits) == 1:
                                return t_hits[0]
                            if t_hits:
                                candidates = t_hits
                        # Stage 2: name-token Jaccard (requires best > 0 AND gap ≥ 0.1)
                        rxn_name = (rxn.name or "").strip()
                        if rxn_name:
                            rxn_tokens = set(_STRIP_RE.sub(" ", rxn_name.lower()).split())
                            rxn_tokens.discard("")
                            if rxn_tokens:
                                scores: dict[str, float] = {}
                                for c in candidates:
                                    desc = mnxr_to_desc.get(c, "")
                                    if not desc:
                                        continue
                                    desc_tokens = set(_STRIP_RE.sub(" ", desc).split())
                                    desc_tokens.discard("")
                                    if desc_tokens:
                                        scores[c] = len(rxn_tokens & desc_tokens) / len(rxn_tokens | desc_tokens)
                                if scores:
                                    sorted_sc = sorted(scores.values(), reverse=True)
                                    best_sc   = sorted_sc[0]
                                    second_sc = sorted_sc[1] if len(sorted_sc) > 1 else 0.0
                                    if best_sc > 0 and (best_sc - second_sc) >= 0.1:
                                        return max(scores, key=lambda m: scores[m])
                        logger.debug(
                            f"B5 ambiguous {source_label} {rxn.id!r} ({rxn.name!r}): "
                            f"mnxm={mnxm_set}  candidates={candidates[:8]}"
                        )
                        return None

                    if len(mnxm_set) == 1:
                        solo_mnxm = next(iter(mnxm_set))
                        b5_cands = list(single_fingerprint_index.get(solo_mnxm, []))
                        hit = _b5_resolve(b5_cands, "single")
                        if hit:
                            mnxr_id = hit
                            strategy = "B5"
                    elif len(mnxm_set) >= 2:
                        b5_cands = list(fingerprint_index.get(mnxm_set, []))
                        hit = _b5_resolve(b5_cands, "multi")
                        if hit:
                            mnxr_id = hit
                            strategy = "B5"

        # ── Strategy C: GPR → EC → MetaNetX reaction ─────────────────────
        if mnxr_id is None and rxn.gene_reaction_rule and rxn.gene_reaction_rule.strip():
            ec_numbers: list[str] = []
            for gene in rxn.genes:
                g_ann = gene.annotation if isinstance(gene.annotation, dict) else {}
                ec_raw = g_ann.get("ec-code", [])
                if isinstance(ec_raw, str):
                    ec_raw = [ec_raw]
                ec_numbers.extend(ec_raw)

            if ec_numbers:
                candidates: list[str] = []
                for ec in set(ec_numbers):
                    candidates.extend(ec_to_mnxr.get(ec, []))
                candidates = list(dict.fromkeys(candidates))

                if len(candidates) == 1:
                    mnxr_id = candidates[0]
                    strategy = "C"
                elif len(candidates) > 1:
                    # Remove ubiquitous metabolites (H+, H2O) before Jaccard — they
                    # appear in nearly all reactions and dilute the similarity score.
                    full_mnxm_set = frozenset(
                        _get_met_mnxm_norm(met)
                        for met in rxn.metabolites
                        if _get_met_mnxm_norm(met)
                    )
                    mnxm_set = full_mnxm_set - _UBIQUITOUS_MNXM
                    if not mnxm_set:
                        mnxm_set = full_mnxm_set   # fall back if all metabolites are ubiquitous
                    resolved = False

                    # Layer 1: exact MNXM fingerprint match (requires ≥2 metabolites)
                    if fingerprint_index and len(mnxm_set) >= 2:
                        fp_candidates = fingerprint_index.get(mnxm_set, [])
                        overlap = [c for c in candidates if c in fp_candidates]
                        if len(overlap) == 1:
                            mnxr_id = overlap[0]
                            strategy = "C"
                            resolved = True

                    # Layer 2: Jaccard similarity (lowered thresholds: 0.3 / 0.1).
                    # Uses mnxr_to_fingerprint for O(N_candidates) lookup.
                    if not resolved and fingerprint_index and len(mnxm_set) >= 1:
                        cand_scores: dict[str, float] = {}
                        for c in candidates:
                            fp = mnxr_to_fingerprint.get(c)
                            if fp:
                                fp_filtered = fp - _UBIQUITOUS_MNXM or fp
                                score = len(mnxm_set & fp_filtered) / len(mnxm_set | fp_filtered)
                                cand_scores[c] = score

                        if cand_scores:
                            sorted_scores = sorted(cand_scores.values(), reverse=True)
                            best_score   = sorted_scores[0]
                            second_score = sorted_scores[1] if len(sorted_scores) > 1 else 0.0
                            if best_score >= 0.3 and (best_score - second_score) >= 0.1:
                                best_mnxr = max(cand_scores, key=lambda m: cand_scores[m])
                                mnxr_id = best_mnxr
                                strategy = "C"
                                resolved = True

                    # Layer 3: reaction name vs MNXR description (token Jaccard ≥ 0.5)
                    if not resolved and (rxn.name or "").strip():
                        rxn_tokens = set(_STRIP_RE.sub(" ", rxn.name.lower()).split())
                        rxn_tokens.discard("")
                        best_name_score = 0.0
                        best_name_mnxr  = None
                        for c in candidates:
                            desc = mnxr_to_desc.get(c, "")
                            if not desc:
                                continue
                            desc_tokens = set(_STRIP_RE.sub(" ", desc).split())
                            desc_tokens.discard("")
                            if not rxn_tokens or not desc_tokens:
                                continue
                            jscore = len(rxn_tokens & desc_tokens) / len(rxn_tokens | desc_tokens)
                            if jscore > best_name_score:
                                best_name_score = jscore
                                best_name_mnxr  = c
                        if best_name_score >= 0.5 and best_name_mnxr:
                            mnxr_id = best_name_mnxr
                            strategy = "C"
                            resolved = True

                    if not resolved:
                        ambig_C += 1

        # ── Strategy D: Jaccard + name + EC joint matching ───────────────────
        # Exact fingerprint matching fails for Y. lipolytica-specific metabolite
        # variants (e.g. chain-length-specific acyl-CoAs) because reac_prop only
        # contains generic-substrate versions.  Jaccard similarity over the
        # inverted index finds the closest MetaNetX entry, then name-token overlap
        # and EC agreement break ties.
        if mnxr_id is None and fingerprint_index:
            full_mnxm_set = frozenset(
                _get_met_mnxm_norm(met) for met in rxn.metabolites if _get_met_mnxm_norm(met)
            )
            # Remove ubiquitous stopwords before similarity scoring
            d_query = full_mnxm_set - _UBIQUITOUS_MNXM
            if not d_query:
                d_query = full_mnxm_set

            if len(d_query) >= 2:
                # Step 1: collect candidate fingerprints via inverted index —
                # only frozensets that share ≥1 MNXM with d_query are considered.
                candidate_fps: set[frozenset] = set()
                for mnxm in d_query:
                    for fp in mnxm_to_fingerprints.get(mnxm, []):
                        candidate_fps.add(fp)

                # Step 2: score each candidate fingerprint by Jaccard ≥ 0.3.
                # Lowered from 0.4; disambiguation in step 3 prevents false positives.
                jacc_hits: dict[str, float] = {}   # mnxr_id → score
                jacc_hits_strict: dict[str, float] = {}   # subset with score ≥ 0.4 (for logging)
                for fp in candidate_fps:
                    fp_filt = fp - _UBIQUITOUS_MNXM or fp
                    score = len(d_query & fp_filt) / len(d_query | fp_filt)
                    if score >= 0.3:
                        for mnxr in fingerprint_index[fp]:
                            if score > jacc_hits.get(mnxr, -1):
                                jacc_hits[mnxr] = score
                            if score >= 0.4 and score > jacc_hits_strict.get(mnxr, -1):
                                jacc_hits_strict[mnxr] = score
                d_new_at_03 += len(set(jacc_hits) - set(jacc_hits_strict))

                if jacc_hits:
                    if len(jacc_hits) == 1:
                        mnxr_id = next(iter(jacc_hits))
                        strategy = "D"
                    else:
                        # Step 3: name-token overlap disambiguation
                        # tokenise reaction name once
                        rxn_tokens = set(_STRIP_RE.sub(" ", (rxn.name or "").lower()).split())
                        rxn_tokens.discard("")

                        name_overlap: dict[str, int] = {}
                        if rxn_tokens:
                            for mnxr in jacc_hits:
                                desc = mnxr_to_desc.get(mnxr, "")
                                if not desc:
                                    continue
                                desc_tokens = set(_STRIP_RE.sub(" ", desc).split())
                                desc_tokens.discard("")
                                name_overlap[mnxr] = len(rxn_tokens & desc_tokens)

                        if name_overlap:
                            sorted_overlap = sorted(name_overlap.values(), reverse=True)
                            best_ov  = sorted_overlap[0]
                            second_ov = sorted_overlap[1] if len(sorted_overlap) > 1 else -1
                            if best_ov >= 1 and (best_ov - second_ov) >= 1:
                                winners = [m for m, v in name_overlap.items() if v == best_ov]
                                if len(winners) == 1:
                                    w = winners[0]
                                    j = jacc_hits[w]
                                    # Compound gate: (J≥0.3 AND overlap≥2) OR (J≥0.5 AND overlap≥1)
                                    if (j >= 0.3 and best_ov >= 2) or (j >= 0.5 and best_ov >= 1):
                                        mnxr_id = w
                                        strategy = "D"

                        # Step 4: EC agreement tiebreak (only when name step is tied)
                        if mnxr_id is None and rxn.gene_reaction_rule:
                            rxn_ecs: set[str] = set()
                            for gene in rxn.genes:
                                g_ann = gene.annotation if isinstance(gene.annotation, dict) else {}
                                ec_raw = g_ann.get("ec-code", [])
                                if isinstance(ec_raw, str):
                                    ec_raw = [ec_raw]
                                rxn_ecs.update(ec_raw)

                            if rxn_ecs:
                                ec_matches = [
                                    m for m in jacc_hits
                                    if mnxr_to_ecs.get(m, set()) & rxn_ecs
                                ]
                                if len(ec_matches) == 1:
                                    mnxr_id = ec_matches[0]
                                    strategy = "D"

                        if mnxr_id is None:
                            ambig_D += 1
                            logger.debug(
                                f"Strategy D ambiguous {rxn.id!r} ({rxn.name!r}): "
                                f"query={sorted(d_query)}  "
                                f"top_hits={sorted(jacc_hits, key=lambda m: jacc_hits[m], reverse=True)[:5]}"
                            )

        if mnxr_id is None:
            no_match += 1
            continue

        _apply_reaction_annotation(rxn, mnxr_id, by_mnxr)
        if strategy == "A":
            hit_A += 1
        elif strategy == "B5":
            hit_B5 += 1
        elif strategy == "B":
            hit_B += 1
        elif strategy == "D":
            hit_D += 1
        else:
            hit_C += 1

    # Diagnostic: log exchange reactions that could not be annotated
    if unmatched_exchanges:
        logger.info(f"  Strategy A: {len(unmatched_exchanges)} exchange reactions unmatched:")
        for rxn_id, rxn_name, met_info in unmatched_exchanges:
            logger.info(f"    {rxn_id} ({rxn_name!r}): {met_info}")

    total = hit_A + hit_B + hit_B5 + hit_C + hit_D
    unannotated_after = sum(1 for rxn in model.reactions if _is_unannotated(rxn))
    logger.info(
        f"annotate_remaining_reactions: {total} newly annotated | "
        f"A(exchange)={hit_A}  B(transport)={hit_B}  B5(transport-fp)={hit_B5}  "
        f"C(EC→MNXR)={hit_C}  D(Jaccard+name+EC)={hit_D}  "
        f"ambiguous-C={ambig_C}  ambiguous-D={ambig_D}  unmatched={no_match}  "
        f"D_extra_candidates_at_0.3={d_new_at_03}"
    )
    logger.info(
        f"  Reactions still unannotated after both passes: {unannotated_after}"
    )
