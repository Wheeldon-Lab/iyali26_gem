"""
annotate_reactions_extended.py — Supplementary reaction annotation strategies
applied after the primary annotate_reactions() pass.

Strategy A  — Exchange reaction auto-annotation via metabolite BiGG ID / MNXM
Strategy B  — Transport reaction name matching via MetaNetX desc_index / BiGG suffixes
Strategy B5 — Transport fingerprint: MNXM frozenset lookup in fingerprint_index,
              filtered to is_transport MNXR when multiple hits
Strategy C  — GPR → gene EC numbers → MetaNetX MNXR, disambiguated by metabolite fingerprint
"""

import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

_STRIP_RE       = re.compile(r"[^a-z0-9]")
_COMP_SUFFIX_RE = re.compile(r"_[a-z]$")   # strips compartment suffix from BiGG IDs

# Longest/most-specific suffixes first; "r" (reverse/ER) and "p" (too ambiguous) excluded
_TRANSPORT_SUFFIXES = ("tex", "tpp", "abc", "t2", "t3", "pp", "t")

# Ubiquitous metabolites that dilute Jaccard similarity in Strategy C
_UBIQUITOUS_MNXM = frozenset({"MNXM1", "WATER"})


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


def annotate_remaining_reactions(model, reac_xref: dict, reac_prop: dict | None = None) -> None:
    """
    Apply three supplementary annotation strategies to reactions that are still
    unannotated after the primary annotate_reactions() pass.

    Parameters
    ----------
    model      : COBRApy model (modified in-place)
    reac_xref  : dict returned by load_reac_xref() — must include ec_to_mnxr key
    reac_prop  : dict returned by load_reac_prop() — used for fingerprint disambiguation
    """
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

    hit_A = hit_B = hit_B5 = hit_C = ambig_C = no_match = 0
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
                    mnxm = _get_met_mnxm(met)
                    if mnxm:
                        mnxr_id = single_met_mnxm_to_mnxr.get(mnxm)
                # Fallback A2: scan bigg_to_mnxr for EX_/DM_/SK_ BiGG IDs whose
                # single-metabolite equation matches this metabolite's MNXM.
                if mnxr_id is None:
                    mnxm = _get_met_mnxm(met)
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
                # Collect all MNXM IDs from participating metabolites, look up
                # fingerprint_index, then filter hits to is_transport MNXR entries.
                if mnxr_id is None and fingerprint_index:
                    mnxm_set = frozenset(
                        _get_met_mnxm(m) for m in rxn.metabolites if _get_met_mnxm(m)
                    )
                    if len(mnxm_set) >= 2:
                        fp_hits = fingerprint_index.get(mnxm_set, [])
                        if len(fp_hits) == 1:
                            mnxr_id = fp_hits[0]
                            strategy = "B5"
                        elif len(fp_hits) > 1 and transport_mnxr:
                            transport_hits = [m for m in fp_hits if m in transport_mnxr]
                            if len(transport_hits) == 1:
                                mnxr_id = transport_hits[0]
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
                        _get_met_mnxm(met)
                        for met in rxn.metabolites
                        if _get_met_mnxm(met)
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
        else:
            hit_C += 1

    # Diagnostic: log exchange reactions that could not be annotated
    if unmatched_exchanges:
        logger.info(f"  Strategy A: {len(unmatched_exchanges)} exchange reactions unmatched:")
        for rxn_id, rxn_name, met_info in unmatched_exchanges:
            logger.info(f"    {rxn_id} ({rxn_name!r}): {met_info}")

    total = hit_A + hit_B + hit_B5 + hit_C
    unannotated_after = sum(1 for rxn in model.reactions if _is_unannotated(rxn))
    logger.info(
        f"annotate_remaining_reactions: {total} newly annotated | "
        f"A(exchange)={hit_A}  B(transport)={hit_B}  B5(transport-fp)={hit_B5}  "
        f"C(EC→MNXR)={hit_C}  ambiguous-C={ambig_C}  unmatched={no_match}"
    )
    logger.info(
        f"  Reactions still unannotated after both passes: {unannotated_after}"
    )
