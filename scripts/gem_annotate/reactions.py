"""
reactions.py — reaction annotation via MetaNetX reac_xref.
"""

import logging
from collections import defaultdict

from .annotate_reactions_extended import _MNXM_NORMALIZE, _resolve_multi_candidates

logger = logging.getLogger(__name__)


def annotate_reactions(model, reac_xref: dict, reac_prop: dict | None = None) -> None:
    """
    Annotate reactions with MNXR IDs and cross-database identifiers.

    Three strategies, applied in order until one succeeds:

    A — bigg.reaction annotation → bigg_to_mnxr lookup  (fast, high-precision)
    C — metabolite-set fingerprint in reac_prop.tsv      (structure-based)
        Requires reac_prop dict from load_reac_prop().
        Uses the frozenset of MNXM IDs annotated on the reaction's metabolites.
        Unique match → high-confidence annotation.
        Multiple matches → resolved by compartment structure / xref count / lex order.
    B — reaction name → MetaNetX description short-name  (string match, last resort)

    Modifies model in-place.
    """
    by_mnxr = reac_xref["by_mnxr"]
    desc_index = reac_xref["desc_index"]
    bigg_to_mnxr = reac_xref["bigg_to_mnxr"]
    fingerprint_index        = (reac_prop or {}).get("fingerprint_index", {})
    mnxr_compartment_type    = (reac_prop or {}).get("mnxr_compartment_type", {})

    # EC→MNXR from reac_prop for use by _resolve_multi_candidates
    ec_from_prop: dict[str, list[str]] = (reac_prop or {}).get("ec_to_mnxr", {})
    ec_to_mnxr: dict[str, list[str]]   = ec_from_prop
    mnxr_to_ecs: dict[str, set[str]] = defaultdict(set)
    for ec, mnxr_list in ec_from_prop.items():
        for mnxr in mnxr_list:
            mnxr_to_ecs[mnxr].add(ec)

    annotated_A = 0
    annotated_B = 0
    annotated_C = 0
    ambiguous_C = 0
    multi_resolved_compartment = 0
    multi_resolved_xref = 0
    no_match = 0

    def _apply_annotation(rxn, mnxr_id: str, extra_mnxrs: list[str] | None = None) -> None:
        new_ann: dict[str, list] = defaultdict(list)
        seen: dict[str, set] = defaultdict(set)
        for source_id in [mnxr_id] + (extra_mnxrs or []):
            for db_prefix, db_id in by_mnxr.get(source_id, []):
                if db_id not in seen[db_prefix]:
                    seen[db_prefix].add(db_id)
                    new_ann[db_prefix].append(db_id)
        new_ann["metanetx.reaction"] = [mnxr_id]
        merged = dict(rxn.annotation)
        for key, val in new_ann.items():
            if key not in merged:
                merged[key] = val
        rxn.annotation = merged

    for rxn in model.reactions:
        mnxr_id = None
        strategy = None
        extra_mnxrs: list[str] | None = None

        # Strategy A: existing bigg.reaction annotation
        bigg_raw = rxn.annotation.get("bigg.reaction")
        if bigg_raw:
            bid = bigg_raw[0] if isinstance(bigg_raw, list) else bigg_raw
            mnxr_id = bigg_to_mnxr.get(bid)
            if mnxr_id:
                strategy = "A"

        # Strategy C: stoichiometric fingerprint via metabolite MNXM set.
        # Apply _MNXM_NORMALIZE so non-standard model IDs map to reac_prop canonical IDs.
        if mnxr_id is None and fingerprint_index:
            mnxm_set = frozenset(
                _MNXM_NORMALIZE.get(raw_id, raw_id)
                for met in rxn.metabolites
                for ann in [met.annotation.get("metanetx.chemical")]
                if ann
                for raw_id in [ann[0] if isinstance(ann, list) else ann]
            )
            if len(mnxm_set) >= 2:
                candidates = fingerprint_index.get(mnxm_set, [])
                if len(candidates) == 1:
                    mnxr_id = candidates[0]
                    strategy = "C"
                elif len(candidates) > 1:
                    winner, equiv, method = _resolve_multi_candidates(
                        list(candidates), rxn, by_mnxr, mnxr_compartment_type,
                        ec_to_mnxr, mnxr_to_ecs,
                        context=f"C {rxn.id!r}",
                    )
                    if method == "compartment":
                        multi_resolved_compartment += 1
                    elif method in ("xref_count", "lex"):
                        multi_resolved_xref += 1
                    elif method == "ambiguous":
                        ambiguous_C += 1
                        logger.debug(
                            f"R:{rxn.id} ({rxn.name!r}): truly ambiguous fingerprint match "
                            f"→ {candidates} — skipping"
                        )
                    if winner:
                        mnxr_id = winner
                        extra_mnxrs = [m for m in (equiv or []) if m != winner]
                        strategy = "C"

        # Strategy B: reaction name → MetaNetX description short-name (last resort)
        if mnxr_id is None:
            name_lower = (rxn.name or "").lower().strip()
            if name_lower:
                mnxr_id = desc_index.get(name_lower)
                if mnxr_id:
                    strategy = "B"

        if mnxr_id is None:
            no_match += 1
            continue

        _apply_annotation(rxn, mnxr_id, extra_mnxrs)
        if strategy == "A":
            annotated_A += 1
        elif strategy == "B":
            annotated_B += 1
        else:
            annotated_C += 1

    total = annotated_A + annotated_B + annotated_C
    logger.info(
        f"Reactions annotated: {total} total  "
        f"(A={annotated_A}, B={annotated_B}, C={annotated_C}), "
        f"C-multi(compartment={multi_resolved_compartment}, xref/lex={multi_resolved_xref}), "
        f"ambiguous-C={ambiguous_C}, unmatched={no_match}"
    )


_BACKFILL_KEYS = ("bigg.reaction", "kegg.reaction", "rhea", "ec-code")


def backfill_reaction_xrefs(model, reac_xref: dict, reac_prop: dict | None = None) -> None:
    """
    Safety-net pass: for every reaction that already has metanetx.reaction,
    fill in any missing cross-references from reac_xref["by_mnxr"] and,
    as a fallback for ec-code, from reac_prop["mnxr_to_ec"].

    Only adds; never overwrites existing annotation values.
    """
    by_mnxr = reac_xref["by_mnxr"]
    mnxr_to_ec: dict[str, list[str]] = (reac_prop or {}).get("mnxr_to_ec", {})

    counts = {k: 0 for k in _BACKFILL_KEYS}

    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        mnxr_raw = ann.get("metanetx.reaction")
        if not mnxr_raw:
            continue
        mnxr = mnxr_raw[0] if isinstance(mnxr_raw, list) else mnxr_raw

        # Step 1: fill from reac_xref cross-references, grouped by prefix
        xref_by_prefix: dict[str, list[str]] = defaultdict(list)
        seen: dict[str, set] = defaultdict(set)
        for db_prefix, db_id in by_mnxr.get(mnxr, []):
            if db_prefix in _BACKFILL_KEYS and db_id not in seen[db_prefix]:
                seen[db_prefix].add(db_id)
                xref_by_prefix[db_prefix].append(db_id)

        for key in _BACKFILL_KEYS:
            if key not in ann and xref_by_prefix.get(key):
                if not isinstance(rxn.annotation, dict):
                    rxn.annotation = {}
                rxn.annotation[key] = xref_by_prefix[key]
                ann = rxn.annotation
                counts[key] += 1

        # Step 2: ec-code fallback from reac_prop classifs if still missing
        if "ec-code" not in ann:
            ec_list = mnxr_to_ec.get(mnxr, [])
            if ec_list:
                if not isinstance(rxn.annotation, dict):
                    rxn.annotation = {}
                rxn.annotation["ec-code"] = ec_list
                counts["ec-code"] += 1

    logger.info(
        f"backfill_reaction_xrefs: "
        f"bigg={counts['bigg.reaction']}, "
        f"kegg={counts['kegg.reaction']}, "
        f"rhea={counts['rhea']}, "
        f"ec-code={counts['ec-code']}"
    )
