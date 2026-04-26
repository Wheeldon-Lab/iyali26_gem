"""
reactions.py — reaction annotation via MetaNetX reac_xref.
"""

import logging
from collections import defaultdict

logger = logging.getLogger(__name__)


def annotate_reactions(model, reac_xref: dict, reac_prop: dict | None = None) -> None:
    """
    Annotate reactions with MNXR IDs and cross-database identifiers.

    Three strategies, applied in order until one succeeds:

    A — bigg.reaction annotation → bigg_to_mnxr lookup  (fast, high-precision)
    B — reaction name → MetaNetX description short-name  (string match)
    C — metabolite-set fingerprint in reac_prop.tsv      (structure-based)
        Requires reac_prop dict from load_reac_prop().
        Uses the frozenset of MNXM IDs annotated on the reaction's metabolites.
        Unique match → high-confidence annotation.
        Multiple matches → logged as ambiguous, skipped.

    Modifies model in-place.
    """
    by_mnxr = reac_xref["by_mnxr"]
    desc_index = reac_xref["desc_index"]
    bigg_to_mnxr = reac_xref["bigg_to_mnxr"]
    fingerprint_index = (reac_prop or {}).get("fingerprint_index", {})

    annotated_A = 0
    annotated_B = 0
    annotated_C = 0
    ambiguous_C = 0
    no_match = 0

    def _apply_annotation(rxn, mnxr_id: str) -> None:
        new_ann: dict[str, list] = defaultdict(list)
        for db_prefix, db_id in by_mnxr.get(mnxr_id, []):
            new_ann[db_prefix].append(db_id)
        new_ann["metanetx.reaction"] = [mnxr_id]
        merged = dict(rxn.annotation)
        merged.update(new_ann)
        rxn.annotation = merged

    for rxn in model.reactions:
        mnxr_id = None
        strategy = None

        # Strategy A: existing bigg.reaction annotation
        bigg_raw = rxn.annotation.get("bigg.reaction")
        if bigg_raw:
            bid = bigg_raw[0] if isinstance(bigg_raw, list) else bigg_raw
            mnxr_id = bigg_to_mnxr.get(bid)
            if mnxr_id:
                strategy = "A"

        # Strategy B: reaction name → MetaNetX description short-name
        if mnxr_id is None:
            mnxr_id = desc_index.get(rxn.name.lower().strip())
            if mnxr_id:
                strategy = "B"

        # Strategy C: stoichiometric fingerprint via metabolite MNXM set
        if mnxr_id is None and fingerprint_index:
            mnxm_set = frozenset(
                ann[0] if isinstance(ann, list) else ann
                for met in rxn.metabolites
                for ann in [met.annotation.get("metanetx.chemical")]
                if ann
            )
            if len(mnxm_set) >= 2:
                candidates = fingerprint_index.get(mnxm_set, [])
                if len(candidates) == 1:
                    mnxr_id = candidates[0]
                    strategy = "C"
                elif len(candidates) > 1:
                    logger.debug(
                        f"R:{rxn.id} ({rxn.name!r}): ambiguous fingerprint match "
                        f"→ {candidates} — skipping"
                    )
                    ambiguous_C += 1

        if mnxr_id is None:
            no_match += 1
            continue

        _apply_annotation(rxn, mnxr_id)
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
        f"ambiguous-C={ambiguous_C}, unmatched={no_match}"
    )
