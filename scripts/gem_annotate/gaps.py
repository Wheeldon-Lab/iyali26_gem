"""
gaps.py — gap analysis, MIS audit, and metabolite merging.
"""

import logging
import re

logger = logging.getLogger(__name__)

# Known duplicate metabolite pairs identified by MIS analysis.
# (keep_id, drop_id) — keep_id survives; drop_id is removed.
DUPLICATE_PAIRS: list[tuple[str, str]] = [
    ("m200",  "m1855"),
    ("m772",  "m2043"),
    ("m878",  "m1963"),
]


def find_gaps(model) -> dict:
    """
    Identify blocked reactions and classify their metabolites as orphan/dead-end.

    A reaction is "blocked" if FVA shows max ≈ 0 AND min ≈ 0
    (i.e. it can carry no flux in any feasible steady state).

    For each blocked reaction's metabolites:
      orphan    — metabolite has NO reaction that can produce it
                  (positive stoichiometric coefficient in any reaction)
      dead_end  — metabolite has NO reaction that can consume it
                  (negative stoichiometric coefficient in any reaction)

    Returns
    -------
    dict with keys:
      "blocked_reactions"    : list of reaction IDs
      "orphan_metabolites"   : list of metabolite IDs
      "dead_end_metabolites" : list of metabolite IDs
      "fva_result"           : cobra FVA DataFrame (full, for inspection)
    """
    from cobra.flux_analysis import flux_variability_analysis

    logger.info("Running FVA to identify blocked reactions …")
    logger.info(f"  Model: {len(model.reactions)} reactions, solver=glpk")

    model.solver = "glpk"

    # Run FVA on all reactions; loopless=False for speed
    # fraction_of_optimum=0 means we don't constrain the objective —
    # we want to know which reactions are STRUCTURALLY blocked (no feasible flux
    # regardless of the objective), so we use the full feasible polytope.
    fva = flux_variability_analysis(
        model,
        reaction_list=model.reactions,
        fraction_of_optimum=0.0,
        loopless=False,
    )

    TOL = 1e-6

    # Precompute per-reaction flux reachability from FVA:
    #   can_carry_forward[rxn_id]  = fva.maximum > TOL   (net positive flux possible)
    #   can_carry_backward[rxn_id] = fva.minimum < -TOL  (net negative flux possible)
    # A reaction is "blocked" when neither direction can carry flux.
    can_forward  = fva["maximum"] >  TOL    # Series indexed by rxn_id
    can_backward = fva["minimum"] < -TOL

    blocked_mask = ~can_forward & ~can_backward
    blocked_ids  = fva.index[blocked_mask].tolist()

    logger.info(f"  Blocked reactions: {len(blocked_ids)} / {len(model.reactions)}")

    # Classify metabolites using FVA-derived reachability, not static graph edges.
    #
    # A metabolite M is a "functional orphan" (no reachable producer) if:
    #   every reaction R that could produce M (positive stoichiometric coefficient in any reaction)
    #   has zero reachable flux in that direction.
    #
    # Concretely, for each reaction R involving M:
    #   - Forward direction (coeff > 0) contributes production  if can_forward[R]
    #   - Reverse direction (coeff < 0) contributes production  if can_backward[R]
    #     (running backward turns a consumer into a producer)
    #
    # Symmetrically for dead-end (no reachable consumer).
    #
    # We only evaluate metabolites that appear in at least one blocked reaction
    # (otherwise they are clearly reachable).

    orphan_mets: set[str] = set()
    dead_end_mets: set[str] = set()

    candidate_mets: set = set()
    for rxn_id in blocked_ids:
        candidate_mets.update(model.reactions.get_by_id(rxn_id).metabolites)

    for met in candidate_mets:
        can_produce = False
        can_consume = False
        for rxn in met.reactions:
            coeff = rxn.metabolites[met]
            fwd = bool(can_forward.get(rxn.id, False))
            rev = bool(can_backward.get(rxn.id, False))
            # Forward run produces M if coeff > 0; reverse run produces M if coeff < 0
            if (coeff > 0 and fwd) or (coeff < 0 and rev):
                can_produce = True
            # Forward run consumes M if coeff < 0; reverse run consumes M if coeff > 0
            if (coeff < 0 and fwd) or (coeff > 0 and rev):
                can_consume = True
            if can_produce and can_consume:
                break   # no need to check further

        if not can_produce:
            orphan_mets.add(met.id)
        if not can_consume:
            dead_end_mets.add(met.id)

    logger.info(f"  Orphan metabolites (no reachable producer):  {len(orphan_mets)}")
    logger.info(f"  Dead-end metabolites (no reachable consumer): {len(dead_end_mets)}")

    return {
        "blocked_reactions": blocked_ids,
        "orphan_metabolites": sorted(orphan_mets),
        "dead_end_metabolites": sorted(dead_end_mets),
        "fva_result": fva,
    }


def report_gaps(gaps: dict) -> None:
    """Log a human-readable gap analysis summary."""
    fva = gaps["fva_result"]
    blocked = gaps["blocked_reactions"]
    orphans = gaps["orphan_metabolites"]
    dead_ends = gaps["dead_end_metabolites"]

    logger.info("─── Gap Analysis Report ───────────────────────────────")
    logger.info(f"  Blocked reactions : {len(blocked)}")
    logger.info(f"  Orphan metabolites (no producer)  : {len(orphans)}")
    logger.info(f"  Dead-end metabolites (no consumer): {len(dead_ends)}")

    # Show worst blocked subsystems if reaction subsystem info available
    # (iYli21 may not have subsystems — fall back to first 20 IDs)
    logger.info(f"  First 20 blocked: {blocked[:20]}")
    if orphans:
        logger.info(f"  Orphans (first 10): {orphans[:10]}")
    if dead_ends:
        logger.info(f"  Dead-ends (first 10): {dead_ends[:10]}")
    logger.info("────────────────────────────────────────────────────────")


def audit_mis_reactions(model, mis_metabolite_ids: list[str]) -> None:
    """
    Targeted stoichiometric audit for metabolites identified by MIS analysis.

    For each metabolite in mis_metabolite_ids, prints every reaction that
    involves it along with:
      - reaction ID, name
      - per-element imbalance (C / H / N / O / P / S / charge)
      - reaction equation string (for manual inspection)

    Use this after `merge_duplicate_metabolites()` to verify that the
    duplicates were the cause of inconsistency, and to find any residual
    stoichiometric errors.
    """
    def _element_balance(rxn) -> dict[str, float]:
        """Return per-element net balance across a reaction (should be 0)."""
        balance: dict[str, float] = {}
        for met, coeff in rxn.metabolites.items():
            formula = met.formula or ""
            for elem, cnt_str in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
                cnt = int(cnt_str) if cnt_str else 1
                balance[elem] = balance.get(elem, 0.0) + coeff * cnt
        return balance

    seen_rxns: set[str] = set()
    any_found = False

    for mid in mis_metabolite_ids:
        try:
            met = model.metabolites.get_by_id(mid)
        except KeyError:
            logger.warning(f"audit_mis_reactions: metabolite '{mid}' not in model")
            continue

        for rxn in met.reactions:
            if rxn.id in seen_rxns:
                continue
            seen_rxns.add(rxn.id)
            any_found = True

            balance = _element_balance(rxn)
            imbalanced_elems = {e: v for e, v in balance.items() if abs(v) > 1e-6}

            # Charge balance
            charge_bal = sum(
                coeff * (met_.charge or 0)
                for met_, coeff in rxn.metabolites.items()
            )

            status = "OK" if (not imbalanced_elems and abs(charge_bal) < 1e-6) else "IMBALANCED"
            imb_str = ", ".join(
                f"{e}:{v:+.3g}" for e, v in sorted(imbalanced_elems.items())
            )
            if abs(charge_bal) > 1e-6:
                imb_str += f", charge:{charge_bal:+.3g}"

            logger.info(
                f"  [{status}] {rxn.id}  ({rxn.name or 'no name'})  "
                f"imbalance={imb_str or 'none'}  "
                f"eq={rxn.reaction}"
            )

    if not any_found:
        logger.info("audit_mis_reactions: no reactions found for the given metabolite IDs")


def merge_duplicate_metabolites(
    model,
    pairs: list[tuple[str, str]],
) -> None:
    """
    Merge pairs of metabolites that represent the same chemical species
    but exist as separate nodes in the model (causing stoichiometric
    inconsistency).

    For each (keep_id, drop_id) pair:
      1. In every reaction that uses drop_id, replace drop_id with keep_id,
         summing coefficients if both appear in the same reaction.
      2. Remove drop_id from the model.
      3. Copy annotations from drop to keep (don't overwrite existing keys).

    pairs : list of (keep_id, drop_id) — keep_id survives, drop_id is removed.

    After merging, run `cobra.util.check_mass_balance()` manually to verify
    the reactions involving these metabolites are now balanced.
    """
    merged_count  = 0
    skipped_count = 0

    for keep_id, drop_id in pairs:
        try:
            keep_met = model.metabolites.get_by_id(keep_id)
        except KeyError:
            logger.warning(f"merge_duplicate_metabolites: keep ID '{keep_id}' not found — skipping pair ({keep_id}, {drop_id})")
            skipped_count += 1
            continue
        try:
            drop_met = model.metabolites.get_by_id(drop_id)
        except KeyError:
            logger.warning(f"merge_duplicate_metabolites: drop ID '{drop_id}' not found — skipping pair ({keep_id}, {drop_id})")
            skipped_count += 1
            continue

        # ── Replace drop_met with keep_met in all reactions ───────────────
        for rxn in list(drop_met.reactions):
            drop_coeff = rxn.metabolites[drop_met]
            # If keep_met already appears in this reaction, sum the coefficients
            existing_coeff = rxn.metabolites.get(keep_met, 0.0)
            new_coeff = existing_coeff + drop_coeff
            # Remove drop_met from reaction
            rxn.subtract_metabolites({drop_met: drop_coeff})
            if new_coeff != 0.0:
                rxn.add_metabolites({keep_met: new_coeff - existing_coeff
                                     if existing_coeff else new_coeff})

        # ── Copy annotations (don't overwrite existing keys in keep) ──────
        for key, val in (drop_met.annotation or {}).items():
            if key not in keep_met.annotation:
                keep_met.annotation[key] = val

        # ── Remove drop_met from model ────────────────────────────────────
        model.remove_metabolites([drop_met])
        logger.info(f"  Merged {drop_id} → {keep_id}")
        merged_count += 1

    logger.info(
        f"merge_duplicate_metabolites: {merged_count} pairs merged, "
        f"{skipped_count} skipped"
    )
