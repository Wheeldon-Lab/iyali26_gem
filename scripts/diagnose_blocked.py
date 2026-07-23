#!/usr/bin/env python3
"""
diagnose_blocked.py — classify blocked reactions by actionable root cause.

READ-ONLY: runs FVA on model.xml, finds blocked reactions, and assigns each a
root-cause category based on metabolite reachability (orphan = no producer,
dead-end = no consumer) and reaction bounds.  Writes a per-reaction CSV and
prints a category breakdown.  Does NOT modify the model.

Root-cause categories (first match wins, in priority order):
  exchange_closed   — exchange/demand/sink reaction whose bounds are closed.
  bound_zero        — non-exchange reaction with lb==ub==0 (hard-disabled).
  isolated_cluster  — every non-currency metabolite is BOTH orphan and dead-end
                      (the reaction sits in a disconnected sub-network).
  missing_producer  — at least one substrate is an orphan (nothing makes it).
  missing_consumer  — at least one product is a dead-end (nothing uses it).
  transport_gap     — transport reaction whose moved metabolite is unreachable
                      on one side.
  other             — blocked but none of the above cleanly applies.

Output: data/blocked_root_cause.csv  (overwrites prior content)

Usage:
    python scripts/diagnose_blocked.py [MODEL]
        MODEL  default: model.xml
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_CSV = REPO_ROOT / "data" / "blocked_root_cause.csv"

# Currency / ubiquitous metabolites are ignored when judging "isolated":
# their presence as orphan/dead-end is not what blocks a reaction.
_CURRENCY_HINTS = (
    "h2o", "h+", "h(+)", "proton", "atp", "adp", "amp", "nad", "nadp",
    "nadh", "nadph", "co2", "phosphate", "diphosphate", "ppi", "pi",
    "coenzyme a", "o2", "oxygen",
)


def _is_currency(met) -> bool:
    n = (met.name or "").lower()
    return any(h in n for h in _CURRENCY_HINTS)


def main():
    model_path = sys.argv[1] if len(sys.argv) > 1 else str(REPO_ROOT / "model.xml")
    print(f"Loading model: {model_path}")
    model = read_sbml_model(model_path)
    model.solver = "glpk"

    print("Running FVA to find blocked reactions (1-3 min)…")
    blocked = set(find_blocked_reactions(model, processes=1))
    print(f"  {len(blocked)} blocked reactions")

    # Reachability: orphan = produced by no (non-blocked... here any) reaction,
    # dead-end = consumed by none.  Use stoichiometric sign across all reactions.
    producers: dict = {}
    consumers: dict = {}
    for met in model.metabolites:
        producers[met.id] = 0
        consumers[met.id] = 0
    for rxn in model.reactions:
        for met, coeff in rxn.metabolites.items():
            # reversible reaction both produces and consumes
            rev = rxn.lower_bound < 0 and rxn.upper_bound > 0
            if coeff > 0 or rev:
                producers[met.id] += 1
            if coeff < 0 or rev:
                consumers[met.id] += 1
    orphan = {mid for mid, p in producers.items() if p == 0}
    deadend = {mid for mid, c in consumers.items() if c == 0}

    exchanges = set(model.exchanges)
    demands = set(model.demands)
    sinks = set(model.sinks)
    boundary = exchanges | demands | sinks

    def classify(rxn) -> tuple[str, str]:
        comps = {m.compartment for m in rxn.metabolites}
        is_transport = len(comps) >= 2
        mets = list(rxn.metabolites)
        non_currency = [m for m in mets if not _is_currency(m)]

        if rxn in boundary:
            if rxn.lower_bound == 0 and rxn.upper_bound == 0:
                return "exchange_closed", "bounds lb==ub==0"
            # an open exchange that is still blocked → its metabolite is isolated
            return "exchange_closed", "boundary, metabolite unreachable"

        if rxn.lower_bound == 0 and rxn.upper_bound == 0:
            return "bound_zero", "lb==ub==0"

        check_mets = non_currency or mets
        all_orphan_and_dead = all(
            (m.id in orphan and m.id in deadend) for m in check_mets
        )
        if check_mets and all_orphan_and_dead:
            return "isolated_cluster", "all non-currency mets orphan+deadend"

        sub_orphan = [m.id for m in mets if m.id in orphan]
        prod_deadend = [m.id for m in mets if m.id in deadend]

        if is_transport and (sub_orphan or prod_deadend):
            tag = (sub_orphan + prod_deadend)[0]
            return "transport_gap", f"unreachable: {tag}"
        if sub_orphan:
            return "missing_producer", f"orphan substrate: {sub_orphan[0]}"
        if prod_deadend:
            return "missing_consumer", f"deadend product: {prod_deadend[0]}"
        return "other", "reachable but blocked (loop/stoichiometry)"

    rows = []
    cat_counts: Counter = Counter()
    for rxn in model.reactions:
        if rxn.id not in blocked:
            continue
        cat, detail = classify(rxn)
        cat_counts[cat] += 1
        comps = ",".join(sorted({m.compartment for m in rxn.metabolites}))
        o_mets = [m.id for m in rxn.metabolites if m.id in orphan]
        d_mets = [m.id for m in rxn.metabolites if m.id in deadend]
        rxn_type = ("exchange" if rxn in exchanges else
                    "demand" if rxn in demands else
                    "sink" if rxn in sinks else
                    "transport" if len({m.compartment for m in rxn.metabolites}) >= 2 else
                    "metabolic")
        rows.append({
            "rxn_id": rxn.id,
            "rxn_name": rxn.name or "",
            "rxn_type": rxn_type,
            "root_cause": cat,
            "detail": detail,
            "has_gpr": bool(rxn.gene_reaction_rule),
            "lb": rxn.lower_bound,
            "ub": rxn.upper_bound,
            "compartments": comps,
            "n_orphan_mets": len(o_mets),
            "n_deadend_mets": len(d_mets),
            "orphan_mets": ";".join(o_mets),
            "deadend_mets": ";".join(d_mets),
            "reaction_str": rxn.reaction,
        })

    rows.sort(key=lambda r: (r["root_cause"], r["rxn_id"]))
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                           ["rxn_id", "root_cause"])
        w.writeheader()
        w.writerows(rows)

    print(f"\nBlocked reactions: {len(blocked)}")
    print(f"Orphan metabolites (no producer):  {len(orphan)}")
    print(f"Dead-end metabolites (no consumer): {len(deadend)}")
    print("\nRoot-cause breakdown:")
    for cat, n in cat_counts.most_common():
        print(f"  {n:5}  {cat}")
    print(f"\nReport written: {OUT_CSV}")


if __name__ == "__main__":
    main()
