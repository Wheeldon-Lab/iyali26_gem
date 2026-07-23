#!/usr/bin/env python3
"""Funnel stage 6 (pre-processing): tier the 90 new_reaction candidates so
manual curation can focus on the few high-value, easy-to-add ones.

Read-only. No network. For each candidate gene (reaction_status=new_reaction)
compute automatable signals and assign a tier:

  T1_clear     : EC -> exactly one MNXR, and that MNXR is balanced
                 (single unambiguous reaction, safe to hand-curate)
  T2_priority  : balanced MNXR(s) AND the gene hits a priority pathway
                 (sphingolipid/lipid — ties to prior ceramide work)
  T3_ambiguous : EC maps to multiple MNXR (must pick the right one by hand)
  T4_no_balanced : no balanced MNXR available (likely reject)

NOTE: an earlier version tried to flag "all substrate/product metabolites
already in the model" as a tier signal.  That proved UNRELIABLE: model
metabolite MNXM annotations contain empty-shell IDs (e.g. CoA = MNXM1094981
with no InChIKey) and MetaNetX assigns multiple co-existing MNXM ids to the
same compound (CoA = MNXM727276 / 727277 / 12 / 1094981), so neither MNXM-id
nor chem_depr nor InChIKey matching resolves them.  Metabolite-coverage is
therefore left to manual inspection of the equation during T1 curation.

Signals come entirely from local data:
  - data/yali1_yali0_map/funnel4_reaction_candidates.csv  (the 90 candidates)
  - data/yali1_yali0_map/s2_metabolic_genes.csv           (pathway hits)
  - data/metanetx/reac_prop.tsv                           (MNXR eqn/balance)

Output: data/yali1_yali0_map/funnel6_curation_tiers.csv
"""
import csv
import os
import re

import cobra

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EC_RE = re.compile(r"\d+\.\d+\.\d+\.\d+")
MNXM_RE = re.compile(r"(MNXM\d+)")

# pathways we consider high-priority (sphingolipid / lipid, tied to prior work)
PRIORITY_PW = {"00600", "00561", "00564", "00565", "01040", "00071", "00062"}


def load_candidates():
    p = os.path.join(ROOT, "data/yali1_yali0_map/funnel4_reaction_candidates.csv")
    return [r for r in csv.DictReader(open(p))
            if r["reaction_status"] == "new_reaction"]


def load_gene_pathways():
    """yali1 -> set(metab pathway codes), from stage-3 CSV."""
    p = os.path.join(ROOT, "data/yali1_yali0_map/s2_metabolic_genes.csv")
    m = {}
    for r in csv.DictReader(open(p)):
        m[r["yali1"]] = {x for x in r["metab_pathways"].split(";") if x}
    return m


def load_ec_to_mnxr_and_props():
    """Return (ec->set(mnxr), mnxr->(set(mnxm), is_balanced_bool))."""
    p = os.path.join(ROOT, "data/metanetx/reac_prop.tsv")
    ec2mnxr = {}
    mnxr_props = {}
    with open(p) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            mnxr, eqn, classifs, balanced = c[0], c[1], c[3], c[4]
            mnxms = set(MNXM_RE.findall(eqn))
            mnxr_props[mnxr] = (mnxms, balanced == "B")
            for ec in EC_RE.findall(classifs):
                ec2mnxr.setdefault(ec, set()).add(mnxr)
    return ec2mnxr, mnxr_props


def main():
    cands = load_candidates()
    gene_pw = load_gene_pathways()
    ec2mnxr, mnxr_props = load_ec_to_mnxr_and_props()
    print(f"new_reaction candidates: {len(cands)}")

    rows = []
    tier_counts = {}
    for r in cands:
        ecs = [e for e in r["ec"].split(";") if EC_RE.fullmatch(e)]
        mnxrs = set()
        for e in ecs:
            mnxrs |= ec2mnxr.get(e, set())
        mnxr_unique = len(mnxrs) == 1
        balanced_mnxrs = sorted(mx for mx in mnxrs
                                if mnxr_props.get(mx, (set(), False))[1])
        rep = (balanced_mnxrs[0] if balanced_mnxrs
               else (sorted(mnxrs)[0] if mnxrs else ""))
        rep_balanced = mnxr_props.get(rep, (set(), False))[1]

        pw = gene_pw.get(r["yali1"], set())
        is_priority = bool(pw & PRIORITY_PW)

        # tier logic — does NOT use metabolite coverage (unreliable; see header)
        if not balanced_mnxrs:
            tier, note = "T4_no_balanced", "no balanced MNXR"
        elif mnxr_unique:
            tier, note = "T1_clear", "exactly one balanced MNXR"
        elif is_priority:
            tier = "T2_priority"
            note = f"priority pathway; EC->{len(mnxrs)} MNXR ({len(balanced_mnxrs)} balanced)"
        else:
            tier = "T3_ambiguous"
            note = f"EC->{len(mnxrs)} MNXR ({len(balanced_mnxrs)} balanced)"

        tier_counts[tier] = tier_counts.get(tier, 0) + 1
        rows.append({
            "yali1": r["yali1"],
            "yali0": r["yali0"],
            "ec": r["ec"],
            "n_mnxr": len(mnxrs),
            "n_balanced_mnxr": len(balanced_mnxrs),
            "rep_mnxr": rep,
            "mnxr_unique": "yes" if mnxr_unique else "no",
            "rep_balanced": "yes" if rep_balanced else "no",
            "priority_pathway": "yes" if is_priority else "no",
            "pathways": ";".join(sorted(pw)),
            "tier": tier,
            "note": note,
        })

    out = os.path.join(ROOT, "data/yali1_yali0_map/funnel6_curation_tiers.csv")
    order = {"T1_clear": 0, "T2_priority": 1, "T3_ambiguous": 2,
             "T4_no_balanced": 3}
    rows.sort(key=lambda x: (order[x["tier"]], x["priority_pathway"] != "yes",
                             x["yali1"]))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("\n=== TIER COUNTS ===")
    for t in ("T1_clear", "T2_priority", "T3_ambiguous", "T4_no_balanced"):
        print(f"  {t}: {tier_counts.get(t, 0)}")
    print(f"CSV written: {out}")


if __name__ == "__main__":
    main()
