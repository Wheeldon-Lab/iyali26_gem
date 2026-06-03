#!/usr/bin/env python3
"""Funnel stage 4: map the 1087 metabolic-but-not-in-model candidates to
reactions and de-duplicate against reactions already in iYli21.

Read-only. No network. Inputs:
  - data/yali1_yali0_map/s2_metabolic_genes.csv  (verdict=metabolic, in_model=no)
  - data/metanetx/reac_prop.tsv                   (classifs col -> EC -> MNXR)
  - model.xml                                      (existing reaction EC/MNXR)

De-duplication judged on EC numbers (most reliable, fully local):
  - reaction availability: candidate EC exists in MetaNetX (has >=1 MNXR)
  - isozyme: candidate EC already among the model's reaction EC set

Per-candidate reaction_status:
  no_reaction  = no EC at all (reaction undefined by this data line)
  isozyme      = has EC, and EC already in model -> only needs GPR, not a new rxn
  new_reaction = has EC, EC maps to a MetaNetX MNXR, and EC NOT in model
  ec_no_mnxr   = has EC, but EC not found in MetaNetX reac_prop (rare; rxn unclear)

Output: data/yali1_yali0_map/funnel4_reaction_candidates.csv
"""
import csv
import os
import re

import cobra

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

EC_RE = re.compile(r"\d+\.\d+\.\d+\.\d+")


def load_candidates():
    """Candidates = metabolic, not in model, AND mapped to a metabolic
    pathway (metab_pathways non-empty).

    The stage-3 verdict used "metabolic pathway OR has EC". Sampling the
    "EC but no metabolic pathway" half showed it is dominated by GEM-boundary
    enzymes (DNA topoisomerases, ubiquitin ligases, proteasome peptidases,
    signaling kinases/phosphatases) -- these have EC numbers but are not
    metabolic reactions. Requiring a metabolic-pathway hit removes that
    contamination. Returns (strict_candidates, n_dropped_ec_only).
    """
    path = os.path.join(ROOT, "data/yali1_yali0_map/s2_metabolic_genes.csv")
    strict, ec_only = [], 0
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["verdict"] == "metabolic" and r["in_model"] == "no":
                if r["metab_pathways"]:
                    strict.append(r)
                else:
                    ec_only += 1
    return strict, ec_only


def load_ec_to_mnxr():
    """Full 4-segment EC -> set(MNXR) from reac_prop classifs column."""
    path = os.path.join(ROOT, "data/metanetx/reac_prop.tsv")
    m = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 or not c[3]:
                continue
            mnxr = c[0]
            for ec in EC_RE.findall(c[3]):
                m.setdefault(ec, set()).add(mnxr)
    return m


def load_model_ec_mnxr():
    """(model_ec set, model_mnxr set) from model.xml annotations."""
    m = cobra.io.read_sbml_model(os.path.join(ROOT, "model.xml"))
    ec, mnxr = set(), set()
    for r in m.reactions:
        a = r.annotation
        e = a.get("ec-code")
        if e:
            for x in (e if isinstance(e, list) else [e]):
                # keep only full 4-segment EC for clean comparison
                if EC_RE.fullmatch(x):
                    ec.add(x)
        mx = a.get("metanetx.reaction")
        if mx:
            for x in (mx if isinstance(mx, list) else [mx]):
                mnxr.add(x)
    return ec, mnxr


def main():
    cands, dropped_ec_only = load_candidates()
    ec2mnxr = load_ec_to_mnxr()
    model_ec, model_mnxr = load_model_ec_mnxr()

    print(f"candidates (metabolic-pathway, not in model): {len(cands)}")
    print(f"dropped (EC but no metabolic pathway = boundary enzymes): "
          f"{dropped_ec_only}")
    print(f"MetaNetX EC->MNXR: {len(ec2mnxr)} unique EC")
    print(f"model EC set: {len(model_ec)}; model MNXR set: {len(model_mnxr)}")

    rows = []
    counts = {"no_reaction": 0, "isozyme": 0, "new_reaction": 0, "ec_no_mnxr": 0}
    for r in cands:
        ec_field = r["ec"]
        ecs = [e for e in ec_field.split(";") if EC_RE.fullmatch(e)] if ec_field else []

        if not ecs:
            status = "no_reaction"
            mnxrs = set()
            overlap = ""
        else:
            # union of MNXRs reachable from this gene's ECs
            mnxrs = set()
            for e in ecs:
                mnxrs |= ec2mnxr.get(e, set())
            ec_in_model = [e for e in ecs if e in model_ec]
            if ec_in_model:
                status = "isozyme"
                overlap = ";".join(sorted(ec_in_model))
            elif mnxrs:
                status = "new_reaction"
                overlap = ""
            else:
                status = "ec_no_mnxr"
                overlap = ""
        counts[status] += 1

        rows.append({
            "yali1": r["yali1"],
            "yali0": r["yali0"],
            "ec": ec_field,
            "candidate_mnxr": ";".join(sorted(mnxrs)) if mnxrs else "",
            "reaction_status": status,
            "ec_already_in_model": overlap,
        })

    assert sum(counts.values()) == len(cands)

    out_path = os.path.join(ROOT, "data/yali1_yali0_map/funnel4_reaction_candidates.csv")
    order = {"new_reaction": 0, "isozyme": 1, "ec_no_mnxr": 2, "no_reaction": 3}
    rows.sort(key=lambda x: (order[x["reaction_status"]], x["yali1"]))
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["yali1", "yali0", "ec",
                                          "candidate_mnxr", "reaction_status",
                                          "ec_already_in_model"])
        w.writeheader()
        w.writerows(rows)

    print("\n=== REACTION STATUS COUNTS ===")
    for k in ("new_reaction", "isozyme", "ec_no_mnxr", "no_reaction"):
        print(f"  {k}: {counts[k]}")

    # how many DISTINCT new reactions (cluster new_reaction genes by EC)
    new_rows = [x for x in rows if x["reaction_status"] == "new_reaction"]
    new_ecs = set()
    for x in new_rows:
        for e in x["ec"].split(";"):
            if EC_RE.fullmatch(e) and e not in model_ec:
                new_ecs.add(e)
    print(f"\n=== NEW-REACTION subset ===")
    print(f"  genes: {len(new_rows)}")
    print(f"  distinct new EC numbers (not in model): {len(new_ecs)}")
    new_no_yali0 = [x for x in new_rows if not x["yali0"]]
    print(f"  CLIB89-specific (no YALI0) among them: {len(new_no_yali0)}")
    print(f"\nCSV written: {out_path}")


if __name__ == "__main__":
    main()
