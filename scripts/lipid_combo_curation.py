#!/usr/bin/env python3
"""Curate chain-length combinations by physiological abundance, deterministically.

Read-only: reads the real acyl-chain weights from the xPOOL_AC_EM pool reaction in
model.xml, enumerates the DISTINCT lipid species per acyl layer (mono/di/tri), collapsing
sn-position permutations into multisets (oleoyl+oleoyl+palmitoyl is one TAG regardless of
position), scores each by physiological probability prob = product(weights), and marks
keep/dead at a probability threshold. Writes data/lipid_combo_curation.csv. Touches no model.

This is the OBJECTIVE input for the lipid-combo-validator agent, which checks the weights
and the keep/dead split against the Y. lipolytica literature.
"""
import csv
import itertools
import os
from collections import defaultdict

import cobra

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")
OUT_CSV = os.path.join(ROOT, "data", "lipid_combo_curation.csv")
MENU_POOL = "xPOOL_AC_EM"
PROB_THRESHOLD = 1e-3   # keep combinations with physiological probability >= this


def read_weights(model):
    """Normalized acyl-chain weights from the pool reaction substrates."""
    pool = model.reactions.get_by_id(MENU_POOL)
    raw = []
    for met, c in pool.metabolites.items():
        if c < 0 and met.formula and "*" not in met.formula:
            raw.append((met.name.split("_")[0], abs(c)))
    tot = sum(w for _, w in raw)
    return sorted([(n, w / tot) for n, w in raw], key=lambda t: -t[1])


def distinct_species(weights, k):
    """Distinct k-acyl species (multisets) with summed probability over permutations."""
    prob = defaultdict(float)
    for combo in itertools.product(weights, repeat=k):
        names = tuple(sorted(x[0] for x in combo))
        p = 1.0
        for _, w in combo:
            p *= w
        prob[names] += p
    return sorted(prob.items(), key=lambda kv: -kv[1])


def main():
    model = cobra.io.read_sbml_model(MODEL)
    weights = read_weights(model)
    print("normalized acyl-chain weights (from %s):" % MENU_POOL)
    for n, w in weights:
        print("  %-14s %.4f" % (n, w))

    rows = []
    layers = [(1, "mono"), (2, "di"), (3, "tri")]
    for k, label in layers:
        species = distinct_species(weights, k)
        cum = 0.0
        for names, prob in species:
            cum += prob
            rows.append({
                "layer": label, "n_chains": k,
                "combination": "+".join(s.replace("-CoA", "") for s in names),
                "member_chains": ";".join(names),
                "prob": "%.6g" % prob,
                "cumulative_coverage": "%.4f" % cum,
                "verdict": "keep" if prob >= PROB_THRESHOLD else "dead",
            })

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "layer", "n_chains", "combination", "member_chains",
            "prob", "cumulative_coverage", "verdict"])
        w.writeheader()
        w.writerows(rows)

    # summary
    print("\nthreshold: keep prob >= %g" % PROB_THRESHOLD)
    for k, label in layers:
        lr = [r for r in rows if r["layer"] == label]
        keep = [r for r in lr if r["verdict"] == "keep"]
        cov = sum(float(r["prob"]) for r in keep)
        print("  %-4s : %3d distinct species -> keep %3d / dead %3d  (kept coverage %.3f%%)" % (
            label, len(lr), len(keep), len(lr) - len(keep), cov * 100))

    print("\nwrote", OUT_CSV)


if __name__ == "__main__":
    main()
