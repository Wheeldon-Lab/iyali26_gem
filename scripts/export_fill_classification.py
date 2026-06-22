#!/usr/bin/env python3
"""Persist the formula-fill-classifier agent's verdicts as a structured CSV.

Read-only over data/fill_candidates.csv; writes data/fill_classification.csv.
The agent's explicit per-ID verdicts (reject_wrong, reject_out_of_scope, and the
fillable / needs_review rows that have an authoritative formula + source) are recorded
below as data. Every other candidate row defaults to needs_review with a 'construct'
note (engine-built chain-resolved lipid whose mass balance is already code-verified).
Touches no model.
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_CSV = os.path.join(ROOT, "data", "fill_candidates.csv")
OUT_CSV = os.path.join(ROOT, "data", "fill_classification.csv")

# metabolite_id -> (verdict, authoritative_formula, source_url, confidence, note)
# from the formula-fill-classifier agent run (254 rows).
VERDICTS = {
    # ---- reject_wrong (7): candidate contradicts the authoritative formula ----
    "m2028[C_er]": ("reject_wrong", "C16H32O2", "https://www.kegg.jp/entry/C00249", "verified",
                    "palmitate back-solved to acetate (C2); name encodes C16"),
    "m1404[C_cy]": ("reject_wrong", "C10H7NO2", "https://www.kegg.jp/entry/C06325", "verified",
                    "quinaldate has 2 O, candidate has 3"),
    "m1956[C_cy]": ("reject_wrong", "C8H8O5", "https://www.kegg.jp/entry/C05580", "verified",
                    "3,4-dihydroxymandelate is 2 H short"),
    "m1891[C_cy]": ("reject_wrong", "C6H10O6", "https://www.kegg.jp/entry/C00198", "verified",
                    "gluconolactone written as open-chain hexose"),
    "m1808[C_cy]": ("reject_wrong", "C10H13N4O7P", "https://pubchem.ncbi.nlm.nih.gov/compound/dimp",
                    "verified", "2'-deoxyinosine-5'-MP is H-deficient by 5"),
    "m1975[C_cy]": ("reject_wrong", "C4H6O3", "https://www.kegg.jp/entry/C06002", "verified",
                    "methylmalonate semialdehyde 2 H short"),
    "m2026[C_cy]": ("reject_wrong", "C8H8O2", "https://www.kegg.jp/entry/C07086", "verified",
                    "phenylacetate 2 H short"),
    # ---- reject_out_of_scope (6): x100 aggregate IPC pool pseudo-formulas ----
    "m1995[C_cy]": ("reject_out_of_scope", "", "", "inferred", "x100 aggregate (N100, ~C4200) IPC pool"),
    "m1996[C_cy]": ("reject_out_of_scope", "", "", "inferred", "x100 aggregate IPC pool"),
    "m1997[C_cy]": ("reject_out_of_scope", "", "", "inferred", "x100 aggregate IPC pool"),
    "m1998[C_cy]": ("reject_out_of_scope", "", "", "inferred", "x100 aggregate IPC pool"),
    "m1999[C_cy]": ("reject_out_of_scope", "", "", "inferred", "x100 aggregate IPC pool"),
    "m2000[C_cy]": ("reject_out_of_scope", "", "", "inferred", "x100 aggregate IPC pool"),
    # ---- fillable with verified authoritative source (named small molecules) ----
    "m474[C_cy]":  ("fillable", "C5H11O7P", "https://www.kegg.jp/entry/C00672", "verified", ""),
    "m1816[C_cy]": ("fillable", "C5H4N4O4", "https://www.kegg.jp/entry/C11821", "verified", ""),
    "m2035[C_cy]": ("fillable", "C5H10N2O3S", "https://www.kegg.jp/entry/C01419", "verified", "Cys-Gly"),
    "m1827[C_cy]": ("fillable", "C4H8O4", "https://www.kegg.jp/entry/C01796", "verified", "erythrose"),
    "m1815[C_cy]": ("fillable", "C3H4O4", "https://www.kegg.jp/entry/C00168", "verified", "hydroxypyruvate"),
    "m2031[C_cy]": ("fillable", "C3H4O3S", "https://www.kegg.jp/entry/C00957", "verified", "mercaptopyruvate"),
    "m1814[C_cy]": ("fillable", "C3H5NO2", "https://www.kegg.jp/entry/C02218", "verified", "dehydroalanine"),
    "m1813[C_cy]": ("fillable", "C12H22O11", "https://www.kegg.jp/entry/C00185", "verified", "cellobiose"),
    "m1789[C_cy]": ("fillable", "C10H22O", "https://pubchem.ncbi.nlm.nih.gov/compound/1-decanol", "verified", "decanol"),
    "m1889[C_go]": ("fillable", "C15H24N2O17P2", "https://www.kegg.jp/entry/C00052", "verified", "UDP-D-galactose"),
    "m2020[C_cy]": ("fillable", "C15H28O7P2", "https://www.kegg.jp/entry/C00448", "verified", "farnesyl-PP"),
    "m1863[C_va]": ("fillable", "C2H5NO2", "https://www.kegg.jp/entry/C00037", "verified", "glycine"),
    "m1809[C_cy]": ("fillable", "C5H7NO3", "https://www.kegg.jp/entry/C01879", "verified", "5-oxoproline"),
    "m1981[C_cy]": ("fillable", "C18H39NO3", "https://www.kegg.jp/entry/C12144", "verified", "phytosphingosine"),
    # ---- needs_review: differs only by +/-H vs a verified DB formula (charge convention) ----
    "m1988[C_er]": ("needs_review", "C15H24N2O17P2", "https://www.kegg.jp/entry/C00029", "verified", "UDP-glucose, dH=-2 (charge)"),
    "m1810[C_cy]": ("needs_review", "C6H13NO5", "https://www.kegg.jp/entry/C00329", "verified", "D-glucosamine, dH=+1"),
    "m2048[C_cy]": ("needs_review", "C3H7O6P", "https://www.kegg.jp/entry/C00111", "verified", "DHAP, dH=-2 (charge)"),
    "m1949[C_cy]": ("needs_review", "C3H6O4", "https://www.kegg.jp/entry/C00258", "verified", "D-glycerate, dH=-1"),
    "m1723[C_cy]": ("needs_review", "C15H22N2O18P2", "https://www.kegg.jp/entry/C00167", "verified", "UDP-glucuronate, dH=-5 over-deprotonated"),
    "m1961[C_er]": ("needs_review", "C16H25N5O16P2", "https://www.kegg.jp/entry/C00096", "verified", "GDP-mannose, dH=-3"),
    # ---- needs_review: name<->formula mismatch suggesting MIS-ANNOTATION (not arithmetic) ----
    "m971[C_ex]":  ("needs_review", "C10H16N2O11P2", "https://www.kegg.jp/entry/C00363", "verified", "named TDP but has S+N4; likely mis-annotation, resolve identity"),
    "m975[C_cy]":  ("needs_review", "C10H16N2O11P2", "https://www.kegg.jp/entry/C00363", "verified", "named TDP but has S+N4; likely mis-annotation"),
    "m1116[C_mi]": ("needs_review", "C10H16N2O11P2", "https://www.kegg.jp/entry/C00363", "verified", "named TDP but has S+N4; likely mis-annotation"),
    "m976[C_cy]":  ("needs_review", "C10H17N2O14P3", "https://www.kegg.jp/entry/C00459", "verified", "named TTP but has S+N4; likely mis-annotation"),
    # ---- needs_review: same-name outlier formula (likely bad back-solve / mislabel) ----
    "m874[C_gm]":  ("needs_review", "", "", "inferred", "PE group is C41-42/O8; this row is C22H48NO5P (lyso-sized) outlier"),
}

LIPID_CONSTRUCT_NOTE = ("engine-constructed chain-resolved lipid; mass balance code-verified, "
                        "no single canonical DB formula")


def main():
    rows = list(csv.DictReader(open(IN_CSV)))
    out = []
    for r in rows:
        mid = r["metabolite_id"]
        if mid in VERDICTS:
            verdict, auth, src, conf, note = VERDICTS[mid]
        elif r["is_lipid"] == "yes":
            verdict, auth, src, conf, note = "needs_review", "", "", "inferred", LIPID_CONSTRUCT_NOTE
        else:
            verdict, auth, src, conf, note = "needs_review", "", "", "inferred", "no opened source; plausible by structure"
        out.append({
            "metabolite_id": mid, "name": r["name"], "candidate_formula": r["formula"],
            "is_lipid": r["is_lipid"], "verdict": verdict, "authoritative_formula": auth,
            "source_url": src, "confidence": conf, "note": note,
        })

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "metabolite_id", "name", "candidate_formula", "is_lipid", "verdict",
            "authoritative_formula", "source_url", "confidence", "note"])
        w.writeheader()
        order = {"reject_wrong": 0, "reject_out_of_scope": 1, "fillable": 2, "needs_review": 3}
        w.writerows(sorted(out, key=lambda x: (order[x["verdict"]], x["name"])))

    from collections import Counter
    c = Counter(x["verdict"] for x in out)
    print("wrote", OUT_CSV, "(%d rows)" % len(out))
    for v in ("reject_wrong", "reject_out_of_scope", "fillable", "needs_review"):
        print("  %-20s %d" % (v, c[v]))


if __name__ == "__main__":
    main()
