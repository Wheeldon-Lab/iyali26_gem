#!/usr/bin/env python3
"""Classify the 8765 S2-Table YALI1 genes by metabolic relevance.

Read-only analysis. Data lineage (no KEGG KO hop needed):

    S2 YALI1 ──(NCBI feature table)──► GeneID ──(KEGG yli)──► pathway / EC

Metabolic-relevance verdict (external evidence, not the model):
    metabolic     = gene maps to >=1 KEGG pathway in the 00xxx/01xxx range
                    (KEGG "Metabolism" super-class)  OR  has an EC number
    non_metabolic = only 03/04/05xxx pathways AND no EC
    no_data       = no pathway and no EC (KEGG didn't cover it; NOT a verdict
                    of "non-metabolic", just "this data line can't tell")

Also reports the subset that is metabolic AND not already in iYli21 -> the
real curation candidate pool.

Outputs: data/yali1_yali0_map/s2_metabolic_genes.csv
"""
import csv
import os
import re

import cobra
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_s2():
    """Return list of (yali1, yali0_or_empty) from S2 table."""
    path = os.path.join(ROOT, "data/yali1_yali0_map/S2_table_YALI1_YALI0_map.xlsx")
    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    # header: YALI0 locus ID | YALI1 Locus ID | ...
    out = []
    for r in rows[1:]:
        yali0 = (r[0] or "").strip()
        yali1 = (r[1] or "").strip()
        # Keep S2's underscored form (YALI1_A00014g): the NCBI feature table
        # uses the same underscored locus_tag, so the bridge matches directly.
        # Only the COBRA model uses the non-underscored form -> normalize there.
        if yali1:
            out.append((yali1, yali0))
    return out


def load_yali1_to_geneid():
    """YALI1 locus_tag -> GeneID, from NCBI feature table (cols 17,16)."""
    path = os.path.join(ROOT, "data/ncbi/clib89_feature_table.txt")
    m = {}
    with open(path) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 17:
                continue
            locus, gid = c[16].strip(), c[15].strip()
            if locus and gid:
                m.setdefault(locus, gid)
    return m


def load_geneid_pathways():
    """GeneID -> set(pathway codes like '00010')."""
    path = os.path.join(ROOT, "data/kegg/yli_pathway.tsv")
    m = {}
    with open(path) as f:
        for line in f:
            gene, pw = line.rstrip("\n").split("\t")
            gid = gene.replace("yli:", "")
            code = pw.replace("path:yli", "")
            m.setdefault(gid, set()).add(code)
    return m


def load_geneid_ec():
    """GeneID -> set(EC numbers)."""
    path = os.path.join(ROOT, "data/kegg/yli_ec.tsv")
    m = {}
    with open(path) as f:
        for line in f:
            gene, ec = line.rstrip("\n").split("\t")
            gid = gene.replace("yli:", "")
            m.setdefault(gid, set()).add(ec.replace("ec:", ""))
    return m


def load_model_genes():
    """Model gene ids, normalized to underscored form to match S2."""
    m = cobra.io.read_sbml_model(os.path.join(ROOT, "model.xml"))
    return {g.id.replace("YALI1", "YALI1_").replace("YALI0", "YALI0_")
            for g in m.genes}


def is_metabolic_pathway(code):
    """KEGG Metabolism super-class = pathway codes 00xxx and 01xxx."""
    return code.startswith("00") or code.startswith("01")


def main():
    s2 = load_s2()
    y2g = load_yali1_to_geneid()
    g2pw = load_geneid_pathways()
    g2ec = load_geneid_ec()
    model_genes = load_model_genes()

    print(f"S2 YALI1 genes: {len(s2)}")
    print(f"YALI1->GeneID bridge entries: {len(y2g)}")
    print(f"GeneID with pathways: {len(g2pw)}; GeneID with EC: {len(g2ec)}")
    print(f"model genes: {len(model_genes)}")

    rows = []
    counts = {"metabolic": 0, "non_metabolic": 0, "no_data": 0}
    no_bridge = 0
    for yali1, yali0 in s2:
        gid = y2g.get(yali1, "")
        pws = g2pw.get(gid, set()) if gid else set()
        ecs = g2ec.get(gid, set()) if gid else set()
        if not gid:
            no_bridge += 1

        metab_pws = {p for p in pws if is_metabolic_pathway(p)}
        has_ec = bool(ecs)

        if metab_pws or has_ec:
            verdict = "metabolic"
        elif pws:  # only non-metabolic pathways, no EC
            verdict = "non_metabolic"
        else:  # no pathway and no EC
            verdict = "no_data"
        counts[verdict] += 1

        rows.append({
            "yali1": yali1,
            "yali0": yali0,
            "geneid": gid,
            "metab_pathways": ";".join(sorted(metab_pws)),
            "ec": ";".join(sorted(ecs)),
            "verdict": verdict,
            "in_model": "yes" if yali1 in model_genes else "no",
        })

    assert sum(counts.values()) == len(s2), "verdict counts must sum to S2 total"

    # write CSV
    out_path = os.path.join(ROOT, "data/yali1_yali0_map/s2_metabolic_genes.csv")
    order = {"metabolic": 0, "non_metabolic": 1, "no_data": 2}
    rows.sort(key=lambda r: (order[r["verdict"]], r["in_model"], r["yali1"]))
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["yali1", "yali0", "geneid",
                                          "metab_pathways", "ec", "verdict",
                                          "in_model"])
        w.writeheader()
        w.writerows(rows)

    # report
    print("\n=== VERDICT COUNTS ===")
    for k in ("metabolic", "non_metabolic", "no_data"):
        print(f"  {k}: {counts[k]}")
    print(f"  (genes with no YALI1->GeneID bridge: {no_bridge})")

    metab = [r for r in rows if r["verdict"] == "metabolic"]
    metab_in = [r for r in metab if r["in_model"] == "yes"]
    metab_out = [r for r in metab if r["in_model"] == "no"]
    print("\n=== METABOLIC subset vs iYli21 ===")
    print(f"  metabolic total: {len(metab)}")
    print(f"  already in iYli21: {len(metab_in)}")
    print(f"  NOT in iYli21 (curation candidate pool): {len(metab_out)}")

    metab_out_no_yali0 = [r for r in metab_out if not r["yali0"]]
    print(f"    of which CLIB89-specific (no YALI0): {len(metab_out_no_yali0)}")
    print(f"\nCSV written: {out_path}")


if __name__ == "__main__":
    main()
