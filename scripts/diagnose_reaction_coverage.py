"""
diagnose_reaction_coverage.py — Comprehensive reaction annotation coverage diagnosis.
Run: python scripts/diagnose_reaction_coverage.py
Output: data/reaction_coverage_diagnosis.csv + terminal summary
"""
import csv
import logging
from collections import Counter
from pathlib import Path

from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

MODEL_PATH = Path(__file__).resolve().parent.parent / "model.xml"
OUTPUT_CSV = Path(__file__).resolve().parent.parent / "data" / "reaction_coverage_diagnosis.csv"


def _rxn_type(rxn, exchanges: set, demands: set) -> str:
    if rxn in exchanges:
        return "exchange"
    if rxn in demands:
        return "demand"
    if len({m.compartment for m in rxn.metabolites}) >= 2:
        return "transport"
    return "metabolic"


def _bool(v) -> str:
    return "True" if v else "False"


def main():
    logging.info(f"Loading model: {MODEL_PATH}")
    model = read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"

    logging.info("Running FVA to identify blocked reactions …")
    blocked_ids: set[str] = set(find_blocked_reactions(model))
    logging.info(f"  {len(blocked_ids)} blocked reactions found")

    exchanges = set(model.exchanges)
    demands   = set(model.demands)

    rows = []
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        ann_keys = sorted(ann.keys())
        real_keys = [k for k in ann_keys if k != "sbo"]
        has_ann       = bool(ann)
        has_real_ann  = bool(real_keys)

        # Gene EC and UniProt info
        gene_ec_all: list[str] = []
        gene_has_uniprot = False
        for gene in rxn.genes:
            g_ann = gene.annotation if isinstance(gene.annotation, dict) else {}
            ec_raw = g_ann.get("ec-code", [])
            if isinstance(ec_raw, str):
                ec_raw = [ec_raw]
            gene_ec_all.extend(ec_raw)
            uni = g_ann.get("uniprot", [])
            if uni and (isinstance(uni, list) and uni or isinstance(uni, str) and uni):
                gene_has_uniprot = True

        gene_ec_dedup = list(dict.fromkeys(gene_ec_all))

        # Metabolite MNXM coverage
        mets = list(rxn.metabolites)
        mnxm_with = 0
        for met in mets:
            m_ann = met.annotation if isinstance(met.annotation, dict) else {}
            raw = m_ann.get("metanetx.chemical", "")
            val = raw[0] if isinstance(raw, list) else str(raw) if raw else ""
            if val:
                mnxm_with += 1
        n_mets = len(mets)
        met_mnxm_cov = round(mnxm_with / n_mets, 4) if n_mets else 0.0
        met_all_mnxm = mnxm_with == n_mets and n_mets > 0

        rows.append({
            "reaction_id":         rxn.id,
            "reaction_name":       rxn.name or "",
            "reaction_type":       _rxn_type(rxn, exchanges, demands),
            "has_annotation":      _bool(has_ann),
            "has_real_annotation": _bool(has_real_ann),
            "annotation_keys":     ",".join(ann_keys),
            "has_bigg":            _bool("bigg.reaction" in ann),
            "has_kegg":            _bool("kegg.reaction" in ann),
            "has_metanetx":        _bool("metanetx.reaction" in ann),
            "has_rhea":            _bool("rhea" in ann),
            "has_ec":              _bool("ec-code" in ann),
            "has_sbo":             _bool("sbo" in ann),
            "has_gpr":             _bool(bool(rxn.gene_reaction_rule and rxn.gene_reaction_rule.strip())),
            "n_genes":             len(rxn.genes),
            "gene_has_ec":         _bool(bool(gene_ec_dedup)),
            "gene_ec_codes":       "|".join(gene_ec_dedup),
            "gene_has_uniprot":    _bool(gene_has_uniprot),
            "n_metabolites":       n_mets,
            "met_mnxm_coverage":   met_mnxm_cov,
            "met_all_have_mnxm":   _bool(met_all_mnxm),
            "compartments":        ",".join(sorted({m.compartment for m in mets})),
            "is_blocked":          _bool(rxn.id in blocked_ids),
        })

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nCSV written → {OUTPUT_CSV}  ({len(rows)} rows)")

    # ── Terminal summary ───────────────────────────────────────────────────────
    total       = len(rows)
    annotated   = [r for r in rows if r["has_annotation"] == "True"]
    unannotated = [r for r in rows if r["has_annotation"] == "False"]
    has_real    = [r for r in rows if r["has_real_annotation"] == "True"]
    sbo_only    = [r for r in annotated if r["has_real_annotation"] == "False"]

    print(f"\n{'='*60}")
    print(f"REACTION ANNOTATION COVERAGE SUMMARY")
    print(f"{'='*60}")
    print(f"  Total reactions        : {total}")
    print(f"  Annotated (any)        : {len(annotated)}  ({100*len(annotated)/total:.1f}%)")
    print(f"  Has real annotation    : {len(has_real)}  ({100*len(has_real)/total:.1f}%)")
    print(f"  Unannotated            : {len(unannotated)}  ({100*len(unannotated)/total:.1f}%)")
    print(f"  Annotated but SBO-only : {len(sbo_only)}")

    print(f"\n── By reaction type ──────────────────────────────────────────")
    type_order = ["metabolic", "transport", "exchange", "demand"]
    type_total  = Counter(r["reaction_type"] for r in rows)
    type_real   = Counter(r["reaction_type"] for r in has_real)
    print(f"  {'Type':<12}  {'Total':>6}  {'Has real ann':>12}  {'Coverage':>8}")
    for t in type_order:
        tot = type_total.get(t, 0)
        ann = type_real.get(t, 0)
        pct = 100 * ann / tot if tot else 0.0
        print(f"  {t:<12}  {tot:>6}  {ann:>12}  {pct:>7.1f}%")

    print(f"\n── Unannotated breakdown  (n={len(unannotated)}) ──────────────────")
    u_gpr      = sum(1 for r in unannotated if r["has_gpr"] == "True")
    u_gene_ec  = sum(1 for r in unannotated if r["gene_has_ec"] == "True")
    u_uniprot  = sum(1 for r in unannotated if r["gene_has_uniprot"] == "True")
    u_all_mnxm = sum(1 for r in unannotated if r["met_all_have_mnxm"] == "True")
    print(f"  With GPR               : {u_gpr}")
    print(f"  Gene has EC number     : {u_gene_ec}")
    print(f"  Gene has UniProt       : {u_uniprot}")
    print(f"  All mets have MNXM     : {u_all_mnxm}")

    print(f"\n── met_mnxm_coverage distribution (unannotated) ─────────────")
    bins = {"100%": 0, "50–99%": 0, "1–49%": 0, "0%": 0}
    for r in unannotated:
        cov = float(r["met_mnxm_coverage"])
        if cov == 1.0:
            bins["100%"] += 1
        elif cov >= 0.5:
            bins["50–99%"] += 1
        elif cov > 0.0:
            bins["1–49%"] += 1
        else:
            bins["0%"] += 1
    for label, count in bins.items():
        pct = 100 * count / len(unannotated) if unannotated else 0.0
        print(f"  {label:<8} : {count:>4}  ({pct:.1f}%)")

    print(f"\n── met_mnxm_coverage distribution (all reactions) ───────────")
    bins_all = {"100%": 0, "50–99%": 0, "1–49%": 0, "0%": 0}
    for r in rows:
        cov = float(r["met_mnxm_coverage"])
        if cov == 1.0:
            bins_all["100%"] += 1
        elif cov >= 0.5:
            bins_all["50–99%"] += 1
        elif cov > 0.0:
            bins_all["1–49%"] += 1
        else:
            bins_all["0%"] += 1
    for label, count in bins_all.items():
        pct = 100 * count / total if total else 0.0
        print(f"  {label:<8} : {count:>4}  ({pct:.1f}%)")

    print(f"\n{'='*60}\n")


if __name__ == "__main__":
    main()
