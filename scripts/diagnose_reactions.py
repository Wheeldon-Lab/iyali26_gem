"""
diagnose_reactions.py — Reaction annotation gap diagnosis.
Run: python scripts/diagnose_reactions.py
Output: data/reaction_diagnosis.csv + terminal summary
"""
import csv
import logging
from collections import Counter
from pathlib import Path

from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
MODEL_PATH = Path(__file__).resolve().parent.parent / "model.xml"
OUTPUT_CSV = Path(__file__).resolve().parent.parent / "data" / "reaction_diagnosis.csv"


def rxn_type(r, model):
    if r in model.exchanges:
        return "exchange"
    if r in model.demands:
        return "demand"
    comps = set(met.compartment for met in r.metabolites)
    if len(comps) >= 2:
        return "transport"
    return "metabolic"


def main():
    m = read_sbml_model(str(MODEL_PATH))
    m.solver = "glpk"
    blocked = set(find_blocked_reactions(m))

    rows = []
    for r in m.reactions:
        ann = r.annotation if isinstance(r.annotation, dict) else {}
        has_ann = bool(ann)
        rt = rxn_type(r, m)

        # Check if any gene has ec-code
        gene_ec_codes = []
        for g in r.genes:
            g_ann = g.annotation if isinstance(g.annotation, dict) else {}
            ec = g_ann.get("ec-code")
            if ec:
                if isinstance(ec, list):
                    gene_ec_codes.extend(ec)
                else:
                    gene_ec_codes.append(str(ec))

        # For exchange reactions, check if metabolite has bigg annotation
        met_bigg = ""
        if rt == "exchange" and len(r.metabolites) == 1:
            met = list(r.metabolites.keys())[0]
            m_ann = met.annotation if isinstance(met.annotation, dict) else {}
            raw = m_ann.get("bigg.metabolite", "")
            met_bigg = raw[0] if isinstance(raw, list) else str(raw) if raw else ""

        rows.append({
            "reaction_id": r.id,
            "reaction_name": r.name or "",
            "reaction_type": rt,
            "gene_reaction_rule": r.gene_reaction_rule or "",
            "n_genes": len(r.genes),
            "n_metabolites": len(r.metabolites),
            "compartments": ",".join(sorted(set(met.compartment for met in r.metabolites))),
            "has_annotation": has_ann,
            "annotation_keys": ",".join(sorted(ann.keys())) if ann else "",
            "has_bigg": "bigg.reaction" in ann,
            "has_kegg": "kegg.reaction" in ann,
            "has_metanetx": "metanetx.reaction" in ann,
            "has_rhea": "rhea" in ann,
            "has_ec": "ec-code" in ann,
            "has_sbo": "sbo" in ann,
            "gene_ec_codes": "|".join(sorted(set(gene_ec_codes))),
            "met_bigg_id": met_bigg,
            "is_blocked": r.id in blocked,
        })

    # Write CSV
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"CSV written to {OUTPUT_CSV} ({len(rows)} rows)")

    # Terminal summary
    unann = [r for r in rows if not r["has_annotation"]]
    ann = [r for r in rows if r["has_annotation"]]
    print(f"\nTotal: {len(rows)}, Annotated: {len(ann)} ({100*len(ann)/len(rows):.1f}%), Unannotated: {len(unann)} ({100*len(unann)/len(rows):.1f}%)")

    print("\nBy type:")
    type_counts = Counter(r["reaction_type"] for r in rows)
    type_unann = Counter(r["reaction_type"] for r in unann)
    for t in ["metabolic", "transport", "exchange", "demand"]:
        total = type_counts.get(t, 0)
        un = type_unann.get(t, 0)
        pct = 100 * un / total if total else 0
        print(f"  {t:12s}: {un}/{total} unannotated ({pct:.0f}%)")

    has_gpr = sum(1 for r in unann if r["gene_reaction_rule"].strip())
    has_ec_gene = sum(1 for r in unann if r["gene_ec_codes"].strip())
    has_name = sum(1 for r in unann if r["reaction_name"].strip())
    print(f"\nUnannotated breakdown:")
    print(f"  With GPR: {has_gpr}")
    print(f"  With gene that has EC number: {has_ec_gene}")
    print(f"  With reaction name: {has_name}")
    print(f"  No GPR and no name: {sum(1 for r in unann if not r['gene_reaction_rule'].strip() and not r['reaction_name'].strip())}")

    # Samples
    print(f"\n=== Unannotated metabolic + GPR + gene has EC (top candidates) ===")
    count = 0
    for r in unann:
        if r["reaction_type"] == "metabolic" and r["gene_ec_codes"]:
            print(f"  {r['reaction_id']:15s}  name={r['reaction_name'][:35]:35s}  EC={r['gene_ec_codes'][:30]}")
            count += 1
            if count >= 15:
                break

    print(f"\n=== Unannotated transport (top 15) ===")
    count = 0
    for r in unann:
        if r["reaction_type"] == "transport":
            print(f"  {r['reaction_id']:15s}  name={r['reaction_name'][:45]:45s}  GPR={r['gene_reaction_rule'][:25]}")
            count += 1
            if count >= 15:
                break

    print(f"\n=== Unannotated exchange with met_bigg_id (top 15) ===")
    count = 0
    for r in unann:
        if r["reaction_type"] == "exchange" and r["met_bigg_id"]:
            print(f"  {r['reaction_id']:15s}  met_bigg={r['met_bigg_id']}")
            count += 1
            if count >= 15:
                break


if __name__ == "__main__":
    main()
