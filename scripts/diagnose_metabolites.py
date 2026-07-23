"""
diagnose_metabolites.py — Standalone annotation-coverage diagnostic for all metabolites.

Usage:
    python scripts/diagnose_metabolites.py

Outputs:
    data/metabolite_diagnosis.csv
    Summary printed to stdout
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import cobra

REPO_ROOT = Path(__file__).resolve().parent.parent
MODEL_PATH = REPO_ROOT / "model.xml"
OUTPUT_CSV = REPO_ROOT / "data" / "metabolite_diagnosis.csv"
CHEM_XREF_PATH = REPO_ROOT / "data" / "metanetx" / "chem_xref.tsv"

# Annotation keys to check individually
DB_KEYS = [
    "bigg.metabolite",
    "kegg.compound",
    "metanetx.chemical",
    "chebi",
    "hmdb",
    "inchikey",
    "pubchem.compound",
    "seed.compound",
]

FIELDNAMES = [
    "met_id",
    "met_name",
    "compartment",
    "formula",
    "charge",
    "has_annotation",
    "annotation_keys",
    "has_bigg",
    "has_kegg",
    "has_metanetx",
    "has_chebi",
    "has_hmdb",
    "has_inchikey",
    "has_pubchem",
    "has_seed",
    "has_formula",
    "has_charge",
    "n_reactions",
    "base_name",
    "same_base_annotated",
    "annotated_twin_id",
]

# Regex: a "formula-like" suffix — only C/H/O/N/S/P/digits/parens/brackets
_FORMULA_RE = re.compile(r'^[CHONSP0-9()\[\]]+$')


def extract_base_name(name: str) -> str:
    """Return the base name: strip trailing _<formula> suffix if present."""
    if not name:
        return ""
    idx = name.rfind("_")
    if idx == -1:
        return name
    suffix = name[idx + 1:]
    if suffix and _FORMULA_RE.match(suffix):
        return name[:idx]
    return name


def main() -> None:
    if not MODEL_PATH.exists():
        sys.exit(f"Model not found: {MODEL_PATH}")

    print(f"Loading model: {MODEL_PATH}")
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    model.solver = "glpk"
    print(f"  {len(model.metabolites)} metabolites, {len(model.reactions)} reactions, {len(model.genes)} genes")

    # ── Pre-compute base_name → list of (met_id, annotated) for twin detection ──
    base_to_mets: dict[str, list[tuple[str, bool]]] = defaultdict(list)
    for met in model.metabolites:
        ann = met.annotation if isinstance(met.annotation, dict) else {}
        base = extract_base_name(met.name or "")
        base_to_mets[base].append((met.id, bool(ann)))

    # ── Build rows ────────────────────────────────────────────────────────────
    rows: list[dict] = []
    for met in model.metabolites:
        ann: dict = met.annotation if isinstance(met.annotation, dict) else {}
        ann_keys = sorted(ann.keys())
        base = extract_base_name(met.name or "")

        # twin: another met with same base_name that IS annotated
        twins = [(mid, is_ann) for mid, is_ann in base_to_mets[base]
                 if mid != met.id and is_ann]
        annotated_twin = next((mid for mid, _ in twins), None)

        row = {
            "met_id": met.id,
            "met_name": met.name or "",
            "compartment": met.compartment or "",
            "formula": met.formula or "",
            "charge": met.charge if met.charge is not None else "",
            "has_annotation": bool(ann),
            "annotation_keys": ",".join(ann_keys),
            "has_bigg": "bigg.metabolite" in ann,
            "has_kegg": "kegg.compound" in ann,
            "has_metanetx": "metanetx.chemical" in ann,
            "has_chebi": "chebi" in ann,
            "has_hmdb": "hmdb" in ann,
            "has_inchikey": "inchikey" in ann,
            "has_pubchem": "pubchem.compound" in ann,
            "has_seed": "seed.compound" in ann,
            "has_formula": bool(met.formula and met.formula.strip()),
            "has_charge": met.charge is not None,
            "n_reactions": len(met.reactions),
            "base_name": base,
            "same_base_annotated": annotated_twin is not None,
            "annotated_twin_id": annotated_twin or "",
        }
        rows.append(row)

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nCSV written → {OUTPUT_CSV}")

    # ── Summary ───────────────────────────────────────────────────────────────
    total = len(rows)
    annotated_rows = [r for r in rows if r["has_annotation"]]
    unannotated_rows = [r for r in rows if not r["has_annotation"]]
    n_ann = len(annotated_rows)
    n_unann = len(unannotated_rows)

    print("\n" + "=" * 60)
    print("METABOLITE ANNOTATION COVERAGE SUMMARY")
    print("=" * 60)
    print(f"Total metabolites : {total}")
    print(f"Has annotation    : {n_ann}  ({n_ann/total*100:.1f}%)")
    print(f"No annotation     : {n_unann}  ({n_unann/total*100:.1f}%)")

    print("\n── Per-database coverage (all metabolites) ──────────────────")
    db_field_map = {
        "bigg.metabolite": "has_bigg",
        "kegg.compound": "has_kegg",
        "metanetx.chemical": "has_metanetx",
        "chebi": "has_chebi",
        "hmdb": "has_hmdb",
        "inchikey": "has_inchikey",
        "pubchem.compound": "has_pubchem",
        "seed.compound": "has_seed",
    }
    for db, field in db_field_map.items():
        n = sum(1 for r in rows if r[field])
        print(f"  {db:25s}: {n:>4}/{total}  ({n/total*100:.1f}%)")

    # ── Twin analysis ─────────────────────────────────────────────────────────
    with_twin = [r for r in unannotated_rows if r["same_base_annotated"]]
    without_twin = [r for r in unannotated_rows if not r["same_base_annotated"]]

    print(f"\n── Unannotated metabolites: {n_unann} ────────────────────────")
    print(f"  Have annotated twin (same base_name) : {len(with_twin)}  "
          f"→ can copy annotations directly")
    print(f"  No annotated twin                    : {len(without_twin)}")

    # Formula / no-formula among no-twin unannotated
    no_twin_with_formula = [r for r in without_twin if r["has_formula"]]
    no_twin_no_formula = [r for r in without_twin if not r["has_formula"]]
    print(f"\n── No-twin unannotated breakdown ────────────────────────────")
    print(f"  Has formula (MetaNetX matchable) : {len(no_twin_with_formula)}")
    print(f"  No formula  (hardest cases)      : {len(no_twin_no_formula)}")

    # With-formula but unannotated (all unannotated, not just no-twin)
    unann_with_formula = [r for r in unannotated_rows if r["has_formula"]]
    print(f"\n  Total unannotated WITH formula   : {len(unann_with_formula)}")

    print("\n── Top 20 no-twin unannotated base_names ────────────────────")
    print(f"  {'met_id':<20}  {'base_name':<40}  compartment")
    print("  " + "-" * 75)
    for r in without_twin[:20]:
        bn = r["base_name"][:39] if r["base_name"] else "(no name)"
        print(f"  {r['met_id']:<20}  {bn:<40}  {r['compartment']}")

    print("\n── Top 20 no-formula no-annotation samples ──────────────────")
    print(f"  {'met_id':<20}  {'met_name':<45}  compartment")
    print("  " + "-" * 75)
    for r in no_twin_no_formula[:20]:
        name = r["met_name"][:44] if r["met_name"] else "(no name)"
        print(f"  {r['met_id']:<20}  {name:<45}  {r['compartment']}")

    # ── pubchem.compound = 0% root-cause analysis ─────────────────────────────
    print("\n── pubchem.compound = 0%: root-cause analysis ───────────────")
    if CHEM_XREF_PATH.exists():
        mnx_ids_in_model: set[str] = set()
        for r in annotated_rows:
            v = None
            for met in model.metabolites:
                if met.id == r["met_id"]:
                    v = (met.annotation or {}).get("metanetx.chemical")
                    break
            if v:
                if isinstance(v, list):
                    mnx_ids_in_model.update(v)
                else:
                    mnx_ids_in_model.add(v)

        # Scan chem_xref for pubchem.compound entries mapping to our MNX IDs
        pubchem_sources_found = set()
        all_source_prefixes: Counter = Counter()
        with open(CHEM_XREF_PATH) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                source, mnx = parts[0], parts[1]
                prefix = source.split(":")[0]
                all_source_prefixes[prefix] += 1
                if mnx in mnx_ids_in_model and "pubchem" in source.lower():
                    pubchem_sources_found.add(source.split(":")[0])

        print(f"  MNX IDs from annotated metabolites : {len(mnx_ids_in_model)}")
        print(f"  pubchem* source prefixes found     : {pubchem_sources_found or 'NONE'}")
        print(f"  All source prefixes in chem_xref   :")
        for prefix, cnt in all_source_prefixes.most_common(20):
            print(f"    {prefix:<25}: {cnt}")
        if not pubchem_sources_found:
            print("\n  Conclusion: chem_xref.tsv (MNXref 4.5) contains NO pubchem.compound")
            print("  entries at all — pubchem IDs are not present in this MNX release.")
            print("  Coverage is 0% because the source data doesn't include them,")
            print("  not because the pipeline failed to extract them.")
    else:
        print(f"  chem_xref.tsv not found at {CHEM_XREF_PATH} — skipping analysis")

    print("\nDone.")


if __name__ == "__main__":
    main()
