#!/usr/bin/env python3
"""Annotate the 15 isozyme genes added by the GPR step (they entered the model
with no annotation, regressing memote gene-SBO / gene-product-annotation).

For each new gene add:
  - sbo: SBO:0000243           (gene)                        [local]
  - ncbigene: <GeneID>         (from NCBI feature table)     [local]
  - kegg.genes: yli:<GeneID>   (verified present in yli)     [local]
  - uniprot: <accession>       (UniProt xref:geneid lookup)  [network, best-effort]

Reads model.xml, annotates only the genes listed in
data/gpr_isozyme_additions.csv (add_gene column), writes model.xml back.
Does not touch reactions, metabolites, or existing genes.
"""
import csv
import os
import urllib.parse
import urllib.request

import cobra

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")


def load_new_genes():
    p = os.path.join(ROOT, "data/gpr_isozyme_additions.csv")
    return sorted({r["add_gene"] for r in csv.DictReader(open(p))})


def load_geneid_map():
    """YALI1 (no underscore) -> GeneID, from NCBI feature table."""
    p = os.path.join(ROOT, "data/ncbi/clib89_feature_table.txt")
    m = {}
    with open(p) as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 17:
                continue
            locus, gid = c[16].strip(), c[15].strip()
            if locus and gid:
                m.setdefault(locus.replace("YALI1_", "YALI1"), gid)
    return m


def load_kegg_geneids():
    """Set of GeneIDs present in KEGG yli."""
    p = os.path.join(ROOT, "data/kegg/yli_genes.tsv")
    s = set()
    with open(p) as f:
        for line in f:
            gid = line.split("\t", 1)[0].replace("yli:", "")
            s.add(gid)
    return s


def fetch_uniprot(geneid):
    """Best-effort UniProt accession from a GeneID xref. None on failure."""
    q = urllib.parse.quote(f"xref:geneid-{geneid}")
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           f"&fields=accession&format=tsv&size=1")
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            lines = resp.read().decode().splitlines()
        if len(lines) >= 2 and lines[1].strip():
            return lines[1].strip()
    except Exception:
        pass
    return None


def main():
    m = cobra.io.read_sbml_model(MODEL)
    n_rxn, n_gene = len(m.reactions), len(m.genes)

    new_genes = load_new_genes()
    y2g = load_geneid_map()
    kegg_ids = load_kegg_geneids()

    print(f"annotating {len(new_genes)} new genes")
    counts = {"sbo": 0, "ncbigene": 0, "kegg.genes": 0, "uniprot": 0}
    for g in new_genes:
        gene = m.genes.get_by_id(g)
        ann = dict(gene.annotation) if gene.annotation else {}
        ann["sbo"] = "SBO:0000243"
        counts["sbo"] += 1

        gid = y2g.get(g)
        if gid:
            ann["ncbigene"] = gid
            counts["ncbigene"] += 1
            if gid in kegg_ids:
                ann["kegg.genes"] = f"yli:{gid}"
                counts["kegg.genes"] += 1
            up = fetch_uniprot(gid)
            if up:
                ann["uniprot"] = up
                counts["uniprot"] += 1
        gene.annotation = ann
        print(f"  {g}: ncbigene={ann.get('ncbigene','-')} "
              f"kegg={ann.get('kegg.genes','-')} uniprot={ann.get('uniprot','-')}")

    assert len(m.reactions) == n_rxn and len(m.genes) == n_gene, "counts changed!"
    cobra.io.write_sbml_model(m, MODEL)
    print(f"\nannotated: {counts}")
    print(f"model written: {MODEL}")


if __name__ == "__main__":
    main()
