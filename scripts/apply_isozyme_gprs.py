#!/usr/bin/env python3
"""Apply ONLY the curated isozyme GPR additions to the current model.xml.

Reads model.xml, runs patches.add_isozyme_gprs (which applies
data/gpr_isozyme_additions.csv), and writes model.xml back. Does not re-run
the full pipeline or touch the network. Only existing reactions' GPR get new
genes via 'or'; stoichiometry and existing genes are untouched.

The same step is wired into the main pipeline (patches.add_isozyme_gprs,
called last in main.py) so a future full rebuild keeps it.
"""
import os

import cobra

from scripts.gem_annotate.patches import add_isozyme_gprs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")


def main():
    m = cobra.io.read_sbml_model(MODEL)
    n_rxn, n_gene = len(m.reactions), len(m.genes)

    added = add_isozyme_gprs(m)
    print(f"GPR additions made: {added}")

    assert len(m.reactions) == n_rxn, "reaction count changed!"
    print(f"genes: {n_gene} -> {len(m.genes)} (+{len(m.genes) - n_gene})")

    cobra.io.write_sbml_model(m, MODEL)
    print(f"model written: {MODEL}")


if __name__ == "__main__":
    main()
