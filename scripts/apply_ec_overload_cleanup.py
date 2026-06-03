#!/usr/bin/env python3
"""Apply ONLY the EC-overload cleanup to the current model.xml in place.

Reads model.xml, runs patches.clean_ec_overload (which applies the
action=clean rows from data/ec_overload_audit.csv), and writes model.xml
back. Does not re-run the full pipeline, does not touch the network, does
not change anything other than the curated reactions' ec-code.

The same cleanup is wired into the main pipeline (patches.clean_ec_overload,
called last in main.py) so a future full rebuild keeps it.
"""
import os

import cobra

from scripts.gem_annotate.patches import clean_ec_overload

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")


def ec_of(rxn):
    e = rxn.annotation.get("ec-code")
    if not e:
        return []
    return e if isinstance(e, list) else [e]


def main():
    m = cobra.io.read_sbml_model(MODEL)
    n_rxn, n_gene = len(m.reactions), len(m.genes)

    # snapshot the reactions the audit will touch, for before/after print
    import csv
    audit = os.path.join(ROOT, "data/ec_overload_audit.csv")
    targets = [r["reaction"] for r in csv.DictReader(open(audit))
               if r.get("action") == "clean"]
    before = {rid: list(ec_of(m.reactions.get_by_id(rid)))
              for rid in targets if rid in {r.id for r in m.reactions}}

    cleaned = clean_ec_overload(m)
    print(f"reactions cleaned: {cleaned}")

    # integrity checks: nothing else changed in counts
    assert len(m.reactions) == n_rxn, "reaction count changed!"
    assert len(m.genes) == n_gene, "gene count changed!"

    print("\n=== before -> after (cleaned reactions) ===")
    for rid in targets:
        if rid not in before:
            continue
        after = list(ec_of(m.reactions.get_by_id(rid)))
        if before[rid] != after:
            print(f"{rid}: {sorted(before[rid])} -> {sorted(after)}")

    cobra.io.write_sbml_model(m, MODEL)
    print(f"\nmodel written: {MODEL}")
    print(f"reactions={n_rxn} genes={n_gene} (unchanged)")


if __name__ == "__main__":
    main()
