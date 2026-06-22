#!/usr/bin/env python3
"""Apply ONLY the C16:1 acyl-CoA pool extension to model.xml and write it back.

Adds palmitoleoyl-CoA (C16:1) to the 3 acyl-CoA pools (xPOOL_AC_EM/LP/MM) at 8.3%,
re-scaling the existing 6 chains so the substrate weight sum (= product coefficient,
0.951) is preserved. palmitoleoyl-CoA already exists in all 3 compartments, so no
metabolite/reaction is created. Guard assertions must all pass before the model is
written; otherwise nothing is saved.

Not wired into the main pipeline (run once, explicitly).
"""
import os

import cobra

from scripts.gem_annotate.patches import extend_acyl_pool_c161, _AC_POOLS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")


def pool_sub_sum(rxn):
    return sum(-c for m, c in rxn.metabolites.items() if c < 0)


def pool_product_coef(rxn):
    return [c for _, c in rxn.metabolites.items() if c > 0][0]


def main():
    m = cobra.io.read_sbml_model(MODEL)
    n_rxn, n_met, n_gene = len(m.reactions), len(m.metabolites), len(m.genes)

    extended = extend_acyl_pool_c161(m)
    print("pools extended: %d" % extended)

    # ── guard assertions (nothing is written unless all pass) ──
    assert len(m.reactions) == n_rxn, "reaction count changed!"
    assert len(m.metabolites) == n_met, "metabolite count changed (no new met expected)!"
    assert len(m.genes) == n_gene, "gene count changed!"
    for rid in _AC_POOLS:
        r = m.reactions.get_by_id(rid)
        ssum, prod = pool_sub_sum(r), pool_product_coef(r)
        assert abs(ssum - prod) < 1e-6, f"{rid}: substrate sum {ssum} != product {prod}"
        c161 = next((-c for mt, c in r.metabolites.items()
                     if c < 0 and "palmitoleoyl-coa" in (mt.name or "").lower()), 0.0)
        frac = c161 / ssum
        assert abs(frac - 0.083) < 0.002, f"{rid}: C16:1 fraction {frac} off"
        print("  %s ok: sub_sum=%.4f product=%.4f C16:1=%.1f%%" % (rid, ssum, prod, frac * 100))

    # FBA still feasible? (use GLPK — the Gurobi trial license is size-limited)
    m.solver = "glpk"
    sol = m.optimize()
    print("FBA objective after extension: %s (%s)" % (round(sol.objective_value or 0, 4), sol.status))
    assert sol.status == "optimal", "model no longer feasible!"

    cobra.io.write_sbml_model(m, MODEL)
    print("\nwrote", MODEL)

    # reload to confirm it persisted
    m2 = cobra.io.read_sbml_model(MODEL)
    for rid in _AC_POOLS:
        r = m2.reactions.get_by_id(rid)
        has = any("palmitoleoyl-coa" in (mt.name or "").lower() for mt in r.metabolites)
        print("  reload check %s: palmitoleoyl-CoA present = %s" % (rid, has))


if __name__ == "__main__":
    main()
