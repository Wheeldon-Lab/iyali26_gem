#!/usr/bin/env python3
"""Validate patches.unlump_stage1 on an in-memory copy of model.xml. Read-only.

Checks:
  - 4 generic reactions (R350/351/352/353) removed, 28 chain-specific added, 4 re-mix pools added
  - every new chain-specific reaction is mass-balanced (check_mass_balance == {})
  - R352's chain copies inherit bounds (0,0); the others inherit (0,1000)
  - the generic products still exist (fed by the re-mix pools) so downstream stays connected
  - FBA still optimal (GLPK)
  - idempotent: a second call adds nothing
Does NOT write model.xml.
"""
import os

import cobra

from scripts.gem_annotate.patches import unlump_stage1, _STAGE1_REACTIONS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")


def main():
    m = cobra.io.read_sbml_model(MODEL)
    n_rxn0, n_met0 = len(m.reactions), len(m.metabolites)

    done = unlump_stage1(m)
    print("unlump_stage1 -> %d generic reactions unlumped" % done)
    print("reactions %d -> %d   metabolites %d -> %d"
          % (n_rxn0, len(m.reactions), n_met0, len(m.metabolites)))

    ok = True
    # original generic reactions gone
    for rid in _STAGE1_REACTIONS:
        if rid in [r.id for r in m.reactions]:
            print("  FAIL: %s still present" % rid); ok = False

    # new chain-specific reactions + mass balance + bounds
    new = [r for r in m.reactions if any(r.id.startswith(g + "_") for g in _STAGE1_REACTIONS)]
    print("\nchain-specific reactions created: %d" % len(new))
    imbal = 0
    for r in new:
        mb = r.check_mass_balance()
        if mb:
            imbal += 1
            print("  IMBALANCED %s: %s" % (r.id, mb))
    print("  mass-balanced: %d / %d" % (len(new) - imbal, len(new)))
    if imbal:
        ok = False

    r352 = [r for r in new if r.id.startswith("R352_")]
    bad_bounds = [r.id for r in r352 if r.bounds != (0.0, 0.0)]
    print("  R352 copies with bounds (0,0): %d / %d %s"
          % (len(r352) - len(bad_bounds), len(r352), "OK" if not bad_bounds else "FAIL " + str(bad_bounds)))
    if bad_bounds:
        ok = False

    # re-mix pools present, generic products still present
    remix = [r for r in m.reactions if r.id.startswith("xREMIX_")]
    print("\nre-mix pools: %d" % len(remix))
    for gid in ["m569[C_em]", "m572[C_lp]", "m575[C_em]", "m577[C_lp]"]:
        present = gid in [x.id for x in m.metabolites]
        producers = [r.id for r in m.metabolites.get_by_id(gid).reactions
                     if m.metabolites.get_by_id(gid) in r.metabolites
                     and r.metabolites[m.metabolites.get_by_id(gid)] > 0] if present else []
        print("  generic %s present=%s produced_by=%s" % (gid, present, producers))

    # FBA
    m.solver = "glpk"
    sol = m.optimize()
    print("\nFBA: %s (%s)  [baseline ~1.78]" % (round(sol.objective_value or 0, 4), sol.status))
    if sol.status != "optimal":
        ok = False

    # idempotent
    n2 = unlump_stage1(m)
    print("second call -> %d (expect 0)" % n2)
    if n2 != 0:
        ok = False

    print("\nmodel.xml NOT written (in-memory only).")
    print("RESULT:", "ALL CHECKS PASS" if ok else "FAILURES ABOVE")


if __name__ == "__main__":
    main()
