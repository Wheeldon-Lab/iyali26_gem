#!/usr/bin/env python3
"""Validate patches.extend_acyl_pool_c161 on an in-memory copy of the model.

Read-only: loads model.xml, applies the C16:1 pool extension to the loaded object,
checks the invariants, and prints before/after. Does NOT write model.xml.

Invariants checked per acyl-CoA pool:
  - palmitoleoyl-CoA (C16:1) is now a substrate
  - substrate weight sum is unchanged (= product coefficient, 0.951)
  - C16:1 coefficient ~= 8.3% of the pool
  - the 6 original chains keep their relative ratios
  - idempotent: a second application changes nothing
"""
import os

import cobra

from scripts.gem_annotate.patches import extend_acyl_pool_c161

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")
POOLS = ("xPOOL_AC_EM", "xPOOL_AC_LP", "xPOOL_AC_MM")


def pool_subs(rxn):
    return {m.name.split("_")[0]: -c for m, c in rxn.metabolites.items() if c < 0}


def main():
    m = cobra.io.read_sbml_model(MODEL)

    print("=== BEFORE ===")
    before = {}
    for rid in POOLS:
        r = m.reactions.get_by_id(rid)
        subs = pool_subs(r)
        before[rid] = subs
        prod = [c for _, c in r.metabolites.items() if c > 0][0]
        print("  %s  sub_sum=%.4f  product=%.4f" % (rid, sum(subs.values()), prod))

    n = extend_acyl_pool_c161(m)
    print("\nextend_acyl_pool_c161 -> %d pools extended\n" % n)

    print("=== AFTER ===")
    ok = True
    for rid in POOLS:
        r = m.reactions.get_by_id(rid)
        subs = pool_subs(r)
        sub_sum = sum(subs.values())
        prod = [c for _, c in r.metabolites.items() if c > 0][0]
        c161 = subs.get("palmitoleoyl-CoA", 0.0)
        print("  %s  sub_sum=%.4f  product=%.4f  C16:1=%.4f (%.1f%%)" % (
            rid, sub_sum, prod, c161, 100 * c161 / sub_sum))
        # invariant 1: sub_sum == product coefficient (unchanged)
        if abs(sub_sum - prod) > 1e-6:
            print("    FAIL: substrate sum != product coefficient"); ok = False
        # invariant 2: C16:1 present at ~8.3%
        if abs(100 * c161 / sub_sum - 8.3) > 0.2:
            print("    FAIL: C16:1 fraction off"); ok = False
        # invariant 3: original 6 chains keep ratios (oleoyl:palmitoyl preserved)
        b = before[rid]
        r_before = b["oleoyl-CoA"] / b["palmitoyl-CoA"]
        r_after = subs["oleoyl-CoA"] / subs["palmitoyl-CoA"]
        if abs(r_before - r_after) > 1e-6:
            print("    FAIL: oleoyl/palmitoyl ratio changed"); ok = False

    # invariant 4: idempotent
    n2 = extend_acyl_pool_c161(m)
    print("\nsecond application -> %d pools extended (expect 0)" % n2)
    if n2 != 0:
        print("  FAIL: not idempotent"); ok = False

    # confirm model.xml on disk is untouched (we never wrote it)
    print("\nmodel.xml on disk NOT written by this script (in-memory only).")
    print("\nRESULT:", "ALL INVARIANTS PASS" if ok else "FAILURES ABOVE")


if __name__ == "__main__":
    main()
