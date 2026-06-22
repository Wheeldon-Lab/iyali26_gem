#!/usr/bin/env python3
"""Deterministic filter over the unlump dry-run output: drop 'formal garbage' so a
downstream agent only spends effort verifying candidates worth checking.

Read-only: reads data/unlump_stage0_plan.csv, writes data/fill_candidates.csv. Touches
no model. Makes NO chemical-correctness judgement (that is the agent's job) — only
deterministic, rule-based rejection.

IDENTITY RULE: a metabolite's identity is its model ID (e.g. m575[C_em]), NEVER its name.
Two metabolites that share a name but have different IDs are DIFFERENT metabolites
(different compartments) and are handled independently. Conflict = the SAME ID solved to
two different non-empty formulas.

Reject categories (all deterministic):
  reject_empty        - back-solved formula is empty
  reject_nonformula   - not a parseable formula / zero carbons / junk token (e.g. 'CAP')
  reject_placeholder  - solved formula still contains a '*' or 'R' (not actually resolved)
  reject_out_of_scope - ion / polymer / biomass aggregate (cannot have a definite formula)
  reject_conflict     - same ID solved to >=2 different non-empty formulas
Everything else -> keep (written to fill_candidates.csv for agent verification).
"""
import csv
import os
import re
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN_CSV = os.path.join(ROOT, "data", "unlump_stage0_plan.csv")
OUT_CSV = os.path.join(ROOT, "data", "fill_candidates.csv")

_ELEM_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
# elements that legitimately appear in metabolite formulas here. Anything outside this
# set means the "formula" is junk that merely happens to parse (e.g. 'CAP' -> Ca + P).
VALID_ELEMENTS = {"C", "H", "N", "O", "P", "S", "Na", "K", "Mg", "Fe", "Cl", "Co", "Zn", "Ca", "Mn", "Cu", "Mo"}
# names that legitimately have no single integer formula (out of unlump scope)
OUT_OF_SCOPE = ["glucan", "pectin", "starch", "amylose", "cellulose", "biomass",
                "iron", "magnesium", "lipids_", "protein_"]
# lipid-identity hints (for the is_lipid flag only; does not affect keep/reject)
LIPID_HINT = ["acyl", "glycer", "phosphatid", "lipid", "ceramide", "sphing",
              "cardiolipin", "glyc", "fatty", "-coa", "coa", "enoyl", "oxo",
              "palmito", "oleo", "stearo", "behenate", "arachidate", "icos",
              "dodec", "tetradec", "hexadec", "docos", "inositol-p-ceramide"]


def parse_elements(f):
    """Parse a formula to element->count, or None if not a definite integer formula."""
    if not f or "*" in f or "R" in f:
        return None
    out = {}
    pos = 0
    for m in _ELEM_RE.finditer(f):
        sym, num = m.group(1), m.group(2)
        if not sym:
            continue
        out[sym] = out.get(sym, 0) + (int(num) if num else 1)
        pos += len(m.group(0))
    # the regex must have consumed the whole string, AND every symbol must be a real
    # element (else it is junk that merely parses, e.g. 'CAP' -> C + A + P).
    if pos != len(f) or not out:
        return None
    if any(sym not in VALID_ELEMENTS for sym in out):
        return None
    return out


def carbons(elems):
    return elems.get("C", 0) if elems else 0


def hint_match(name, hints):
    nm = (name or "").lower()
    return any(h in nm for h in hints)


def main():
    rows = list(csv.DictReader(open(IN_CSV)))

    # group every solved row by metabolite IDENTITY = ID (never name)
    by_id = defaultdict(list)
    for r in rows:
        by_id[r["solved_met"]].append(r)

    keep, rejects = [], []
    for mid, recs in by_id.items():
        name = recs[0]["solved_name"]
        formulas = [r["representative_formula"].strip() for r in recs]
        nonempty = sorted({f for f in formulas if f})
        src_rxns = sorted({r["reaction_id"] for r in recs})

        def rej(cat, detail):
            rejects.append((cat, mid, name, detail))

        # 1. out of scope by identity (ion/polymer/biomass) — before formula checks
        if hint_match(name, OUT_OF_SCOPE):
            rej("reject_out_of_scope", "name implies ion/polymer/aggregate")
            continue
        # 2. empty
        if not nonempty:
            rej("reject_empty", "no non-empty formula solved")
            continue
        # 3. conflict: same ID -> 2+ different non-empty formulas
        if len(nonempty) > 1:
            rej("reject_conflict", "same ID solved to %s" % nonempty)
            continue
        formula = nonempty[0]
        # 4. placeholder still present
        if "*" in formula or "R" in formula:
            rej("reject_placeholder", "formula still generic: %s" % formula)
            continue
        # 5. not a real formula / zero carbon
        elems = parse_elements(formula)
        if elems is None:
            rej("reject_nonformula", "unparseable: %r" % formula)
            continue
        if carbons(elems) == 0:
            rej("reject_nonformula", "zero carbon: %s" % formula)
            continue
        # keep
        keep.append({
            "metabolite_id": mid,
            "name": name,
            "formula": formula,
            "is_lipid": "yes" if hint_match(name, LIPID_HINT) else "no",
            "n_reactions_solved": len(src_rxns),
            "source_reactions": ";".join(src_rxns),
        })

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "metabolite_id", "name", "formula", "is_lipid",
            "n_reactions_solved", "source_reactions"])
        w.writeheader()
        w.writerows(sorted(keep, key=lambda r: (r["is_lipid"], r["name"])))

    # ---- summary (never silently drop: list every reject) ----
    total = len(by_id)
    print("fill-candidate filter (identity = metabolite ID, never name). Read-only.")
    print("distinct metabolite IDs solved: %d" % total)
    print("  keep              : %d" % len(keep))
    rej_by_cat = defaultdict(list)
    for cat, mid, name, detail in rejects:
        rej_by_cat[cat].append((mid, name, detail))
    for cat in sorted(rej_by_cat):
        print("  %-18s: %d" % (cat, len(rej_by_cat[cat])))

    print("\n=== kept candidates: %d (%d lipid, %d non-lipid) ===" % (
        len(keep), sum(k["is_lipid"] == "yes" for k in keep),
        sum(k["is_lipid"] == "no" for k in keep)))

    for cat in sorted(rej_by_cat):
        print("\n--- %s (%d) ---" % (cat, len(rej_by_cat[cat])))
        for mid, name, detail in sorted(rej_by_cat[cat]):
            print("  %-14s %-40s %s" % (mid, name[:40], detail))

    print("\nwrote", OUT_CSV)


if __name__ == "__main__":
    main()
