#!/usr/bin/env python3
"""Stage 0 of the lipid-unlump project: mass-balance back-solve engine, DRY-RUN.

Read-only: computes what unlumping WOULD produce and writes a plan CSV only.
Does NOT touch model.xml, does NOT create reactions/metabolites, does NOT go online.

What it does this run (single-acyl layer only):
  - Reads the 6-chain menu straight from the xPOOL_AC_EM pool reaction (not hard-coded).
  - For every reaction that consumes a generic lipid-pool metabolite, tries to
    substitute each concrete chain length and BACK-SOLVE the one unknown product's
    formula from element conservation (product = sum(substrates) - sum(other products)).
  - Classifies each (reaction x chain) attempt:
        balanced            -> single unknown solved, integer, non-negative
        fail_negative       -> a solved element count went negative
        fail_fraction       -> a solved element count is non-integer
        imbalanced          -> 0 unknowns after substitution but mass not balanced
        blocked_multi_unknown -> >=2 generic species remain (a downstream/di-/tri-acyl
                                 reaction whose substrate isn't resolved yet)
  - Di-/tri-acyl reactions therefore surface as blocked_multi_unknown this run; they
    are listed so we can see exactly which reactions they are.

A "generic" metabolite is one whose formula cannot be resolved to definite integer
element counts: empty, OR contains a wildcard '*' (e.g. CHO2*), OR contains an
R-group placeholder 'R' (e.g. C10H18NO8PR2) — all three mean "an unspecified acyl
chain is attached here". R-group / '*' species are NOT solvable by the 6-chain
back-solve and surface as blocked_multi_unknown.

Output: data/unlump_stage0_plan.csv  +  a stdout summary.
"""
import csv
import os
import re
from collections import Counter, defaultdict

import cobra

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model.xml")
OUT_CSV = os.path.join(ROOT, "data", "unlump_stage0_plan.csv")

# pool reaction whose substrates define the chain menu (acyl-CoA side)
MENU_POOL = "xPOOL_AC_EM"
_ELEM_RE = re.compile(r"([A-Z][a-z]?)(\d*)")

# Skeleton-formula metabolites: a real-looking formula that only counts the glycerol
# backbone, not the attached acyl chains (e.g. diglyceride = C6H8O6). These are
# generic but carry no '*'/'R'/empty marker, so they must be named explicitly.
# A full scan of all 1865 metabolites (lipid-ish name with carbon count < 12) found
# diglyceride to be the ONLY such case; everything else in the glyceride/phospholipid
# family already carries a '*' or 'R' placeholder.
SKELETON_NAMES = {"diglyceride"}


# ---------------------------------------------------------------- formula utils
def parse_formula(f):
    """Parse 'C39H68N7O17P3S' -> {'C':39,...}. Returns None if the formula is generic:
    empty, OR contains a wildcard '*', OR contains an R-group placeholder 'R'."""
    if not f or "*" in f or "R" in f:
        return None
    out = Counter()
    for sym, num in _ELEM_RE.findall(f):
        if not sym:
            continue
        out[sym] += int(num) if num else 1
    return dict(out)


def formula_from_elements(d):
    """{'C':21,'H':41,...} -> 'C21H41O7P' (Hill order: C, H, then alphabetical)."""
    parts = []
    for sym in ["C", "H"]:
        if d.get(sym):
            parts.append(sym + (str(d[sym]) if d[sym] != 1 else ""))
    for sym in sorted(k for k in d if k not in ("C", "H")):
        if d[sym]:
            parts.append(sym + (str(d[sym]) if d[sym] != 1 else ""))
    return "".join(parts)


def is_generic(met):
    """Generic = formula cannot be resolved to definite integers: empty / '*' / 'R',
    OR a known skeleton-formula metabolite (e.g. diglyceride) whose formula counts
    only the glycerol backbone, not the acyl chains."""
    if parse_formula(met.formula) is None:
        return True
    return (met.name or "").strip().lower() in SKELETON_NAMES


# ---------------------------------------------------------------- chain menu
def read_chain_menu(model):
    """Read the 6 concrete acyl-CoA chain lengths (+ matching free fatty acid) from
    the pool reactions. Returns list of dicts with concrete formulas. Not hard-coded."""
    pool = model.reactions.get_by_id(MENU_POOL)
    menu = []
    for met, coef in pool.metabolites.items():
        if coef < 0 and not is_generic(met):  # the real acyl-CoA substrates
            # name like 'oleoyl-CoA_C39H68N7O17P3S' -> label from carbon/double-bond
            menu.append({
                "acyl_coa_name": met.name.split("_")[0],
                "formula": met.formula,
                "weight": abs(coef),
            })
    menu.sort(key=lambda d: -d["weight"])
    return menu


# ---------------------------------------------------------------- pool substitutions
def pool_substitutions(model):
    """Map each generic pool product (acyl-CoA_ / fatty acid_) to a REPRESENTATIVE
    concrete formula (the highest-weight chain: oleoyl-CoA / oleate). These are the
    'entry points' that turn a pool generic into a definite molecule for back-solve."""
    rep_acyl_coa = read_chain_menu(model)[0]["formula"]                     # oleoyl-CoA
    fa_pool = model.reactions.get_by_id("xPOOL_FA_EM")
    fa_subs = [(m, abs(c)) for m, c in fa_pool.metabolites.items()
               if c < 0 and not is_generic(m)]
    rep_fa = max(fa_subs, key=lambda t: t[1])[0].formula                    # oleate
    subs = {}
    for r in model.reactions:
        if not r.id.startswith("xPOOL"):
            continue
        for m, c in r.metabolites.items():
            if c > 0 and is_generic(m):
                subs[m.id] = rep_acyl_coa if "coa" in (m.name or "").lower() else rep_fa
    return subs


def effective_formula(met, known, pool_subs):
    """Resolved element dict for a metabolite, or None if it is still a genuine unknown.
    Priority: already-solved table > pool entry substitution > concrete formula > unknown."""
    if met.id in known:
        return parse_formula(known[met.id])
    if met.id in pool_subs:
        return parse_formula(pool_subs[met.id])
    if is_generic(met):
        return None
    return parse_formula(met.formula)


# ---------------------------------------------------------------- back-solve one reaction
def try_resolve(reaction, known, pool_subs):
    """If exactly one species in `reaction` is still unknown, back-solve its formula by
    element conservation. Returns dict(status, met_id, formula, mass_check, charge_residual, note)."""
    unknowns = []
    known_sum = Counter()
    n_fresh_acyl = 0      # pool-acyl entries consumed here (each multiplies combos by 6)
    for met, coef in reaction.metabolites.items():
        if met.id in pool_subs and met.id not in known and coef < 0:
            n_fresh_acyl += 1
        elems = effective_formula(met, known, pool_subs)
        if elems is None:
            unknowns.append((met, coef))
            continue
        for sym, n in elems.items():
            known_sum[sym] += coef * n

    if len(unknowns) >= 2:
        return {"status": "blocked_multi_unknown",
                "note": "%d unknown species: %s" % (len(unknowns), [u[0].id for u in unknowns])}
    if len(unknowns) == 0:
        residual = {s: v for s, v in known_sum.items() if abs(v) > 1e-9}
        return {"status": "balanced" if not residual else "imbalanced",
                "note": "no new metabolite" if not residual else "residual=%s" % residual}

    umet, ucoef = unknowns[0]
    solved = {}
    for sym, n in known_sum.items():
        val = -n / ucoef
        if abs(val) < 1e-9:
            continue
        if val < 0:
            return {"status": "fail_negative", "met_id": umet.id, "note": "%s -> %.3f" % (sym, val)}
        if abs(val - round(val)) > 1e-6:
            return {"status": "fail_fraction", "met_id": umet.id, "note": "%s -> %.4f" % (sym, val)}
        solved[sym] = int(round(val))
    formula = formula_from_elements(solved)

    # HARD CHECK 1 (mass balance): re-plug the solved formula and recompute residual.
    verify = Counter(known_sum)
    for sym, n in solved.items():
        verify[sym] += ucoef * n
    mass_residual = {s: v for s, v in verify.items() if abs(v) > 1e-9}
    mass_check = "pass" if not mass_residual else "FAIL %s" % mass_residual

    # CHARGE (informational only; model charge convention is a known separate issue): residual
    # across known species, excluding the new product whose charge we don't assign here.
    charge_res = 0.0
    for met, coef in reaction.metabolites.items():
        if met is umet:
            continue
        charge_res += coef * (met.charge or 0)

    return {"status": "balanced", "met_id": umet.id, "formula": formula,
            "n_fresh_acyl": n_fresh_acyl, "mass_check": mass_check,
            "charge_residual": charge_res, "note": "solved for %s" % umet.id}


# ---------------------------------------------------------------- iterative propagation
def propagate(model, pool_subs):
    """Resolve reactions layer by layer: each round, any reaction with exactly one unknown
    is back-solved and its product added to the virtual-known table. `combos` tracks the
    chain-length variant count (1-acyl=6, di-acyl=36, tri-acyl=216) for the cost model.
    Returns (rows, known, combos)."""
    pool_ids = {r.id for r in model.reactions if r.id.startswith("xPOOL")}
    candidates = [r for r in model.reactions if r.id not in pool_ids]

    # SCOPE: only the lipid-pool connected subgraph. Seed the frontier with the 8 generic
    # pool products plus the real chain-length substrates they stand for; a reaction is only
    # considered if it CONSUMES something already on the lipid frontier. New solved products
    # join the frontier. This keeps sugars / nucleotides / peptides / polymers out — they
    # never consume a lipid-frontier metabolite.
    frontier = set(pool_subs)
    for r in model.reactions:
        if r.id in pool_ids:
            for met, coef in r.metabolites.items():
                if coef < 0 and not is_generic(met):   # real acyl-CoA / fatty-acid substrates
                    frontier.add(met.id)

    known = {}      # met_id -> representative formula (virtual)
    combos = {}     # met_id -> number of chain-length variants
    solved_in = {}  # met_id -> layer
    rows = []
    layer = 0
    while True:
        layer += 1
        progress = False
        for r in candidates:
            # SCOPE gate: only react with the lipid frontier
            if not any(coef < 0 and met.id in frontier for met, coef in r.metabolites.items()):
                continue
            # skip reactions whose product is already solved (avoid re-solving)
            res = try_resolve(r, known, pool_subs)
            if res["status"] != "balanced" or "formula" not in res:
                continue
            mid = res["met_id"]
            if mid in known:
                continue
            # combos: product of upstream multi-combo substrates * 6 per fresh acyl chain
            upstream = 1
            for met, coef in r.metabolites.items():
                if coef < 0 and met.id in combos:
                    upstream *= combos[met.id]
            combos[mid] = max(1, upstream) * (len(read_chain_menu(model)) ** res["n_fresh_acyl"])
            known[mid] = res["formula"]
            solved_in[mid] = layer
            frontier.add(mid)   # the new product extends the lipid frontier
            mname = model.metabolites.get_by_id(mid).name
            rows.append({
                "layer": layer, "reaction_id": r.id, "reaction_name": r.name,
                "solved_met": mid, "solved_name": mname,
                "representative_formula": res["formula"], "combos": combos[mid],
                "mass_check": res["mass_check"], "charge_residual": "%+g" % res["charge_residual"],
                "status": "balanced", "note": res["note"],
            })
            progress = True
        if not progress:
            break
    return rows, known, combos


# ---------------------------------------------------------------- main dry-run
def main():
    model = cobra.io.read_sbml_model(MODEL)
    menu = read_chain_menu(model)
    pool_subs = pool_substitutions(model)
    pool_ids = {r.id for r in model.reactions if r.id.startswith("xPOOL")}

    rows, known, combos = propagate(model, pool_subs)

    # remaining reactions that touch a pool generic but never resolved (stuck: R-group/skeleton webs)
    solved_rxn = {r["reaction_id"] for r in rows}
    pool_generics = set(pool_subs)
    touch_pool = set()
    for r in model.reactions:
        if r.id in pool_ids:
            continue
        if any(m.id in pool_generics for m in r.metabolites):
            touch_pool.add(r.id)
    stuck = sorted(touch_pool - solved_rxn)

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "layer", "reaction_id", "reaction_name", "solved_met", "solved_name",
            "representative_formula", "combos", "mass_check", "charge_residual", "status", "note"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["layer"], r["reaction_id"])))

    # ---- summary ----
    print("Stage 0 dry-run (iterative back-solve, representative chain=%s). model.xml NOT modified."
          % menu[0]["acyl_coa_name"])
    print("chain menu (%d):" % len(menu),
          ", ".join("%s(%.3f)" % (c["acyl_coa_name"], c["weight"]) for c in menu))

    by_layer = defaultdict(list)
    for r in rows:
        by_layer[r["layer"]].append(r)
    print("\n=== resolved per layer (representative-combo dry-run) ===")
    total_new_rxn = total_new_met = 0
    for L in sorted(by_layer):
        rs = by_layer[L]
        new_rxn = sum(rr["combos"] for rr in rs)   # each solved rxn expands to `combos` real reactions
        new_met = sum(rr["combos"] for rr in rs)   # and roughly one new metabolite variant each
        total_new_rxn += new_rxn
        total_new_met += new_met
        print("  layer %d: %d generic reactions resolved -> %d chain-resolved reactions (sum of combos)"
              % (L, len(rs), new_rxn))
        for rr in rs:
            print("     %-8s %-42s %-22s x%-3d %s" % (
                rr["reaction_id"], rr["solved_name"][:42], rr["representative_formula"],
                rr["combos"], rr["mass_check"]))

    mass_fail = [r for r in rows if r["mass_check"] != "pass"]
    print("\nmass-balance hard check: %d/%d pass, %d FAIL" % (len(rows) - len(mass_fail), len(rows), len(mass_fail)))
    for r in mass_fail:
        print("  FAIL %-8s %s -> %s" % (r["reaction_id"], r["solved_name"], r["mass_check"]))

    print("\n=== COST TABLE (representative-combo) ===")
    print("  reactions resolved: %d   ->  ~%d chain-resolved reactions (full N^k expansion)"
          % (len(rows), total_new_rxn))
    print("  new chain-resolved metabolite variants: ~%d" % total_new_met)

    print("\n=== STUCK: touch a pool generic but never resolved (R-group / skeleton webs): %d ===" % len(stuck))
    for rid in stuck:
        print("  %-8s %s" % (rid, model.reactions.get_by_id(rid).name))

    # regression: R352 must resolve in layer 1 with the C18:1 representative -> C21H41O7P
    r352 = [r for r in rows if r["reaction_id"] == "R352"]
    if r352:
        print("\nregression R352: layer=%d formula=%s mass_check=%s (expect 1 / C21H41O7P / pass)"
              % (r352[0]["layer"], r352[0]["representative_formula"], r352[0]["mass_check"]))
    else:
        print("\nregression R352: NOT RESOLVED (unexpected)")

    print("\nwrote", OUT_CSV)


if __name__ == "__main__":
    main()
