#!/usr/bin/env python3
"""Audit reactions with EC-code overload (>=5 EC numbers) and compute a clean
EC set from the authoritative KEGG reaction-level ENZYME field.

Read-only: writes a before/after comparison CSV only, does NOT touch model.xml.

Cleaning rule (subset-only, never invents EC):
    new_ec = current_ec  INTERSECT  KEGG_reaction_ENZYME
Edge cases are left untouched and flagged for manual review:
    - no kegg.reaction id            -> action = skip_no_kegg
    - empty intersection             -> action = skip_empty_intersection
    - KEGG fetch failed              -> action = skip_fetch_error

Output: data/ec_overload_audit.csv
"""
import csv
import os
import re
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EC_RE = re.compile(r"\d+\.\d+\.\d+\.\d+")

import cobra


def fetch_kegg_enzyme(kegg_rxn_id):
    """Return set of EC numbers from KEGG reaction ENZYME field, or None on error."""
    url = f"https://rest.kegg.jp/get/rn:{kegg_rxn_id}"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as resp:
                text = resp.read().decode("utf-8", "replace")
            ecs = set()
            for line in text.splitlines():
                if line.startswith("ENZYME"):
                    ecs.update(EC_RE.findall(line))
            return ecs
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None


def ec_list(ann):
    e = ann.get("ec-code")
    if not e:
        return []
    return [x for x in (e if isinstance(e, list) else [e]) if EC_RE.fullmatch(x)]


def kegg_ids(ann):
    k = ann.get("kegg.reaction")
    if not k:
        return []
    out = []
    for x in (k if isinstance(k, list) else [k]):
        if x and x != "None" and x.startswith("R"):
            out.append(x)
    return out


def main():
    m = cobra.io.read_sbml_model(os.path.join(ROOT, "model.xml"))
    overload = [r for r in m.reactions if len(ec_list(r.annotation)) >= 5]
    print(f"reactions with >=5 EC: {len(overload)}")

    # collect unique KEGG ids, fetch once each
    all_kegg = set()
    for r in overload:
        all_kegg.update(kegg_ids(r.annotation))
    print(f"unique KEGG reaction ids to fetch: {len(all_kegg)}")

    kegg_ec = {}
    for i, kid in enumerate(sorted(all_kegg), 1):
        kegg_ec[kid] = fetch_kegg_enzyme(kid)
        n = "ERR" if kegg_ec[kid] is None else len(kegg_ec[kid])
        print(f"  [{i}/{len(all_kegg)}] {kid}: {n} EC")

    rows = []
    counts = {"clean": 0, "skip_no_kegg": 0, "skip_empty_intersection": 0,
              "skip_fetch_error": 0}
    for r in overload:
        cur = set(ec_list(r.annotation))
        kids = kegg_ids(r.annotation)
        if not kids:
            action, keep, drop = "skip_no_kegg", cur, set()
        else:
            kegg_union = set()
            fetch_failed = False
            for kid in kids:
                ke = kegg_ec.get(kid)
                if ke is None:
                    fetch_failed = True
                else:
                    kegg_union |= ke
            inter = cur & kegg_union
            if fetch_failed and not inter:
                action, keep, drop = "skip_fetch_error", cur, set()
            elif not inter:
                action, keep, drop = "skip_empty_intersection", cur, set()
            else:
                action, keep, drop = "clean", inter, cur - inter
        counts[action] += 1
        rows.append({
            "reaction": r.id,
            "name": r.name,
            "kegg": ";".join(kids),
            "n_current_ec": len(cur),
            "current_ec": ";".join(sorted(cur)),
            "kegg_ec": ";".join(sorted(set().union(
                *[kegg_ec.get(k) or set() for k in kids]))) if kids else "",
            "action": action,
            "keep_ec": ";".join(sorted(keep)),
            "drop_ec": ";".join(sorted(drop)),
        })

    out = os.path.join(ROOT, "data/ec_overload_audit.csv")
    order = {"clean": 0, "skip_empty_intersection": 1, "skip_fetch_error": 2,
             "skip_no_kegg": 3}
    rows.sort(key=lambda x: (order[x["action"]], x["reaction"]))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction", "name", "kegg",
                                          "n_current_ec", "current_ec",
                                          "kegg_ec", "action", "keep_ec",
                                          "drop_ec"])
        w.writeheader()
        w.writerows(rows)

    print("\n=== ACTION COUNTS ===")
    for k in ("clean", "skip_empty_intersection", "skip_fetch_error",
              "skip_no_kegg"):
        print(f"  {k}: {counts[k]}")
    n_drop = sum(len(r["drop_ec"].split(";")) if r["drop_ec"] else 0 for r in rows)
    print(f"\ntotal EC numbers that would be dropped: {n_drop}")
    print(f"CSV written: {out}")


if __name__ == "__main__":
    main()
