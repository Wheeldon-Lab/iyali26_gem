#!/usr/bin/env python3
"""Audit candidate metabolites that lack a formula but have a definite chemical
identity, by querying PubChem for the authoritative molecular formula.

Read-only: writes a comparison CSV only; does NOT touch model.xml.

Candidate set = metabolites with empty formula, in the "other" name category,
not generic (no _-suffix / R-group / polymer / (n)), and with NO usable
metanetx.chemical (the MNXM ids on these are empty-shell ids that carry no
formula in chem_prop, so PubChem-by-name is the route).

For each unique name:
  - normalize the name (model uses underscores like 'Gly_Glu',
    'Indole_3_acetonitrile' that PubChem won't resolve)
  - query PubChem REST: MolecularFormula + Charge
  - record the model's own charge for comparison (neutral PubChem formula vs
    the model's ionic state — flagged, NOT auto-applied)

Output: data/missing_formula_audit.csv
Columns: name, model_charge, n_copies, pubchem_formula, pubchem_charge,
         source_url, status
status: ok | charge_differs | not_found
"""
import csv
import os
import re
import urllib.parse
import urllib.request

import cobra

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def candidate_names(model):
    """unique name -> (list of met ids, model charge or None)."""
    def cat(x):
        n = (x.name or "").lower()
        if "trna" in n or "trna" in x.id.lower():
            return "tRNA"
        if any(k in n for k in ["glycogen", "starch", "chitin", "mannan",
                                "glucan", "dolichol", "dolichyl"]):
            return "polymer"
        if any(k in n for k in ["protein", "peptide", "apo-", "holo-",
                                "acyl-carrier", "-acp", "acp_"]):
            return "carrier"
        if (any(k in n for k in ["phosphatidate", "acyl-coa_",
                                 "cdp-diacylglycerol", "diacylglycerol",
                                 "lipid", "phosphatidyl"])
                and "coa" not in n.split("_")[0]):
            return "generic_lipid"
        return "other"

    def is_generic(x):
        n = x.name or ""
        return (n.endswith("_")
                or bool(re.search(r"\bR\d*\b|\(n\)|polymer|repeating", n))
                or bool(re.match(r"^(a|an) ", n.lower())))

    out = {}
    for x in model.metabolites:
        if x.formula:
            continue
        if cat(x) != "other" or is_generic(x):
            continue
        if x.annotation.get("metanetx.chemical"):
            continue
        nm = x.name or x.id
        if nm == x.id or nm.startswith("m2"):  # empty-shell (name == id)
            continue
        e = out.setdefault(nm, {"ids": [], "charge": None})
        e["ids"].append(x.id)
        if x.charge is not None:
            e["charge"] = x.charge
    return out


def normalize_for_pubchem(name):
    """Convert model name to a PubChem-resolvable form."""
    s = name
    # model uses underscores where standard names use hyphens/spaces
    s = s.replace("_", "-")
    # strip a trailing 'mitochondrial' style qualifier already handled by name
    return s.strip()


def pubchem_lookup(name):
    """Return (formula, charge, url) or (None, None, url)."""
    enc = urllib.parse.quote(name)
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc}"
           f"/property/MolecularFormula,Charge/JSON")
    human = f"https://pubchem.ncbi.nlm.nih.gov/#query={enc}"
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            import json
            d = json.load(resp)
        p = d["PropertyTable"]["Properties"][0]
        return p.get("MolecularFormula"), p.get("Charge"), human
    except Exception:
        return None, None, human


def main():
    m = cobra.io.read_sbml_model(os.path.join(ROOT, "model.xml"))
    cands = candidate_names(m)
    print(f"candidate unique names: {len(cands)}")

    rows = []
    counts = {"ok": 0, "charge_differs": 0, "not_found": 0}
    for nm, info in cands.items():
        query = normalize_for_pubchem(nm)
        formula, pcharge, url = pubchem_lookup(query)
        mcharge = info["charge"]
        if formula is None:
            status = "not_found"
        elif mcharge is not None and pcharge is not None and mcharge != pcharge:
            status = "charge_differs"
        else:
            status = "ok"
        counts[status] += 1
        rows.append({
            "name": nm,
            "query": query,
            "model_charge": mcharge if mcharge is not None else "",
            "n_copies": len(info["ids"]),
            "pubchem_formula": formula or "",
            "pubchem_charge": pcharge if pcharge is not None else "",
            "source_url": url,
            "status": status,
        })
        print(f"  [{status:13}] {nm[:40]:40} -> {formula or '-'} "
              f"(q={pcharge}, model={mcharge})")

    out = os.path.join(ROOT, "data/missing_formula_audit.csv")
    order = {"ok": 0, "charge_differs": 1, "not_found": 2}
    rows.sort(key=lambda r: (order[r["status"]], r["name"]))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "query", "model_charge",
                                          "n_copies", "pubchem_formula",
                                          "pubchem_charge", "source_url",
                                          "status"])
        w.writeheader()
        w.writerows(rows)

    print("\n=== STATUS COUNTS ===")
    for k in ("ok", "charge_differs", "not_found"):
        print(f"  {k}: {counts[k]}")
    print(f"CSV written: {out}")


if __name__ == "__main__":
    main()
