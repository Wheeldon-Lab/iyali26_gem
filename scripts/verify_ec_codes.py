#!/usr/bin/env python3
"""
verify_ec_codes.py — cross-check every reaction's EC annotation against the
authoritative ExPASy ENZYME database (enzyme.dat).  READ-ONLY: produces a CSV
report, never modifies the model.

Catches the "auto-annotation matched the wrong EC" class of bug — the EC
analogue of the ceramide name-collision problem.  For each (reaction, EC):

  ok            — full EC, active, and the official name/aliases share at least
                  one word with the reaction name.
  transferred   — enzyme.dat says "Transferred entry: X" → EC is obsolete,
                  give the replacement.
  deleted       — enzyme.dat says "Deleted entry".
  not_found     — full 4-level EC but absent from enzyme.dat (illegal/outdated).
  preliminary   — preliminary EC (contains 'n', e.g. 1.1.1.n12); not in
                  enzyme.dat by design — informational, not a bug.
  partial       — partial EC (contains '-', e.g. 2.5.1.-); resolved to its
                  class name via enzclass.txt — informational, not a bug.
  name_mismatch — EC exists & active, but its official name/aliases share NO
                  word with the reaction name → SUSPECT (heuristic; e.g. IPC
                  synthase reaction annotated with a desaturase EC).

`suspect` = True for transferred / deleted / not_found / name_mismatch.

Usage:
    python scripts/verify_ec_codes.py [MODEL] [--out CSV]
        MODEL  default: model.xml
        --out  default: data/ec_verification.csv
"""
from __future__ import annotations

import argparse
import csv
import logging
import re
from collections import Counter
from pathlib import Path

from cobra.io import read_sbml_model

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
ENZYME_DAT = REPO_ROOT / "data" / "expasy" / "enzyme.dat"
ENZCLASS = REPO_ROOT / "data" / "expasy" / "enzclass.txt"

# Words that carry no discriminative meaning for the name-overlap heuristic.
_STOPWORDS = {
    "the", "a", "an", "of", "and", "or", "to", "with", "for", "in", "on",
    "ec", "enzyme", "protein", "putative", "probable", "reaction", "rxn",
    "yeast", "specific", "cytosolic", "mitochondrial", "c", "n", "type",
}
_TRANSFERRED_RE = re.compile(r"Transferred entry:\s*(.+?)\.?\s*$")


def _tokenize(text: str) -> set[str]:
    """Lowercase word tokens, minus stopwords and pure numbers."""
    toks = re.split(r"[^a-z0-9]+", text.lower())
    return {t for t in toks if t and t not in _STOPWORDS and not t.isdigit()}


def load_enzyme_dat(path: Path) -> dict[str, dict]:
    """
    Parse enzyme.dat into {EC: {names: set[str], status: str, transferred_to: str|None}}.

    status: 'active' | 'transferred' | 'deleted'
    """
    db: dict[str, dict] = {}
    ec = None
    names: list[str] = []
    status = "active"
    transferred_to = None

    def flush():
        if ec is not None:
            db[ec] = {
                "names": set(names),
                "status": status,
                "transferred_to": transferred_to,
            }

    with path.open(encoding="latin-1") as fh:
        for line in fh:
            tag = line[:2]
            body = line[5:].rstrip("\n").strip()
            if tag == "ID":
                flush()
                ec = body.strip()
                names, status, transferred_to = [], "active", None
            elif tag in ("DE", "AN") and ec is not None:
                if body.startswith("Transferred entry"):
                    status = "transferred"
                    m = _TRANSFERRED_RE.search(body)
                    if m:
                        transferred_to = m.group(1).strip()
                elif body.startswith("Deleted entry"):
                    status = "deleted"
                else:
                    names.append(body.rstrip("."))
    flush()
    return db


def load_enzclass(path: Path) -> dict[str, str]:
    """Parse enzclass.txt → {partial_EC: class_name}.  Normalises '1. 1. 1.-' → '1.1.1.-'."""
    classes: dict[str, str] = {}
    line_re = re.compile(r"^([\d ]+\.[\d \-]+\.[\d \-]+\.[\d \-]+)\s+(.+?)\.?\s*$")
    with path.open(encoding="latin-1") as fh:
        for line in fh:
            m = line_re.match(line)
            if m:
                key = m.group(1).replace(" ", "")
                classes[key] = m.group(2).strip()
    return classes


def normalize_ec(raw: str) -> str:
    """Strip 'ec-code/', 'EC ' prefixes and surrounding whitespace."""
    s = raw.strip()
    s = re.sub(r"^ec-code[:/]", "", s, flags=re.IGNORECASE)
    s = re.sub(r"^EC[: ]+", "", s, flags=re.IGNORECASE)
    return s.strip()


def classify(ec: str, rxn_name: str, db: dict, classes: dict) -> dict:
    """Return a classification dict for one (EC, reaction-name) pair."""
    out = {
        "ec": ec,
        "status": "",
        "enzyme_name": "",
        "note": "",
        "suspect": False,
    }
    # preliminary EC (contains a letter 'n' segment)
    if re.search(r"\bn\d+\b|\.n\d+", ec) or re.search(r"[a-zA-Z]", ec):
        out["status"] = "preliminary"
        out["note"] = "preliminary EC (not in enzyme.dat by design)"
        return out
    # partial EC (contains '-')
    if "-" in ec:
        out["status"] = "partial"
        out["enzyme_name"] = classes.get(ec, "")
        out["note"] = "partial EC → class name" if out["enzyme_name"] else "partial EC (no class match)"
        return out
    # full EC
    rec = db.get(ec)
    if rec is None:
        out["status"] = "not_found"
        out["note"] = "full EC absent from enzyme.dat (illegal/outdated)"
        out["suspect"] = True
        return out
    if rec["status"] == "transferred":
        out["status"] = "transferred"
        out["enzyme_name"] = f"→ {rec['transferred_to']}"
        out["note"] = f"obsolete; transferred to {rec['transferred_to']}"
        out["suspect"] = True
        return out
    if rec["status"] == "deleted":
        out["status"] = "deleted"
        out["note"] = "obsolete; deleted from ENZYME"
        out["suspect"] = True
        return out
    # active — name-overlap heuristic
    official = "; ".join(sorted(rec["names"]))
    out["enzyme_name"] = official
    ec_tokens = set()
    for nm in rec["names"]:
        ec_tokens |= _tokenize(nm)
    rxn_tokens = _tokenize(rxn_name)
    overlap = ec_tokens & rxn_tokens
    if rxn_tokens and ec_tokens and not overlap:
        out["status"] = "name_mismatch"
        out["note"] = "EC name shares no word with reaction name (heuristic — verify)"
        out["suspect"] = True
    else:
        out["status"] = "ok"
        out["note"] = f"matched on: {', '.join(sorted(overlap))}" if overlap else "no name to compare"
    return out


def get_ec_list(rxn) -> list[str]:
    """Extract ec-code annotation as a list of normalized EC strings."""
    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    raw = ann.get("ec-code", [])
    if isinstance(raw, str):
        raw = [raw]
    return [normalize_ec(x) for x in raw if x]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("model", nargs="?", default=str(REPO_ROOT / "model.xml"))
    ap.add_argument("--out", default=str(REPO_ROOT / "data" / "ec_verification.csv"))
    args = ap.parse_args()

    for p in (ENZYME_DAT, ENZCLASS):
        if not p.exists():
            logger.error(f"Missing required file: {p}")
            return

    logger.info(f"Loading enzyme.dat …")
    db = load_enzyme_dat(ENZYME_DAT)
    classes = load_enzclass(ENZCLASS)
    logger.info(f"  {len(db)} EC records, {len(classes)} class nodes")

    logger.info(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)

    rows = []
    counts: Counter = Counter()
    rxns_with_ec = 0
    for rxn in model.reactions:
        ecs = get_ec_list(rxn)
        if not ecs:
            continue
        rxns_with_ec += 1
        for ec in ecs:
            c = classify(ec, rxn.name or "", db, classes)
            counts[c["status"]] += 1
            rows.append({
                "reaction_id": rxn.id,
                "reaction_name": rxn.name or "",
                "ec_code": ec,
                "status": c["status"],
                "enzyme_dat_name": c["enzyme_name"],
                "note": c["note"],
                "suspect": c["suspect"],
            })

    rows.sort(key=lambda r: (not r["suspect"], r["status"], r["reaction_id"]))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "reaction_id", "reaction_name", "ec_code", "status",
            "enzyme_dat_name", "note", "suspect",
        ])
        w.writeheader()
        w.writerows(rows)

    n_suspect = sum(1 for r in rows if r["suspect"])
    logger.info("")
    logger.info(f"Reactions with EC: {rxns_with_ec} | (EC, reaction) pairs: {len(rows)}")
    logger.info("Status breakdown:")
    for status in ("ok", "partial", "preliminary", "name_mismatch",
                   "transferred", "deleted", "not_found"):
        if counts.get(status):
            mark = "  ⚠" if status in ("name_mismatch", "transferred", "deleted", "not_found") else "   "
            logger.info(f"{mark} {status:14} {counts[status]}")
    logger.info(f"\nSUSPECT (need review): {n_suspect}")
    logger.info(f"Report written: {out_path}")


if __name__ == "__main__":
    main()
