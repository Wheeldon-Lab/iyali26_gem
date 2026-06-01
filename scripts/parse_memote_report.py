"""
parse_memote_report.py — extract memote test results into CSVs for triage.

memote HTML reports embed a JSON object with one entry per test
(test_reaction_mass_balance, test_blocked_reactions, …).
Each test has: title, summary, metric, result, data (list of offender IDs).

This script produces two CSVs:

  data/memote_summary.csv
    One row per test. Columns: section, test_id, title, result, metric,
    metric_direction (lower_better / higher_better), n_offenders, message.

  data/memote_offenders.csv
    One row per (test, offender). Columns: section, test_id, title,
    offender_kind (reaction / metabolite / gene), offender_id,
    offender_name, offender_formula, offender_charge, offender_compartment.

Usage:
    # latest report by mtime
    python scripts/parse_memote_report.py

    # specific report
    python scripts/parse_memote_report.py results/20260528_162040.html

    # compare two reports → produces _diff.csv with delta (new - old)
    # IMPORTANT: pass OLD report first, then NEW. The first report is also
    # what summary/offenders CSVs are written from.
    python scripts/parse_memote_report.py results/20260528_162040.html \\
                                          --compare-against results/20260526_230202.html
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MODEL_PATH = REPO_ROOT / "model.xml"


# ─────────────────────────────────────────────────────────────
# Test classification: section + metric direction
# ─────────────────────────────────────────────────────────────

# Hand-curated mapping based on memote test_consistency.py / test_annotation.py / test_sbo.py / test_basic.py
TEST_SECTIONS = {
    # consistency
    "test_stoichiometric_consistency":              "consistency",
    "test_reaction_mass_balance":                   "consistency",
    "test_reaction_charge_balance":                 "consistency",
    "test_find_disconnected":                       "consistency",
    "test_find_metabolites_not_consumed_with_open_bounds": "consistency",
    "test_find_metabolites_not_produced_with_open_bounds": "consistency",
    "test_find_reactions_unbounded_flux_default_condition":"consistency",
    "test_find_stoichiometrically_balanced_cycles": "consistency",
    "test_blocked_reactions":                       "consistency",
    "test_metabolites_charge_presence":             "consistency",
    "test_metabolites_formula_presence":            "consistency",
    "test_unconserved_metabolites":                 "consistency",
    "test_inconsistent_min_stoichiometry":          "consistency",
    "test_find_orphans":                            "consistency",
    "test_find_deadends":                           "consistency",

    # annotation: metabolites
    "test_metabolite_annotation_presence":          "annotation_met",
    "test_metabolite_annotation_overview":          "annotation_met",
    "test_metabolite_annotation_conformity":        "annotation_met",
    "test_metabolite_id_namespace_consistency":     "annotation_met",

    # annotation: reactions
    "test_reaction_annotation_presence":            "annotation_rxn",
    "test_reaction_annotation_overview":            "annotation_rxn",
    "test_reaction_annotation_conformity":          "annotation_rxn",
    "test_reaction_id_namespace_consistency":       "annotation_rxn",

    # annotation: genes
    "test_gene_product_annotation_presence":        "annotation_gene",
    "test_gene_product_annotation_overview":        "annotation_gene",
    "test_gene_product_annotation_conformity":      "annotation_gene",

    # SBO
    "test_metabolite_sbo_presence":                 "annotation_sbo",
    "test_metabolite_specific_sbo_presence":        "annotation_sbo",
    "test_reaction_sbo_presence":                   "annotation_sbo",
    "test_metabolic_reaction_specific_sbo_presence":"annotation_sbo",
    "test_transport_reaction_specific_sbo_presence":"annotation_sbo",
    "test_exchange_specific_sbo_presence":          "annotation_sbo",
    "test_demand_specific_sbo_presence":            "annotation_sbo",
    "test_sink_specific_sbo_presence":              "annotation_sbo",
    "test_biomass_specific_sbo_presence":           "annotation_sbo",
    "test_gene_sbo_presence":                       "annotation_sbo",
    "test_gene_specific_sbo_presence":              "annotation_sbo",

    # basic / metadata
    "test_biomass_presence":                        "basic",
    "test_ngam_presence":                           "basic",
    "test_compartments_presence":                   "basic",
    "test_fbc_presence":                            "basic",
    "test_sbml_level":                              "basic",
    "test_model_id_presence":                       "basic",
    "test_metabolites_presence":                    "basic",
    "test_reactions_presence":                      "basic",
    "test_genes_presence":                          "basic",
    "test_gene_protein_reaction_rule_presence":     "basic",
    "test_protein_complex_presence":                "basic",
    "test_transport_reaction_gpr_presence":         "basic",
    "test_metabolic_coverage":                      "basic",
    "test_find_pure_metabolic_reactions":           "basic",
    "test_find_transport_reactions":                "basic",
    "test_find_unique_metabolites":                 "basic",
    "test_find_duplicate_metabolites_in_compartments": "basic",
    "test_find_duplicate_reactions":                "basic",
    "test_find_reactions_with_identical_genes":     "basic",
    "test_find_reactions_with_partially_identical_annotations": "basic",
    "test_find_constrained_pure_metabolic_reactions":  "basic",
    "test_find_constrained_transport_reactions":       "basic",
    "test_find_reversible_oxygen_reactions":        "basic",
    "test_find_medium_metabolites":                 "basic",
    "test_matrix_rank":                             "basic",
    "test_number_independent_conservation_relations":"basic",
    "test_degrees_of_freedom":                      "basic",
    "test_absolute_extreme_coefficient_ratio":      "basic",
}

# Tests where lower metric = better
LOWER_BETTER = {
    "test_reaction_mass_balance",
    "test_reaction_charge_balance",
    "test_find_disconnected",
    "test_find_metabolites_not_consumed_with_open_bounds",
    "test_find_metabolites_not_produced_with_open_bounds",
    "test_find_reactions_unbounded_flux_default_condition",
    "test_find_stoichiometrically_balanced_cycles",
    "test_blocked_reactions",
    "test_metabolites_charge_presence",
    "test_metabolites_formula_presence",
    "test_unconserved_metabolites",
    "test_find_duplicate_metabolites_in_compartments",
    "test_metabolite_annotation_presence",
    "test_metabolite_sbo_presence",
    "test_metabolite_specific_sbo_presence",
    "test_reaction_annotation_presence",
    "test_reaction_sbo_presence",
    "test_gene_product_annotation_presence",
    "test_gene_sbo_presence",
    "test_gene_specific_sbo_presence",
    "test_biomass_presence",
    "test_biomass_specific_sbo_presence",
    "test_demand_specific_sbo_presence",
    "test_exchange_specific_sbo_presence",
    "test_reaction_id_namespace_consistency",
    "test_metabolite_id_namespace_consistency",
    "test_find_orphans",
    "test_find_deadends",
}


# ─────────────────────────────────────────────────────────────
# Extract test blocks from HTML
# ─────────────────────────────────────────────────────────────

def _find_block(content: str, start: int) -> tuple[int, str]:
    """Return (end_index, block_string) for a JSON object starting at content[start]=='{'."""
    depth = 0
    for i in range(start, len(content)):
        c = content[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i + 1, content[start:i + 1]
    return start, ""


def extract_test_results(report_path: Path) -> dict[str, dict]:
    """
    Walk the memote HTML and extract every "test_<name>": { … } block,
    parse as JSON, return {test_name: dict}.
    """
    content = report_path.read_text()
    tests: dict[str, dict] = {}

    # Use a careful regex that won't false-match suite IDs.
    # memote nests tests at depth — we accept any "test_<word>": { occurrence.
    for match in re.finditer(r'"(test_\w+)"\s*:\s*\{', content):
        name = match.group(1)
        if name in tests:
            continue   # first occurrence wins
        start = match.end() - 1
        _, block = _find_block(content, start)
        if not block:
            continue
        try:
            parsed = json.loads(block)
        except json.JSONDecodeError:
            continue
        # heuristic: only keep blocks that look like real tests
        if not any(k in parsed for k in ("metric", "result", "data")):
            continue
        tests[name] = parsed

    return tests


# ─────────────────────────────────────────────────────────────
# Look up model context for offenders (reaction/metabolite name, etc)
# ─────────────────────────────────────────────────────────────

def load_model_context(model_path: Path) -> dict:
    """Load model and build O(1) lookup for reactions / metabolites / genes."""
    import warnings
    warnings.filterwarnings("ignore")
    from cobra.io import read_sbml_model

    model = read_sbml_model(str(model_path))
    rxn = {r.id: r for r in model.reactions}
    met = {m.id: m for m in model.metabolites}
    gene = {g.id: g for g in model.genes}
    # Also accept M_/R_ prefixed lookups defensively
    for prefix, table in (("M_", met), ("R_", rxn), ("G_", gene)):
        extra = {}
        for k in table:
            extra[prefix + k] = table[k]
        table.update(extra)
    return {"rxn": rxn, "met": met, "gene": gene}


def classify_offender(test_id: str) -> str:
    """Return 'reaction', 'metabolite', 'gene' based on what the test reports."""
    n = test_id.lower()
    if "reaction" in n or "blocked" in n or "stoichiometric" in n or "cycles" in n \
            or "rule" in n or "transport_reaction" in n:
        return "reaction"
    if "metabolite" in n or "orphan" in n or "deadend" in n or "consumed" in n or "produced" in n \
            or "disconnected" in n or "duplicate_metabolites" in n or "unconserved" in n:
        return "metabolite"
    if "gene" in n:
        return "gene"
    return "unknown"


# ─────────────────────────────────────────────────────────────
# Write CSVs
# ─────────────────────────────────────────────────────────────

def write_summary(tests: dict[str, dict], out_path: Path) -> None:
    """One row per test."""
    rows = []
    for test_id, payload in tests.items():
        section = TEST_SECTIONS.get(test_id, "uncategorized")
        direction = "lower_better" if test_id in LOWER_BETTER else "higher_better"
        metric = payload.get("metric")
        # metric may be a number, a dict (per-db breakdown), or missing
        if isinstance(metric, dict):
            metric_repr = json.dumps(metric, separators=(",", ":"))
        else:
            metric_repr = metric

        data = payload.get("data")
        if isinstance(data, list):
            n_off = len(data)
        elif isinstance(data, dict):
            n_off = sum(len(v) for v in data.values() if isinstance(v, list))
        else:
            n_off = ""

        # truncate message (may be str, dict, or missing)
        raw_msg = payload.get("message") or ""
        if isinstance(raw_msg, dict):
            raw_msg = json.dumps(raw_msg, separators=(",", ":"))
        msg = str(raw_msg).replace("\n", " ").strip()
        if len(msg) > 200:
            msg = msg[:200] + "…"

        rows.append({
            "section":          section,
            "test_id":          test_id,
            "title":            payload.get("title", ""),
            "result":           payload.get("result", ""),
            "metric":           metric_repr,
            "metric_direction": direction,
            "n_offenders":      n_off,
            "message":          msg,
        })

    rows.sort(key=lambda r: (r["section"], r["test_id"]))
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows):>4} rows → {out_path}")


def write_offenders(tests: dict[str, dict], out_path: Path,
                    model_ctx: dict | None) -> None:
    """One row per (test, offender_id)."""
    rows = []
    for test_id, payload in tests.items():
        data = payload.get("data")
        if not data:
            continue
        section = TEST_SECTIONS.get(test_id, "uncategorized")
        title = payload.get("title", "")
        kind = classify_offender(test_id)

        # `data` may be: list of IDs, list of [id, value] pairs, or dict
        if isinstance(data, dict):
            # e.g. per-database annotation breakdowns: {"kegg": [...], ...}
            for db, sub in data.items():
                if not isinstance(sub, list):
                    continue
                for item in sub:
                    rows.append(_make_offender_row(
                        section, test_id, title, kind,
                        item, model_ctx, sub_category=db,
                    ))
        elif isinstance(data, list):
            for item in data:
                rows.append(_make_offender_row(
                    section, test_id, title, kind,
                    item, model_ctx, sub_category="",
                ))

    if not rows:
        print(f"No offenders to write — skipping {out_path}")
        return

    fieldnames = [
        "section", "test_id", "title", "offender_kind", "sub_category",
        "offender_id", "offender_name", "offender_compartment",
        "offender_formula", "offender_charge", "offender_subsystem",
        "extra_value",
    ]
    rows.sort(key=lambda r: (r["section"], r["test_id"], r["offender_id"]))
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows):>4} rows → {out_path}")


def _make_offender_row(section, test_id, title, kind, item,
                       model_ctx, sub_category) -> dict:
    """Build one offender row, looking up name/formula/etc from the model when possible."""
    # item may be a string ID, or a [id, value] pair, or a dict
    offender_id = ""
    extra = ""
    if isinstance(item, str):
        offender_id = item
    elif isinstance(item, (list, tuple)) and len(item) >= 1:
        offender_id = str(item[0])
        if len(item) >= 2:
            extra = str(item[1])
    elif isinstance(item, dict):
        offender_id = item.get("id") or item.get("name") or json.dumps(item, separators=(",", ":"))
        # stash anything else
        rest = {k: v for k, v in item.items() if k not in ("id", "name")}
        if rest:
            extra = json.dumps(rest, separators=(",", ":"))
    else:
        offender_id = str(item)

    row = {
        "section": section,
        "test_id": test_id,
        "title": title,
        "offender_kind": kind,
        "sub_category": sub_category,
        "offender_id": offender_id,
        "offender_name": "",
        "offender_compartment": "",
        "offender_formula": "",
        "offender_charge": "",
        "offender_subsystem": "",
        "extra_value": extra,
    }

    if model_ctx is None:
        return row

    # Look up in the appropriate table
    if kind == "reaction":
        obj = model_ctx["rxn"].get(offender_id)
        if obj:
            row["offender_name"] = obj.name or ""
            row["offender_compartment"] = ",".join(sorted({m.compartment for m in obj.metabolites}))
            row["offender_subsystem"] = obj.subsystem or ""
    elif kind == "metabolite":
        obj = model_ctx["met"].get(offender_id)
        if obj:
            row["offender_name"] = obj.name or ""
            row["offender_compartment"] = obj.compartment or ""
            row["offender_formula"] = obj.formula or ""
            row["offender_charge"] = obj.charge if obj.charge is not None else ""
    elif kind == "gene":
        obj = model_ctx["gene"].get(offender_id)
        if obj:
            row["offender_name"] = obj.name or ""
    else:
        # try all three
        for k in ("rxn", "met", "gene"):
            obj = model_ctx[k].get(offender_id)
            if obj:
                row["offender_name"] = obj.name or ""
                break
    return row


def write_diff(tests_old: dict, tests_new: dict, out_path: Path) -> None:
    """Per-test metric delta between two reports."""
    rows = []
    all_keys = sorted(set(tests_old) | set(tests_new))
    for k in all_keys:
        old_metric = tests_old.get(k, {}).get("metric")
        new_metric = tests_new.get(k, {}).get("metric")
        if isinstance(old_metric, dict) or isinstance(new_metric, dict):
            continue   # skip multi-value metrics for diff
        try:
            o = float(old_metric); n = float(new_metric)
            delta = n - o
        except (TypeError, ValueError):
            o = n = delta = ""

        # Offender count change
        def _n_off(t):
            d = t.get("data") if t else None
            if isinstance(d, list): return len(d)
            return ""
        n_off_old = _n_off(tests_old.get(k, {}))
        n_off_new = _n_off(tests_new.get(k, {}))

        section = TEST_SECTIONS.get(k, "uncategorized")
        direction = "lower_better" if k in LOWER_BETTER else "higher_better"
        if delta == "" or not isinstance(delta, float):
            verdict = "n/a"
        elif abs(delta) < 1e-6:
            verdict = "unchanged"
        elif (direction == "lower_better" and delta < 0) or \
             (direction == "higher_better" and delta > 0):
            verdict = "improved"
        else:
            verdict = "regressed"

        rows.append({
            "section":   section,
            "test_id":   k,
            "direction": direction,
            "old_metric": o,
            "new_metric": n,
            "delta":      delta,
            "old_n_offenders": n_off_old,
            "new_n_offenders": n_off_new,
            "verdict":   verdict,
        })

    rows.sort(key=lambda r: (r["section"], r["test_id"]))
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows):>4} rows → {out_path}")


# ─────────────────────────────────────────────────────────────
# Driver
# ─────────────────────────────────────────────────────────────

def latest_report() -> Path:
    reports = sorted((REPO_ROOT / "results").glob("*.html"),
                     key=lambda p: p.stat().st_mtime, reverse=True)
    if not reports:
        raise FileNotFoundError("No reports found in results/")
    return reports[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", nargs="?", type=Path,
                        help="memote HTML report to analyze (default: latest in results/). "
                             "Summary and offenders CSVs are built from this report.")
    parser.add_argument("--compare-against", type=Path, default=None,
                        dest="baseline",
                        help="optional OLDER report to diff against. Produces memote_diff.csv "
                             "with delta = current - baseline.")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL_PATH,
                        help="model.xml to enrich offender rows with name/formula etc")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data",
                        help="output directory (default: data/)")
    args = parser.parse_args()

    report = args.report or latest_report()
    print(f"Parsing report: {report.name}")
    tests = extract_test_results(report)
    print(f"  Found {len(tests)} tests")

    model_ctx = None
    if args.model.exists():
        print(f"Loading model context: {args.model.name}")
        model_ctx = load_model_context(args.model)
    else:
        print(f"Warning: {args.model} not found — offender rows will lack model context")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_path   = args.out_dir / "memote_summary.csv"
    offenders_path = args.out_dir / "memote_offenders.csv"

    write_summary(tests, summary_path)
    write_offenders(tests, offenders_path, model_ctx)

    if args.baseline:
        baseline = args.baseline
        print(f"\nParsing baseline report for diff: {baseline.name}")
        tests_baseline = extract_test_results(baseline)
        diff_path = args.out_dir / "memote_diff.csv"
        # tests = current/new, tests_baseline = older
        write_diff(tests_baseline, tests, diff_path)


if __name__ == "__main__":
    main()
