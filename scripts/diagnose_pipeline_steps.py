"""
diagnose_pipeline_steps.py — trace what each annotation-pipeline step does to a single reaction.

Loads the raw iYli21 model from data/iyli21.xml, then runs the pipeline step by
step, taking a snapshot of the chosen reaction (annotation, equation, metabolite
formulas, GPR genes) before and after each step.  At the end, prints a
per-step diff showing exactly what changed.

Usage
-----
    python scripts/diagnose_pipeline_steps.py R36
    python scripts/diagnose_pipeline_steps.py R36 --output results/r36_trace.txt
    python scripts/diagnose_pipeline_steps.py R36 --steps metabolites,reactions_4a,gene_ec
    python scripts/diagnose_pipeline_steps.py --list-steps
"""

from __future__ import annotations

import argparse
import logging
import sys
from copy import deepcopy
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from gem_annotate.config import CACHE_DIR, MNX_DIR, STARTING_MODEL_PATH

logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")


# Ordered list of pipeline steps mirroring main.py.
# Each entry: (short_name, human_label, callable(model, ctx) -> None)
# ctx is a dict carrying loaded MetaNetX tables / flags between steps.
PIPELINE_STEPS: list[tuple[str, str]] = [
    ("metabolites",        "Priority 1+2a: annotate_metabolites"),
    ("ho_balance_1",       "Priority 2b: H+/H2O balance (first pass)"),
    ("reactions_4a",       "Priority 4a: annotate_reactions"),
    ("merge_duplicates",   "Stoichiometric: merge duplicate metabolites"),
    ("exchange_bounds",    "Priority 2c: set_exchange_bounds"),
    ("medium",             "Priority 2c+: configure_medium"),
    ("biomass",            "Priority 3: fix_biomass_reaction"),
    ("genes",              "Priority 4b: annotate_genes"),
    ("idmapping",          "Priority 4c: ncbigene → uniprot id-mapping"),
    ("gene_ec",            "Priority 4d: gene EC enrichment"),
    ("reactions_4e",       "Priority 4e: annotate_remaining_reactions"),
    ("ec_backfill",        "EC backfill: gene ec-code → reaction"),
    ("xref_backfill",      "Reaction xref backfill (MNXR → bigg/kegg/rhea/ec)"),
    ("ho_balance_2",       "Priority 2b: H+/H2O balance (second pass)"),
    ("sbo",                "SBO term assignment"),
    ("normalize",          "Final: normalize_all_annotations"),
]


# ─────────────────────────────────────────────────────────────
# Snapshot helpers
# ─────────────────────────────────────────────────────────────

def snapshot(model, rxn_id: str) -> dict:
    """Take a snapshot of one reaction's state.  Returns None if rxn missing."""
    try:
        rxn = model.reactions.get_by_id(rxn_id)
    except KeyError:
        return None

    ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
    met_info = {}
    for met in rxn.metabolites:
        mann = met.annotation if isinstance(met.annotation, dict) else {}
        met_info[met.id] = {
            "coeff":   rxn.metabolites[met],
            "name":    met.name or "",
            "formula": met.formula or "",
            "charge":  met.charge,
            "bigg":    mann.get("bigg.metabolite", ""),
            "mnxm":    mann.get("metanetx.chemical", ""),
        }

    gene_info = {}
    for gene in rxn.genes:
        gann = gene.annotation if isinstance(gene.annotation, dict) else {}
        gene_info[gene.id] = {
            "uniprot": gann.get("uniprot", ""),
            "ec-code": gann.get("ec-code", ""),
            "ncbigene": gann.get("ncbigene", ""),
        }

    return {
        "rxn_id":     rxn.id,
        "rxn_name":   rxn.name,
        "equation":   rxn.reaction,
        "lower_bound": rxn.lower_bound,
        "upper_bound": rxn.upper_bound,
        "gpr":        rxn.gene_reaction_rule,
        "annotation": deepcopy(ann),
        "metabolites": met_info,
        "genes":      gene_info,
    }


def _norm_val(v):
    """Normalise an annotation value for comparison: list → sorted tuple, else str."""
    if v is None:
        return None
    if isinstance(v, list):
        return tuple(sorted(str(x) for x in v))
    return str(v)


def diff_snapshots(before: dict | None, after: dict | None) -> list[str]:
    """Return a list of human-readable lines describing what changed.
    Empty list = no changes."""
    if before is None and after is None:
        return ["  (reaction not in model in either snapshot)"]
    if before is None:
        return ["  ← reaction did NOT exist before this step"]
    if after is None:
        return ["  ← reaction REMOVED in this step"]

    lines = []

    # equation
    if before["equation"] != after["equation"]:
        lines.append(f"  equation:")
        lines.append(f"    before: {before['equation']}")
        lines.append(f"    after:  {after['equation']}")

    # bounds
    if before["lower_bound"] != after["lower_bound"] or before["upper_bound"] != after["upper_bound"]:
        lines.append(
            f"  bounds: [{before['lower_bound']}, {before['upper_bound']}] "
            f"→ [{after['lower_bound']}, {after['upper_bound']}]"
        )

    # GPR
    if before["gpr"] != after["gpr"]:
        lines.append(f"  GPR: {before['gpr']!r} → {after['gpr']!r}")

    # annotation diff
    b_ann = before["annotation"]
    a_ann = after["annotation"]
    b_keys = set(b_ann.keys())
    a_keys = set(a_ann.keys())

    added_keys = sorted(a_keys - b_keys)
    removed_keys = sorted(b_keys - a_keys)
    common = sorted(a_keys & b_keys)

    if added_keys:
        lines.append(f"  annotation: + {len(added_keys)} new key(s):")
        for k in added_keys:
            v = a_ann[k]
            disp = _trim_value(v)
            lines.append(f"    + {k} = {disp}")

    if removed_keys:
        lines.append(f"  annotation: - {len(removed_keys)} removed key(s):")
        for k in removed_keys:
            lines.append(f"    - {k} (was {_trim_value(b_ann[k])})")

    changed_keys = []
    for k in common:
        if _norm_val(b_ann[k]) != _norm_val(a_ann[k]):
            changed_keys.append(k)
    if changed_keys:
        lines.append(f"  annotation: ~ {len(changed_keys)} changed key(s):")
        for k in changed_keys:
            lines.append(f"    ~ {k}:")
            lines.append(f"        before: {_trim_value(b_ann[k])}")
            lines.append(f"        after:  {_trim_value(a_ann[k])}")

    # metabolite changes
    b_mets = set(before["metabolites"].keys())
    a_mets = set(after["metabolites"].keys())
    met_added = sorted(a_mets - b_mets)
    met_removed = sorted(b_mets - a_mets)
    met_common = sorted(a_mets & b_mets)

    if met_added:
        lines.append(f"  metabolites: + {len(met_added)} new:")
        for mid in met_added:
            mi = after["metabolites"][mid]
            lines.append(
                f"    + {mid} coeff={mi['coeff']:+g} formula={mi['formula'] or '?'} "
                f"bigg={mi['bigg'] or '-'}"
            )
    if met_removed:
        lines.append(f"  metabolites: - {len(met_removed)} removed:")
        for mid in met_removed:
            mi = before["metabolites"][mid]
            lines.append(f"    - {mid} (was coeff={mi['coeff']:+g})")

    # for common metabolites: did formula/charge/annotation upgrades happen?
    met_upgrade_lines = []
    for mid in met_common:
        b = before["metabolites"][mid]
        a = after["metabolites"][mid]
        sub_changes = []
        if b["formula"] != a["formula"]:
            sub_changes.append(f"formula: {b['formula'] or '?'!r} → {a['formula']!r}")
        if b["charge"] != a["charge"]:
            sub_changes.append(f"charge: {b['charge']} → {a['charge']}")
        if _norm_val(b["bigg"]) != _norm_val(a["bigg"]):
            sub_changes.append(f"bigg.metabolite: {_trim_value(b['bigg']) or '-'} → {_trim_value(a['bigg'])}")
        if _norm_val(b["mnxm"]) != _norm_val(a["mnxm"]):
            sub_changes.append(f"metanetx.chemical: {_trim_value(b['mnxm']) or '-'} → {_trim_value(a['mnxm'])}")
        if b["coeff"] != a["coeff"]:
            sub_changes.append(f"coeff: {b['coeff']:+g} → {a['coeff']:+g}")
        if sub_changes:
            met_upgrade_lines.append(f"    {mid}:")
            for sc in sub_changes:
                met_upgrade_lines.append(f"      {sc}")
    if met_upgrade_lines:
        lines.append(f"  metabolites: ~ upgrades:")
        lines.extend(met_upgrade_lines)

    # gene changes
    b_genes = set(before["genes"].keys())
    a_genes = set(after["genes"].keys())
    gene_added = sorted(a_genes - b_genes)
    gene_common = sorted(a_genes & b_genes)

    if gene_added:
        lines.append(f"  genes: + {len(gene_added)} added:")
        for gid in gene_added:
            lines.append(f"    + {gid}  {after['genes'][gid]}")

    gene_upgrade_lines = []
    for gid in gene_common:
        b = before["genes"][gid]
        a = after["genes"][gid]
        sub_changes = []
        for key in ("uniprot", "ec-code", "ncbigene"):
            if _norm_val(b[key]) != _norm_val(a[key]):
                sub_changes.append(
                    f"{key}: {_trim_value(b[key]) or '-'} → {_trim_value(a[key])}"
                )
        if sub_changes:
            gene_upgrade_lines.append(f"    {gid}:")
            for sc in sub_changes:
                gene_upgrade_lines.append(f"      {sc}")
    if gene_upgrade_lines:
        lines.append(f"  genes: ~ upgrades:")
        lines.extend(gene_upgrade_lines)

    return lines


def _trim_value(v, maxlen: int = 90) -> str:
    """Render an annotation value compactly; trim long lists."""
    if isinstance(v, list):
        s = ", ".join(str(x) for x in v[:4])
        if len(v) > 4:
            s += f", ... ({len(v)} total)"
        return f"[{s}]"
    s = str(v)
    if len(s) > maxlen:
        s = s[:maxlen] + "…"
    return s


# ─────────────────────────────────────────────────────────────
# Pipeline step runners (each takes model + ctx, modifies in place)
# ─────────────────────────────────────────────────────────────

def _ensure_mnx_loaded(ctx):
    """Lazily load MetaNetX tables on first use."""
    if "mnx_loaded" in ctx:
        return
    from gem_annotate.io import (
        load_chem_prop, load_chem_xref, load_mnxm_depr,
        load_reac_prop, load_reac_xref,
    )
    if not (MNX_DIR.exists() and (MNX_DIR / "chem_xref.tsv").exists()):
        ctx["mnx_loaded"] = False
        return
    ctx["chem_xref"]      = load_chem_xref(MNX_DIR / "chem_xref.tsv")
    ctx["chem_prop_data"] = load_chem_prop(MNX_DIR / "chem_prop.tsv")
    ctx["reac_xref"]      = load_reac_xref(MNX_DIR / "reac_xref.tsv")
    reac_prop_path = MNX_DIR / "reac_prop.tsv"
    ctx["reac_prop"] = (
        load_reac_prop(reac_prop_path) if reac_prop_path.exists() else None
    )
    ctx["mnxm_depr"] = load_mnxm_depr(MNX_DIR / "chem_depr.tsv")
    ctx["mnx_loaded"] = True


STEP_FUNCS = {}


def step(name):
    def deco(fn):
        STEP_FUNCS[name] = fn
        return fn
    return deco


@step("metabolites")
def _step_metabolites(model, ctx):
    _ensure_mnx_loaded(ctx)
    if not ctx.get("mnx_loaded"):
        print("    [SKIPPED — MetaNetX not available]")
        return
    from gem_annotate.metabolites import annotate_metabolites
    annotate_metabolites(model, ctx["chem_xref"], ctx["chem_prop_data"])


@step("ho_balance_1")
def _step_ho1(model, ctx):
    from gem_annotate.metabolites import fix_proton_water_balance
    fix_proton_water_balance(model)


@step("reactions_4a")
def _step_4a(model, ctx):
    _ensure_mnx_loaded(ctx)
    if not ctx.get("mnx_loaded"):
        print("    [SKIPPED — MetaNetX not available]")
        return
    from gem_annotate.reactions import annotate_reactions
    annotate_reactions(model, ctx["reac_xref"], ctx["reac_prop"])


@step("merge_duplicates")
def _step_merge(model, ctx):
    from gem_annotate.gaps import DUPLICATE_PAIRS, merge_duplicate_metabolites
    merge_duplicate_metabolites(model, DUPLICATE_PAIRS)


@step("exchange_bounds")
def _step_ex(model, ctx):
    from gem_annotate.exchange import set_exchange_bounds
    set_exchange_bounds(model)


@step("medium")
def _step_medium(model, ctx):
    from gem_annotate.exchange import configure_medium
    configure_medium(model)


@step("biomass")
def _step_biomass(model, ctx):
    from gem_annotate.biomass import fix_biomass_reaction
    fix_biomass_reaction(model)


@step("genes")
def _step_genes(model, ctx):
    from gem_annotate.genes import annotate_genes
    annotate_genes(model)


@step("idmapping")
def _step_idmap(model, ctx):
    from gem_annotate.idmapping import _enrich_via_idmapping
    _enrich_via_idmapping(model)


@step("gene_ec")
def _step_gene_ec(model, ctx):
    from gem_annotate.ec_annotation import enrich_genes_with_ec
    enrich_genes_with_ec(model)


@step("reactions_4e")
def _step_4e(model, ctx):
    _ensure_mnx_loaded(ctx)
    if not ctx.get("mnx_loaded"):
        print("    [SKIPPED — MetaNetX not available]")
        return
    from gem_annotate.annotate_reactions_extended import annotate_remaining_reactions
    annotate_remaining_reactions(model, ctx["reac_xref"], ctx["reac_prop"],
                                  mnxm_depr=ctx.get("mnxm_depr"))


@step("ec_backfill")
def _step_ec_backfill(model, ctx):
    """Inline EC backfill: copy gene EC numbers to reaction.ec-code (only when missing)."""
    count = 0
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        if "ec-code" in ann:
            continue
        ec_set = set()
        for gene in rxn.genes:
            g_ann = gene.annotation if isinstance(gene.annotation, dict) else {}
            ec_raw = g_ann.get("ec-code", [])
            if isinstance(ec_raw, str):
                ec_raw = [ec_raw]
            for ec in ec_raw:
                ec = ec.strip()
                if ec:
                    ec_set.add(ec)
        if ec_set:
            if not isinstance(rxn.annotation, dict):
                rxn.annotation = {}
            rxn.annotation["ec-code"] = sorted(ec_set)
            count += 1
    print(f"    ec_backfill: filled {count} reactions")


@step("xref_backfill")
def _step_xref_backfill(model, ctx):
    _ensure_mnx_loaded(ctx)
    if not ctx.get("mnx_loaded"):
        print("    [SKIPPED — MetaNetX not available]")
        return
    from gem_annotate.reactions import backfill_reaction_xrefs
    backfill_reaction_xrefs(model, ctx["reac_xref"], ctx["reac_prop"])


@step("ho_balance_2")
def _step_ho2(model, ctx):
    from gem_annotate.metabolites import fix_proton_water_balance
    fix_proton_water_balance(model)


@step("sbo")
def _step_sbo(model, ctx):
    """Inline SBO assignment mirroring main.py."""
    import re as _re
    _BIOMASS_RE     = _re.compile(r"BIOMASS|biomass|newBiom|R1372")
    _MAINTENANCE_RE = _re.compile(r"MAINTENANCE|ATPM")

    def _set_sbo(obj, term):
        ann = obj.annotation if isinstance(obj.annotation, dict) else {}
        if "sbo" in ann:
            return False
        if not isinstance(obj.annotation, dict):
            obj.annotation = {}
        obj.annotation["sbo"] = term
        return True

    for m in model.metabolites:
        _set_sbo(m, "SBO:0000247")
    for g in model.genes:
        _set_sbo(g, "SBO:0000243")
    _ex = set(model.exchanges); _dm = set(model.demands); _sk = set(model.sinks)
    for r in model.reactions:
        if r in _ex:   term = "SBO:0000627"
        elif r in _dm: term = "SBO:0000628"
        elif r in _sk: term = "SBO:0000632"
        elif _BIOMASS_RE.search(r.id):     term = "SBO:0000629"
        elif _MAINTENANCE_RE.search(r.id): term = "SBO:0000630"
        elif r.id.startswith(("xLIPID", "xAMINOACID", "xPOOL_")): term = "SBO:0000395"
        elif len({mt.compartment for mt in r.metabolites}) >= 2:  term = "SBO:0000185"
        else: term = "SBO:0000176"
        _set_sbo(r, term)


@step("normalize")
def _step_normalize(model, ctx):
    from gem_annotate.metabolites import normalize_all_annotations
    normalize_all_annotations(model)


# ─────────────────────────────────────────────────────────────
# Main driver
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rxn_id", nargs="?", help="Reaction ID to trace (e.g. R36)")
    parser.add_argument("--output", "-o", type=Path, default=None,
                        help="Write the trace to this file in addition to stdout")
    parser.add_argument("--steps", default=None,
                        help="Comma-separated step short-names to run (default: all). "
                             "Use --list-steps to see options.")
    parser.add_argument("--list-steps", action="store_true",
                        help="List available step short-names and exit")
    parser.add_argument("--model", default=str(STARTING_MODEL_PATH),
                        help="Starting SBML model (default: data/iyli21.xml)")
    parser.add_argument("--show-empty", action="store_true",
                        help="Print '(no changes)' for steps that didn't touch the reaction "
                             "(default: skip those steps in output)")
    args = parser.parse_args()

    if args.list_steps:
        print("Available pipeline steps (in execution order):")
        for sn, label in PIPELINE_STEPS:
            mark = " " if sn in STEP_FUNCS else "X"
            print(f"  [{mark}] {sn:18s}  {label}")
        return

    if not args.rxn_id:
        parser.error("rxn_id is required (use --list-steps to inspect available steps)")

    # Decide which steps to run
    if args.steps:
        chosen = [s.strip() for s in args.steps.split(",") if s.strip()]
        steps_to_run = [(sn, lbl) for sn, lbl in PIPELINE_STEPS if sn in chosen]
        missing = [s for s in chosen if s not in {sn for sn, _ in PIPELINE_STEPS}]
        if missing:
            print(f"WARNING: unknown step names ignored: {missing}", file=sys.stderr)
    else:
        steps_to_run = list(PIPELINE_STEPS)

    # Load model
    from cobra.io import read_sbml_model
    print(f"Loading model: {args.model}")
    model = read_sbml_model(args.model)
    print(f"  {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, "
          f"{len(model.genes)} genes")

    # Confirm reaction exists
    try:
        model.reactions.get_by_id(args.rxn_id)
    except KeyError:
        print(f"\nERROR: reaction {args.rxn_id!r} not found in {args.model}")
        sys.exit(1)

    # Collect output as a list of strings so we can dump to file AND stdout
    out_lines: list[str] = []

    def emit(line: str = ""):
        out_lines.append(line)
        print(line)

    # Initial snapshot
    ctx: dict = {}
    current = snapshot(model, args.rxn_id)
    emit("=" * 80)
    emit(f"TRACING: {args.rxn_id}  ({current['rxn_name']!r})")
    emit("=" * 80)
    emit("")
    emit("INITIAL STATE (raw iyli21.xml):")
    emit(f"  equation : {current['equation']}")
    emit(f"  bounds   : [{current['lower_bound']}, {current['upper_bound']}]")
    emit(f"  GPR      : {current['gpr'] or '(none)'}")
    emit(f"  annotation keys ({len(current['annotation'])}): "
         f"{sorted(current['annotation'].keys()) or '(empty)'}")
    emit(f"  metabolites ({len(current['metabolites'])}):")
    for mid, mi in current["metabolites"].items():
        emit(f"    {mid:14s} coeff={mi['coeff']:+g} formula={mi['formula'] or '?':12s} "
             f"bigg={_trim_value(mi['bigg']) or '-':14s} mnxm={_trim_value(mi['mnxm']) or '-'}")
    if current["genes"]:
        emit(f"  genes ({len(current['genes'])}):")
        for gid, gi in current["genes"].items():
            emit(f"    {gid:18s} uniprot={_trim_value(gi['uniprot']) or '-':12s} "
                 f"ec={_trim_value(gi['ec-code']) or '-'}")
    emit("")

    # Run each step, capture diff
    for idx, (sn, label) in enumerate(steps_to_run, 1):
        fn = STEP_FUNCS.get(sn)
        emit("-" * 80)
        emit(f"STEP {idx}/{len(steps_to_run)}: {sn}  —  {label}")
        emit("-" * 80)
        if fn is None:
            emit("  [not implemented]")
            continue

        before = current
        try:
            fn(model, ctx)
        except Exception as e:
            emit(f"  [ERROR running step: {type(e).__name__}: {e}]")
            continue

        after = snapshot(model, args.rxn_id)
        diff_lines = diff_snapshots(before, after)
        if diff_lines:
            for L in diff_lines:
                emit(L)
        else:
            if args.show_empty:
                emit("  (no changes to this reaction)")
            else:
                emit("  (no changes — skipped)")
        current = after
        emit("")

    # Final snapshot summary
    emit("=" * 80)
    emit("FINAL STATE:")
    emit("=" * 80)
    emit(f"  equation : {current['equation']}")
    emit(f"  GPR      : {current['gpr'] or '(none)'}")
    emit(f"  annotation keys ({len(current['annotation'])}):")
    for k in sorted(current["annotation"].keys()):
        emit(f"    {k:25s} = {_trim_value(current['annotation'][k])}")
    emit(f"  metabolites:")
    for mid, mi in current["metabolites"].items():
        emit(f"    {mid:14s} coeff={mi['coeff']:+g} formula={mi['formula'] or '?'}")
    if current["genes"]:
        emit(f"  genes:")
        for gid, gi in current["genes"].items():
            emit(f"    {gid:18s} uniprot={_trim_value(gi['uniprot']) or '-':12s} "
                 f"ec={_trim_value(gi['ec-code']) or '-'}")

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text("\n".join(out_lines))
        print(f"\n(trace also written to {args.output})")


if __name__ == "__main__":
    main()
