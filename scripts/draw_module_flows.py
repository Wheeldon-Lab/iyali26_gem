#!/usr/bin/env python3
"""Draw a per-module data-flow diagram for each gem_annotate module.

Each diagram shows, in English: what flows IN (model state / external tables),
what the module DOES (key functions / steps), what external data it consults,
and how the model is mutated on the way OUT.

Output: docs/module_flows/<module>.png
"""
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "docs", "module_flows")

# colors per node role
C_IN = "#cfe8ff"      # input
C_STEP = "#fff2cc"    # processing step
C_EXT = "#e2f0d9"     # external data source
C_OUT = "#f8cbad"     # output / model mutation


def _box(ax, x, y, w, h, text, color, fontsize=8.5):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.1, edgecolor="#555555", facecolor=color, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, zorder=3, wrap=True)


def _arrow(ax, x1, y1, x2, y2, style="-|>", color="#444444", ls="-"):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle=style, mutation_scale=14,
        linewidth=1.2, color=color, linestyle=ls, zorder=1,
        shrinkA=2, shrinkB=2))


def draw_metabolites():
    """METABOLITE layer: name/BiGG/formula -> MNXM annotation + formula/charge."""
    fig, ax = plt.subplots(figsize=(12, 7.5))
    ax.set_xlim(0, 12); ax.set_ylim(0, 7.5); ax.axis("off")
    ax.set_title("gem_annotate / metabolites.py  —  metabolite annotation & H/H2O balance",
                 fontsize=12, fontweight="bold", loc="left")

    # INPUT
    _box(ax, 0.2, 5.6, 2.4, 1.1,
         "IN: model metabolites\n(name, BiGG id,\nraw formula in name)", C_IN)
    # EXTERNAL
    _box(ax, 0.2, 0.4, 2.4, 1.6,
         "EXTERNAL (loaded by io.py):\nchem_xref  (id -> MNXM)\nchem_prop  (MNXM -> formula,\ncharge, xrefs)\n_DIRECT_MNXM_TABLE", C_EXT)

    # STEPS
    _box(ax, 3.2, 5.6, 3.0, 1.1,
         "annotate_metabolites()\nmatch each met to an MNXM,\ntry strategies in order", C_STEP)
    _box(ax, 6.7, 4.3, 5.0, 2.6,
         "Match strategy order:\n"
         "A  = BiGG id\n"
         "BD = _DIRECT_MNXM_TABLE (curated)\n"
         "B / B0 / B1 / B2a / B2b = name\n"
         "       (exact / excel-fix / synonym /\n"
         "        normalized / prefix)\n"
         "C  = formula match\n"
         "+ carbon-count guard rejects\n"
         "  name matches whose C-count\n"
         "  disagrees (name-collision guard)", C_STEP, fontsize=8)
    _box(ax, 3.2, 3.5, 3.0, 1.1,
         "_apply_mnxm()\npull formula / charge /\ncross-refs from chem_prop", C_STEP)
    _box(ax, 3.2, 1.9, 3.0, 1.1,
         "fix_proton_water_balance()\nadd/adjust H+ / H2O so\nreactions mass-balance", C_STEP)
    _box(ax, 3.2, 0.4, 3.0, 1.0,
         "normalize_all_annotations()\nfinal annotation cleanup\n(run last in pipeline)", C_STEP)

    # OUTPUT
    _box(ax, 9.4, 1.6, 2.4, 1.6,
         "OUT: model metabolites\nnow carry\nmetanetx.chemical,\nformula, charge,\ncross-ref annotations", C_OUT)

    # arrows
    _arrow(ax, 2.6, 6.15, 3.2, 6.15)                 # IN -> annotate
    _arrow(ax, 6.2, 6.15, 6.7, 5.6)                  # annotate -> strategies
    _arrow(ax, 2.6, 1.2, 2.6, 5.6, ls="--")          # external -> (consulted)
    _arrow(ax, 2.6, 5.9, 3.2, 5.9, color="#888")     # external feeds annotate
    _arrow(ax, 4.7, 5.6, 4.7, 4.6)                   # annotate -> apply_mnxm
    _arrow(ax, 4.7, 3.5, 4.7, 3.0)                   # apply_mnxm -> balance
    _arrow(ax, 4.7, 1.9, 4.7, 1.4)                   # balance -> normalize
    _arrow(ax, 6.2, 4.05, 9.4, 2.6)                  # apply_mnxm -> OUT

    # legend row at bottom-right
    lx = 6.9
    for c, lbl in [(C_IN, "input"), (C_STEP, "step"),
                   (C_EXT, "external"), (C_OUT, "output")]:
        _box(ax, lx, 0.45, 0.3, 0.3, "", c, fontsize=6)
        ax.text(lx + 0.4, 0.6, lbl, fontsize=7.5, va="center")
        lx += 1.25

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "metabolites.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def _flow(module, title, in_text, steps, ext_text, out_text):
    """Generic vertical flow: IN -> step1 -> step2 ... -> OUT, with one
    EXTERNAL box feeding the steps.  `steps` is a list of multi-line strings."""
    n = len(steps)
    fig_h = max(6.0, 2.2 + 1.55 * n)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    ax.set_xlim(0, 12); ax.set_ylim(0, fig_h); ax.axis("off")
    ax.set_title(f"gem_annotate / {module}.py  —  {title}",
                 fontsize=12, fontweight="bold", loc="left")

    # leave a clear band at the bottom for the legend
    legend_y = 0.15
    base = 1.2
    top = fig_h - 1.2
    # INPUT
    _box(ax, 0.2, top - 0.55, 2.6, 1.1, in_text, C_IN)
    # EXTERNAL (bottom-left) if any
    if ext_text:
        _box(ax, 0.2, base + 0.1, 2.6, 1.7, ext_text, C_EXT, fontsize=8)

    # STEPS column (center) — kept above the legend band
    sx, sw, sh = 3.4, 4.2, 1.05
    gap = (top - 0.55 - (base + 0.2)) / max(n, 1)
    ys = []
    for i, s in enumerate(steps):
        sy = top - 0.55 - (i + 1) * gap
        ys.append(sy)
        _box(ax, sx, sy, sw, sh, s, C_STEP, fontsize=8.5)
    # OUTPUT (bottom-right, above the legend band)
    _box(ax, 8.6, base + 0.1, 3.0, 1.7, out_text, C_OUT)

    # arrows
    _arrow(ax, 2.8, top, sx, ys[0] + sh / 2)             # IN -> step1
    for i in range(n - 1):
        _arrow(ax, sx + sw / 2, ys[i], sx + sw / 2, ys[i + 1] + sh)
    _arrow(ax, sx + sw / 2, ys[-1], 9.8, base + 1.0)     # last step -> OUT
    if ext_text:
        _arrow(ax, 1.5, base + 1.8, 1.5, top - 0.55, ls="--")  # external consulted
        _arrow(ax, 2.8, top - 0.2, sx, ys[0] + sh - 0.2, color="#999")

    # legend
    lx = 6.6
    for c, lbl in [(C_IN, "input"), (C_STEP, "step"),
                   (C_EXT, "external"), (C_OUT, "output")]:
        _box(ax, lx, 0.15, 0.3, 0.3, "", c, fontsize=6)
        ax.text(lx + 0.4, 0.3, lbl, fontsize=7.5, va="center")
        lx += 1.25

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"{module}.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def draw_io():
    _flow("io", "load MetaNetX reference tables once",
          "IN: file paths to\nMetaNetX TSV dumps",
          ["load_chem_xref()\nbuild id -> [(MNXM, source)]",
           "load_chem_prop()\nMNXM -> {formula, charge, xrefs}",
           "load_reac_xref() / load_reac_prop()\nreaction id -> MNXR, MNXR props",
           "load_mnxm_depr()\ndeprecated MNXM -> current MNXM"],
          "EXTERNAL (on disk):\ndata/metanetx/\n  chem_xref.tsv\n  chem_prop.tsv\n  reac_xref.tsv\n  reac_prop.tsv\n  chem_depr.tsv",
          "OUT: in-memory dicts\nreused by every\ndownstream step\n(no model mutation)")


def draw_reactions():
    _flow("reactions", "MNXR fingerprint annotation + xref backfill",
          "IN: model reactions\n(metabolites already\ncarry MNXM)",
          ["annotate_reactions()\nfingerprint each reaction by its\nMNXM substrate/product set\n-> match a MNXR",
           "backfill_reaction_xrefs()\nfrom matched MNXR fill missing\nbigg / kegg / rhea / ec-code"],
          "EXTERNAL (io.py):\nreac_xref (fingerprint\n  -> MNXR)\nreac_prop (MNXR ->\n  bigg/kegg/rhea/ec)",
          "OUT: reactions carry\nmetanetx.reaction +\nbigg/kegg/rhea/\nec-code annotations")


def draw_annotate_reactions_extended():
    _flow("annotate_reactions_extended", "annotate reactions still unmatched",
          "IN: reactions left\nunannotated after\nreactions.py",
          ["_is_unannotated()\nfind reactions with no MNXR yet",
           "annotate_remaining_reactions()\nclassify: exchange / transport /\nEC->MNXR fallback",
           "_resolve_multi_candidates()\npick best MNXR when several\nfingerprints match"],
          "EXTERNAL (io.py):\nreac_xref / reac_prop\n(EC -> MNXR,\n transport patterns)",
          "OUT: previously\nunannotated reactions\nnow carry MNXR /\ntransport / EC tags")


def draw_exchange():
    _flow("exchange", "exchange bounds + minimal medium",
          "IN: model exchange\n& transport reactions",
          ["set_exchange_bounds()\nclose all uptake, then open\nonly the defined medium",
           "configure_medium()\napply minimal-medium uptake\nrates + required vitamins"],
          "EXTERNAL (config.py):\nmedium_bigg dict\n(allowed uptake\n metabolites + rates)",
          "OUT: exchange lower/\nupper bounds set so\nFBA grows on the\nintended medium")


def draw_biomass():
    _flow("biomass", "fix the biomass pseudo-reaction",
          "IN: model biomass\nreaction (R1372)\n+ its precursors",
          ["_formula_mw()\ncompute MW of each\nbiomass component",
           "fix_biomass_reaction()\nrescale coefficients so the\nbiomass MW is ~1 g/mmol\nand FBA is well-posed"],
          "",
          "OUT: biomass reaction\nstoichiometry rescaled;\nmodel.optimize() gives\na meaningful growth rate")


def draw_genes():
    _flow("genes", "annotate genes from the UniProt proteome",
          "IN: model genes\n(YALI locus tags)",
          ["_fetch_proteome()\npull the W29 proteome\nfrom UniProt",
           "_tier_a / _tier_b / _tier_ncbi()\nmatch locus tags to UniProt\nby decreasing confidence",
           "_enrich_kegg_genes()\nadd kegg.genes via\nUniProt<->KEGG index",
           "annotate_genes()\nwrite uniprot / ncbigene /\nkegg.genes / refseq onto genes"],
          "EXTERNAL (network):\nUniProt proteome &\nsearch REST API\nKEGG gene index",
          "OUT: genes carry\nuniprot / ncbigene /\nkegg.genes / refseq\ncross-references")


def draw_idmapping():
    _flow("idmapping", "second-pass gene mapping via GeneID",
          "IN: genes still\nunmapped after\ngenes.py",
          ["_search_uniprot_by_geneids()\nquery UniProt by\nxref:geneid-<id>",
           "_enrich_via_idmapping()\nmerge any newly found\nUniProt accessions onto genes"],
          "EXTERNAL (network):\nUniProt REST\n(xref:geneid search)",
          "OUT: extra genes gain\nuniprot / cross-ref\nannotations that the\nproteome pass missed")


def draw_ec_annotation():
    _flow("ec_annotation", "attach EC numbers to genes",
          "IN: genes with\nUniProt accessions",
          ["_fetch_ec_for_accessions()\nstream EC numbers for the\naccessions from UniProt",
           "enrich_genes_with_ec()\nwrite ec-code onto each gene"],
          "EXTERNAL (network):\nUniProt stream API\n(accession -> EC)",
          "OUT: genes carry\nec-code; later copied\nto reactions in the\nEC backfill step")


def draw_gaps():
    _flow("gaps", "gap analysis (FVA): blocked / dead-end / orphan",
          "IN: fully annotated,\nmedium-configured\nmodel",
          ["find_gaps()\nFVA: a reaction is blocked if\nmax flux approx 0 AND min approx 0;\nfind dead-end / orphan metabolites",
           "merge_duplicate_metabolites()\nmerge known duplicate met pairs\nfor stoichiometric consistency",
           "add_gap_fill_reactions()\ninsert prioritized gap-fill\nreactions from a CSV",
           "report_gaps()\nsummarize counts to the log"],
          "EXTERNAL:\ndata/ gap-fill CSV\n(prioritized\n reactions to add)",
          "OUT: duplicate mets\nmerged; P0 gap-fill\nreactions added; gap\ncounts reported")


def draw_gap_fill_prioritize():
    _flow("gap_fill_prioritize", "rank candidate gap-fill reactions",
          "IN: list of missing /\nblocked metabolites\n& candidate reactions",
          ["has_formula() / mnxm_in_model()\ncheck each candidate's mets\nfor formula + presence in model",
           "recount_missing()\nhow many mets are still\nmissing if this rxn is added",
           "assign_priority()\nrank P0/P1/... by metabolite\ncoverage + bigg/kegg support"],
          "EXTERNAL (io.py):\nchem_prop (has formula?)\nreac_xref (bigg/kegg?)",
          "OUT: prioritized\ngap-fill candidate\ntable (CSV) consumed\nby gaps.add_gap_fill")


_DRAWERS = [
    draw_metabolites, draw_io, draw_reactions,
    draw_annotate_reactions_extended, draw_exchange, draw_biomass,
    draw_genes, draw_idmapping, draw_ec_annotation, draw_gaps,
    draw_gap_fill_prioritize,
]


if __name__ == "__main__":
    for d in _DRAWERS:
        d()
    print(f"\n{len(_DRAWERS)} diagrams written to {OUTDIR}")
