"""
patches.py — known data-bug patches for the iYli21 model.

Each patch fixes a discrete, documented error in the source SBML.  Patches
are applied early in the pipeline, immediately after metabolite annotation,
so that downstream mass-balance, FBA, and gap analyses operate on a
chemically correct model.

Currently applied patches:

  1. NADP+ formula fix
     iYli21 stores NADP+ as C21H28N7O17P3 (one H short of the KEGG C00006
     reference C21H29N7O17P3).  This affects 6 compartmental copies and
     ~126 reactions that use NADP+/NADPH, which had cascading effects on
     ceramide synthesis (see patch 2).

  2. Ceramide formula corrections
     iYli21 stores ceramide-1-(C24/C26) with formula C63H125NO6, which is
     inconsistent with the LIPID MAPS / Yeast-GEM reference value
     C42H85NO3 (C24) / C44H89NO3 (C26).  ceramide-2/2'/3/4 are simply
     missing formulas.  These were left wrong/missing because the upstream
     NADP+ error made the relevant reactions unbalanceable — fixing patch 1
     now makes patch 2 work.

     Naming map (verified against LIPID MAPS):
       ceramide-1  = Cer(d18:0/24:0)         = C42H85NO3  / C44H89NO3
       ceramide-2  = Cer(t18:0/24:0)         = C42H85NO4  / C44H89NO4
       ceramide-2' = Cer(d18:0/24:0(2OH))    = C42H85NO4  / C44H89NO4
       ceramide-3  = Cer(t18:0/24:0(2OH))    = C42H85NO5  / C44H89NO5
       ceramide-4  = (extrapolated +1 O from ceramide-3)
                                              C42H85NO6  / C44H89NO6
"""

from __future__ import annotations

import logging
import re

logger = logging.getLogger(__name__)


# ── Patch 1: NADP+ formula fix ────────────────────────────────────────────

# iYli21 NADP+ copies (all compartments). Identified by metabolite name
# prefix "NADP(+)_" which is the iYli21 naming convention.
# Formula correction: C21H28N7O17P3 → C21H29N7O17P3 (KEGG C00006).
_NADP_PLUS_OLD_FORMULA = "C21H28N7O17P3"
_NADP_PLUS_NEW_FORMULA = "C21H29N7O17P3"


def _is_nadp_plus(met) -> bool:
    """Identify NADP+ by name (handles iYli21 'NADP(+)_C21H28N7O17P3' naming)."""
    name = (met.name or "")
    return name.startswith("NADP(+)") or name.lower().startswith("nadp(+)")


def fix_nadp_plus_formula(model) -> int:
    """
    Fix the NADP+ formula bug across all compartments.

    Returns the number of metabolites patched.
    """
    fixed = 0
    for met in model.metabolites:
        if not _is_nadp_plus(met):
            continue
        if met.formula == _NADP_PLUS_OLD_FORMULA:
            met.formula = _NADP_PLUS_NEW_FORMULA
            fixed += 1
            logger.debug(f"  NADP+ patch: {met.id} formula → {_NADP_PLUS_NEW_FORMULA}")
        elif met.formula == _NADP_PLUS_NEW_FORMULA:
            pass   # already correct
        else:
            logger.warning(
                f"  NADP+ patch: {met.id} has unexpected formula {met.formula!r} "
                f"— skipped (expected {_NADP_PLUS_OLD_FORMULA})"
            )
    return fixed


# ── Patch 2: Ceramide formula corrections ─────────────────────────────────

# Map: lowercased base name (compartment-independent) → (formula_C24, formula_C26)
# Base name is what comes BEFORE the trailing "_" or "_FORMULA" in iYli21 names.
# E.g. "ceramide-1-(C24)_" matches base "ceramide-1-(c24)".
_CERAMIDE_FORMULAS_C24 = {
    "ceramide-1-(c24)":  "C42H85NO3",   # Cer(d18:0/24:0)        LMSP02020012
    "ceramide-2-(c24)":  "C42H85NO4",   # Cer(t18:0/24:0)        LMSP02030004
    "ceramide-2'-(c24)": "C42H85NO4",   # Cer(d18:0/24:0(2OH))   (α-OH dihydroceramide)
    "ceramide-3-(c24)":  "C42H85NO5",   # Cer(t18:0/24:0(2OH))   LMSP02030002
    "ceramide-4-(c24)":  "C42H85NO6",   # (extrapolated, R208 product)
}
_CERAMIDE_FORMULAS_C26 = {
    "ceramide-1-(c26)":  "C44H89NO3",
    "ceramide-2-(c26)":  "C44H89NO4",
    "ceramide-2'-(c26)": "C44H89NO4",
    "ceramide-3-(c26)":  "C44H89NO5",
    "ceramide-4-(c26)":  "C44H89NO6",
}
_CERAMIDE_FORMULAS = {**_CERAMIDE_FORMULAS_C24, **_CERAMIDE_FORMULAS_C26}

# Old (wrong) formula on ceramide-1 — explicitly remove this value when we
# overwrite, so the log warns if it's anything else (defensive).
_CERAMIDE_1_KNOWN_BAD_FORMULA = "C63H125NO6"


def _base_name(met) -> str:
    """Return the iYli21 base name (lowercased, strip trailing _ and _FORMULA)."""
    n = (met.name or "").lower().rstrip("_").strip()
    # Strip trailing "_<FORMULA>" if present (e.g. "ceramide-1-(c24)_c63h125no6")
    if "_" in n:
        parts = n.rsplit("_", 1)
        last = parts[-1]
        # heuristic: looks like a chemical formula → strip
        if last and last[0] in "cChH" and any(c.isdigit() for c in last):
            n = parts[0].rstrip("_")
    return n


def fix_ceramide_formulas(model) -> int:
    """
    Fix ceramide-1 wrong formula + fill ceramide-2/2'/3/4 missing formulas.

    Authoritative source: LIPID MAPS LMSD (verified manually) cross-checked
    against Yeast-GEM (https://github.com/SysBioChalmers/yeast-GEM).

    Applies to all compartmental copies (C_er, C_mi, C_go, etc).
    Returns number of metabolites patched.
    """
    fixed = 0
    seen_bases = set()
    for met in model.metabolites:
        base = _base_name(met)
        if base not in _CERAMIDE_FORMULAS:
            continue
        new_formula = _CERAMIDE_FORMULAS[base]
        old_formula = met.formula or ""
        if old_formula == new_formula:
            continue   # already correct
        if old_formula and old_formula != _CERAMIDE_1_KNOWN_BAD_FORMULA:
            logger.warning(
                f"  Ceramide patch: {met.id} ({base!r}) has unexpected existing "
                f"formula {old_formula!r}, overwriting with {new_formula!r}"
            )
        met.formula = new_formula
        fixed += 1
        seen_bases.add(base)
        logger.debug(
            f"  Ceramide patch: {met.id} [{met.compartment}] {old_formula or '(none)'} "
            f"→ {new_formula}  ({base})"
        )
    # Report which base names did not appear in the model
    missing = set(_CERAMIDE_FORMULAS) - seen_bases
    if missing:
        logger.info(
            f"  Ceramide patch: {len(missing)} mapped base name(s) not found in model "
            f"(possibly C24 or C26 variants absent): {sorted(missing)[:5]}…"
        )
    return fixed


# ── Patch 3: free CoA charge fix ──────────────────────────────────────────

# iYli21 stores free coenzyme A (KEGG C00010) with charge 0, but the correct
# physiological charge is -4.  Only the 9 free-CoA copies are affected; the
# acyl-CoA species (tetracosanoyl-CoA, acetyl-CoA, …) already carry charge -4.
# Identified by: name starts with "coenzyme A", formula == C21H36N7O16P3S, charge == 0.
_COA_FREE_FORMULA   = "C21H36N7O16P3S"
_COA_CORRECT_CHARGE = -4


def fix_coa_charge(model) -> int:
    """
    Fix the free-CoA charge bug (KEGG C00010: charge should be -4, not 0).

    Patches only free coenzyme A; acyl-CoA species are left untouched because
    they already carry the correct charge.  Returns the number patched.
    """
    fixed = 0
    for met in model.metabolites:
        name = met.name or ""
        if not name.startswith("coenzyme A"):
            continue
        if met.formula == _COA_FREE_FORMULA and met.charge == 0:
            met.charge = _COA_CORRECT_CHARGE
            fixed += 1
            logger.debug(f"  CoA patch: {met.id} charge → {_COA_CORRECT_CHARGE}")
    return fixed


# ── Patch 6: cation formula/charge self-consistency ───────────────────────

# A set of protonated amines/cations carry a NEUTRAL formula but a positive
# charge — the two are mutually inconsistent (a +1 ammonium must have one more
# H than its neutral form).  This unbalances every reaction they participate
# in (e.g. sphinganine's 8 sphingolipid reactions).  Fix the formula to the
# protonated (ionic) form, matching the MetaNetX chem_prop value for the
# metabolite's own MNXM.
#
# Whitelist keyed by (current_formula, charge) → target_formula, so a same-
# formula species with a different charge is never touched.  Each target is
# the MetaNetX formula for that metabolite's MNXM (verified 2026-06).
_CATION_FORMULA_FIX: dict[tuple[str, int], str] = {
    ("C11H18N4O4", 1): "C11H19N4O4",   # MNXM735139  2-[3-carboxy-3-(methylammonio)propyl]-L-his
    ("C15H22N6O5S", 1): "C15H23N6O5S", # MNXM1363767 S-adenosyl-L-methionine
    ("C18H39NO2", 1): "C18H40NO2",     # MNXM733692  sphinganine
    ("C2H7NO", 1): "C2H8NO",           # MNXM218     ethanolamine
    ("C3H7NO", 1): "C3H8NO",           # MNXM736082  3-aminopropanal
    ("C4H12N2", 2): "C4H14N2",         # MNXM118     putrescine (+2)
    ("C4H9NO", 1): "C4H10NO",          # MNXM422     4-aminobutanal
    ("C5H12N4O", 1): "C5H13N4O",       # MNXM2617    4-guanidinobutanamide
    ("C6H11N3O", 1): "C6H12N3O",       # MNXM1281    L-histidinol
    ("C6H13NO2S", 1): "C6H14NO2S",     # MNXM681265  S-methyl-L-methionine
    ("C6H14N2O2", 1): "C6H15N2O2",     # MNXM1364268 L-lysine
    ("C6H14N4O2", 1): "C6H15N4O2",     # MNXM739527  L-arginine
    ("C7H19N3", 3): "C7H22N3",         # MNXM124     spermidine (+3)
    ("C8H12N2O2", 1): "C8H13N2O2",     # MNXM548     pyridoxamine
}


def fix_cation_formula_consistency(model) -> int:
    """
    Fix protonated cations whose formula is the neutral form but charge is
    positive.  Only metabolites matching both (formula, charge) in the
    whitelist are changed.  Returns the number of metabolites patched.
    """
    fixed = 0
    for met in model.metabolites:
        if met.formula is None or met.charge is None:
            continue
        key = (met.formula.strip(), int(met.charge))
        target = _CATION_FORMULA_FIX.get(key)
        if target:
            met.formula = target
            fixed += 1
            logger.debug(f"  cation formula patch: {met.id} {key[0]} → {target}")
    return fixed


# ── Patch 0: clean Excel-corrupted "ActiveX VT_ERROR" names ───────────────

# iYli21 was exported from Excel; some metabolite names had their trailing
# chemical-formula token replaced by the literal string "ActiveX VT_ERROR:",
# e.g. "butyrate_ActiveX VT_ERROR:".  The corrupted name prevents
# annotate_metabolites from matching the metabolite to MetaNetX (no MNXM, no
# formula).  Stripping the garbage suffix restores a clean name that the
# annotation step can match.
#
# MUST run BEFORE annotate_metabolites so the clean name is used for matching
# — so it is invoked from main() before annotation, not from apply_all_patches.
_ACTIVEX_RE = re.compile(r"_?ActiveX VT_ERROR:?\s*$")


def fix_activex_names(model) -> int:
    """
    Strip the Excel-corruption suffix '_ActiveX VT_ERROR:' from metabolite
    names.  Returns the number of names cleaned.
    """
    fixed = 0
    for met in model.metabolites:
        name = met.name or ""
        if "ActiveX" in name or "VT_ERROR" in name:
            cleaned = _ACTIVEX_RE.sub("", name).rstrip("_ ").strip()
            if cleaned and cleaned != name:
                met.name = cleaned
                fixed += 1
                logger.debug(f"  ActiveX name patch: {met.id} → {cleaned!r}")
    return fixed


# ── Patch 5: move misfiled TCDB numbers out of ec-code ────────────────────

# Some transport reactions carry a TCDB (Transporter Classification Database)
# number in their 'ec-code' annotation, e.g. "2.A.1.44.1" or "9.A.6.1.1".
# These are not EC numbers (second segment is a letter) and fail Memote's
# ec-code conformity check.  They are valid identifiers, just in the wrong
# field — move them to the identifiers.org-registered 'tcdb' annotation.
#
# NOTE: this must run after all EC codes are populated, so it is invoked from
# main() (before fix_ec_code_format), not from apply_all_patches.
_TCDB_RE = re.compile(r"^\d+\.[A-Z]\.")   # TCDB class id, e.g. 2.A.1.44.1


def move_tcdb_out_of_ec(model) -> int:
    """
    Move TCDB transporter numbers misfiled in 'ec-code' to the 'tcdb'
    annotation field.  Returns the number of TCDB ids moved.
    """
    moved = 0
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        ecs = ann.get("ec-code")
        if not ecs:
            continue
        if isinstance(ecs, str):
            ecs = [ecs]
        tcdb = [e for e in ecs if _TCDB_RE.match(e)]
        if not tcdb:
            continue
        kept = [e for e in ecs if not _TCDB_RE.match(e)]
        existing = ann.get("tcdb", [])
        if isinstance(existing, str):
            existing = [existing]
        ann["tcdb"] = sorted(set(existing) | set(tcdb))
        if kept:
            ann["ec-code"] = kept
        else:
            ann.pop("ec-code", None)
        rxn.annotation = ann
        moved += len(tcdb)
        logger.debug(f"  TCDB cleanup: {rxn.id} moved {tcdb} ec-code → tcdb")
    return moved


# ── Patch 4: EC-code format compliance ───────────────────────────────────

# Some reactions carry EC codes missing the final level (exactly three
# numeric segments, e.g. "3.1.3" or "1.2.1").  These violate the
# identifiers.org EC regex and are flagged by Memote's
# test_reaction_annotation_wrong_ids (ec-code conformity).  Padding them with
# a trailing ".-" makes them valid partial EC codes (e.g. "3.1.3.-") without
# inventing a more specific level we cannot verify.
#
# Only "exactly three pure-numeric segments" are touched.  EC strings with
# letters (preliminary IDs like "2.7.1.M29", malformed "7.4.2.i") or
# misfiled non-EC identifiers (e.g. TCDB number "2.A.29.8.3") are NOT padded
# — those need case-by-case handling and are left untouched (logged).
_EC_THREE_SEGMENT_RE = re.compile(r"^\d+\.\d+\.\d+$")


def fix_ec_code_format(model) -> int:
    """
    Pad three-segment numeric EC codes with a trailing '.-' for
    identifiers.org compliance (e.g. '3.1.3' -> '3.1.3.-').

    Only exact three-segment pure-numeric EC codes are modified.  EC strings
    containing letters or extra/fewer segments are left untouched.  Returns
    the number of individual EC strings padded.
    """
    fixed = 0
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        ecs = ann.get("ec-code")
        if not ecs:
            continue
        if isinstance(ecs, str):
            ecs = [ecs]
        new = []
        changed = False
        for ec in ecs:
            if _EC_THREE_SEGMENT_RE.match(ec):
                new.append(ec + ".-")
                fixed += 1
                changed = True
                logger.debug(f"  EC-format patch: {rxn.id} {ec} → {ec}.-")
            else:
                new.append(ec)
        if changed:
            rxn.annotation["ec-code"] = new
    return fixed


# ── Patch 7: EC-overload cleanup ──────────────────────────────────────────

# Some reactions accumulated 5+ EC numbers from MetaNetX MNXR `classifs`
# back-fill — a single, well-defined reaction tagged with EC numbers spanning
# multiple enzyme classes (e.g. R1893 mannitol dehydrogenase carried 11 EC
# including glutathione transferase 2.5.1.18).  This "EC soup" pollutes
# EC-based analysis and inflates Memote's EC inconsistency.
#
# The authoritative fix uses KEGG reaction-level ENZYME data: keep only the
# intersection of the reaction's current EC set with the EC numbers KEGG
# assigns to its kegg.reaction id.  That curation is done offline by
# scripts/audit_ec_overload.py, which writes data/ec_overload_audit.csv with
# a per-reaction keep/drop decision.  This patch *applies* the rows marked
# action=clean — it never invents EC numbers (keep_ec is always a subset of
# the reaction's current ec-code, double-checked below).
#
# Runs LAST in the pipeline (after all EC back-fill and format steps) so the
# back-fill cannot re-introduce the dropped EC numbers.
_EC_OVERLOAD_AUDIT_CSV = "data/ec_overload_audit.csv"


def clean_ec_overload(model, audit_csv: str | None = None) -> int:
    """
    Apply the KEGG-curated EC-overload cleanup from the audit CSV.

    For each audit row with action=='clean', replace the reaction's 'ec-code'
    with the curated keep_ec set — but ONLY if keep_ec is a subset of the
    reaction's current ec-code (guards against the CSV drifting out of sync
    with the model).  Returns the number of reactions cleaned.
    """
    import csv
    import os

    if audit_csv is None:
        # patches.py is at scripts/gem_annotate/patches.py -> repo root is ../../
        root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        audit_csv = os.path.join(root, _EC_OVERLOAD_AUDIT_CSV)

    if not os.path.exists(audit_csv):
        logger.warning(f"  EC-overload audit CSV not found: {audit_csv} — skipping")
        return 0

    cleaned = 0
    with open(audit_csv) as f:
        for row in csv.DictReader(f):
            if row.get("action") != "clean":
                continue
            rid = row["reaction"]
            keep = {e for e in row["keep_ec"].split(";") if e}
            if not keep:
                continue
            try:
                rxn = model.reactions.get_by_id(rid)
            except KeyError:
                logger.warning(f"  EC-overload: reaction {rid} not in model — skipping")
                continue
            ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
            cur = ann.get("ec-code")
            if not cur:
                continue
            if isinstance(cur, str):
                cur = [cur]
            cur_set = set(cur)
            # subset guard: only proceed if keep_ec really is a subset of current
            if not keep <= cur_set:
                logger.warning(
                    f"  EC-overload: {rid} keep_ec {sorted(keep)} not subset of "
                    f"current {sorted(cur_set)} — skipping (CSV out of sync)")
                continue
            if cur_set == keep:
                continue  # already clean, nothing to drop
            rxn.annotation["ec-code"] = sorted(keep)
            cleaned += 1
            logger.debug(
                f"  EC-overload: {rid} {sorted(cur_set)} → {sorted(keep)}")
    return cleaned


# ── Patch 8: isozyme GPR additions ────────────────────────────────────────

# CLIB89 expansion funnel (S2 Table -> KEGG metabolic relevance -> MetaNetX
# reaction mapping -> de-dup against model EC) identified genes whose EC is
# already carried by an existing model reaction: these are isozymes that
# should be added to that reaction's GPR (NOT new reactions).
#
# This patch applies only the SAFE subset curated in
# data/gpr_isozyme_additions.csv: each gene maps to <=3 reactions and none of
# those reactions has an 'and' (multi-subunit complex) GPR — so adding the
# gene with 'or' is unambiguous.  Genes hitting broad ECs (many reactions) or
# complex GPRs are excluded and left for manual review.
#
# Only 'or' additions are made; existing genes and reaction stoichiometry are
# never touched.  Idempotent: a gene already in the rule is skipped.
_GPR_ADDITIONS_CSV = "data/gpr_isozyme_additions.csv"


def add_isozyme_gprs(model, additions_csv: str | None = None) -> int:
    """
    Add curated isozyme genes to existing reactions' GPR via 'or'.

    Reads data/gpr_isozyme_additions.csv (columns: reaction, add_gene, ...).
    For each (reaction, gene): if the gene is not already in the reaction's
    gene_reaction_rule, append it with 'or' (or set it as the sole rule if the
    reaction had no GPR).  Returns the number of (reaction, gene) additions made.
    """
    import csv
    import os

    if additions_csv is None:
        root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        additions_csv = os.path.join(root, _GPR_ADDITIONS_CSV)

    if not os.path.exists(additions_csv):
        logger.warning(f"  GPR additions CSV not found: {additions_csv} — skipping")
        return 0

    added = 0
    with open(additions_csv) as f:
        for row in csv.DictReader(f):
            rid = row["reaction"]
            gene = row["add_gene"].strip()
            if not gene:
                continue
            try:
                rxn = model.reactions.get_by_id(rid)
            except KeyError:
                logger.warning(f"  GPR add: reaction {rid} not in model — skipping")
                continue
            rule = rxn.gene_reaction_rule.strip()
            # idempotent: skip if gene already present as a token
            existing = set(re.findall(r"[A-Za-z0-9_]+", rule))
            if gene in existing:
                continue
            new_rule = gene if not rule else f"{rule} or {gene}"
            rxn.gene_reaction_rule = new_rule
            added += 1
            logger.debug(f"  GPR add: {rid} += {gene}")
    return added


# ── Patch 8b: remove mis-annotated genes from reaction GPRs ────────────────
#
# Two genes were assigned to reactions whose enzyme they do NOT catalyse — a
# reaction-mis-annotation (the gene's real protein is a different enzyme). Both
# are OR-isozyme partners, so removing them leaves the true catalyst in place and
# does NOT change any FBA/growth result — it only corrects GPR biology and removes
# a false-negative contaminant from essential-gene screening.
#
# Verified this session against UniProt (opened):
#   YALI1E07744g -> UniProt Q6C6P1 = glycoside hydrolase family 65 (trehalase),
#     NOT transketolase (EC 2.2.1.1). Wrongly placed in R765/R766; real
#     transketolase is the partner YALI1D02625g.
#     https://rest.uniprot.org/uniprotkb/Q6C6P1.json
#   YALI1E11370g -> similar to PET112/GatB, UniProt P33893 = glutamyl-tRNA(Gln)
#     amidotransferase subunit B (EC 6.3.5.-), NOT prephenate dehydrogenase.
#     Wrongly placed in R671; real prephenate dehydrogenase is the partner
#     YALI1F23441g.  https://rest.uniprot.org/uniprotkb/P33893.json
#
# (R671 also carries a self-contradictory EC 6.3.5.7 vs its "prephenate
# dehydrogenase" name — flagged for separate curation, NOT changed here.)
_GPR_MISANNOTATION_REMOVALS = {
    "R765": "YALI1E07744g",
    "R766": "YALI1E07744g",
    "R671": "YALI1E11370g",
}


def remove_misannotated_gprs(model) -> int:
    """Remove mis-annotated isozyme genes from reaction GPRs (safe: all 'or' cases,
    true partner retained). Returns the number of (reaction, gene) removals made.
    Idempotent; refuses to empty a GPR."""
    removed = 0
    for rid, gene in _GPR_MISANNOTATION_REMOVALS.items():
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            logger.warning(f"  GPR removal: reaction {rid} not in model — skipping")
            continue
        toks = set(re.findall(r"YALI1[A-Za-z0-9]+", rxn.gene_reaction_rule))
        if gene not in toks:
            continue  # idempotent: already removed
        remaining = sorted(t for t in toks if t != gene)
        if not remaining:
            logger.warning(f"  GPR removal: {rid} -= {gene} would empty GPR — skipping")
            continue
        rxn.gene_reaction_rule = " or ".join(remaining) if len(remaining) > 1 else remaining[0]
        removed += 1
        logger.info(f"  GPR removal (mis-annotation): {rid} -= {gene} -> '{rxn.gene_reaction_rule}'")
    return removed


# ── Patch 9: annotate the isozyme genes added by patch 8 ───────────────────

# The genes added by add_isozyme_gprs enter the model with no annotation.
# The pipeline's main gene-annotation step (genes.annotate_genes) and SBO step
# both run BEFORE the GPR additions, so they never reach these new genes —
# leaving them with empty annotation (regresses memote gene-SBO /
# gene-product-annotation, and they re-appear empty on every full rebuild).
#
# This patch annotates exactly those genes: sbo (SBO:0000243), ncbigene +
# kegg.genes (from NCBI feature table / KEGG yli, local), and uniprot
# (best-effort network lookup via xref:geneid).  Runs right after
# add_isozyme_gprs.  Idempotent: genes already carrying an sbo are skipped.
_FEATURE_TABLE = "data/ncbi/clib89_feature_table.txt"
_KEGG_GENES = "data/kegg/yli_genes.tsv"


def _fetch_uniprot_for_geneid(geneid: str):
    """Best-effort UniProt accession from a GeneID xref. None on failure."""
    import urllib.parse
    import urllib.request
    q = urllib.parse.quote(f"xref:geneid-{geneid}")
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           f"&fields=accession&format=tsv&size=1")
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            lines = resp.read().decode().splitlines()
        if len(lines) >= 2 and lines[1].strip():
            return lines[1].strip()
    except Exception:
        pass
    return None


def annotate_isozyme_genes(model, additions_csv: str | None = None,
                           network: bool = True) -> int:
    """
    Annotate the isozyme genes added by add_isozyme_gprs.

    Adds sbo / ncbigene / kegg.genes (local) and uniprot (network, optional).
    Idempotent: a gene that already has an 'sbo' annotation is skipped.
    Returns the number of genes annotated.
    """
    import csv
    import os

    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    if additions_csv is None:
        additions_csv = os.path.join(root, _GPR_ADDITIONS_CSV)
    if not os.path.exists(additions_csv):
        logger.warning(f"  gene annotate: {additions_csv} not found — skipping")
        return 0

    # YALI1 (no underscore) -> GeneID from NCBI feature table
    y2g = {}
    ft = os.path.join(root, _FEATURE_TABLE)
    if os.path.exists(ft):
        with open(ft) as f:
            next(f)
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 17:
                    continue
                locus, gid = c[16].strip(), c[15].strip()
                if locus and gid:
                    y2g.setdefault(locus.replace("YALI1_", "YALI1"), gid)

    # GeneIDs present in KEGG yli
    kegg_ids = set()
    kg = os.path.join(root, _KEGG_GENES)
    if os.path.exists(kg):
        with open(kg) as f:
            for line in f:
                kegg_ids.add(line.split("\t", 1)[0].replace("yli:", ""))

    genes = sorted({r["add_gene"] for r in csv.DictReader(open(additions_csv))})
    annotated = 0
    for g in genes:
        try:
            gene = model.genes.get_by_id(g)
        except KeyError:
            continue
        ann = dict(gene.annotation) if gene.annotation else {}
        if ann.get("sbo") == "SBO:0000243":
            continue  # idempotent
        ann["sbo"] = "SBO:0000243"
        gid = y2g.get(g)
        if gid:
            ann["ncbigene"] = gid
            if gid in kegg_ids:
                ann["kegg.genes"] = f"yli:{gid}"
            if network:
                up = _fetch_uniprot_for_geneid(gid)
                if up:
                    ann["uniprot"] = up
        gene.annotation = ann
        annotated += 1
        logger.debug(f"  gene annotate: {g} ncbigene={ann.get('ncbigene','-')} "
                     f"uniprot={ann.get('uniprot','-')}")
    return annotated


# ── Patch 10: fill formulas for definite neutral metabolites ──────────────

# A set of metabolites carried no formula but have a definite chemical
# identity.  Their MetaNetX (metanetx.chemical) ids turned out to be empty
# shells with no formula in chem_prop, so the formula was looked up by name in
# PubChem instead (scripts/audit_missing_formula.py ->
# data/missing_formula_fill.csv).
#
# This patch applies ONLY the safe subset: metabolites that are neutral in
# their physiological state (terpenes, esters, amines, nitriles, alcohols,
# aldehydes, peroxide) where charge = 0 is unambiguous.  Charged species
# (carboxylates, CoA-thioesters, phosphate esters, dipeptides) are deliberately
# excluded — the model's own charge convention for those is internally
# inconsistent (e.g. UDP-glucose charge 0 vs GDP-mannose -2), so there is no
# single correct value to fill and getting it wrong would unbalance reactions.
#
# Model formula convention is neutral-H + separate charge (H_model = H_neutral),
# so the PubChem neutral formula is exactly what the model wants.  Idempotent:
# metabolites that already have a formula are skipped.
_NEUTRAL_FILL_CSV = "data/missing_formula_fill.csv"


def fill_neutral_formulas(model, fill_csv: str | None = None) -> int:
    """
    Fill formula (and charge 0) for definite-neutral metabolites listed in the
    fill CSV with status=ready and charge=0, excluding dipeptides.

    Matches by metabolite name; only fills metabolites whose formula is
    currently empty.  Returns the number of metabolite copies filled.
    """
    import csv
    import os
    import re as _re

    if fill_csv is None:
        root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        fill_csv = os.path.join(root, _NEUTRAL_FILL_CSV)
    if not os.path.exists(fill_csv):
        logger.warning(f"  neutral fill: {fill_csv} not found — skipping")
        return 0

    dipeptide = _re.compile(
        r"^(gly|ala|cys|pro|asp|glu|ser|thr|val|leu|ile|phe|tyr|trp|his|lys|"
        r"arg|asn|gln|met)[-_]", _re.I)

    targets = {}  # name -> formula
    with open(fill_csv) as f:
        for row in csv.DictReader(f):
            if row.get("status") != "ready" or row.get("charge") != "0":
                continue
            if dipeptide.match(row["name"]):
                continue
            if row["formula"] and row["formula"] != "(查不到)":
                targets[row["name"]] = row["formula"]

    filled = 0
    for met in model.metabolites:
        if met.formula:
            continue
        f = targets.get(met.name)
        if not f:
            continue
        met.formula = f
        met.charge = 0
        filled += 1
        logger.debug(f"  neutral fill: {met.id} ({met.name}) <- {f} charge=0")
    return filled


# ── Lipid chain-menu extension: add C16:1 palmitoleoyl-CoA to acyl-CoA pools ──

# Y. lipolytica W29 makes ~8% palmitoleate (C16:1) but the iYli21 acyl-CoA pools
# omit it (Carsanba 2020, Table 3: C16:1 = 8.3% on glucose, day 2 —
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7409262/, verified 2026-06).
# We add palmitoleoyl-CoA to the 3 acyl-CoA pools (xPOOL_AC_EM/LP/MM), giving it
# 8.3% of the pool and scaling the existing 6 chains by (1 - 0.083) so the substrate
# weight sum stays 0.951 (= the unchanged product coefficient). palmitoleoyl-CoA
# already exists in all three compartments, so no new metabolite is created.
# Fatty-acid pools (xPOOL_FA_*) are intentionally NOT touched here.
_AC_POOL_C161_FRACTION = 0.083                       # C16:1 share within the pool
_AC_POOLS = ("xPOOL_AC_EM", "xPOOL_AC_LP", "xPOOL_AC_MM")
# palmitoleoyl-CoA id per acyl-CoA-pool compartment (verified present in model)
_PALMITOLEOYL_COA = {"C_em": "m243[C_em]", "C_lp": "m1486[C_lp]", "C_mm": "m1624[C_mm]"}


def extend_acyl_pool_c161(model) -> int:
    """
    Add C16:1 palmitoleoyl-CoA to the 3 acyl-CoA pools and re-scale the existing
    6 chains so the substrate weight sum (= product coefficient, 0.951) is unchanged.

    Idempotent: a pool that already contains its palmitoleoyl-CoA is skipped.
    Returns the number of pools extended. Does not create metabolites and does not
    touch the fatty-acid pools.
    """
    extended = 0
    for rid in _AC_POOLS:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        comp = next(iter({m.compartment for m in rxn.metabolites}))
        c161_id = _PALMITOLEOYL_COA.get(comp)
        if c161_id is None:
            continue
        c161 = model.metabolites.get_by_id(c161_id)
        if c161 in rxn.metabolites:          # already extended → idempotent skip
            continue

        # substrate weight sum (must equal the product coefficient, preserved)
        sub_sum = sum(-c for m, c in rxn.metabolites.items() if c < 0)
        scale = 1.0 - _AC_POOL_C161_FRACTION
        # scale the 6 existing acyl-CoA substrates in place
        delta = {}
        for met, coef in rxn.metabolites.items():
            if coef < 0:
                delta[met] = coef * scale - coef     # bring coef to coef*scale
        rxn.add_metabolites(delta)
        # add palmitoleoyl-CoA at its share of the (unchanged) total
        rxn.add_metabolites({c161: -sub_sum * _AC_POOL_C161_FRACTION})
        extended += 1
        logger.info("  C16:1 pool extension: %s += %s (%.4f), 6 chains x%.3f"
                    % (rid, c161_id, sub_sum * _AC_POOL_C161_FRACTION, scale))
    return extended


# ── Charge fix, Stage 1: free monocarboxylate/sulfonate anions left at charge 0 ──
#
# These 4 metabolites store the DEPROTONATED (anion) formula but were left at
# charge 0 — a missing-charge bug, not a wrong formula. Setting charge = -1 makes
# each metabolite formula/charge self-consistent and balances the reactions they
# touch, WITHOUT any formula edit. This is the "provably-safe" subset from the
# charge-convention strategy review: each touches only reactions that become
# balanced (or were already mass-imbalanced for other reasons), breaking none.
#
# Verified anchors (pH ~7.3 major microspecies):
#   m965 taurocholic acid C26H44NO7S -> -1: ChEBI:36257 taurocholate(1-), formula
#        matches the model exactly (verified this session)
#        https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:36257
#   m1403 kynurenic acid C10H7NO3 -> -1: monocarboxylic acid, KEGG C01717 (verified)
#        https://www.kegg.jp/entry/C01717
#   m286 ureidoglycolic acid C3H6N2O4 -> -1: carboxylate, KEGG C00603 (verified)
#        https://www.kegg.jp/entry/C00603
#   m753 1-pyrroline-3-hydroxy-5-carboxylate C5H6NO3 -> -1: name is the carboxylate;
#        formula already deprotonated; in-model R1194 balances and R500 improves
#        (charge gap -3 -> -2), consistent with a -1 anion (inferred, model-self-consistent)
#
# Note: m1403 = -1 also adds a charge term to R1375 (spontaneous kynurenic ->
# quinaldic), but R1375 is ALREADY mass-imbalanced (drops C10H7NO3 from nothing) —
# a pre-existing modeling error, flagged separately for curation. The -1 charge is
# the chemically correct value; we do not freeze it to flatter a broken reaction.
_CHARGE_STAGE1 = {
    "m286[C_cy]":  -1,   # ureidoglycolic acid
    "m753[C_mi]":  -1,   # 1-pyrroline-3-hydroxy-5-carboxylate
    "m965[C_cy]":  -1,   # taurocholic acid
    "m1403[C_cy]": -1,   # kynurenic acid
}


def fix_charge_stage1(model) -> int:
    """Set 4 free anionic metabolites from charge 0 to -1 (formula unchanged).

    Idempotent: a metabolite already at the target charge is skipped. Returns the
    number of metabolites changed. See _CHARGE_STAGE1 for per-metabolite sources.
    """
    changed = 0
    for mid, target in _CHARGE_STAGE1.items():
        try:
            met = model.metabolites.get_by_id(mid)
        except KeyError:
            continue
        if met.charge == target:
            continue
        old = met.charge
        met.charge = target
        changed += 1
        logger.info("  Charge Stage 1: %s (%s) charge %s -> %d"
                    % (mid, met.name, old, target))
    return changed


# ── Charge fix, Stage 2: anion-stored metabolites identified by InChI dH ──────
#
# Found by the InChI-dH discriminator (full-model scan): a metabolite whose stored
# formula has FEWER H than its own embedded-InChI neutral formula (dH < 0) is already
# the deprotonated anion, so setting charge = dH leaves the (formula, charge) pair a
# real protonation microspecies WITHOUT any formula edit. This avoids the FMN/FAD
# trap (those are dH == 0 = neutral-stored, and are deliberately NOT touched here).
#
# All 7 below were reaction-safety-gated on a cold model reload: balanced 914 -> 921
# (+7: R313/R322/R344/R567/R764/R768/R769), 0 reactions broken (full mass+charge).
#
# The 3 farnesyl-diphosphate copies MUST be changed together: R1116 is an FPP
# cross-membrane transport (m512[C_cy] <=> m610[C_mi]) that only balances when both
# ends carry the same charge.
#
# Verified charges (pH ~7.3 major microspecies, opened this session):
#   glyceraldehyde-3-phosphate(2-) C3H5O6P  -> -2  ChEBI:59776 (verified, formula matches)
#       https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:59776
#   farnesyl diphosphate(3-)       C15H25O7P2 -> -3  ChEBI:175763 (verified, formula matches)
#       https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:175763
#   ADP-ribose                     C15H21N5O14P2 -> -2  (InChI dH = -2, model-embedded InChI)
#   1-pyrroline-3-OH-5-carboxylate C5H6NO3       -> -1  (InChI dH = -1; sibling m753 already -1)
#
# NOT included: m966[C_va] vacuolar taurocholate — chemically -1 but R746 needs an
# H+ inserted first (see the r746-taurocholate-m966-needs-hplus note); deferred.
_CHARGE_STAGE2 = {
    "m536[C_cy]":  -2,   # D-glyceraldehyde 3-phosphate
    "m818[C_nu]":  -2,   # ADP-ribose (nucleus)
    "m1725[C_cy]": -2,   # ADP-ribose (cytosol)
    "m512[C_cy]":  -3,   # farnesyl diphosphate  ┐
    "m610[C_mi]":  -3,   # farnesyl diphosphate  ├ change together (R1116 symmetry)
    "m411[C_lp]":  -3,   # farnesyl diphosphate  ┘
    "m764[C_cy]":  -1,   # 1-pyrroline-3-hydroxy-5-carboxylate
}


def fix_charge_stage2(model) -> int:
    """Set 7 anion-stored metabolites (InChI dH < 0) to their dH charge (formula unchanged).

    Idempotent. Returns the number of metabolites changed. See _CHARGE_STAGE2 for sources.
    """
    changed = 0
    for mid, target in _CHARGE_STAGE2.items():
        try:
            met = model.metabolites.get_by_id(mid)
        except KeyError:
            continue
        if met.charge == target:
            continue
        old = met.charge
        met.charge = target
        changed += 1
        logger.info("  Charge Stage 2: %s (%s) charge %s -> %d"
                    % (mid, met.name, old, target))
    return changed


# ── Top-level driver ──────────────────────────────────────────────────────

def apply_all_patches(model) -> dict:
    """
    Apply all known model patches. Returns a dict with per-patch counts.
    Safe to call multiple times (idempotent — already-correct values are skipped).
    """
    logger.info("Applying iYli21 known-bug patches …")
    # NOTE: fix_ec_code_format is intentionally NOT called here.  Reaction EC
    # codes are populated later in the pipeline (gene EC enrichment, EC
    # backfill, reaction xref backfill), so EC formatting must run after those
    # steps — it is invoked separately near the end of main().
    # NOTE: fix_coa_charge is intentionally NOT called.  Setting free CoA to
    # charge -4 in isolation unbalances every reaction that produces/consumes
    # it (the other side is not adjusted), which regressed Memote's
    # reaction_charge_balance (+10 offenders).  The function is kept for
    # reference but disabled until a charge fix is applied consistently across
    # whole reactions rather than to the metabolite alone.
    counts = {
        "nadp_plus_fixed":   fix_nadp_plus_formula(model),
        "ceramide_fixed":    fix_ceramide_formulas(model),
        "cation_formula_fixed": fix_cation_formula_consistency(model),
    }
    logger.info(
        f"  patches applied: NADP+={counts['nadp_plus_fixed']} copies, "
        f"ceramide={counts['ceramide_fixed']} copies, "
        f"cation-formula={counts['cation_formula_fixed']} copies"
    )
    return counts
