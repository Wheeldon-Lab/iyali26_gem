"""
validate_essential_genes.py — Gene essentiality validation for iYli21 GEM.

Compares model-predicted essential genes (via single-gene deletion FBA) against
experimental data, then reports a confusion matrix and ranked FN/FP lists that
point directly to gaps needing manual curation.

Usage
-----
    python scripts/validate_essential_genes.py \\
        --experimental data/essentiality/essential_genes.csv \\
        --model model.xml \\
        [--output results/essentiality_report.tsv] \\
        [--growth-cutoff 0.01] \\
        [--solver glpk]

Experimental CSV format
-----------------------
Two required columns (additional columns are ignored):

    gene_id,essential
    YALI1C08548g,1
    YALI1D11769g,0
    ...

    essential: 1 / True / yes  → gene is experimentally essential
               0 / False / no  → gene is experimentally non-essential

Output
------
Prints to stdout:
  - Wild-type growth rate
  - Confusion matrix (TP / FP / TN / FN) with MCC and accuracy
  - False Negatives: model predicts viable but experiment says essential
    (each line: gene_id, affected reactions — shows WHERE the model is wrong)
  - False Positives: model predicts lethal but experiment says non-essential
    (each line: gene_id, affected reactions — possible redundancy not modelled)

Writes TSV with full per-gene predictions if --output is specified.
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import single_gene_deletion

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent


# ─────────────────────────────────────────────
# Parsing
# ─────────────────────────────────────────────

def load_experimental(path: Path) -> pd.DataFrame:
    """
    Load experimental essentiality CSV.
    Returns DataFrame with columns [gene_id, essential (bool)].
    Accepts 1/0, True/False, yes/no (case-insensitive) in the essential column.
    """
    df = pd.read_csv(path, dtype=str).rename(columns=str.strip)
    required = {"gene_id", "essential"}
    missing = required - set(df.columns.str.lower())
    if missing:
        raise ValueError(
            f"Experimental CSV missing required columns: {missing}. "
            f"Found: {list(df.columns)}"
        )
    # normalise column names
    df.columns = df.columns.str.lower().str.strip()

    truthy  = {"1", "true", "yes", "essential", "y"}
    falsy   = {"0", "false", "no", "non-essential", "nonessential", "n"}

    def parse_bool(v: str) -> bool:
        v = str(v).strip().lower()
        if v in truthy:
            return True
        if v in falsy:
            return False
        raise ValueError(f"Unrecognised essentiality value: {v!r}")

    df["essential"] = df["essential"].map(parse_bool)
    df["gene_id"]   = df["gene_id"].str.strip()
    return df[["gene_id", "essential"]].drop_duplicates("gene_id")


# ─────────────────────────────────────────────
# Prediction
# ─────────────────────────────────────────────

def predict_essential(model, gene_ids: list[str], growth_cutoff: float,
                      solver: str) -> pd.DataFrame:
    """
    Run single-gene deletion FBA for each gene in gene_ids.

    Returns DataFrame with columns:
      gene_id, wt_growth, ko_growth, predicted_essential
    """
    model.solver = solver

    wt_sol = model.optimize()
    if wt_sol.status != "optimal":
        raise RuntimeError(
            f"Wild-type FBA infeasible (status={wt_sol.status}). "
            "Check model bounds and objective."
        )
    wt_growth = wt_sol.objective_value
    logger.info(f"Wild-type growth rate: {wt_growth:.4f}")

    # Filter to genes actually in the model
    model_gene_ids = {g.id for g in model.genes}
    valid_ids = [g for g in gene_ids if g in model_gene_ids]
    missing   = [g for g in gene_ids if g not in model_gene_ids]
    if missing:
        logger.warning(f"{len(missing)} experimental genes not found in model: {missing[:10]}{'…' if len(missing)>10 else ''}")

    logger.info(f"Running single_gene_deletion for {len(valid_ids)} genes …")
    deletion_results = single_gene_deletion(
        model,
        gene_list=valid_ids,
    )
    # deletion_results index is frozensets; flatten to gene_id
    rows = []
    for idx, row in deletion_results.iterrows():
        # idx is a frozenset with one gene id
        gid = next(iter(idx))
        ko_growth = row["growth"] if row["status"] == "optimal" else 0.0
        predicted_essential = (ko_growth is None) or (ko_growth < growth_cutoff * wt_growth)
        rows.append({
            "gene_id": gid,
            "wt_growth": wt_growth,
            "ko_growth": ko_growth if ko_growth is not None else 0.0,
            "predicted_essential": bool(predicted_essential),
        })

    return pd.DataFrame(rows), wt_growth


# ─────────────────────────────────────────────
# Confusion matrix + metrics
# ─────────────────────────────────────────────

def confusion_matrix(merged: pd.DataFrame) -> dict:
    """
    merged must have columns: essential (bool), predicted_essential (bool).
    Returns dict with TP, FP, TN, FN, accuracy, MCC.
    """
    tp = int(( merged["essential"] &  merged["predicted_essential"]).sum())
    fp = int((~merged["essential"] &  merged["predicted_essential"]).sum())
    tn = int((~merged["essential"] & ~merged["predicted_essential"]).sum())
    fn = int(( merged["essential"] & ~merged["predicted_essential"]).sum())

    accuracy = (tp + tn) / len(merged) if len(merged) else 0.0

    denom = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc   = (tp * tn - fp * fn) / denom if denom else 0.0

    return {"TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "accuracy": accuracy, "MCC": mcc}


# ─────────────────────────────────────────────
# Affected-reaction helper
# ─────────────────────────────────────────────

def _affected_reactions(model, gene_id: str) -> list[str]:
    """Return reaction IDs catalysed solely by this gene (GPR has no alternatives)."""
    gene = model.genes.get_by_id(gene_id)
    sole = []
    for rxn in gene.reactions:
        # knock out this gene; if reaction becomes inactive → it's solely dependent
        other_genes = rxn.genes - {gene}
        # simple check: no other gene associated → solely dependent
        if not other_genes:
            sole.append(rxn.id)
        # TODO: full GPR evaluation via rxn.gene_reaction_rule for AND/OR logic
    return sole


# ─────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────

def print_report(merged: pd.DataFrame, model, cm: dict, wt_growth: float,
                 growth_cutoff: float) -> None:
    sep = "─" * 60

    print(f"\n{sep}")
    print("  Gene Essentiality Validation Report")
    print(sep)
    print(f"  Wild-type growth rate : {wt_growth:.4f}")
    print(f"  Growth cutoff (% WT)  : {growth_cutoff*100:.1f}%")
    print(f"  Genes evaluated       : {len(merged)}")
    print()
    print("  Confusion Matrix")
    print(f"                       Predicted")
    print(f"                  Essential  Non-ess")
    print(f"  Exp Essential      {cm['TP']:5d}     {cm['FN']:5d}   ← FN = missed essentials")
    print(f"  Exp Non-essential  {cm['FP']:5d}     {cm['TN']:5d}   ← FP = spurious lethals")
    print()
    print(f"  Accuracy : {cm['accuracy']:.3f}")
    print(f"  MCC      : {cm['MCC']:.3f}   (−1 worst → 0 random → +1 perfect)")
    print(sep)

    # ── False Negatives ────────────────────────────────────────────────────
    fn_df = merged[merged["essential"] & ~merged["predicted_essential"]].copy()
    print(f"\n  FALSE NEGATIVES ({len(fn_df)})  — experimentally essential, model predicts viable")
    print("  These indicate missing reactions / blocked pathways in the model.")
    print(f"  {'Gene':<20} {'KO growth':>10}  Solely-dependent reactions")
    print(f"  {'─'*20} {'─'*10}  {'─'*30}")
    for _, row in fn_df.sort_values("ko_growth", ascending=False).iterrows():
        gid = row["gene_id"]
        affected = _affected_reactions(model, gid)
        rxn_str = ", ".join(affected[:5]) + ("…" if len(affected) > 5 else "")
        print(f"  {gid:<20} {row['ko_growth']:>10.4f}  {rxn_str or '(shared GPR)'}")

    # ── False Positives ────────────────────────────────────────────────────
    fp_df = merged[~merged["essential"] & merged["predicted_essential"]].copy()
    print(f"\n  FALSE POSITIVES ({len(fp_df)})  — model predicts lethal, experiment says viable")
    print("  These suggest missing isozymes / redundant pathways not yet in the model.")
    print(f"  {'Gene':<20} {'KO growth':>10}  Solely-dependent reactions")
    print(f"  {'─'*20} {'─'*10}  {'─'*30}")
    for _, row in fp_df.sort_values("ko_growth").iterrows():
        gid = row["gene_id"]
        affected = _affected_reactions(model, gid)
        rxn_str = ", ".join(affected[:5]) + ("…" if len(affected) > 5 else "")
        print(f"  {gid:<20} {row['ko_growth']:>10.4f}  {rxn_str or '(shared GPR)'}")

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

def validate_essential_genes(
    experimental_path: Path,
    model_path: Path,
    output_path: Path | None,
    growth_cutoff: float,
    solver: str,
) -> None:
    # Load experimental data
    logger.info(f"Loading experimental data: {experimental_path}")
    exp_df = load_experimental(experimental_path)
    logger.info(
        f"  {len(exp_df)} genes: "
        f"{exp_df['essential'].sum()} essential, "
        f"{(~exp_df['essential']).sum()} non-essential"
    )

    # Load model
    logger.info(f"Loading model: {model_path}")
    model = read_sbml_model(str(model_path))
    model.solver = solver

    # Run single-gene deletion
    pred_df, wt_growth = predict_essential(
        model,
        gene_ids=exp_df["gene_id"].tolist(),
        growth_cutoff=growth_cutoff,
        solver=solver,
    )

    # Merge experimental + predicted
    merged = exp_df.merge(pred_df, on="gene_id", how="inner")
    if len(merged) == 0:
        logger.error("No gene IDs matched between experimental data and model.")
        sys.exit(1)
    logger.info(f"Matched {len(merged)} genes for evaluation")

    # Confusion matrix
    cm = confusion_matrix(merged)

    # Print report
    print_report(merged, model, cm, wt_growth, growth_cutoff)

    # Optional TSV output
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_df = merged[["gene_id", "essential", "predicted_essential",
                          "wt_growth", "ko_growth"]].copy()
        out_df["TP"] = out_df["essential"] &  out_df["predicted_essential"]
        out_df["FP"] = ~out_df["essential"] &  out_df["predicted_essential"]
        out_df["TN"] = ~out_df["essential"] & ~out_df["predicted_essential"]
        out_df["FN"] =  out_df["essential"] & ~out_df["predicted_essential"]
        out_df.to_csv(output_path, sep="\t", index=False)
        logger.info(f"Full results written to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Validate GEM gene essentiality predictions against experimental data."
    )
    parser.add_argument(
        "--experimental", "-e", required=True, type=Path,
        help="CSV file with columns: gene_id, essential (1/0)"
    )
    parser.add_argument(
        "--model", "-m", type=Path,
        default=REPO_ROOT / "model.xml",
        help="SBML model file (default: model.xml)"
    )
    parser.add_argument(
        "--output", "-o", type=Path, default=None,
        help="Optional TSV output path for full per-gene predictions"
    )
    parser.add_argument(
        "--growth-cutoff", type=float, default=0.01,
        help="KO growth / WT growth threshold below which a gene is predicted essential (default: 0.01)"
    )
    parser.add_argument(
        "--solver", type=str, default="glpk",
        help="LP solver to use (default: glpk)"
    )
    args = parser.parse_args()

    if not args.experimental.exists():
        logger.error(f"Experimental file not found: {args.experimental}")
        sys.exit(1)
    if not args.model.exists():
        logger.error(f"Model file not found: {args.model}")
        sys.exit(1)

    validate_essential_genes(
        experimental_path=args.experimental,
        model_path=args.model,
        output_path=args.output,
        growth_cutoff=args.growth_cutoff,
        solver=args.solver,
    )


if __name__ == "__main__":
    main()
