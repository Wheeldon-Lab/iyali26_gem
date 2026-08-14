"""Matched B-group/PO1f attribution for uracil-salvage essentiality calls.

All models are disposable in-memory copies.  The script writes reports only;
it never updates model.xml, reaction bounds/GPRs in the release, or the durable
FN ledger/dossiers.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Callable

import pandas as pd
from cobra.flux_analysis import pfba
from cobra.io import read_sbml_model
from cobra.manipulation.delete import knock_out_model_genes

from .config import load_project_paths
from .essentiality_simulation_context import (
    load_effective_simulation_context,
    sha256_payload,
)
from .trna_biomass_pipeline import restore_free_amino_acid_biomass
from .validate_essential_genes import (
    DEFAULT_CUTOFFS,
    PRIMARY_CUTOFF,
    build_per_gene_table,
    load_experimental,
    run_single_gene_deletions,
)


REFERENCE_SCENARIO = "B_PO1f"
SCREEN_SCENARIOS = (
    ("A_W29", False, True, "free_amino_acid_biomass", "W29", ()),
    ("B_W29", False, False, "split_aa_trna_biomass", "W29", ()),
    ("A_PO1f", True, True, "free_amino_acid_biomass", "PO1f", ()),
    ("B_PO1f", True, False, "split_aa_trna_biomass", "PO1f", ()),
    (
        "B_PO1f_no_R935_GPR",
        True,
        False,
        "split_aa_trna_biomass",
        "PO1f",
        ("detach_R935_GPR",),
    ),
    (
        "B_PO1f_R1308_forward_only",
        True,
        False,
        "split_aa_trna_biomass",
        "PO1f",
        ("R1308_forward_only",),
    ),
    (
        "B_PO1f_R1927_closed",
        True,
        False,
        "split_aa_trna_biomass",
        "PO1f",
        ("close_R1927",),
    ),
)

REACTION_AUDIT = {
    "YALI1D07232g": {
        "established_name": None,
        "protein_function": "predicted NCS1-family nucleobase/solute:cation symporter",
        "evidence_status": "uncharacterized; family-level computational annotation only",
        "crosswalk": {
            "YALI1": "YALI1_D07232g",
            "YALI0": "YALI0D05621g",
            "YALI2": "YALI2_D01012g",
            "NCBI_Gene": "2911032",
            "RefSeq": "XP_502451.1",
            "UniProt": "Q6CA61",
        },
        "model_roles": ["R819 allantoin transport", "R935 uracil:H+ symport", "R937 uridine:H+ symport"],
        "decision": "retain provisional GPR; do not rename, delete, or expand",
        "references": [
            "https://www.uniprot.org/uniprotkb/Q6CA61/entry",
            "https://www.ncbi.nlm.nih.gov/gene/2911032",
            "https://www.kegg.jp/entry/yli:2911032",
            "https://publication-theses.unistra.fr/public/theses_doctorat/2008/FRITSCH_Emilie_2008.pdf",
            "https://www.nature.com/articles/s42003-023-04996-8",
        ],
    },
    "R935": {
        "equation": "H+[C_ex] + uracil[C_ex] -> H+[C_cy] + uracil[C_cy]",
        "bounds": [0.0, 1000.0],
        "compartments": ["C_ex", "C_cy"],
        "gpr": "YALI1D07232g",
        "assessment": "direction and compartments are plausible; exact substrate specificity and proton stoichiometry remain provisional",
        "references": ["https://www.rhea-db.org/rhea/29239"],
    },
    "R1308": {
        "equation": "uridine[C_cy] + phosphate[C_cy] <=> uracil[C_cy] + alpha-D-ribose-1-phosphate[C_cy]",
        "gpr": "",
        "assessment": "generic reversible chemistry, but no Yarrowia enzyme identity or SD-Leu reverse-direction evidence",
        "references": ["https://www.rhea-db.org/rhea/24388"],
    },
    "R1927": {
        "gpr": "YALI1E36477g",
        "gene_identity": "YALI1E36477g — no established gene name — uridine kinase EC 2.7.1.48 (curated annotation)",
        "assessment": "one of nine interchangeable model reactions assigned to the same gene; not a unique bypass step",
    },
}


@dataclass
class ScenarioResult:
    name: str
    per_gene: pd.DataFrame
    wt_growth: float
    context: dict[str, object]
    operations: tuple[str, ...]
    split_inverse_audit: list[dict[str, object]]


def _apply_operations(model, operations: tuple[str, ...]) -> None:
    for operation in operations:
        if operation == "detach_R935_GPR":
            model.reactions.get_by_id("R935").gene_reaction_rule = ""
        elif operation == "R1308_forward_only":
            model.reactions.get_by_id("R1308").lower_bound = 0.0
        elif operation == "close_R1927":
            model.reactions.get_by_id("R1927").bounds = (0.0, 0.0)
        else:
            raise ValueError(f"Unknown attribution operation: {operation}")


def _run_scenario(
    *,
    name: str,
    use_profile: bool,
    restore_a: bool,
    operations: tuple[str, ...],
    model_path: Path,
    media_path: Path,
    profile_path: Path,
    experimental: pd.DataFrame,
    solver: str,
) -> ScenarioResult:
    context = load_effective_simulation_context(
        model_path=model_path,
        media_path=media_path,
        strain_profile_path=profile_path if use_profile else None,
    )
    model = context.model
    inverse_audit = restore_free_amino_acid_biomass(model) if restore_a else []
    _apply_operations(model, operations)
    predictions, wt_growth = run_single_gene_deletions(
        model, solver, excluded_gene_ids=context.excluded_runtime_genes
    )
    per_gene = build_per_gene_table(
        experimental, predictions, DEFAULT_CUTOFFS, PRIMARY_CUTOFF
    )
    provenance = {
        **context.provenance(),
        "canonical_model_sha256": context.canonical_model_sha256,
        "medium_sha256": context.medium_sha256,
        "scenario_fingerprint": sha256_payload(
            {
                "base_simulation_context_fingerprint": context.simulation_context_fingerprint,
                "scenario": name,
                "restore_free_amino_acid_biomass": restore_a,
                "operations": list(operations),
            }
        ),
    }
    return ScenarioResult(name, per_gene, wt_growth, provenance, operations, inverse_audit)


def _call(row: pd.Series, cutoff: float) -> str:
    if not bool(row["experimental_essential"]):
        return "outside_positive_reference"
    return "TP" if float(row["ko_growth_ratio"]) < cutoff else "FN"


def _threshold_calls(results: dict[str, ScenarioResult]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for scenario, result in results.items():
        positives = result.per_gene[
            result.per_gene["experimental_essential"].eq(True)
            & result.per_gene["in_model"].eq(True)
        ]
        for cutoff in DEFAULT_CUTOFFS:
            for row in positives.itertuples(index=False):
                rows.append(
                    {
                        "scenario": scenario,
                        "scenario_fingerprint": result.context["scenario_fingerprint"],
                        "cutoff_fraction_of_wt": cutoff,
                        "gene_id": row.gene_id,
                        "source_gene_id": row.source_gene_id,
                        "function": row.function,
                        "in_model": bool(row.in_model),
                        "ko_growth_ratio": row.ko_growth_ratio,
                        "call": "TP" if row.ko_growth_ratio < cutoff else "FN",
                        "exploratory_only": scenario != REFERENCE_SCENARIO,
                    }
                )
    return pd.DataFrame(rows)


def _scenario_deltas(results: dict[str, ScenarioResult]) -> pd.DataFrame:
    contrasts = {
        "B_group_in_W29": ("A_W29", "B_W29", "B_group"),
        "B_group_in_PO1f": ("A_PO1f", "B_PO1f", "B_group"),
        "PO1f_in_A": ("A_W29", "A_PO1f", "PO1f_background"),
        "PO1f_in_B": ("B_W29", "B_PO1f", "PO1f_background"),
        "R935_GPR_dependency": ("B_PO1f", "B_PO1f_no_R935_GPR", "suspect_GPR"),
        "R1308_reverse_dependency": ("B_PO1f", "B_PO1f_R1308_forward_only", "GPR_less_bypass"),
        "R1927_specificity_control": ("B_PO1f", "B_PO1f_R1927_closed", "solver_basis_control"),
    }
    rows: list[dict[str, object]] = []
    for contrast, (source, target, mechanism) in contrasts.items():
        source_table = results[source].per_gene.set_index("gene_id")
        target_table = results[target].per_gene.set_index("gene_id")
        for cutoff in DEFAULT_CUTOFFS:
            for gene_id in sorted(set(source_table.index) & set(target_table.index)):
                left = source_table.loc[gene_id]
                right = target_table.loc[gene_id]
                if not bool(left["experimental_essential"]):
                    continue
                before, after = _call(left, cutoff), _call(right, cutoff)
                rows.append(
                    {
                        "contrast": contrast,
                        "mechanism": mechanism,
                        "from_scenario": source,
                        "to_scenario": target,
                        "cutoff_fraction_of_wt": cutoff,
                        "gene_id": gene_id,
                        "function": left["function"],
                        "from_ko_growth_ratio": left["ko_growth_ratio"],
                        "to_ko_growth_ratio": right["ko_growth_ratio"],
                        "from_call": before,
                        "to_call": after,
                        "transition": f"{before}->{after}",
                        "changed": before != after,
                    }
                )
    return pd.DataFrame(rows)


def _attribution_table(results: dict[str, ScenarioResult]) -> pd.DataFrame:
    tables = {name: result.per_gene.set_index("gene_id") for name, result in results.items()}
    rows: list[dict[str, object]] = []
    for cutoff in DEFAULT_CUTOFFS:
        for gene_id in sorted(tables[REFERENCE_SCENARIO].index):
            reference = tables[REFERENCE_SCENARIO].loc[gene_id]
            if not bool(reference["experimental_essential"]) or not bool(reference["in_model"]):
                continue
            calls = {name: _call(table.loc[gene_id], cutoff) for name, table in tables.items()}
            rows.append(
                {
                    "gene_id": gene_id,
                    "source_gene_id": reference["source_gene_id"],
                    "function": reference["function"],
                    "cutoff_fraction_of_wt": cutoff,
                    **{f"call_{name}": call for name, call in calls.items()},
                    "changed_by_B_group": (
                        calls["A_W29"] != calls["B_W29"]
                        or calls["A_PO1f"] != calls["B_PO1f"]
                    ),
                    "changed_by_PO1f": (
                        calls["A_W29"] != calls["A_PO1f"]
                        or calls["B_W29"] != calls["B_PO1f"]
                    ),
                    "depends_on_R935_GPR": calls["B_PO1f"] != calls["B_PO1f_no_R935_GPR"],
                    "depends_on_R1308_reverse": calls["B_PO1f"] != calls["B_PO1f_R1308_forward_only"],
                    "R1927_only_effect": calls["B_PO1f"] != calls["B_PO1f_R1927_closed"],
                    "B_by_PO1f_interaction": (
                        (calls["B_PO1f"] != calls["B_W29"])
                        != (calls["A_PO1f"] != calls["A_W29"])
                    ),
                }
            )
    return pd.DataFrame(rows)


def _mechanism_audit(model_path: Path, media_path: Path, profile_path: Path, solver: str) -> pd.DataFrame:
    context = load_effective_simulation_context(
        model_path=model_path, media_path=media_path, strain_profile_path=profile_path
    )
    model = context.model
    model.solver = solver
    growth = float(model.slim_optimize(error_value=0.0) or 0.0)
    solution = pfba(model)
    rows = [
        {
            "analysis": "baseline_pfba_flux",
            "perturbation": "none",
            "growth": growth,
            "reaction_id": reaction_id,
            "flux": float(solution.fluxes.get(reaction_id, 0.0)),
            "interpretation": "baseline PO1f effective-model flux",
        }
        for reaction_id in ("R935", "R1308", "R783", "R2085", "R786", "R1921", "R1927")
    ]
    for direct_closed in (False, True):
        for reverse_blocked in (False, True):
            probe = model.copy()
            if direct_closed:
                for reaction_id in ("R783", "R2085"):
                    probe.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)
            if reverse_blocked:
                probe.reactions.get_by_id("R1308").lower_bound = 0.0
            growth = float(probe.slim_optimize(error_value=0.0) or 0.0)
            rows.append(
                {
                    "analysis": "salvage_2x2",
                    "perturbation": f"direct_UPRT_closed={direct_closed};R1308_reverse_blocked={reverse_blocked}",
                    "growth": growth,
                    "reaction_id": "",
                    "flux": None,
                    "interpretation": "direct UPRT arm versus GPR-less reverse-R1308 arm",
                }
            )
    for label, genes, reverse_blocked in (
        ("FUR1_KO", ("YALI1F38986g",), False),
        ("UPRT_double_KO", ("YALI1F38986g", "YALI1E24479g"), False),
        ("UPRT_double_KO_R1308_forward_only", ("YALI1F38986g", "YALI1E24479g"), True),
        ("UPRT_double_KO_uridine_kinase_KO", ("YALI1F38986g", "YALI1E24479g", "YALI1E36477g"), False),
    ):
        probe = model.copy()
        if reverse_blocked:
            probe.reactions.get_by_id("R1308").lower_bound = 0.0
        knock_out_model_genes(probe, [probe.genes.get_by_id(gene_id) for gene_id in genes])
        rows.append(
            {
                "analysis": "gene_epistasis",
                "perturbation": label,
                "growth": float(probe.slim_optimize(error_value=0.0) or 0.0),
                "reaction_id": "",
                "flux": None,
                "interpretation": "FUR1 rescue requires either direct UPRT alternative or R1308 reverse plus uridine kinase",
            }
        )
    return pd.DataFrame(rows)


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _materialize_reaction_audit(model_path: Path) -> dict[str, object]:
    """Attach actual current-model reaction fields to curated evidence notes."""
    model = read_sbml_model(str(model_path))
    audit = json.loads(json.dumps(REACTION_AUDIT))
    for reaction_id in ("R935", "R1308", "R1927"):
        reaction = model.reactions.get_by_id(reaction_id)
        audit[reaction_id].update(
            {
                "current_model_equation": reaction.reaction,
                "current_model_bounds": [
                    float(reaction.lower_bound),
                    float(reaction.upper_bound),
                ],
                "current_model_compartments": sorted(
                    {metabolite.compartment for metabolite in reaction.metabolites}
                ),
                "current_model_gpr": reaction.gene_reaction_rule,
            }
        )
    return audit


def run_attribution(
    *, model_path: Path, experimental_path: Path, media_path: Path, profile_path: Path,
    output_dir: Path, solver: str,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=False)
    experimental = load_experimental(experimental_path, positive_only=True)
    results: dict[str, ScenarioResult] = {}
    for name, use_profile, restore_a, _biomass, _strain, operations in SCREEN_SCENARIOS:
        results[name] = _run_scenario(
            name=name, use_profile=use_profile, restore_a=restore_a, operations=operations,
            model_path=model_path, media_path=media_path, profile_path=profile_path,
            experimental=experimental, solver=solver,
        )
    threshold_calls = _threshold_calls(results)
    deltas = _scenario_deltas(results)
    attribution = _attribution_table(results)
    mechanism = _mechanism_audit(model_path, media_path, profile_path, solver)
    threshold_calls.to_csv(output_dir / "essentiality_threshold_calls.tsv", sep="\t", index=False)
    deltas.to_csv(output_dir / "essentiality_scenario_deltas.tsv", sep="\t", index=False)
    attribution.to_csv(output_dir / "essentiality_attribution.tsv", sep="\t", index=False)
    mechanism.to_csv(output_dir / "essentiality_mechanism_audit.tsv", sep="\t", index=False)
    fn_sets: dict[str, object] = {}
    for cutoff in DEFAULT_CUTOFFS:
        label = f"{int(cutoff * 100)}pct"
        rows = threshold_calls[(threshold_calls["cutoff_fraction_of_wt"] == cutoff) & (threshold_calls["call"] == "FN")]
        rows.to_csv(output_dir / f"essentiality_fn_{label}.tsv", sep="\t", index=False)
        fn_sets[label] = {
            scenario: sorted(group["gene_id"].tolist())
            for scenario, group in rows.groupby("scenario")
        }
        for scenario, genes in fn_sets[label].items():
            fn_sets[label][scenario] = {"gene_ids": genes, "set_sha256": sha256_payload(genes)}
    _write_json(output_dir / "essentiality_fn_sets.json", fn_sets)
    _write_json(output_dir / "gene_reaction_audit.json", _materialize_reaction_audit(model_path))
    manifest = {
        "schema_version": "1.0",
        "analysis": "PO1f_uracil_salvage_FN_attribution",
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "model_xml_written": False,
        "durable_ledger_written": False,
        "reference_scenario": REFERENCE_SCENARIO,
        "durable_case_eligible_scenario": REFERENCE_SCENARIO,
        "counterfactuals_exploratory_only": True,
        "inputs": {
            "model": str(model_path.resolve()),
            "experimental": str(experimental_path.resolve()),
            "medium": str(media_path.resolve()),
            "PO1f_profile": str(profile_path.resolve()),
        },
        "scenarios": {
            name: {
                "wt_growth": result.wt_growth,
                "operations": list(result.operations),
                "split_inverse_audit": result.split_inverse_audit,
                "simulation_context": result.context,
                "exploratory_only": name != REFERENCE_SCENARIO,
            }
            for name, result in results.items()
        },
        "outputs": {
            "threshold_calls": "essentiality_threshold_calls.tsv",
            "scenario_deltas": "essentiality_scenario_deltas.tsv",
            "attribution": "essentiality_attribution.tsv",
            "mechanism_audit": "essentiality_mechanism_audit.tsv",
            "fn_sets": "essentiality_fn_sets.json",
            "gene_reaction_audit": "gene_reaction_audit.json",
        },
    }
    _write_json(output_dir / "attribution_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--solver", default="gurobi")
    args = parser.parse_args()
    paths = load_project_paths(args.research_root, required=True)
    paths.require(paths.output_model, paths.essentiality, paths.media, paths.strain_profiles)
    output_dir = args.output_dir or (
        paths.results / "essentiality" / "po1f_uracil_salvage_attribution" /
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    )
    manifest = run_attribution(
        model_path=paths.output_model,
        experimental_path=paths.essentiality / "consensus_essential_genes.csv",
        media_path=paths.media / "sd_leu.csv",
        profile_path=paths.strain_profiles / "po1f_sd_leu.json",
        output_dir=output_dir,
        solver=args.solver,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
