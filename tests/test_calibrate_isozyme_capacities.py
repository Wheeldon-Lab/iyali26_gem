import csv
import json
from pathlib import Path

import pandas as pd
import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model, write_sbml_model

from scripts.gem_annotate import calibrate_isozyme_capacities as calibration
from scripts.gem_annotate.essentiality_evidence import sha256_file, target_fingerprint
from scripts.gem_annotate.provisional_capacity import load_provisional_capacity_table
from scripts.gem_annotate.validate_essential_genes import load_assay_fitness


REPO_ROOT = Path(__file__).resolve().parents[1]
GRID_PATH = (
    REPO_ROOT
    / "data"
    / "essentiality"
    / "scenarios"
    / "provisional_capacity_joint_grid.csv"
)


def _toy_model() -> Model:
    model = Model("toy_joint_capacity")
    a = Metabolite("a_c", compartment="c")
    x = Metabolite("x_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    c = Metabolite("c_c", compartment="c")
    waste = Metabolite("waste_c", compartment="c")

    ex_a = Reaction("EX_a")
    ex_a.bounds = (-1000.0, 1000.0)
    ex_a.add_metabolites({a: -1.0})

    ex_x = Reaction("EX_x")
    ex_x.bounds = (-1000.0, 1000.0)
    ex_x.add_metabolites({x: -1.0})

    r4 = Reaction("R4")
    r4.bounds = (0.0, 1000.0)
    r4.add_metabolites({a: -1.0, b: 1.0})
    r4.gene_reaction_rule = "g_r4_backup or g_r4_main"

    r1846 = Reaction("R1846")
    r1846.bounds = (0.0, 1000.0)
    r1846.add_metabolites({x: -1.0, c: 1.0})
    r1846.gene_reaction_rule = "g_r1846_backup_a or g_r1846_backup_b or g_r1846_main"

    # A fixed, small alternate contribution keeps the 5% R1846 grid point
    # comfortably inside the inclusive 5% gate despite LP feasibility tolerance.
    r1846_bypass = Reaction("R1846_TEST_BYPASS")
    r1846_bypass.bounds = (0.01, 0.01)
    r1846_bypass.add_metabolites({x: -1.0, c: 1.0})

    r1843 = Reaction("R1843")
    r1843.bounds = (0.0, 0.0)
    r1843.add_metabolites({b: -1.0, c: 1.0})
    r1843.gene_reaction_rule = "g_blocked"

    proxy = Reaction("R_PROXY")
    proxy.bounds = (0.0, 1000.0)
    proxy.add_metabolites({a: -1.0, waste: 1.0})
    proxy.gene_reaction_rule = "g_proxy"

    waste_drain = Reaction("DM_waste")
    waste_drain.bounds = (0.0, 1000.0)
    waste_drain.add_metabolites({waste: -1.0})

    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1000.0)
    biomass.add_metabolites({b: -1.0, c: -1.0})

    model.add_reactions(
        [
            ex_a,
            ex_x,
            r4,
            r1846,
            r1846_bypass,
            r1843,
            proxy,
            waste_drain,
            biomass,
        ]
    )
    model.objective = biomass
    return model


def _reaction_fingerprint(model: Model, reaction_id: str) -> str:
    reaction = model.reactions.get_by_id(reaction_id)
    return target_fingerprint(
        [
            {
                "reaction_id": reaction.id,
                "stoichiometry": {
                    metabolite.id: float(coefficient)
                    for metabolite, coefficient in reaction.metabolites.items()
                },
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
                "gpr": reaction.gene_reaction_rule,
            }
        ]
    )


def _write_profile(path: Path, model_path: Path) -> None:
    model = read_sbml_model(str(model_path))
    fieldnames = [
        "capacity_id",
        "status",
        "source_reaction_id",
        "expected_gpr",
        "primary_gpr",
        "backup_gpr",
        "backup_reaction_id",
        "provisional_upper_bound",
        "units",
        "parameter_basis",
        "case_id",
        "validated_model_sha256",
        "target_fingerprint",
        "requires_protein_abundance",
        "requires_kcat",
        "replacement_formula",
        "rationale",
    ]
    rows = [
        {
            "capacity_id": "PCAP-R4-TOY",
            "status": "active_exploratory",
            "source_reaction_id": "R4",
            "expected_gpr": "g_r4_backup or g_r4_main",
            "primary_gpr": "g_r4_main",
            "backup_gpr": "g_r4_backup",
            "backup_reaction_id": "R4__PCAP_BACKUP",
            "provisional_upper_bound": "0.1",
            "units": "mmol_gDW_h",
            "parameter_basis": "toy_not_measured",
            "case_id": "EGC-r4-toy",
            "validated_model_sha256": sha256_file(model_path),
            "target_fingerprint": _reaction_fingerprint(model, "R4"),
            "requires_protein_abundance": "true",
            "requires_kcat": "true",
            "replacement_formula": "sum_i(kcat_i*E_i*3600)",
            "rationale": "test only",
        },
        {
            "capacity_id": "PCAP-R1846-TOY",
            "status": "active_exploratory",
            "source_reaction_id": "R1846",
            "expected_gpr": ("g_r1846_backup_a or g_r1846_backup_b or g_r1846_main"),
            "primary_gpr": "g_r1846_main",
            "backup_gpr": "g_r1846_backup_a or g_r1846_backup_b",
            "backup_reaction_id": "R1846__PCAP_BACKUP",
            "provisional_upper_bound": "0.1",
            "units": "mmol_gDW_h",
            "parameter_basis": "toy_not_measured",
            "case_id": "EGC-r1846-toy",
            "validated_model_sha256": sha256_file(model_path),
            "target_fingerprint": _reaction_fingerprint(model, "R1846"),
            "requires_protein_abundance": "true",
            "requires_kcat": "true",
            "replacement_formula": "sum_i(kcat_i*E_i*3600)",
            "rationale": "test only",
        },
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_assay_fitness(path: Path, gene_ids: list[str]) -> None:
    fieldnames = [
        "gene_id",
        "source_gene_id",
        "assay",
        "fitness_score",
        "raw_p_value",
        "q_value",
        "experimental_call",
        "source_sheet",
        "source_row",
        "source_sha256",
    ]
    target_genes = {"g_r4_main", "g_r1846_main"}
    rows = []
    for assay in ("Cas9", "Cas12a"):
        for source_row, gene_id in enumerate(sorted(gene_ids), start=2):
            essential = gene_id in target_genes
            rows.append(
                {
                    "gene_id": gene_id,
                    "source_gene_id": gene_id,
                    "assay": assay,
                    "fitness_score": "-2" if essential else "0",
                    "raw_p_value": "0.001" if essential else "0.5",
                    "q_value": "0.01" if essential else "0.8",
                    "experimental_call": "essential" if essential else "nonessential",
                    "source_sheet": assay,
                    "source_row": source_row,
                    "source_sha256": "a" * 64,
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def toy_inputs(tmp_path: Path) -> dict[str, Path]:
    model_path = tmp_path / "model.xml"
    write_sbml_model(_toy_model(), str(model_path))
    profile_path = tmp_path / "profile.csv"
    _write_profile(profile_path, model_path)
    model = read_sbml_model(str(model_path))

    experimental_path = tmp_path / "experimental.csv"
    experimental_path.write_text("gene_id\ng_r4_main\ng_r1846_main\n")

    assay_path = tmp_path / "assay.csv"
    _write_assay_fitness(assay_path, sorted(gene.id for gene in model.genes))

    media_path = tmp_path / "media.csv"
    media_path.write_text("exchange,uptake\nEX_a,1\nEX_x,1\n")
    return {
        "model_path": model_path,
        "profile_path": profile_path,
        "grid_path": GRID_PATH,
        "experimental_path": experimental_path,
        "assay_fitness_path": assay_path,
        "media_path": media_path,
        "output_dir": tmp_path / "results",
    }


def test_approved_grid_is_exact_and_cli_defaults_to_gurobi() -> None:
    grid = calibration.load_joint_capacity_grid(GRID_PATH)

    assert len(grid) == 9
    assert set(grid["r4_capacity_fraction_of_wt_flux"]) == {0.025, 0.075, 0.15}
    assert set(grid["r1846_capacity_fraction_of_wt_flux"]) == {0.01, 0.025, 0.05}
    args = calibration.build_parser().parse_args([])
    explicit = calibration.build_parser().parse_args(
        ["run", "--growth-cutoffs", "0.01,0.05,0.10,0.15"]
    )
    assert args.solver == "gurobi"
    assert args.workers == 1
    assert explicit.command == "run"
    assert explicit.growth_cutoffs == (0.01, 0.05, 0.10, 0.15)


def test_ranking_uses_feasibility_margin_change_then_id() -> None:
    common = {
        "r4_capacity_fraction_of_wt_flux": 0.075,
        "r1846_capacity_fraction_of_wt_flux": 0.05,
    }
    ranked = calibration.rank_scenarios(
        [
            {
                **common,
                "scenario_id": "infeasible",
                "feasible": False,
                "boundary_margin": 0.04,
                "non_target_ratio_change_score": 0.0,
            },
            {
                **common,
                "scenario_id": "b",
                "feasible": True,
                "boundary_margin": 0.02,
                "non_target_ratio_change_score": 0.1,
            },
            {
                **common,
                "scenario_id": "a",
                "feasible": True,
                "boundary_margin": 0.02,
                "non_target_ratio_change_score": 0.1,
            },
            {
                **common,
                "scenario_id": "narrower",
                "feasible": True,
                "boundary_margin": 0.01,
                "non_target_ratio_change_score": 0.0,
            },
        ]
    )

    assert ranked["scenario_id"].tolist() == ["a", "b", "narrower", "infeasible"]
    assert ranked["rank"].tolist() == [1, 2, 3, 4]


def test_per_gene_outputs_ignore_sub_tolerance_solver_noise() -> None:
    baseline = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "ko_status": ["optimal", "optimal"],
            "ko_growth": [1.0, 0.5],
            "ko_growth_ratio": [1.0, 0.5],
        }
    )
    first = baseline.copy()
    second = baseline.copy()
    first.loc[0, ["ko_growth", "ko_growth_ratio"]] += 2e-12
    second.loc[0, ["ko_growth", "ko_growth_ratio"]] -= 3e-12

    tables = [
        calibration._build_scenario_per_gene(
            baseline,
            scenario,
            target_genes={"g2"},
            backup_genes={"g1"},
            proxy_genes={"g1"},
            experimental_positive_genes={"g2"},
            cutoffs=(0.01, 0.05, 0.1, 0.15),
        )
        for scenario in (first, second)
    ]

    pd.testing.assert_frame_equal(tables[0], tables[1])
    assert tables[0]["abs_ko_growth_ratio_delta"].max() == 0.0


def test_checkpoint_resume_gates_and_candidate_are_provenance_locked(
    toy_inputs: dict[str, Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        calibration,
        "_load_assay_fitness",
        lambda path: load_assay_fitness(path, expected_rows=None),
    )
    model_before = sha256_file(toy_inputs["model_path"])
    profile_before = sha256_file(toy_inputs["profile_path"])
    original_evaluate = calibration._evaluate_scenario
    calls = 0

    def interrupt_after_one(**kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("simulated interruption")
        return original_evaluate(**kwargs)

    monkeypatch.setattr(calibration, "_evaluate_scenario", interrupt_after_one)
    with pytest.raises(RuntimeError, match="simulated interruption"):
        calibration.calibrate_isozyme_capacities(solver="glpk", **toy_inputs)

    first_id = "R4_0p025__R1846_0p010"
    first_dir = toy_inputs["output_dir"] / "scenarios" / first_id
    assert (first_dir / calibration.SCENARIO_TABLE_NAME).is_file()
    assert (first_dir / calibration.SCENARIO_SUMMARY_NAME).is_file()
    assert not (toy_inputs["output_dir"] / calibration.GRID_SUMMARY_NAME).exists()
    assert not (toy_inputs["output_dir"] / calibration.CANDIDATE_PROFILE_NAME).exists()
    running = json.loads(
        (toy_inputs["output_dir"] / calibration.MANIFEST_NAME).read_text()
    )
    assert running["status"] == "running"

    resumed_calls = 0

    def count_evaluations(**kwargs):
        nonlocal resumed_calls
        resumed_calls += 1
        return original_evaluate(**kwargs)

    monkeypatch.setattr(calibration, "_evaluate_scenario", count_evaluations)
    manifest = calibration.calibrate_isozyme_capacities(
        solver="glpk", resume=True, **toy_inputs
    )

    assert resumed_calls == 8
    assert manifest["status"] == "complete"
    assert manifest["scenario_count"] == 9
    assert manifest["feasible_scenario_count"] == 1
    assert manifest["selected_scenario_id"] == "R4_0p075__R1846_0p050"
    assert manifest["configuration"]["workers"] == 1
    assert sha256_file(toy_inputs["model_path"]) == model_before
    assert sha256_file(toy_inputs["profile_path"]) == profile_before

    output_dir = toy_inputs["output_dir"]
    for filename in (
        calibration.GRID_SUMMARY_NAME,
        calibration.COLLATERAL_NAME,
        calibration.ASSAY_METRICS_NAME,
        calibration.MANIFEST_NAME,
        calibration.CANDIDATE_PROFILE_NAME,
    ):
        assert (output_dir / filename).is_file()
    for scenario_id in calibration.load_joint_capacity_grid(GRID_PATH)["scenario_id"]:
        scenario_dir = output_dir / "scenarios" / scenario_id
        assert (scenario_dir / calibration.SCENARIO_TABLE_NAME).is_file()
        summary_path = scenario_dir / calibration.SCENARIO_SUMMARY_NAME
        assert summary_path.is_file()
        summary = json.loads(summary_path.read_text())
        assert summary["workers"] == 1
        assert summary["gates"]["wt_growth_gate_pass"] is True
        assert summary["gates"]["r1843_bounds_gate_pass"] is True
        assert summary["gates"]["reaction_family_total_upper_bound_gate_pass"] is True
        assert summary["gates"]["no_non_target_call_flips_gate_pass"] is True
        assert summary["gates"]["backup_gene_ko_ratio_gate_pass"] is True
        assert summary["gates"]["concordant_nonessential_proxy_gate_pass"] is True

    candidate = load_provisional_capacity_table(
        output_dir / calibration.CANDIDATE_PROFILE_NAME
    )
    candidate_bounds = {
        row["source_reaction_id"]: row["provisional_upper_bound"] for row in candidate
    }
    assert candidate_bounds["R4"] == pytest.approx(0.075)
    assert candidate_bounds["R1846"] == pytest.approx(0.0495)

    # Every complete checkpoint is reused; no scientific scenario is rerun.
    monkeypatch.setattr(
        calibration,
        "_evaluate_scenario",
        lambda **kwargs: (_ for _ in ()).throw(AssertionError("unexpected rerun")),
    )
    calibration.calibrate_isozyme_capacities(solver="glpk", resume=True, **toy_inputs)

    # Changing any one input invalidates resume before a checkpoint is trusted.
    with toy_inputs["media_path"].open("a", encoding="utf-8") as handle:
        handle.write("\n")
    with pytest.raises(RuntimeError, match="configuration or input SHA differs"):
        calibration.calibrate_isozyme_capacities(
            solver="glpk", resume=True, **toy_inputs
        )


def test_workers_other_than_one_are_rejected(toy_inputs: dict[str, Path]) -> None:
    with pytest.raises(ValueError, match="workers=1"):
        calibration.calibrate_isozyme_capacities(solver="glpk", workers=2, **toy_inputs)
