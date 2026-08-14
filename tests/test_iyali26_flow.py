from __future__ import annotations

import csv
import json
from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model, write_sbml_model

import iyali26_flow.cli as flow_cli
import iyali26_flow.experiment as experiment_module
import iyali26_flow.simulator as simulator_module
from iyali26_flow import (
    ExperimentRunner,
    ParameterPoint,
    R4R1846CapacitySimulator,
    SimulationResult,
    load_experiment_config,
)
from iyali26_flow.analysis import (
    AffineDiagnostics,
    make_phase1_decision,
    sobol_parameter_points,
)
from iyali26_flow.experiment import ExperimentOutcome
from iyali26_flow.fmpe import (
    FMPE_EXTRA_ERROR,
    _decision_allows_training,
    bounded_logit,
    bounded_sigmoid,
    require_fmpe_dependencies,
)
from iyali26_flow.provenance import sha256_file
from scripts.gem_annotate.essentiality_evidence import target_fingerprint
from scripts.gem_annotate.validate_essential_genes import load_assay_fitness


REPO_ROOT = Path(__file__).resolve().parents[1]
REAL_CONFIG = REPO_ROOT / "experiments" / "flow_matching" / "r4_r1846_phase1.json"


def _anchor_rows() -> list[dict[str, object]]:
    rows = []
    for r4 in (0.025, 0.075, 0.15):
        for r1846 in (0.01, 0.025, 0.05):
            scenario_id = (
                f"R4_{r4:.3f}__R1846_{r1846:.3f}".replace(".", "p")
            )
            rows.append(
                {
                    "scenario_id": scenario_id,
                    "r4_capacity_fraction_of_wt_flux": r4,
                    "r1846_capacity_fraction_of_wt_flux": r1846,
                    "r4_applied_backup_upper_bound": r4,
                    "r1846_applied_backup_upper_bound": r1846,
                    "r4_target_ko_ratio": r4,
                    "r1846_target_ko_ratio": r1846,
                    "feasible": 0.05 <= r4 < 0.1 and 0.05 <= r1846 < 0.1,
                    "target_window_pass": (
                        0.05 <= r4 < 0.1 and 0.05 <= r1846 < 0.1
                    ),
                    "wt_growth_delta": 0.0,
                    "wt_growth_gate_pass": True,
                    "r1843_bounds_gate_pass": True,
                    "reaction_family_total_upper_bound_gate_pass": True,
                    "non_target_call_flip_count": 0,
                    "no_non_target_call_flips_gate_pass": True,
                    "backup_gene_ko_ratio_gate_pass": True,
                    "proxy_new_essential_count": 0,
                    "concordant_nonessential_proxy_gate_pass": True,
                }
            )
    return rows


def _write_anchor_inputs(root: Path) -> tuple[Path, Path]:
    rows = _anchor_rows()
    grid_path = root / "anchor_grid.csv"
    with grid_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scenario_id",
                "r4_capacity_fraction_of_wt_flux",
                "r1846_capacity_fraction_of_wt_flux",
            ],
        )
        writer.writeheader()
        writer.writerows(
            {
                key: row[key]
                for key in writer.fieldnames
            }
            for row in rows
        )
    reference_path = root / "anchor_reference.tsv"
    with reference_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(rows[0]),
        )
        writer.writeheader()
        writer.writerows(rows)
    return grid_path, reference_path


def _write_config(
    root: Path,
    *,
    inputs: dict[str, Path] | None = None,
    solver: str = "glpk",
    checkpoint_every: int = 2,
) -> Path:
    if inputs is None:
        inputs = {}
        for name in ("model", "profile", "media", "experimental", "assay_fitness"):
            path = root / f"{name}.txt"
            path.write_text(name, encoding="utf-8")
            inputs[name] = path
    anchor_grid, anchor_reference = _write_anchor_inputs(root)
    inputs = {
        **inputs,
        "anchor_grid": anchor_grid,
        "anchor_reference": anchor_reference,
    }
    raw = {
        "schema_version": 1,
        "experiment_id": "toy_flow",
        "seed": 20260720,
        "solver": solver,
        "inputs": {
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in inputs.items()
        },
        "parameter_space": {
            "r4_capacity_fraction": [0.025, 0.15],
            "r1846_capacity_fraction": [0.01, 0.05],
        },
        "design": {
            "sobol_samples": 4,
            "holdout_samples": 4,
            "reference_grid_size": 3,
            "synthetic_truths": 2,
            "full_feasible_cap": 2,
            "full_nonfeasible_count": 2,
            "checkpoint_every": checkpoint_every,
        },
        "observation_noise_sigma": [0.005, 0.02],
        "growth_cutoffs": [0.01, 0.05, 0.1, 0.15],
        "thresholds": {
            "anchor_tolerance": 1e-9,
            "affine_max_normalized_error": 0.01,
            "multimodal_fraction": 0.1,
            "component_minimum_mass": 0.01,
            "coverage_lower": 0.8,
            "coverage_upper": 0.98,
        },
        "soft_time_budget_seconds": 60,
        "hard_timeout_seconds": 120,
        "fmpe": {
            "simulations": 8,
            "hidden_features": [64, 64, 64, 64, 64],
            "epochs": 2,
            "patience": 1,
            "batch_size": 4,
            "learning_rate": 0.001,
            "gradient_clip": 1.0,
            "seeds": [1, 2, 3],
        },
    }
    path = root / "config.json"
    path.write_text(json.dumps(raw, indent=2), encoding="utf-8")
    return path


def _result(
    point: ParameterPoint,
    sample_id: str,
    fidelity: str = "targeted",
) -> SimulationResult:
    feasible = 0.05 <= point.r4_fraction < 0.10 and 0.05 <= point.r1846_fraction < 0.10
    return SimulationResult(
        sample_id=sample_id,
        fidelity=fidelity,
        point=point,
        r4_applied_upper_bound=point.r4_fraction,
        r1846_applied_upper_bound=point.r1846_fraction,
        r4_target_ko_ratio=point.r4_fraction,
        r1846_target_ko_ratio=point.r1846_fraction,
        backup_ko_ratios=(("backup_a", 1.0), ("backup_b", 1.0), ("backup_c", 1.0)),
        baseline_wt_growth=1.0,
        scenario_wt_growth=1.0,
        wt_growth_delta=0.0,
        target_window_pass=feasible,
        backup_gate_pass=True,
        r1843_gate_pass=True,
        family_total_gate_pass=True,
        targeted_feasible=feasible,
        full_feasible=feasible if fidelity == "full" else None,
        non_target_call_flip_count=0 if fidelity == "full" else None,
        proxy_new_essential_count=0 if fidelity == "full" else None,
        no_non_target_call_flips_gate_pass=True if fidelity == "full" else None,
        proxy_gate_pass=True if fidelity == "full" else None,
        gene_count=8 if fidelity == "full" else 5,
        elapsed_seconds=0.001,
    )


class _FakeSimulator:
    calls = 0
    fail_after: int | None = None
    callback = None

    def __init__(self, config) -> None:
        self.config = config
        self.model_path = config.input_map["model"].path
        self.model_sha256 = sha256_file(self.model_path)
        self.target_fingerprints = {"R4": "toy-r4", "R1846": "toy-r1846"}
        self.solver_metadata = {"requested": config.solver, "version": "toy"}
        self.baseline_invariants = {"model_sha256": self.model_sha256}
        self.backup_genes = ("backup_a", "backup_b", "backup_c")

    def simulate(self, point, *, sample_id, fidelity="targeted"):
        if type(self).fail_after is not None and type(self).calls >= type(self).fail_after:
            raise RuntimeError("simulated interruption")
        type(self).calls += 1
        if type(self).callback is not None:
            type(self).callback()
        return _result(point, sample_id, fidelity)

    def assert_input_model_unchanged(self) -> None:
        assert sha256_file(self.model_path) == self.model_sha256


def test_parameter_bounds_logit_roundtrip_and_immutability() -> None:
    point = ParameterPoint(0.075, 0.025)
    with pytest.raises(FrozenInstanceError):
        point.r4_fraction = 0.1
    with pytest.raises(ValueError, match="positive"):
        ParameterPoint(0.0, 0.025)

    values = np.asarray([[0.026, 0.011], [0.075, 0.025], [0.149, 0.049]])
    lower = np.asarray([0.025, 0.01])
    upper = np.asarray([0.15, 0.05])
    recovered = bounded_sigmoid(bounded_logit(values, lower, upper), lower, upper)
    np.testing.assert_allclose(recovered, values, atol=1e-14, rtol=0)
    with pytest.raises(ValueError, match="strictly inside"):
        bounded_logit(lower, lower, upper)


def test_sobol_reproducibility_and_configuration_validation(tmp_path: Path) -> None:
    config_path = _write_config(tmp_path)
    config = load_experiment_config(config_path, repo_root=tmp_path)
    first = sobol_parameter_points(config, count=4, seed=7)
    second = sobol_parameter_points(config, count=4, seed=7)
    different = sobol_parameter_points(config, count=4, seed=8)

    assert first == second
    assert first != different
    assert all(config.r4_bounds.lower <= point.r4_fraction <= config.r4_bounds.upper for point in first)

    raw = json.loads(config_path.read_text(encoding="utf-8"))
    raw["unknown"] = True
    config_path.write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown"):
        load_experiment_config(config_path, repo_root=tmp_path)


def test_input_sha_mismatch_is_rejected(tmp_path: Path) -> None:
    config_path = _write_config(tmp_path)
    config = load_experiment_config(config_path, repo_root=tmp_path)
    config.input_map["media"].path.write_text("changed", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA is stale"):
        config.verify_inputs()


def test_phase1_decision_gates() -> None:
    config = load_experiment_config(REAL_CONFIG, repo_root=REPO_ROOT)
    affine = AffineDiagnostics(((0.0, 0.0),), 0.001, 0.0001, (1.0, 1.0))
    posterior = {
        "coverage_by_sigma": {0.005: 0.9, 0.02: 0.9},
        "multimodal_fraction": 0.0,
    }
    simple = make_phase1_decision(
        config,
        invalid_baseline_reasons=[],
        affine=affine,
        posterior=posterior,
    )
    assert simple.status == "stop_simple_baseline_sufficient"

    complex_posterior = {**posterior, "multimodal_fraction": 0.1}
    proceed = make_phase1_decision(
        config,
        invalid_baseline_reasons=[],
        affine=affine,
        posterior=complex_posterior,
    )
    assert proceed.status == "proceed_fmpe"

    blocked = make_phase1_decision(
        config,
        invalid_baseline_reasons=[],
        affine=affine,
        posterior={**posterior, "coverage_by_sigma": {0.005: 0.7, 0.02: 0.9}},
    )
    assert blocked.status == "blocked_invalid_baseline"


def test_checkpoint_resume_and_partial_timeout_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = load_experiment_config(_write_config(tmp_path), repo_root=tmp_path)
    monkeypatch.setattr(
        experiment_module,
        "R4R1846CapacitySimulator",
        _FakeSimulator,
    )
    points = [
        (f"point_{index}", ParameterPoint(0.03 + index * 0.01, 0.02))
        for index in range(5)
    ]
    _FakeSimulator.calls = 0
    _FakeSimulator.fail_after = 2
    _FakeSimulator.callback = None
    output = tmp_path / "resume_results"
    runner = ExperimentRunner(config, output)
    with pytest.raises(RuntimeError, match="simulated interruption"):
        runner._run_stage("resume_test", points, fidelity="targeted")
    assert len(_read_rows(output / "resume_test.tsv")) == 2

    _FakeSimulator.calls = 0
    _FakeSimulator.fail_after = None
    resumed = ExperimentRunner(config, output)
    results = resumed._run_stage("resume_test", points, fidelity="targeted")
    assert len(results) == 5
    assert _FakeSimulator.calls == 3

    timeout_output = tmp_path / "timeout_results"
    _FakeSimulator.calls = 0
    _FakeSimulator.callback = None
    timeout_runner = ExperimentRunner(config, timeout_output)
    _FakeSimulator.callback = lambda: setattr(timeout_runner, "_soft_deadline", 0.0)
    outcome = timeout_runner.run_phase1()
    assert outcome.status == "partial_timeout"
    manifest = json.loads((timeout_output / "run_manifest.json").read_text())
    assert manifest["status"] == "partial_timeout"
    assert len(_read_rows(timeout_output / "anchors.tsv")) == 1
    checkpoint = json.loads(
        (timeout_output / "checkpoints" / "anchors.json").read_text()
    )
    assert checkpoint["status"] == "running"
    _FakeSimulator.callback = None


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_missing_fmpe_extra_and_decision_gate(tmp_path: Path) -> None:
    def missing(_name: str):
        raise ModuleNotFoundError

    with pytest.raises(RuntimeError, match="uv sync --extra fmpe") as error:
        require_fmpe_dependencies(missing)
    assert str(error.value) == FMPE_EXTRA_ERROR

    (tmp_path / "phase1_decision.json").write_text(
        json.dumps({"status": "stop_simple_baseline_sufficient"}),
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="gated off"):
        _decision_allows_training(tmp_path)


def test_cli_phase1_analyze_and_conditional_train(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = load_experiment_config(_write_config(tmp_path), repo_root=tmp_path)
    output = tmp_path / "cli_results"

    class FakeRunner:
        def __init__(
            self,
            _config,
            output_dir,
            *,
            resume=True,
            research_root=None,
            force_rerun=False,
            reproduction_reason=None,
        ):
            self.output_dir = Path(output_dir)
            self.resume = resume

        def run_phase1(self):
            return ExperimentOutcome(
                "partial_timeout",
                self.output_dir,
                self.output_dir / "run_manifest.json",
                None,
            )

        def analyze_only(self):
            return ExperimentOutcome(
                "complete",
                self.output_dir,
                self.output_dir / "run_manifest.json",
                None,
            )

    monkeypatch.setattr(flow_cli, "load_experiment_config", lambda *args, **kwargs: config)
    monkeypatch.setattr(flow_cli, "ExperimentRunner", FakeRunner)
    monkeypatch.setattr(flow_cli, "_install_alarm", lambda _seconds: None)
    assert flow_cli.main(["phase1", "--config", "x", "--output", str(output)]) == 124
    assert flow_cli.main(["analyze", "--config", "x", "--output", str(output)]) == 0

    monkeypatch.setattr(
        flow_cli,
        "train_fmpe",
        lambda _config, _output: {"acceptance_pass": True},
    )
    assert flow_cli.main(["train-fmpe", "--config", "x", "--output", str(output)]) == 0


def _toy_model() -> Model:
    model = Model("toy_iyali26_flow")
    a = Metabolite("a_c", compartment="c")
    x = Metabolite("x_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    c = Metabolite("c_c", compartment="c")

    ex_a = Reaction("EX_a")
    ex_a.bounds = (-10.0, 1000.0)
    ex_a.add_metabolites({a: -1.0})
    ex_x = Reaction("EX_x")
    ex_x.bounds = (-10.0, 1000.0)
    ex_x.add_metabolites({x: -1.0})

    r4 = Reaction("R4")
    r4.bounds = (0.0, 1000.0)
    r4.add_metabolites({a: -1.0, b: 1.0})
    r4.gene_reaction_rule = "g_r4_backup or g_r4_main"
    r1846 = Reaction("R1846")
    r1846.bounds = (0.0, 1000.0)
    r1846.add_metabolites({x: -1.0, c: 1.0})
    r1846.gene_reaction_rule = (
        "g_r1846_backup_a or g_r1846_backup_b or g_r1846_main"
    )
    r1843 = Reaction("R1843")
    r1843.bounds = (0.0, 0.0)
    r1843.add_metabolites({b: -1.0, c: 1.0})
    r1843.gene_reaction_rule = "g_locked"
    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1000.0)
    biomass.add_metabolites({b: -1.0, c: -1.0})
    model.add_reactions([ex_a, ex_x, r4, r1846, r1843, biomass])
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


def _write_toy_profile(path: Path, model_path: Path) -> None:
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
    rows = []
    for reaction_id, expected, primary, backup in (
        ("R4", "g_r4_backup or g_r4_main", "g_r4_main", "g_r4_backup"),
        (
            "R1846",
            "g_r1846_backup_a or g_r1846_backup_b or g_r1846_main",
            "g_r1846_main",
            "g_r1846_backup_a or g_r1846_backup_b",
        ),
    ):
        rows.append(
            {
                "capacity_id": f"PCAP-{reaction_id}-TOY",
                "status": "active_exploratory",
                "source_reaction_id": reaction_id,
                "expected_gpr": expected,
                "primary_gpr": primary,
                "backup_gpr": backup,
                "backup_reaction_id": f"{reaction_id}__PCAP_BACKUP",
                "provisional_upper_bound": "0.1",
                "units": "mmol_gDW_h",
                "parameter_basis": "toy_not_measured",
                "case_id": f"EGC-{reaction_id}-toy",
                "validated_model_sha256": sha256_file(model_path),
                "target_fingerprint": _reaction_fingerprint(model, reaction_id),
                "requires_protein_abundance": "true",
                "requires_kcat": "true",
                "replacement_formula": "sum_i(kcat_i*E_i*3600)",
                "rationale": "test only",
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_toy_assay(path: Path, gene_ids: list[str]) -> None:
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
    rows = []
    target_genes = {"g_r4_main", "g_r1846_main"}
    for assay in ("Cas9", "Cas12a"):
        for row_number, gene_id in enumerate(gene_ids, start=2):
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
                    "source_row": row_number,
                    "source_sha256": "a" * 64,
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _toy_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    model_path = tmp_path / "model.xml"
    write_sbml_model(_toy_model(), str(model_path))
    profile_path = tmp_path / "profile.csv"
    _write_toy_profile(profile_path, model_path)
    loaded = read_sbml_model(str(model_path))
    assay_path = tmp_path / "assay.csv"
    _write_toy_assay(assay_path, sorted(gene.id for gene in loaded.genes))
    experimental_path = tmp_path / "experimental.csv"
    experimental_path.write_text("gene_id\ng_r4_main\ng_r1846_main\n", encoding="utf-8")
    media_path = tmp_path / "media.csv"
    media_path.write_text("exchange,uptake\nEX_a,1\nEX_x,1\n", encoding="utf-8")
    monkeypatch.setattr(
        simulator_module,
        "_load_assay_fitness",
        lambda path: load_assay_fitness(path, expected_rows=None),
    )
    config_path = _write_config(
        tmp_path,
        inputs={
            "model": model_path,
            "profile": profile_path,
            "media": media_path,
            "experimental": experimental_path,
            "assay_fitness": assay_path,
        },
        solver="glpk",
    )
    return load_experiment_config(config_path, repo_root=tmp_path), model_path


def test_toy_targeted_full_fidelity_and_no_input_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config, model_path = _toy_config(tmp_path, monkeypatch)
    before = sha256_file(model_path)
    simulator = R4R1846CapacitySimulator(config)
    point = ParameterPoint(0.075, 0.025)
    targeted = simulator.simulate(point, sample_id="toy_targeted", fidelity="targeted")
    full = simulator.simulate(point, sample_id="toy_full", fidelity="full")

    assert targeted.r4_target_ko_ratio == pytest.approx(full.r4_target_ko_ratio, abs=1e-9)
    assert targeted.r1846_target_ko_ratio == pytest.approx(
        full.r1846_target_ko_ratio, abs=1e-9
    )
    assert dict(targeted.backup_ko_ratios) == pytest.approx(dict(full.backup_ko_ratios), abs=1e-9)
    assert targeted.scenario_wt_growth == pytest.approx(
        simulator.baseline_wt_growth, abs=1e-9
    )
    assert targeted.family_total_gate_pass
    assert full.no_non_target_call_flips_gate_pass
    assert full.proxy_gate_pass
    assert sha256_file(model_path) == before


@pytest.mark.integration
def test_real_model_rejects_any_stale_configured_input() -> None:
    config = load_experiment_config(REAL_CONFIG, repo_root=REPO_ROOT)
    model_path = config.input_map["model"].path
    before = sha256_file(model_path)
    configured_sha = config.input_map["model"].sha256

    assert configured_sha == (
        "39f4cae11c3f270400c8a227c78b6af3ed412e85b1ade6cb604b0f85c3d8b1d9"
    )
    assert before == (
        "3b0369f25e9d3727642507e35684f3cf036bdc9fcedf290a921121e956da71bf"
    )
    with pytest.raises(ValueError, match=r"configured (?:model|media) SHA is stale"):
        R4R1846CapacitySimulator(config)
    assert sha256_file(model_path) == before
