import json
import math
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model
from cobra.manipulation.delete import knock_out_model_genes

from scripts.gem_annotate import quinone_dfba_essentiality as dfba
from scripts.gem_annotate.summarize_quinone_dfba import (
    FIVE_COQ_CONTROLS,
    NINE_COQ_CANDIDATES,
    summarize_calls,
    write_hypothesis_matrix,
)
from scripts.gem_annotate.coq9_dilution import (
    add_runtime_coq9_pool_source,
    apply_runtime_coq9_dilution,
    runtime_coq9_dilution_record,
)
from scripts.gem_annotate.quinone_dfba_essentiality import _optimize_minimal_pool


def _toy_model():
    model = Model("q9_reserve")
    nutrient, uracil, q9 = (
        Metabolite(name, compartment="e") for name in ("a", "uracil", "m468[C_mi]")
    )
    uptake = Reaction("R1070", lower_bound=-2000, upper_bound=0)
    uptake.add_metabolites({nutrient: -1})
    uracil_uptake = Reaction("R1354", lower_bound=-1000, upper_bound=0)
    uracil_uptake.add_metabolites({uracil: -1})
    synthase = Reaction("Q9_SYNTHASE", lower_bound=0, upper_bound=1000)
    synthase.add_metabolites({nutrient: -1, q9: 1})
    synthase.gene_reaction_rule = "qgene"
    growth = Reaction("biomass_C", lower_bound=0, upper_bound=1000)
    growth.add_metabolites({nutrient: -1, uracil: -1})
    model.add_reactions([uptake, uracil_uptake, synthase, growth])
    return model


def _toy_model_with_maintenance():
    model = _toy_model()
    maintenance = Reaction("xMAINTENANCE", lower_bound=1, upper_bound=1000)
    maintenance.add_metabolites({model.metabolites.get_by_id("a"): -1})
    model.add_reactions([maintenance])
    return model


def test_runtime_q9_reserve_is_used_only_after_q9_synthesis_knockout():
    model = _toy_model()
    with model:
        source = add_runtime_coq9_pool_source(model)
        apply_runtime_coq9_dilution(model, 0.1)
        source.upper_bound = 1.0
        _, solution = _optimize_minimal_pool(model, source)
        assert solution.fluxes[source.id] == 0
        knock_out_model_genes(model, ["qgene"])
        _, solution = _optimize_minimal_pool(model, source)
        assert solution.fluxes[source.id] > 0


def test_seeded_alpha_samples_are_fixed_before_gene_simulations():
    parser = dfba._parser()
    args = parser.parse_args([
        "--research-root", "research", "--alpha-seed", "20260826", "--alpha-replicates", "3",
    ])
    samples = dfba._alpha_samples(args)

    assert [item["replicate_id"] for item in samples] == [0, 1, 2]
    assert samples == dfba._alpha_samples(args)
    assert all(1e-6 <= item["alpha_mmol_gDW"] <= 1e-3 for item in samples)

    incompatible = parser.parse_args([
        "--research-root", "research", "--alpha-seed", "20260826", "--alphas", "1e-4",
    ])
    with pytest.raises(ValueError, match="cannot be combined"):
        dfba._alpha_samples(incompatible)


def test_pool_multiplier_rejects_negative_or_nonfinite_values():
    for value in (-1, math.nan):
        with pytest.raises(ValueError, match="pool multipliers"):
            dfba.simulate_gene(
                _toy_model(), gene_id=None, alpha=0, pool_multiplier=value,
                hours=0.5, dt=0.5, initial_biomass=1,
            )
    with pytest.raises(ValueError, match="empty gene IDs"):
        dfba._gene_ids(_toy_model(), "", (), 0, 1)
    with pytest.raises(ValueError, match="must be unique"):
        dfba._gene_ids(_toy_model(), "qgene,qgene", (), 0, 1)


def test_merge_uses_only_manifest_bound_calls_and_checks_their_hashes(tmp_path):
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    for index, gene_id in enumerate(("g1", "g2")):
        calls_path = output_dir / f"chunk_{index:03d}_calls.tsv"
        trajectory_path = output_dir / f"chunk_{index:03d}_trajectory.tsv"
        pd.DataFrame([{"gene_id": gene_id}]).to_csv(calls_path, sep="\t", index=False)
        pd.DataFrame([{"gene_id": gene_id, "time_h": 0.0}]).to_csv(
            trajectory_path, sep="\t", index=False,
        )
        manifest = {
            "workflow": dfba.WORKFLOW,
            "schema_version": dfba.SCHEMA_VERSION,
            "run_id": "merge-test",
            "chunk_index": index,
            "chunk_count": 2,
            "genes": [gene_id],
            "alpha_sampling": [{"alpha_mmol_gDW": 1e-4}],
            "coq9_dilution": [runtime_coq9_dilution_record(1e-4, pool_source_enabled=True)],
            "coq9_dilution_tool_sha256": "tool",
            "q9_reserve_policy": dfba.Q9_RESERVE_POLICY,
            "runtime_topology": {"enabled": False},
            "output_files": {
                "calls": {"filename": calls_path.name, "sha256": dfba.sha256_file(calls_path)},
                "trajectory": {"filename": trajectory_path.name, "sha256": dfba.sha256_file(trajectory_path)},
            },
            "created_at_utc": f"2026-08-26T00:00:0{index}Z",
        }
        manifest["fingerprint"] = dfba.sha256_payload(manifest)
        (output_dir / f"chunk_{index:03d}_manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8",
        )

    pd.DataFrame([{"gene_id": "stale"}]).to_csv(
        output_dir / "chunk_999_calls.tsv", sep="\t", index=False,
    )
    dfba._merge(SimpleNamespace(output_dir=str(output_dir)))
    assert set(pd.read_csv(output_dir / "essentiality_dynamic_calls.tsv", sep="\t")["gene_id"]) == {"g1", "g2"}

    manifest_path = output_dir / "chunk_001_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    calls_sha = manifest["output_files"]["calls"]["sha256"]
    manifest["output_files"]["calls"]["sha256"] = "bad"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="manifest fingerprint mismatch"):
        dfba._merge(SimpleNamespace(output_dir=str(output_dir)))
    manifest["output_files"]["calls"]["sha256"] = calls_sha
    manifest["fingerprint"] = dfba.sha256_payload({
        key: value for key, value in manifest.items() if key != "fingerprint"
    })
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    pd.DataFrame([{"gene_id": "tampered"}]).to_csv(
        output_dir / "chunk_001_calls.tsv", sep="\t", index=False,
    )
    with pytest.raises(ValueError, match="output hash mismatch"):
        dfba._merge(SimpleNamespace(output_dir=str(output_dir)))


def test_po1f_nonlimiting_uracil_keeps_base_bound_and_has_no_finite_pool(monkeypatch):
    monkeypatch.setattr(dfba, "INITIAL_POOLS_MMOL_L", {"R1070": 1e9, "R1354": 0.1})
    finite, _ = dfba.simulate_gene(
        _toy_model(), gene_id=None, alpha=0, pool_multiplier=0, hours=0.5, dt=0.5,
        initial_biomass=1, uracil_mode="finite_batch",
    )
    nonlimiting, trace = dfba.simulate_gene(
        _toy_model(), gene_id=None, alpha=0, pool_multiplier=0, hours=0.5, dt=0.5,
        initial_biomass=1, uracil_mode="po1f_nonlimiting",
    )

    assert nonlimiting["dynamic_doublings"] > finite["dynamic_doublings"]
    assert trace[0]["growth_h-1"] == 1000
    assert math.isnan(nonlimiting["final_uracil_mmol_L"])
    assert all(row["uracil_mode"] == "po1f_nonlimiting" for row in trace)
    assert all(math.isnan(row["uracil_mmol_L"]) for row in trace)


def test_nonoptimal_step_is_recorded_once_without_reading_invalid_fluxes(monkeypatch):
    monkeypatch.setattr(dfba, "INITIAL_POOLS_MMOL_L", {"R1070": 0.5, "R1354": 1e9})
    result, trace = dfba.simulate_gene(
        _toy_model_with_maintenance(), gene_id=None, alpha=0, pool_multiplier=0,
        hours=1, dt=0.5, initial_biomass=1,
    )

    assert [row["status"] for row in trace] == ["optimal", "infeasible"]
    assert math.isnan(trace[-1]["growth_h-1"])
    assert math.isnan(trace[-1]["q9_source_flux_mmol_gDW_h"])
    assert result["termination_status"] == "infeasible"
    assert result["termination_time_h"] == 0.5


def test_gurobi_feasibility_tolerance_is_bounded_and_explicit():
    parser = dfba._parser()
    args = parser.parse_args([
        "--research-root", "research", "--solver", "gurobi", "--feasibility-tol", "1e-9",
    ])

    assert args.feasibility_tol == 1e-9
    with pytest.raises(SystemExit):
        parser.parse_args(["--research-root", "research", "--feasibility-tol", "1e-10"])
    with pytest.raises(ValueError, match="requires --solver gurobi"):
        dfba._configure_solver(_toy_model(), "glpk", 1e-9)


def test_r39_r19_runtime_topology_is_balanced_and_leaves_canonical_model_unchanged():
    model = read_sbml_model(Path(__file__).parents[1] / "model.xml")
    reaction_ids = ("R39", "R969", "R808", "R19")
    before = {
        reaction_id: (
            dfba._reaction_stoichiometry(model.reactions.get_by_id(reaction_id)),
            model.reactions.get_by_id(reaction_id).bounds,
            model.reactions.get_by_id(reaction_id).gene_reaction_rule,
        )
        for reaction_id in reaction_ids
    }
    counts = (len(model.reactions), len(model.metabolites), len(model.genes))

    with model:
        audit = dfba._apply_r39_r19_runtime_topology(model)
        assert audit["mapping_sha256"] == dfba.R39_R19_RUNTIME_MAPPING_SHA256
        assert audit["mass_balance"] == {
            "R39": {}, "DFBA_R19_HYDROXYLATION": {}, "DFBA_R19_FORMAL_OXIDATION": {},
        }
        assert model.reactions.R969.bounds == (0.0, 0.0)
        assert model.reactions.R808.bounds == (0.0, 0.0)
        assert model.reactions.R19.bounds == (0.0, 0.0)
        assert model.reactions.get_by_id("DFBA_R19_HYDROXYLATION").gene_reaction_rule == ""
        assert model.reactions.get_by_id("DFBA_R19_FORMAL_OXIDATION").gene_reaction_rule == ""

    after = {
        reaction_id: (
            dfba._reaction_stoichiometry(model.reactions.get_by_id(reaction_id)),
            model.reactions.get_by_id(reaction_id).bounds,
            model.reactions.get_by_id(reaction_id).gene_reaction_rule,
        )
        for reaction_id in reaction_ids
    }
    assert counts == (len(model.reactions), len(model.metabolites), len(model.genes))
    assert before == after


def test_preregistered_runtime_gpr_scenarios_are_atomic_and_restore_the_model():
    model = read_sbml_model(Path(__file__).parents[1] / "model.xml")
    orphan_candidates = ("YALI1A08781g", "YALI1B20527g", "YALI1F34675g")

    def signature():
        return {
            reaction.id: (
                dfba._reaction_stoichiometry(reaction),
                reaction.bounds,
                reaction.gene_reaction_rule,
            )
            for reaction in model.reactions
        }

    routewide = {"R39", "R715", "R40", "DFBA_R19_HYDROXYLATION", "R18", "R695", "R385"}
    expected = {
        "baseline_topology": ({}, {gene_id: set() for gene_id in orphan_candidates}),
        "coq6_r39": (
            {"R39": "YALI1A08781g"}, {"YALI1A08781g": {"R39"}},
        ),
        "coq6_r19_hydroxylation": (
            {"DFBA_R19_HYDROXYLATION": "YALI1A08781g"},
            {"YALI1A08781g": {"DFBA_R19_HYDROXYLATION"}},
        ),
        "coq6_both_hydroxylations": (
            {"R39": "YALI1A08781g", "DFBA_R19_HYDROXYLATION": "YALI1A08781g"},
            {"YALI1A08781g": {"R39", "DFBA_R19_HYDROXYLATION"}},
        ),
        "coq8_routewide_absolute": (
            {
                "R39": "YALI1B20527g",
                "R715": "YALI1B20527g and YALI1B20835g",
                "R40": "YALI1B20527g and YALI1F34625g",
                "DFBA_R19_HYDROXYLATION": "YALI1B20527g",
                "R18": "YALI1B20527g and YALI1C25352g",
                "R695": "YALI1B20527g and YALI1E18269g",
                "R385": "YALI1B20527g and YALI1B20835g",
            },
            {"YALI1B20527g": routewide},
        ),
        "coq9_r695_absolute": (
            {"R695": "YALI1F34675g and YALI1E18269g"},
            {"YALI1F34675g": {"R695"}},
        ),
        "coq9_routewide_absolute": (
            {
                "R39": "YALI1F34675g",
                "R715": "YALI1F34675g and YALI1B20835g",
                "R40": "YALI1F34675g and YALI1F34625g",
                "DFBA_R19_HYDROXYLATION": "YALI1F34675g",
                "R18": "YALI1F34675g and YALI1C25352g",
                "R695": "YALI1F34675g and YALI1E18269g",
                "R385": "YALI1F34675g and YALI1B20835g",
            },
            {"YALI1F34675g": routewide},
        ),
    }
    base_ids = (
        tuple(reaction.id for reaction in model.reactions),
        tuple(metabolite.id for metabolite in model.metabolites),
        tuple(gene.id for gene in model.genes),
    )
    base_signature = signature()

    for scenario in expected:
        with pytest.raises(ValueError, match="require the R39/R19 runtime topology"):
            dfba._apply_runtime_gpr_scenario(model, scenario)

    for scenario, (updates, candidate_kos) in expected.items():
        with model:
            dfba._apply_r39_r19_runtime_topology(model)
            topology_signature = signature()
            wt_growth = model.slim_optimize()
            audit = dfba._apply_runtime_gpr_scenario(model, scenario)
            mapped_signature = signature()

            assert audit["scenario_id"] == scenario
            assert dfba.RUNTIME_GPR_SCENARIOS[scenario]["gpr_updates"] == updates
            assert audit["mapping_sha256"] == dfba.sha256_payload(
                dfba._runtime_gpr_mapping_payload(scenario)
            )
            assert {
                reaction_id
                for reaction_id in mapped_signature
                if mapped_signature[reaction_id] != topology_signature[reaction_id]
            } == set(updates), scenario
            assert all(
                mapped_signature[reaction_id][:2] == topology_signature[reaction_id][:2]
                for reaction_id in model.reactions.list_attr("id")
            ), scenario
            assert model.slim_optimize() == pytest.approx(wt_growth, abs=1e-9)
            assert model.reactions.R19.gene_reaction_rule == ""
            assert model.reactions.R19.bounds == (0.0, 0.0)
            assert model.reactions.get_by_id("DFBA_R19_FORMAL_OXIDATION").gene_reaction_rule == ""

            for candidate_gene, expected_closed in candidate_kos.items():
                assert {reaction.id for reaction in model.genes.get_by_id(candidate_gene).reactions} == set(updates)
                mapped_bounds = {reaction.id: reaction.bounds for reaction in model.reactions}
                with model:
                    knock_out_model_genes(model, [candidate_gene])
                    assert {
                        reaction.id for reaction in model.reactions
                        if reaction.bounds != mapped_bounds[reaction.id]
                    } == expected_closed, scenario
                    assert all(model.reactions.get_by_id(reaction_id).bounds == (0.0, 0.0) for reaction_id in expected_closed)
                    assert "R19" not in expected_closed
                    assert "DFBA_R19_FORMAL_OXIDATION" not in expected_closed
                assert signature() == mapped_signature

        assert signature() == base_signature, scenario
        assert base_ids == (
            tuple(reaction.id for reaction in model.reactions),
            tuple(metabolite.id for metabolite in model.metabolites),
            tuple(gene.id for gene in model.genes),
        ), scenario
        assert all(not model.genes.get_by_id(gene_id).reactions for gene_id in orphan_candidates)

    with model:
        dfba._apply_r39_r19_runtime_topology(model)
        model.reactions.R715.gene_reaction_rule = "unexpected_gene"
        unexpected_signature = signature()
        with pytest.raises(ValueError, match="baseline mismatch"):
            dfba._apply_runtime_gpr_scenario(model, "coq6_r39")
        assert signature() == unexpected_signature


def test_dfba_calls_summary_recomputes_grid_controls_and_monotonicity():
    rows = []
    for alpha in (1e-6, 1e-4, 1e-3):
        for pool in (0.0, 0.5, 1.0, 2.0):
            doublings = math.log2(1 + pool)
            for control in FIVE_COQ_CONTROLS:
                ratio = doublings / 10
                rows.append({
                    "gene_id": control["gene_id"], "alpha_mmol_gDW": alpha,
                    "pool_multiplier": pool, "dynamic_doublings": doublings,
                    "dynamic_growth_ratio": ratio,
                    "q9_source_total_mmol_L": alpha * 0.01 * pool,
                    "experimental_essential": False,
                    **{f"essential_at_{cutoff * 100:g}pct": ratio < cutoff for cutoff in (0.01, 0.05, 0.1, 0.15)},
                })

    tables, summary = summarize_calls(pd.DataFrame(rows), initial_biomass=0.01)

    assert len(tables["grid_summary"]) == 12
    assert len(tables["five_gene_pool_summary"]) == 60
    assert summary["pool_monotonicity_pass"]
    assert summary["q9_source_calls_bound_pass"]
    assert summary["five_control_max_abs_theory_error_doublings"] == 0


def test_hypothesis_matrix_validates_and_joins_the_seven_runtime_arms(tmp_path):
    arm_dirs = []
    for scenario, spec in dfba.RUNTIME_GPR_SCENARIOS.items():
        arm_dir = tmp_path / scenario
        arm_dir.mkdir()
        arm_dirs.append(arm_dir)
        mapping = dfba._runtime_gpr_mapping_payload(scenario)
        mapping_sha = dfba.sha256_payload(mapping)
        mapping_audit = {
            "enabled": True, **mapping, "mapping_sha256": mapping_sha,
            "runtime_target_fingerprint_before_mapping": "sha256:target",
        }
        rows = []
        for alpha in dfba.DEFAULT_ALPHAS:
            for pool in dfba.DEFAULT_POOL_MULTIPLIERS:
                doublings = math.log2(1 + pool)
                ratio = doublings / 10
                for gene_id in NINE_COQ_CANDIDATES:
                    rows.append({
                        "gene_id": gene_id,
                        "dynamic_doublings": doublings,
                        "dynamic_growth_ratio": ratio,
                        "q9_source_total_mmol_L": alpha * 0.01 * pool,
                        "termination_status": "completed",
                        "termination_time_h": 24.0,
                        "alpha_mmol_gDW": alpha,
                        "pool_multiplier": pool,
                        "runtime_gpr_scenario": scenario,
                        "runtime_gpr_mapping_sha256": mapping_sha,
                        "wt_dynamic_doublings": 10.0,
                        "wt_termination_status": "completed",
                        "wt_termination_time_h": 24.0,
                        "experimental_essential": False,
                        **{
                            f"essential_at_{cutoff * 100:g}pct": ratio < cutoff
                            for cutoff in (0.01, 0.05, 0.1, 0.15)
                        },
                    })
        pd.DataFrame(rows).to_csv(arm_dir / "essentiality_dynamic_calls.tsv", sep="\t", index=False)
        manifest = {
            "workflow": dfba.WORKFLOW,
            "schema_version": dfba.SCHEMA_VERSION,
            "run_id": f"run_{scenario}",
            "chunk_index": 0,
            "chunk_count": 1,
            "genes": list(NINE_COQ_CANDIDATES),
            "solver": "gurobi",
            "runtime_versions": {"python": "3.11", "cobra": "0.30", "gurobipy": "10.0.3"},
            "solver_feasibility_tolerance": 1e-9,
            "script_sha256": "script",
            "coq9_dilution_tool_sha256": "coq9-tool",
            "hours": 24.0,
            "dt_h": 0.0625,
            "initial_biomass_gDW_L": 0.01,
            "alphas_mmol_gDW": list(dfba.DEFAULT_ALPHAS),
            "alpha_sampling": [
                {
                    "sampler_id": "explicit_alpha_v1", "base_seed": None,
                    "replicate_id": None, "distribution": "explicit",
                    "low_mmol_gDW": None, "mode_mmol_gDW": None,
                    "high_mmol_gDW": None, "alpha_mmol_gDW": alpha,
                }
                for alpha in dfba.DEFAULT_ALPHAS
            ],
            "pool_multipliers": list(dfba.DEFAULT_POOL_MULTIPLIERS),
            "uracil_mode": "po1f_nonlimiting",
            "optimizer": "pfba",
            "nonoptimal_policy": "terminal_record_and_stop",
            "calibration_status": "sensitivity_only_not_calibrated",
            "q9_reserve_definition": "alpha * initial_biomass * pool_multiplier mmol/L",
            "q9_reserve_policy": dfba.Q9_RESERVE_POLICY,
            "coq9_dilution": [
                runtime_coq9_dilution_record(alpha, pool_source_enabled=True)
                for alpha in dfba.DEFAULT_ALPHAS
            ],
            "runtime_topology": {
                "enabled": True, "mapping_sha256": dfba.R39_R19_RUNTIME_MAPPING_SHA256,
            },
            "runtime_gpr_scenario": mapping_audit,
            "simulation_context": {"simulation_context_fingerprint": "context"},
            "input_sha256": {
                "model": "model", "medium": "medium", "profile": "profile",
                "experimental": "experimental",
            },
        }
        manifest["fingerprint"] = dfba.sha256_payload(manifest)
        (arm_dir / "chunk_000_manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8",
        )
        (arm_dir / "merge_manifest.json").write_text(json.dumps({
            "chunks": 1,
            "calls": 108,
            "chunk_fingerprints": [manifest["fingerprint"]],
            "runtime_gpr_scenario": mapping_audit,
        }), encoding="utf-8")

    evidence_path = tmp_path / "evidence.tsv"
    pd.DataFrame([{
        "gene_id": gene_id,
        "cas9_fitness": -1.0,
        "cas9_call": "essential",
        "cas12a_fitness": -0.5,
        "cas12a_call": "nonessential",
    } for gene_id in NINE_COQ_CANDIDATES]).to_csv(evidence_path, sep="\t", index=False)

    output_dir = write_hypothesis_matrix(arm_dirs, evidence_path, tmp_path / "summary")
    matrix = pd.read_csv(output_dir / "hypothesis_matrix.tsv", sep="\t")
    result_manifest = json.loads((output_dir / "hypothesis_matrix_manifest.json").read_text())
    assert len(matrix) == 756
    assert not matrix.duplicated(["mapping_id", "gene_id", "alpha_mmol_gDW", "pool_multiplier"]).any()
    assert result_manifest["checks"] == {
        "exact_seven_scenarios": True,
        "mapping_payloads_and_sha": True,
        "runtime_target_fingerprint_uniform": True,
        "frozen_context_uniform": True,
        "wt_endpoints_unchanged": True,
        "non_target_results_unchanged": True,
        "nested_absolute_constraints_nonincreasing": True,
        "strict_cutoffs_recomputed": True,
        "calls_level_pool_source_and_five_gene_controls": True,
        "evidence_exact_nine_gene_coverage": True,
    }

    calls_path = arm_dirs[-1] / "essentiality_dynamic_calls.tsv"
    valid_calls = calls_path.read_text()
    bad_calls = pd.read_csv(calls_path, sep="\t")
    bad_calls.loc[1, "q9_source_total_mmol_L"] = 999
    bad_calls.to_csv(calls_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="Calls-level QC failed"):
        write_hypothesis_matrix(arm_dirs, evidence_path, tmp_path / "bad_source_summary")
    calls_path.write_text(valid_calls)

    bad_context_manifest_path = arm_dirs[0] / "chunk_000_manifest.json"
    bad_context_manifest = json.loads(bad_context_manifest_path.read_text())
    bad_context_manifest["coq9_dilution_tool_sha256"] = "different-tool"
    bad_context_manifest["fingerprint"] = dfba.sha256_payload({
        key: value for key, value in bad_context_manifest.items() if key != "fingerprint"
    })
    bad_context_manifest_path.write_text(json.dumps(bad_context_manifest), encoding="utf-8")
    bad_context_merge_path = arm_dirs[0] / "merge_manifest.json"
    bad_context_merge = json.loads(bad_context_merge_path.read_text())
    bad_context_merge["chunk_fingerprints"] = [bad_context_manifest["fingerprint"]]
    bad_context_merge_path.write_text(json.dumps(bad_context_merge), encoding="utf-8")
    with pytest.raises(ValueError, match="Runtime context differs"):
        write_hypothesis_matrix(arm_dirs, evidence_path, tmp_path / "bad_context_summary")
    bad_context_manifest["coq9_dilution_tool_sha256"] = "coq9-tool"
    bad_context_manifest["fingerprint"] = dfba.sha256_payload({
        key: value for key, value in bad_context_manifest.items() if key != "fingerprint"
    })
    bad_context_manifest_path.write_text(json.dumps(bad_context_manifest), encoding="utf-8")
    bad_context_merge["chunk_fingerprints"] = [bad_context_manifest["fingerprint"]]
    bad_context_merge_path.write_text(json.dumps(bad_context_merge), encoding="utf-8")

    bad_manifest_path = arm_dirs[0] / "chunk_000_manifest.json"
    bad_manifest = json.loads(bad_manifest_path.read_text())
    bad_manifest["runtime_gpr_scenario"]["mapping_sha256"] = "bad"
    bad_manifest["fingerprint"] = dfba.sha256_payload({
        key: value for key, value in bad_manifest.items() if key != "fingerprint"
    })
    bad_manifest_path.write_text(json.dumps(bad_manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="mapping SHA mismatch"):
        write_hypothesis_matrix(arm_dirs, evidence_path, tmp_path / "bad_summary")
