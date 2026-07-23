from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest
from cobra.io import read_sbml_model

from scripts.gem_annotate.microspecies import load_curated_microspecies
from scripts.gem_annotate.reaction_chemistry import (
    DEFAULT_REACTION_CHEMISTRY_TABLE,
    attach_reference_chemistry_audit,
    audit_global_reference_component,
    audit_reference_reaction_chemistry,
    load_reaction_chemistry_proposals,
    write_reference_chemistry_audits,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = REPO_ROOT / "model.xml"
FRONTIER_EVIDENCE_PATH = (
    REPO_ROOT
    / "data"
    / "essentiality"
    / "evidence"
    / "chemistry"
    / "frontier_microspecies_sources.json"
)
COMPONENT_HISTORY_PATH = FRONTIER_EVIDENCE_PATH.with_name(
    "component_migration_history.json"
)
GLOBAL_COMPONENT_AUDIT_PATH = FRONTIER_EVIDENCE_PATH.with_name(
    "global_component_migration.json"
)
DEFERRED_FORMULA_EVIDENCE_PATH = FRONTIER_EVIDENCE_PATH.with_name(
    "deferred_formula_candidates.json"
)
REACTION_IDENTITY_FLAGS_PATH = FRONTIER_EVIDENCE_PATH.with_name(
    "reaction_identity_flags.json"
)


def _chemistry_snapshot(model) -> tuple:
    metabolite_state = tuple(
        sorted(
            (metabolite.id, metabolite.formula, metabolite.charge)
            for metabolite in model.metabolites
        )
    )
    reaction_state = tuple(
        sorted(
            (
                reaction.id,
                tuple(
                    sorted(
                        (metabolite.id, float(coefficient))
                        for metabolite, coefficient in reaction.metabolites.items()
                    )
                ),
            )
            for reaction in model.reactions
        )
    )
    return metabolite_state, reaction_state


def test_reference_chemistry_table_covers_open_cases_and_nad_component():
    proposals = load_reaction_chemistry_proposals()
    observed = {(proposal.case_id, proposal.reaction_id) for proposal in proposals}

    assert {
        ("EGC-1fd6c310af7f", "R60"),
        ("EGC-1fd6c310af7f", "R545"),
        ("EGC-324e120cbc71", "R2081"),
        ("EGC-324e120cbc71", "R334"),
        ("EGC-41703ef278c1", "R540"),
        ("EGC-549ae7942847", "R344"),
        ("EGC-7aa18c1a4b82", "R568"),
        ("EGC-7aa18c1a4b82", "R569"),
    } <= observed
    assert all(proposal.status == "component_review" for proposal in proposals)


def test_frontier_microspecies_evidence_matches_deferred_curation_and_model():
    evidence = json.loads(FRONTIER_EVIDENCE_PATH.read_text(encoding="utf-8"))
    rows = {row.family_id: row for row in load_curated_microspecies()}
    model = read_sbml_model(str(MODEL_PATH))

    for record in evidence["records"]:
        row = rows[record["family_id"]]
        assert record["status"] == row.status
        assert row.status in {"component_review", "verified_current"}
        assert record["source"] == row.source_url
        assert record["review_pair"] == {
            "formula": row.target_formula,
            "charge": row.target_charge,
        }
        observed_model_pairs = set()
        for metabolite_id in row.expected_metabolite_ids:
            metabolite = model.metabolites.get_by_id(metabolite_id)
            observed_model_pairs.add((metabolite.formula, metabolite.charge))
        documented_model_pairs = record.get(
            "model_pairs", [record.get("model_pair")]
        )
        assert observed_model_pairs == {
            (pair["formula"], pair["charge"])
            for pair in documented_model_pairs
            if isinstance(pair, dict)
        }


def test_deferred_formula_candidate_is_inert_and_precleanup_blocked():
    evidence = json.loads(
        DEFERRED_FORMULA_EVIDENCE_PATH.read_text(encoding="utf-8")
    )
    source_model = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))
    curated_families = {row.family_id for row in load_curated_microspecies()}

    record = evidence["records"][0]
    assert record["metabolite_id"] == "m1720[C_cy]"
    assert record["status"] == "blocked_precleanup_selector"
    assert record["activation_allowed"] is False
    source_metabolite = source_model.metabolites.get_by_id(record["metabolite_id"])
    assert "ActiveX VT_ERROR" in source_metabolite.name
    assert source_metabolite.formula is None
    assert source_metabolite.charge == 0
    assert "nicotinate_d_ribonucleoside" not in curated_families


def test_component_migration_history_latest_entry_matches_current_audit():
    history = json.loads(COMPONENT_HISTORY_PATH.read_text(encoding="utf-8"))
    audit = json.loads(GLOBAL_COMPONENT_AUDIT_PATH.read_text(encoding="utf-8"))
    latest = history["milestones"][-1]

    assert history["model_sha256"] == hashlib.sha256(
        MODEL_PATH.read_bytes()
    ).hexdigest()
    assert latest["microspecies_table_sha256"] == hashlib.sha256(
        (REPO_ROOT / "data" / "metabolite_microspecies.csv").read_bytes()
    ).hexdigest()
    assert latest["audit_file_sha256"] == hashlib.sha256(
        GLOBAL_COMPONENT_AUDIT_PATH.read_bytes()
    ).hexdigest()
    assert latest["audit_fingerprint"] == audit["audit_fingerprint"]
    assert latest["fixed_reactions"] == len(audit["fixed_reaction_ids"])
    assert latest["regressed_reactions"] == len(audit["regressed_reaction_ids"])
    assert latest["ready_for_activation"] is audit["ready_for_activation"] is False


@pytest.mark.parametrize(
    "case_id,target_reaction_ids",
    [
        (
            "EGC-1fd6c310af7f",
            {"R60", "R65", "R309", "R310", "R355", "R545", "R749"},
        ),
        ("EGC-324e120cbc71", {"R2081", "R334"}),
        ("EGC-41703ef278c1", {"R540"}),
        ("EGC-549ae7942847", {"R344"}),
        ("EGC-7aa18c1a4b82", {"R568", "R569"}),
    ],
)
def test_reference_equations_balance_but_remain_component_review(
    case_id, target_reaction_ids
):
    model = read_sbml_model(str(MODEL_PATH))
    before = _chemistry_snapshot(model)

    report = audit_reference_reaction_chemistry(model, case_id)

    assert report["reference_equations_balanced"] is True
    assert report["unbalanced_target_reaction_ids"] == []
    assert set(report["target_reaction_balances"]) == target_reaction_ids
    assert all(
        record["after"]["status"] == "balanced"
        for record in report["target_reaction_balances"].values()
    )
    assert report["ready_for_activation"] is False
    assert report["regressed_reaction_ids"]
    assert _chemistry_snapshot(model) == before


def test_r60_reference_changes_two_protons_to_one_without_mutating_model():
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R60")
    proton = model.metabolites.get_by_id("m10[C_cy]")
    assert reaction.metabolites[proton] == -2.0

    report = audit_reference_reaction_chemistry(model, "EGC-1fd6c310af7f")

    r60_change = next(
        change
        for change in report["coefficient_changes"]
        if change["reaction_id"] == "R60"
    )
    assert r60_change == {
        "reaction_id": "R60",
        "metabolite_id": "m10[C_cy]",
        "before_coefficient": -2.0,
        "after_coefficient": -1.0,
    }
    assert reaction.metabolites[proton] == -2.0


def test_r60_component_rejects_hydrogen_only_shortcuts():
    model = read_sbml_model(str(MODEL_PATH))

    report = audit_reference_reaction_chemistry(model, "EGC-1fd6c310af7f")

    assert report["preliminary_regressed_reaction_ids"]
    assert report["microspecies_bookkeeping_changes"] == []
    assert report["microspecies_bookkeeping_rejected"]
    assert {
        row["reason"] for row in report["microspecies_bookkeeping_rejected"]
    } <= {"proton_gate", "multi_compartment"}
    assert report["regressed_reaction_ids"]
    assert report["ready_for_activation"] is False


def test_nad_component_uses_explicit_reference_proton_coefficients():
    model = read_sbml_model(str(MODEL_PATH))
    original = {
        ("R344", "m10[C_cy]"): 2.0,
        ("R568", "m10[C_cy]"): 0.0,
        ("R569", "m627[C_nu]"): 0.0,
    }

    reports = [
        audit_reference_reaction_chemistry(model, "EGC-549ae7942847"),
        audit_reference_reaction_chemistry(model, "EGC-7aa18c1a4b82"),
    ]
    assert reports[0]["reference_reaction_ids"] == {"R344": "RHEA:10300"}
    assert reports[1]["reference_reaction_ids"] == {
        "R568": "RHEA:24385",
        "R569": "RHEA:24385",
    }

    observed = {
        (change["reaction_id"], change["metabolite_id"]): (
            change["before_coefficient"],
            change["after_coefficient"],
        )
        for report in reports
        for change in report["coefficient_changes"]
    }
    assert observed == {
        key: (before, 1.0) for key, before in original.items()
    }
    assert all(
        record["after"]["status"] == "balanced"
        for report in reports
        for record in report["target_reaction_balances"].values()
    )
    assert all(report["ready_for_activation"] is False for report in reports)
    for (reaction_id, metabolite_id), coefficient in original.items():
        reaction = model.reactions.get_by_id(reaction_id)
        metabolite = model.metabolites.get_by_id(metabolite_id)
        assert float(reaction.metabolites.get(metabolite, 0.0)) == coefficient


def test_global_component_audit_ranks_next_frontier_without_mutating_model():
    model = read_sbml_model(str(MODEL_PATH))
    before = _chemistry_snapshot(model)

    report = audit_global_reference_component(model)

    assert report["reference_equations_balanced"] is True
    assert report["unbalanced_proposal_reaction_ids"] == []
    assert len(report["fixed_reaction_ids"]) == 168
    assert len(report["regressed_reaction_ids"]) == 66
    assert len(report["unresolved_reaction_ids"]) == 584
    newly_closed = {
        "R533",
        "R534",
        "R535",
        "R538",
        "R574",
        "R575",
        "R580",
        "R581",
        "R582",
    }
    assert newly_closed <= set(report["fixed_reaction_ids"])
    assert report["reaction_balances"]["R576"]["after"]["status"] == "balanced"
    r490_closure = {"R138", "R216", "R490", "R785", "R1556"}
    assert all(
        report["reaction_balances"][reaction_id]["after"]["status"] == "balanced"
        for reaction_id in r490_closure
    )
    assert all(
        report["reaction_balances"][reaction_id]["after"]["status"]
        == "balanced"
        for reaction_id in newly_closed
    )
    assert {
        change["reaction_id"]
        for change in report["microspecies_bookkeeping_changes"]
    } == {"R567"}
    assert all(
        change["status"] == "provisional_requires_reference_review"
        and change["ready_for_activation"] is False
        for change in report["microspecies_bookkeeping_changes"]
    )
    assert {"R567"} <= set(report["frontier_regressed_reaction_ids"])
    assert len(report["microspecies_bookkeeping_rejected"]) == 66
    assert report["ready_for_activation"] is False
    assert report["actionable_frontier_ranking"][0]["candidate_family"] == (
        "gmp"
    )
    assert report["actionable_frontier_ranking"][0][
        "regressed_reaction_count"
    ] == 3
    assert {
        "isopentenyl_diphosphate",
        "two_oxoglutarate",
        "l_aspartate",
        "alpha_d_ribose_1_phosphate",
        "amp",
        "pyruvate",
        "oxaloacetate",
        "gdp",
        "nad",
        "nadh",
        "deamido_nad",
        "malate",
        "fumarate",
        "succinate",
        "s_dihydroorotate",
        "orotate",
        "nicotinate_d_ribonucleotide",
        "nicotinate",
        "quinolinate",
        "isocitrate",
        "glyoxylate",
        "cis_aconitate",
        "citrate",
        "s_ureidoglycolate",
        "allantoate",
        "ump",
    } <= set(report["family_ids"])
    assert _chemistry_snapshot(model) == before


def test_r538_balance_does_not_hide_the_unresolved_cofactor_identity():
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R538")

    assert {"m27[C_mi]", "m30[C_mi]"} <= {
        metabolite.id for metabolite in reaction.metabolites
    }
    assert {"m176[C_mi]", "m178[C_mi]"}.isdisjoint(
        metabolite.id for metabolite in reaction.metabolites
    )
    assert "NADP" in reaction.name
    assert "18253" in set(reaction.annotation["rhea"])

    evidence = json.loads(REACTION_IDENTITY_FLAGS_PATH.read_text(encoding="utf-8"))
    record = next(
        item for item in evidence["records"] if item["reaction_id"] == "R538"
    )
    assert record["current_equation_cofactor"] == "NAD+/NADH"
    assert record["verdict"] == "supported_metadata_correction_candidate"
    assert record["status"] == "awaiting_pipeline_gate"
    assert record["activation_allowed"] is False
    assert record["proposed_metadata_correction"]["kegg.reaction"] == "R00214"


def test_r580_gene_identity_flag_does_not_authorize_gpr_removal():
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R580")
    evidence = json.loads(REACTION_IDENTITY_FLAGS_PATH.read_text(encoding="utf-8"))
    record = next(
        item for item in evidence["records"] if item["reaction_id"] == "R580"
    )

    assert {gene.id for gene in reaction.genes} == {
        "YALI1A20611g",
        "YALI1E30400g",
    }
    assert " and " not in reaction.gene_reaction_rule
    assert record["verdict"] == "needs_more_evidence"
    assert record["status"] == "needs_more_evidence"
    assert record["activation_allowed"] is False
    assert "not complex subunits" in record["gpr_recommendation"]


def test_r1708_identity_evidence_is_metadata_only_and_keeps_both_gprs():
    model = read_sbml_model(str(MODEL_PATH))
    r1708 = model.reactions.get_by_id("R1708")
    r89 = model.reactions.get_by_id("R89")
    evidence = json.loads(REACTION_IDENTITY_FLAGS_PATH.read_text(encoding="utf-8"))
    record = next(
        item for item in evidence["records"] if item["reaction_id"] == "R1708"
    )

    assert r1708.annotation["ec-code"] == "3.1.2.1"
    assert record["proposed_metadata_correction"]["R1708.ec-code"] == "2.8.3.18"
    assert {gene.id for gene in r1708.genes} == {"YALI1E36437g"}
    assert {gene.id for gene in r89.genes} == {"YALI1E36437g"}
    assert record["verdict"] == "supported_metadata_correction_candidate"
    assert record["status"] == "awaiting_pipeline_gate"
    assert record["activation_allowed"] is False


def test_unknown_microspecies_family_fails_closed(tmp_path):
    with DEFAULT_REACTION_CHEMISTRY_TABLE.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        row = next(reader)
    row["required_microspecies_families"] = "invented_family"
    path = tmp_path / "bad_reaction_chemistry.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)

    with pytest.raises(ValueError, match="unknown microspecies families"):
        load_reaction_chemistry_proposals(path)


def test_nonfinite_reference_coefficient_fails_closed(tmp_path):
    with DEFAULT_REACTION_CHEMISTRY_TABLE.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        row = next(reader)
    row["coefficient_updates_json"] = json.dumps({"m10[C_cy]": float("inf")})
    path = tmp_path / "bad_coefficient.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)

    with pytest.raises(ValueError, match="finite number"):
        load_reaction_chemistry_proposals(path)


def test_durable_audit_writer_records_model_and_input_provenance(tmp_path):
    paths = write_reference_chemistry_audits(
        MODEL_PATH,
        tmp_path,
        case_ids=["EGC-1fd6c310af7f"],
    )

    assert paths == [tmp_path / "EGC-1fd6c310af7f.json"]
    payload = json.loads(paths[0].read_text(encoding="utf-8"))
    assert payload["case_id"] == "EGC-1fd6c310af7f"
    assert payload["model_sha256"] == (
        "b3f60933aa9503ab63ab5d8bca58a9525b8c81534d3eac72cf66d2769ff44f48"
    )
    assert payload["reaction_chemistry_table_sha256"]
    assert payload["microspecies_table_sha256"]
    assert payload["reference_equations_balanced"] is True
    assert payload["ready_for_activation"] is False


def test_guarded_audit_attachment_records_summary_without_promoting_balance(tmp_path):
    audit_path = write_reference_chemistry_audits(
        MODEL_PATH,
        tmp_path / "chemistry",
        case_ids=["EGC-1fd6c310af7f"],
    )[0]
    dossier_path = tmp_path / "EGC-1fd6c310af7f.json"
    dossier_path.write_text(
        json.dumps(
            {
                "case_id": "EGC-1fd6c310af7f",
                "model_sha256": (
                    "b3f60933aa9503ab63ab5d8bca58a9525b8c81534d3eac72cf66d2769ff44f48"
                ),
                "chemistry_review": {"status": "imbalanced"},
            }
        ),
        encoding="utf-8",
    )

    attached = attach_reference_chemistry_audit(dossier_path, audit_path)

    review = attached["chemistry_review"]
    assert review["status"] == "blocked_component_migration"
    assert review["reference_chemistry_audit"]["reference_equations_balanced"] is True
    assert review["reference_chemistry_audit"]["ready_for_activation"] is False
    assert review["reference_chemistry_audit"]["final_regressed_reactions"] > 0


def test_guarded_audit_attachment_rejects_stale_table_hash(tmp_path):
    audit_path = write_reference_chemistry_audits(
        MODEL_PATH,
        tmp_path / "chemistry",
        case_ids=["EGC-1fd6c310af7f"],
    )[0]
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["reaction_chemistry_table_sha256"] = "stale"
    audit_path.write_text(json.dumps(audit), encoding="utf-8")
    dossier_path = tmp_path / "EGC-1fd6c310af7f.json"
    dossier_path.write_text(
        json.dumps(
            {
                "case_id": "EGC-1fd6c310af7f",
                "model_sha256": (
                    "b3f60933aa9503ab63ab5d8bca58a9525b8c81534d3eac72cf66d2769ff44f48"
                ),
                "chemistry_review": {},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="stale chemistry audit"):
        attach_reference_chemistry_audit(dossier_path, audit_path)
