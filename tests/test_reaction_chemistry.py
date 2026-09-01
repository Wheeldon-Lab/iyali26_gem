from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest
from cobra.io import read_sbml_model

from scripts.gem_annotate.microspecies import load_curated_microspecies
from scripts.gem_annotate.reaction_chemistry import (
    attach_reference_chemistry_audit,
    audit_global_reference_component,
    audit_reference_reaction_chemistry,
    load_reaction_chemistry_proposals,
    write_reference_chemistry_audits,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = REPO_ROOT / "model.xml"
CHEMISTRY_EVIDENCE_DIR = Path("data/essentiality/evidence/chemistry")


@pytest.fixture
def chemistry_tables(external_data_file):
    return (
        external_data_file("data/essentiality/reaction_chemistry_curation.csv"),
        external_data_file("data/metabolite_microspecies.csv"),
    )


def _write_minimal_chemistry_tables(tmp_path: Path) -> tuple[Path, Path]:
    microspecies_path = tmp_path / "microspecies.csv"
    microspecies_path.write_text(
        "schema_version,status,family_id,selector_type,selector_value,"
        "target_formula,target_charge,reference_ph,chebi_id,rhea_id,source_url,"
        "min_matches,expected_metabolite_ids,allowed_current_pairs,rationale\n"
        "1,active,proton,base_name,H+,H,1,7.3,CHEBI:15378,,"
        "https://example.org/proton,1,h_c,H|1,test\n",
        encoding="utf-8",
    )
    reaction_path = tmp_path / "reaction_chemistry.csv"
    with reaction_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "schema_version",
                "status",
                "case_id",
                "reaction_id",
                "reference_reaction_id",
                "coefficient_updates_json",
                "required_microspecies_families",
                "evidence_url",
                "rationale",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "schema_version": 1,
                "status": "component_review",
                "case_id": "EGC-test",
                "reaction_id": "R_TEST",
                "reference_reaction_id": "RHEA:12345",
                "coefficient_updates_json": "{}",
                "required_microspecies_families": "proton",
                "evidence_url": "https://example.org/reaction",
                "rationale": "test",
            }
        )
    return reaction_path, microspecies_path


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


@pytest.mark.external_data
@pytest.mark.integration
def test_reference_chemistry_table_covers_open_cases_and_nad_component(
    chemistry_tables,
):
    table_path, microspecies_table_path = chemistry_tables
    proposals = load_reaction_chemistry_proposals(
        table_path, microspecies_table_path=microspecies_table_path
    )
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


@pytest.mark.external_data
@pytest.mark.integration
def test_frontier_microspecies_evidence_matches_deferred_curation_and_model(
    external_data_file, chemistry_tables
):
    evidence_path = external_data_file(
        CHEMISTRY_EVIDENCE_DIR / "frontier_microspecies_sources.json"
    )
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    _, microspecies_table_path = chemistry_tables
    rows = {
        row.family_id: row
        for row in load_curated_microspecies(microspecies_table_path)
    }
    model = read_sbml_model(str(MODEL_PATH))

    model_pair_mismatches = {}
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
        documented_model_pairs = {
            (pair["formula"], pair["charge"])
            for pair in documented_model_pairs
            if isinstance(pair, dict)
        }
        if observed_model_pairs != documented_model_pairs:
            model_pair_mismatches[record["family_id"]] = {
                "observed": observed_model_pairs,
                "documented": documented_model_pairs,
            }

    current_model_sha = hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest()
    assert {
        "model_sha256": evidence["model_sha256"],
        "model_pair_mismatches": model_pair_mismatches,
    } == {
        "model_sha256": current_model_sha,
        "model_pair_mismatches": {},
    }


def test_deferred_formula_source_candidate_is_inert():
    source_model = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))
    source_metabolite = source_model.metabolites.get_by_id("m1720[C_cy]")
    assert "ActiveX VT_ERROR" in source_metabolite.name
    assert source_metabolite.formula is None
    assert source_metabolite.charge == 0


@pytest.mark.external_data
@pytest.mark.integration
def test_deferred_formula_evidence_remains_blocked(
    external_data_file, chemistry_tables
):
    evidence_path = external_data_file(
        CHEMISTRY_EVIDENCE_DIR / "deferred_formula_candidates.json"
    )
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    _, microspecies_table_path = chemistry_tables
    curated_families = {
        row.family_id for row in load_curated_microspecies(microspecies_table_path)
    }
    record = evidence["records"][0]

    assert record["metabolite_id"] == "m1720[C_cy]"
    assert record["status"] == "blocked_precleanup_selector"
    assert record["activation_allowed"] is False
    assert "nicotinate_d_ribonucleoside" not in curated_families


@pytest.mark.external_data
@pytest.mark.integration
def test_component_migration_history_is_current_or_fails_closed(
    external_data_file, chemistry_tables
):
    history_path = external_data_file(
        CHEMISTRY_EVIDENCE_DIR / "component_migration_history.json"
    )
    audit_path = external_data_file(
        CHEMISTRY_EVIDENCE_DIR / "global_component_migration.json"
    )
    history = json.loads(history_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    latest = history["milestones"][-1]
    reaction_table_path, microspecies_table_path = chemistry_tables
    model = read_sbml_model(str(MODEL_PATH))
    live_audit = audit_global_reference_component(
        model,
        table_path=reaction_table_path,
        microspecies_table_path=microspecies_table_path,
    )

    current = {
        "model_sha256": hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest(),
        "reaction_chemistry_table_sha256": hashlib.sha256(
            reaction_table_path.read_bytes()
        ).hexdigest(),
        "microspecies_table_sha256": hashlib.sha256(
            microspecies_table_path.read_bytes()
        ).hexdigest(),
        "audit_fingerprint": live_audit["audit_fingerprint"],
        "fixed_reactions": len(live_audit["fixed_reaction_ids"]),
        "regressed_reactions": len(live_audit["regressed_reaction_ids"]),
        "unresolved_reactions": len(live_audit["unresolved_reaction_ids"]),
        "ready_for_activation": live_audit["ready_for_activation"],
    }
    history_state = {
        "model_sha256": history["model_sha256"],
        "reaction_chemistry_table_sha256": history[
            "reaction_chemistry_table_sha256"
        ],
        **{
            field: latest[field]
            for field in (
                "microspecies_table_sha256",
                "audit_fingerprint",
                "fixed_reactions",
                "regressed_reactions",
                "unresolved_reactions",
                "ready_for_activation",
            )
        },
    }
    stored_audit_state = {
        field: audit[field]
        for field in (
            "model_sha256",
            "reaction_chemistry_table_sha256",
            "microspecies_table_sha256",
            "audit_fingerprint",
            "ready_for_activation",
        )
    }
    stored_audit_state.update(
        fixed_reactions=len(audit["fixed_reaction_ids"]),
        regressed_reactions=len(audit["regressed_reaction_ids"]),
        unresolved_reactions=len(audit["unresolved_reaction_ids"]),
    )

    assert {"history": history_state, "stored_audit": stored_audit_state} == {
        "history": current,
        "stored_audit": current,
    }
    assert latest["audit_file_sha256"] == hashlib.sha256(
        audit_path.read_bytes()
    ).hexdigest()
    assert current["ready_for_activation"] is False


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
@pytest.mark.external_data
@pytest.mark.integration
def test_reference_equations_balance_but_remain_component_review(
    case_id, target_reaction_ids, chemistry_tables
):
    model = read_sbml_model(str(MODEL_PATH))
    before = _chemistry_snapshot(model)
    table_path, microspecies_table_path = chemistry_tables

    report = audit_reference_reaction_chemistry(
        model,
        case_id,
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )

    assert report["reference_equations_balanced"] is True
    assert report["unbalanced_target_reaction_ids"] == []
    assert set(report["target_reaction_balances"]) == target_reaction_ids
    assert all(
        record["after"]["status"] == "balanced"
        for record in report["target_reaction_balances"].values()
    )
    assert report["ready_for_activation"] is False
    assert (
        report["unbalanced_target_reaction_ids"]
        or report["regressed_reaction_ids"]
        or report["microspecies_bookkeeping_changes"]
    )
    assert _chemistry_snapshot(model) == before


@pytest.mark.external_data
@pytest.mark.integration
def test_r60_reference_changes_two_protons_to_one_without_mutating_model(
    chemistry_tables,
):
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R60")
    proton = model.metabolites.get_by_id("m10[C_cy]")
    assert reaction.metabolites[proton] == -2.0

    table_path, microspecies_table_path = chemistry_tables
    report = audit_reference_reaction_chemistry(
        model,
        "EGC-1fd6c310af7f",
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )

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


@pytest.mark.external_data
@pytest.mark.integration
def test_r60_component_rejects_hydrogen_only_shortcuts(chemistry_tables):
    model = read_sbml_model(str(MODEL_PATH))
    before = _chemistry_snapshot(model)
    table_path, microspecies_table_path = chemistry_tables

    report = audit_reference_reaction_chemistry(
        model,
        "EGC-1fd6c310af7f",
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )

    assert report["preliminary_regressed_reaction_ids"]
    changes = report["microspecies_bookkeeping_changes"]
    assert changes
    assert len({row["reaction_id"] for row in changes}) == len(changes)
    assert all(
        row["status"] == "provisional_requires_reference_review"
        and row["ready_for_activation"] is False
        and all(
            "H+" in model.metabolites.get_by_id(metabolite_id).name
            or model.metabolites.get_by_id(metabolite_id).name.startswith("H2O_")
            for metabolite_id in row["additions"]
        )
        for row in changes
    )
    assert report["microspecies_bookkeeping_rejected"]
    assert {
        row["reason"] for row in report["microspecies_bookkeeping_rejected"]
    } <= {"proton_gate", "multi_compartment"}
    assert report["regressed_reaction_ids"]
    assert report["ready_for_activation"] is False
    assert _chemistry_snapshot(model) == before


@pytest.mark.external_data
@pytest.mark.integration
def test_nad_component_uses_explicit_reference_proton_coefficients(
    chemistry_tables,
):
    model = read_sbml_model(str(MODEL_PATH))
    original = {
        ("R344", "m10[C_cy]"): 2.0,
        ("R568", "m10[C_cy]"): 0.0,
        ("R569", "m627[C_nu]"): 0.0,
    }

    table_path, microspecies_table_path = chemistry_tables
    reports = [
        audit_reference_reaction_chemistry(
            model,
            case_id,
            table_path=table_path,
            microspecies_table_path=microspecies_table_path,
        )
        for case_id in ("EGC-549ae7942847", "EGC-7aa18c1a4b82")
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


@pytest.mark.external_data
@pytest.mark.integration
def test_global_component_audit_is_safe_and_deterministic(chemistry_tables):
    model = read_sbml_model(str(MODEL_PATH))
    before = _chemistry_snapshot(model)
    table_path, microspecies_table_path = chemistry_tables

    report = audit_global_reference_component(
        model,
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )
    repeated = audit_global_reference_component(
        model,
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )

    assert repeated == report
    assert report["reference_equations_balanced"] is True
    assert report["unbalanced_proposal_reaction_ids"] == []
    model_reaction_ids = {reaction.id for reaction in model.reactions}
    audited_sets = [
        report[key]
        for key in (
            "fixed_reaction_ids",
            "regressed_reaction_ids",
            "unresolved_reaction_ids",
        )
    ]
    assert all(len(values) == len(set(values)) for values in audited_sets)
    assert set().union(*(set(values) for values in audited_sets)) <= model_reaction_ids
    assert set(report["fixed_reaction_ids"]).isdisjoint(
        report["unresolved_reaction_ids"]
    )
    assert all(
        change["status"] == "provisional_requires_reference_review"
        and change["ready_for_activation"] is False
        for change in report["microspecies_bookkeeping_changes"]
    )
    assert report["ready_for_activation"] is False
    assert _chemistry_snapshot(model) == before


def test_r538_current_cofactor_identity_is_preserved():
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


@pytest.mark.external_data
@pytest.mark.integration
def test_r538_identity_evidence_remains_inactive(external_data_file):
    evidence_path = external_data_file(
        CHEMISTRY_EVIDENCE_DIR / "reaction_identity_flags.json"
    )
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    record = next(
        item for item in evidence["records"] if item["reaction_id"] == "R538"
    )
    assert record["current_equation_cofactor"] == "NAD+/NADH"
    assert record["verdict"] == "supported_metadata_correction_candidate"
    assert record["status"] == "awaiting_pipeline_gate"
    assert record["activation_allowed"] is False
    assert record["proposed_metadata_correction"]["kegg.reaction"] == "R00214"


def test_r580_gene_associations_are_preserved():
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R580")

    assert {gene.id for gene in reaction.genes} == {
        "YALI1A20611g",
        "YALI1E30400g",
    }
    assert " and " not in reaction.gene_reaction_rule


@pytest.mark.external_data
@pytest.mark.integration
def test_r580_gene_identity_evidence_does_not_authorize_removal(
    external_data_file,
):
    evidence_path = external_data_file(
        CHEMISTRY_EVIDENCE_DIR / "reaction_identity_flags.json"
    )
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    record = next(
        item for item in evidence["records"] if item["reaction_id"] == "R580"
    )

    assert record["verdict"] == "needs_more_evidence"
    assert record["status"] == "needs_more_evidence"
    assert record["activation_allowed"] is False
    assert "not complex subunits" in record["gpr_recommendation"]


def test_r1708_and_r89_gene_associations_are_preserved():
    model = read_sbml_model(str(MODEL_PATH))
    r1708 = model.reactions.get_by_id("R1708")
    r89 = model.reactions.get_by_id("R89")

    assert r1708.annotation["ec-code"] == "3.1.2.1"
    assert {gene.id for gene in r1708.genes} == {"YALI1E36437g"}
    assert {gene.id for gene in r89.genes} == {"YALI1E36437g"}


@pytest.mark.external_data
@pytest.mark.integration
def test_r1708_identity_evidence_is_metadata_only(external_data_file):
    evidence_path = external_data_file(
        CHEMISTRY_EVIDENCE_DIR / "reaction_identity_flags.json"
    )
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    record = next(
        item for item in evidence["records"] if item["reaction_id"] == "R1708"
    )

    assert record["proposed_metadata_correction"]["R1708.ec-code"] == "2.8.3.18"
    assert record["verdict"] == "supported_metadata_correction_candidate"
    assert record["status"] == "awaiting_pipeline_gate"
    assert record["activation_allowed"] is False


def test_unknown_microspecies_family_fails_closed(tmp_path):
    source_path, microspecies_path = _write_minimal_chemistry_tables(tmp_path)
    with source_path.open(newline="", encoding="utf-8") as handle:
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
        load_reaction_chemistry_proposals(
            path, microspecies_table_path=microspecies_path
        )


def test_nonfinite_reference_coefficient_fails_closed(tmp_path):
    source_path, microspecies_path = _write_minimal_chemistry_tables(tmp_path)
    with source_path.open(newline="", encoding="utf-8") as handle:
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
        load_reaction_chemistry_proposals(
            path, microspecies_table_path=microspecies_path
        )


@pytest.mark.external_data
@pytest.mark.integration
def test_durable_audit_writer_records_model_and_input_provenance(
    tmp_path, chemistry_tables
):
    table_path, microspecies_table_path = chemistry_tables
    paths = write_reference_chemistry_audits(
        MODEL_PATH,
        tmp_path,
        case_ids=["EGC-1fd6c310af7f"],
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )

    assert paths == [tmp_path / "EGC-1fd6c310af7f.json"]
    payload = json.loads(paths[0].read_text(encoding="utf-8"))
    assert payload["case_id"] == "EGC-1fd6c310af7f"
    assert payload["model_sha256"] == (
        "bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee"
    )
    assert payload["reaction_chemistry_table_sha256"]
    assert payload["microspecies_table_sha256"]
    assert payload["reference_equations_balanced"] is True
    assert payload["ready_for_activation"] is False


@pytest.mark.external_data
@pytest.mark.integration
def test_guarded_audit_attachment_records_summary_without_promoting_balance(
    tmp_path, chemistry_tables
):
    table_path, microspecies_table_path = chemistry_tables
    audit_path = write_reference_chemistry_audits(
        MODEL_PATH,
        tmp_path / "chemistry",
        case_ids=["EGC-1fd6c310af7f"],
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )[0]
    dossier_path = tmp_path / "EGC-1fd6c310af7f.json"
    dossier_path.write_text(
        json.dumps(
            {
                "case_id": "EGC-1fd6c310af7f",
                "model_sha256": (
                    "bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee"
                ),
                "chemistry_review": {"status": "imbalanced"},
            }
        ),
        encoding="utf-8",
    )

    attached = attach_reference_chemistry_audit(
        dossier_path,
        audit_path,
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
    )

    review = attached["chemistry_review"]
    assert review["status"] == "blocked_component_migration"
    assert review["reference_chemistry_audit"]["reference_equations_balanced"] is True
    assert review["reference_chemistry_audit"]["ready_for_activation"] is False
    assert review["reference_chemistry_audit"]["final_regressed_reactions"] > 0


@pytest.mark.external_data
@pytest.mark.integration
def test_guarded_audit_attachment_rejects_stale_table_hash(
    tmp_path, chemistry_tables
):
    table_path, microspecies_table_path = chemistry_tables
    audit_path = write_reference_chemistry_audits(
        MODEL_PATH,
        tmp_path / "chemistry",
        case_ids=["EGC-1fd6c310af7f"],
        table_path=table_path,
        microspecies_table_path=microspecies_table_path,
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
                    "bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee"
                ),
                "chemistry_review": {},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="stale chemistry audit"):
        attach_reference_chemistry_audit(
            dossier_path,
            audit_path,
            table_path=table_path,
            microspecies_table_path=microspecies_table_path,
        )
