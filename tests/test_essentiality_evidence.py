import json

import pandas as pd
import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model, write_sbml_model

from scripts.gem_annotate import essentiality_evidence as evidence_module
from scripts.gem_annotate.essentiality_evidence import (
    LEDGER_FIELDS,
    assert_transition,
    chemistry_fingerprint,
    import_research_batch_results,
    merge_detected_cases,
    read_ledger,
    record_human_decision,
    sha256_file,
    stable_case_id,
    stable_source_id,
    target_fingerprint,
    transition_selected_batch_to_researching,
    transition_case_status,
    validate_evidence_dossier,
    write_ledger,
)
from scripts.gem_annotate.validate_essential_genes import (
    build_agent_cases,
    prepare_agent_case_files,
)


def _reaction_context(
    reaction_id: str = "R1", gpr: str = "g1 or g2"
) -> dict:
    return {
        "reaction_id": reaction_id,
        "stoichiometry": {"a_c": -1.0, "b_c": 1.0},
        "metabolite_chemistry": {
            "a_c": {"formula": "C", "charge": 0, "compartment": "c"},
            "b_c": {"formula": "C", "charge": 0, "compartment": "c"},
        },
        "lower_bound": 0.0,
        "upper_bound": 1000.0,
        "gpr": gpr,
        "gpr_gene_ids": sorted(
            token for token in gpr.replace("(", "").replace(")", "").split()
            if token not in {"and", "or"}
        ),
    }


def _identity_crosscheck(gene_id: str) -> dict:
    return {
        "database": "UniProt",
        "gene_id": gene_id,
        "accession": f"P-{gene_id}",
        "url": f"https://www.uniprot.org/uniprotkb/P-{gene_id}/entry",
        "status": "match",
    }


def _verified_identity_review(model_sha: str, gene_ids: list[str]) -> dict:
    return {
        "status": "verified",
        "model_sha256": model_sha,
        "reviewed_gene_ids": gene_ids,
        "genes": [
            {
                "gene_id": gene_id,
                "identity_status": "verified",
                "function_status": "verified",
                "functional_role": "experimentally supported enzyme role",
                "evidence_refs": [f"uniprot:{gene_id}", "paper:figure-2"],
            }
            for gene_id in gene_ids
        ],
    }


def _verified_chemistry_review(
    model_sha: str, fingerprint: str, reaction_ids: list[str]
) -> dict:
    return {
        "status": "verified_balanced",
        "model_sha256": model_sha,
        "chemistry_fingerprint": fingerprint,
        "ready_for_activation": True,
        "audited_reaction_ids": reaction_ids,
        "residuals_by_reaction": {
            reaction_id: {} for reaction_id in reaction_ids
        },
        "audit_path": "data/essentiality/evidence/chemistry/audit.json",
        "audit_sha256": "a" * 64,
    }


def _source_record(**updates) -> dict:
    source = {
        "url": "https://example.org/primary-paper",
        "doi": "10.0000/example",
        "title": "Direct characterization in Yarrowia lipolytica",
        "year": 2020,
        "species": "Yarrowia lipolytica",
        "strain": "PO1f",
        "culture_conditions": "30 C, defined glucose medium",
        "source_type": "primary_research",
        "evidence_type": "direct_experiment",
        "stance": "supports",
        "claim": "The proteins form one required complex.",
        "location": "Figure 2",
        "evidence_tags": ["complex_membership"],
        "genes": ["g1", "g2"],
        "reactions": ["R1"],
        "pathways": ["test pathway"],
        "methods": "Targeted deletion and complementation",
        "result": "Both proteins were required for the measured activity.",
        "condition_match": "partial",
        "condition_mismatch_reason": "Glucose defined medium was not SD-Leu.",
        "relevance": "Directly tests the represented enzyme function.",
        "confidence": "high",
    }
    source.update(updates)
    source["source_id"] = stable_source_id(source)
    return source


def _search_audit(sources: list[dict], *, direct_found: bool) -> dict:
    return {
        "searched_at": "2026-07-21T12:00:00+00:00",
        "databases": ["PubMed", "Crossref", "UniProt"],
        "queries": ["Yarrowia lipolytica g1 complex"],
        "inclusion_criteria": ["Direct or contextual evidence for the case claim"],
        "exclusion_criteria": ["No stable source record or unrelated pathway"],
        "screened_sources": [
            {"source_id": source["source_id"], "disposition": "included"}
            for source in sources
        ],
        "excluded_sources": [],
        "direct_evidence_found": direct_found,
        "direct_evidence_absence_note": (
            "not_applicable"
            if direct_found
            else "No direct Y. lipolytica experiment was found in the recorded search."
        ),
    }


def _valid_dossier() -> dict:
    context = _reaction_context()
    chemistry_hash = chemistry_fingerprint([context])
    gene_ids = ["g1", "g2"]
    source = _source_record()
    return {
        "case_id": "EGC-123456789abc",
        "model_sha256": "model-sha",
        "experimental_sha256": "experimental-sha",
        "media_sha256": "media-sha",
        "target_fingerprint": "sha256:target",
        "chemistry_fingerprint": chemistry_hash,
        "claim_under_review": "g1 and g2 may be complex subunits",
        "model_context": {"reactions": [context]},
        "primary_sources": [source],
        "identity_crosschecks": [
            _identity_crosscheck(gene_id) for gene_id in gene_ids
        ],
        "contradictions": [],
        "chemistry_review": _verified_chemistry_review(
            "model-sha", chemistry_hash, ["R1"]
        ),
        "identity_review": _verified_identity_review("model-sha", gene_ids),
        "verdict": "supported_patch_candidate",
        "proposed_operation": {
            "operation": "set_gpr",
            "target_id": "R1",
            "value": "g1 and g2",
        },
        "confidence": "high",
        "literature_review": {
            "search_audit": _search_audit([source], direct_found=True),
            "reasoning": "Direct experiment and identity checks support the claim.",
        },
        "adversarial_review": {
            "status": "complete",
            "verdict": "pass",
            "findings": [],
        },
        "human_decision": {
            "decision": "pending",
            "approved_by": "",
            "approved_at": "",
        },
    }


def _write_live_acceptance_inputs(tmp_path, dossier: dict):
    model = Model("acceptance-live-model")
    a = Metabolite("a_c", compartment="c", formula="C", charge=0)
    b = Metabolite("b_c", compartment="c", formula="C", charge=0)
    reaction = Reaction("R1")
    reaction.bounds = (0.0, 1000.0)
    reaction.add_metabolites({a: -1.0, b: 1.0})
    reaction.gene_reaction_rule = "g1 or g2"
    model.add_reactions([reaction])
    model_path = tmp_path / "model.xml"
    write_sbml_model(model, str(model_path))
    experimental_path = tmp_path / "experimental.csv"
    experimental_path.write_text("gene_id\ng1\n", encoding="utf-8")
    media_path = tmp_path / "media.csv"
    media_path.write_text("exchange,uptake\n", encoding="utf-8")
    context = {
        "reaction_id": reaction.id,
        "stoichiometry": {
            metabolite.id: float(coefficient)
            for metabolite, coefficient in reaction.metabolites.items()
        },
        "lower_bound": float(reaction.lower_bound),
        "upper_bound": float(reaction.upper_bound),
        "gpr": reaction.gene_reaction_rule,
        "gpr_gene_ids": sorted(gene.id for gene in reaction.genes),
        "metabolite_chemistry": {
            metabolite.id: {
                "formula": metabolite.formula,
                "charge": metabolite.charge,
                "compartment": metabolite.compartment,
            }
            for metabolite in reaction.metabolites
        },
    }
    chemistry_hash = chemistry_fingerprint([context])
    dossier.update(
        {
            "model_sha256": sha256_file(model_path),
            "experimental_sha256": sha256_file(experimental_path),
            "media_sha256": sha256_file(media_path),
            "target_fingerprint": target_fingerprint([context]),
            "chemistry_fingerprint": chemistry_hash,
            "model_context": {"reactions": [context]},
        }
    )
    dossier["chemistry_review"] = _verified_chemistry_review(
        dossier["model_sha256"], chemistry_hash, ["R1"]
    )
    dossier["identity_review"] = _verified_identity_review(
        dossier["model_sha256"], ["g1", "g2"]
    )
    return model_path, experimental_path, media_path


def _case_packet(
    case_id: str,
    verdict_index: int = 0,
    *,
    category: str = "isozyme_redundancy",
    priority: int = 1,
    ko_growth_ratio: float = 0.5,
) -> dict:
    gene_id = f"g{verdict_index}"
    reaction_id = f"R{verdict_index}"
    context = _reaction_context(reaction_id, gene_id)
    return {
        "schema_version": "2.0",
        "case_id": case_id,
        "category": category,
        "priority": priority,
        "ranking_reason": "synthetic test ranking",
        "gene_ids": [gene_id],
        "reaction_ids": [reaction_id],
        "model_sha256": "batch-model-sha",
        "experimental_sha256": "batch-experimental-sha",
        "media_sha256": "batch-media-sha",
        "target_fingerprint": target_fingerprint([context]),
        "chemistry_fingerprint": chemistry_fingerprint([context]),
        "claim_under_review": f"claim {verdict_index}",
        "model_context": {
            "reactions": [context],
            "diagnostics": [{"ko_growth_ratio": ko_growth_ratio}],
        },
    }


def _reviewer_result(packet: dict, *, supported: bool = False) -> dict:
    result = {
        "case_id": packet["case_id"],
        "claim_under_review": packet["claim_under_review"],
        "search_audit": _search_audit([], direct_found=False),
        "primary_sources": [],
        "identity_crosschecks": [],
        "contradictions": [],
        "verdict": "needs_more_evidence",
        "proposed_operation": {},
        "confidence": "medium",
        "reasoning": "The recorded search did not find direct case-specific evidence.",
        "unresolved_questions": ["more direct evidence required"],
    }
    if supported:
        valid = _valid_dossier()
        result.update(
            {
                "primary_sources": valid["primary_sources"],
                "search_audit": _search_audit(
                    valid["primary_sources"], direct_found=True
                ),
                "identity_crosschecks": [
                    _identity_crosscheck(gene_id) for gene_id in packet["gene_ids"]
                ],
                "verdict": "supported_patch_candidate",
                "proposed_operation": {
                    "operation": "set_gpr",
                    "target_id": packet["reaction_ids"][0],
                    "value": packet["gene_ids"][0],
                },
                "confidence": "high",
                "reasoning": "Direct case-specific evidence supports the operation.",
                "unresolved_questions": [],
            }
        )
    return result


def _skeptic_batch(packets: list[dict]) -> dict:
    return {
        "batch_case_ids": [packet["case_id"] for packet in packets],
        "results": [
            {
                "case_id": packet["case_id"],
                "status": "complete",
                "verdict": "pass",
                "findings": [],
                "unresolved_contradictions": [],
                "corrected_verdict": "",
                "confidence": "high",
            }
            for packet in packets
        ],
    }


def _prepare_batch_state(tmp_path):
    packets = [_case_packet(f"EGC-batch{i:08d}", i) for i in range(3)]
    extra = _case_packet("EGC-extra000000", 9)
    ledger_path = tmp_path / "ledger.csv"
    evidence_dir = tmp_path / "evidence"
    merge_detected_cases(
        [*packets, extra], ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    return packets, extra, ledger_path, evidence_dir


def test_input_sha_refresh_reuses_complete_literature_review(tmp_path) -> None:
    case = _case_packet("EGC-reuse000001", category="metabolic_bypass")
    ledger_path = tmp_path / "ledger.csv"
    evidence_dir = tmp_path / "evidence"
    merge_detected_cases(
        [case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    dossier_path = evidence_dir / f"{case['case_id']}.json"
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    source = _source_record(stance="contextual")
    review = _reviewer_result(case)
    review["primary_sources"] = [source]
    review["search_audit"] = _search_audit([source], direct_found=False)
    review["reasoning"] = "Preserve this raw reviewer reasoning and audit."
    dossier.update(
        {
            "primary_sources": [source],
            "identity_crosschecks": [_identity_crosscheck(case["gene_ids"][0])],
            "contradictions": [],
            "literature_review": review,
            "workflow_status": "reviewed",
        }
    )
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")

    refreshed_case = json.loads(json.dumps(case))
    refreshed_case["model_sha256"] = "new-model-sha"
    merge_detected_cases(
        [refreshed_case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    refreshed = json.loads(dossier_path.read_text(encoding="utf-8"))
    assert refreshed["primary_sources"] == [source]
    assert refreshed["literature_review"] == review
    assert refreshed["literature_review"]["search_audit"] == review["search_audit"]
    assert refreshed["reused_literature_provenance"]["previous_model_sha256"] == (
        case["model_sha256"]
    )
    assert refreshed.get("workflow_status") != "reviewed"
    assert not (evidence_dir / "archive").exists()


def test_case_id_and_target_fingerprint_are_deterministic() -> None:
    assert stable_case_id(
        "isozyme_redundancy", ["g2", "g1"], ["R2", "R1"]
    ) == stable_case_id("isozyme_redundancy", ["g1", "g2"], ["R1", "R2"])
    assert stable_case_id("isozyme_redundancy", ["g1"], ["R1"]).startswith("EGC-")

    first = {
        "reaction_id": "R1",
        "stoichiometry": {"b_c": 1.0, "a_c": -1.0},
        "lower_bound": 0.0,
        "upper_bound": 1000.0,
        "gpr": "g1 or g2",
    }
    reordered = {
        "gpr": "g1 or g2",
        "upper_bound": 1000.0,
        "stoichiometry": {"a_c": -1.0, "b_c": 1.0},
        "reaction_id": "R1",
        "lower_bound": 0.0,
    }
    assert target_fingerprint([first]) == target_fingerprint([reordered])
    changed = dict(first, upper_bound=10.0)
    assert target_fingerprint([first]) != target_fingerprint([changed])

    chemistry_context = _reaction_context()
    reordered_chemistry = dict(chemistry_context)
    reordered_chemistry["stoichiometry"] = {"b_c": 1.0, "a_c": -1.0}
    reordered_chemistry["metabolite_chemistry"] = {
        "b_c": chemistry_context["metabolite_chemistry"]["b_c"],
        "a_c": chemistry_context["metabolite_chemistry"]["a_c"],
    }
    assert chemistry_fingerprint([chemistry_context]) == chemistry_fingerprint(
        [reordered_chemistry]
    )
    changed_formula = json.loads(json.dumps(chemistry_context))
    changed_formula["metabolite_chemistry"]["a_c"]["formula"] = "C2"
    assert chemistry_fingerprint([chemistry_context]) != chemistry_fingerprint(
        [changed_formula]
    )
    changed_charge = json.loads(json.dumps(chemistry_context))
    changed_charge["metabolite_chemistry"]["a_c"]["charge"] = -1
    assert chemistry_fingerprint([chemistry_context]) != chemistry_fingerprint(
        [changed_charge]
    )


@pytest.mark.parametrize(
    "mutator,expected",
    [
        (lambda dossier: dossier.update(primary_sources=[]), "no supporting direct"),
        (
            lambda dossier: dossier["primary_sources"][0].update(
                species="Saccharomyces cerevisiae"
            ),
            "no supporting direct",
        ),
        (
            lambda dossier: dossier.update(
                adversarial_review={"status": "not_run", "verdict": ""}
            ),
            "adversarial review",
        ),
        (
            lambda dossier: dossier.update(
                contradictions=[
                    {
                        "claim": "direct counterevidence",
                        "resolution_status": "unresolved",
                    }
                ]
            ),
            "unresolved direct contradiction",
        ),
        (
            lambda dossier: dossier["chemistry_review"].update(
                status="blocked_component_migration"
            ),
            "reaction chemistry review is not verified_balanced",
        ),
        (
            lambda dossier: dossier["identity_review"].update(status="conflict"),
            "gene identity review is not verified",
        ),
        (
            lambda dossier: dossier["chemistry_review"].update(
                model_sha256="stale-model"
            ),
            "chemistry review model SHA",
        ),
        (
            lambda dossier: dossier["chemistry_review"].update(
                component_migration_audit={"ready_for_activation": False}
            ),
            "component_migration_audit is not ready",
        ),
        (
            lambda dossier: dossier["chemistry_review"].update(
                reference_chemistry_audit={"ready_for_activation": False}
            ),
            "reference_chemistry_audit is not ready",
        ),
        (
            lambda dossier: dossier["chemistry_review"].update(
                nested_component_audit=False
            ),
            "nested_component_audit is not ready",
        ),
        (
            lambda dossier: dossier["identity_review"].update(
                model_sha256="stale-model"
            ),
            "identity review model SHA",
        ),
    ],
)
def test_supported_patch_rejects_weak_or_unreviewed_evidence(mutator, expected) -> None:
    dossier = _valid_dossier()
    mutator(dossier)
    errors = validate_evidence_dossier(dossier, require_supported_patch=True)
    assert any(expected in error for error in errors)


def test_valid_direct_evidence_passes_before_human_decision() -> None:
    assert (
        validate_evidence_dossier(_valid_dossier(), require_supported_patch=True) == []
    )


def test_supported_patch_requires_stable_source_id_and_search_audit() -> None:
    dossier = _valid_dossier()
    dossier["primary_sources"][0]["source_id"] = "SRC-wrong"
    dossier["literature_review"].pop("search_audit")

    errors = validate_evidence_dossier(dossier, require_supported_patch=True)

    assert any("source_id must be" in error for error in errors)
    assert any("search_audit" in error for error in errors)


def test_reviewer_search_audit_is_persisted_verbatim(tmp_path) -> None:
    packets, _extra, ledger_path, evidence_dir = _prepare_batch_state(tmp_path)
    transition_selected_batch_to_researching(
        packets, ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    reviewers = [_reviewer_result(packet) for packet in packets]

    import_research_batch_results(
        packets,
        reviewers,
        _skeptic_batch(packets),
        ledger_path=ledger_path,
        evidence_dir=evidence_dir,
    )

    dossier = json.loads(
        (evidence_dir / f"{packets[0]['case_id']}.json").read_text()
    )
    assert dossier["literature_review"]["search_audit"] == reviewers[0][
        "search_audit"
    ]
    assert dossier["literature_review"]["reasoning"] == reviewers[0]["reasoning"]


@pytest.mark.parametrize(
    "missing_field",
    [
        "audit_path",
        "audit_sha256",
        "audited_reaction_ids",
        "residuals_by_reaction",
        "ready_for_activation",
    ],
)
def test_supported_patch_rejects_incomplete_chemistry_audit(missing_field) -> None:
    dossier = _valid_dossier()
    dossier["chemistry_review"].pop(missing_field)
    errors = validate_evidence_dossier(dossier, require_supported_patch=True)
    assert any("chemistry" in error for error in errors)


def test_supported_patch_rejects_unverified_gene_function() -> None:
    dossier = _valid_dossier()
    dossier["identity_review"]["genes"][1]["function_status"] = "inferred"
    errors = validate_evidence_dossier(dossier, require_supported_patch=True)
    assert any("g2 function_status is not verified" in error for error in errors)


def test_supported_patch_rejects_missing_gpr_partner_identity() -> None:
    dossier = _valid_dossier()
    dossier["identity_review"]["reviewed_gene_ids"] = ["g1"]
    dossier["identity_review"]["genes"] = dossier["identity_review"]["genes"][:1]
    dossier["identity_crosschecks"] = dossier["identity_crosschecks"][:1]
    errors = validate_evidence_dossier(dossier, require_supported_patch=True)
    assert any("required genes: ['g2']" in error for error in errors)


@pytest.mark.parametrize("field,value", [("formula", "C2"), ("charge", -1)])
def test_supported_patch_rejects_stale_context_chemistry(field, value) -> None:
    dossier = _valid_dossier()
    dossier["model_context"]["reactions"][0]["metabolite_chemistry"]["a_c"][
        field
    ] = value
    errors = validate_evidence_dossier(dossier, require_supported_patch=True)
    assert any(
        "chemistry_fingerprint does not match model_context chemistry" in error
        for error in errors
    )


def test_translation_patch_requires_composition_and_carrier_evidence() -> None:
    dossier = _valid_dossier()
    dossier["proposed_operation"] = {"operation": "couple_trna_biomass"}
    dossier["primary_sources"][0]["evidence_tags"] = ["biomass_composition"]
    errors = validate_evidence_dossier(dossier, require_supported_patch=True)
    assert any("carrier_conservation" in error for error in errors)


@pytest.mark.parametrize(
    "case_id,verdict",
    [
        ("EGC-1a767f7a2547", "experimental_conflict"),
        ("EGC-8c9e74291a8f", "outside_metabolic_scope"),
    ],
)
def test_known_nonpatch_verdicts_do_not_require_patch_evidence(case_id, verdict) -> None:
    dossier = _valid_dossier()
    dossier.update(
        {
            "case_id": case_id,
            "primary_sources": [],
            "identity_crosschecks": [],
            "verdict": verdict,
            "proposed_operation": {},
        }
    )
    assert validate_evidence_dossier(dossier) == []


def test_human_acceptance_cannot_skip_awaiting_human(tmp_path) -> None:
    dossier = _valid_dossier()
    model_path, experimental_path, media_path = _write_live_acceptance_inputs(
        tmp_path, dossier
    )
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    (evidence_dir / f"{dossier['case_id']}.json").write_text(
        json.dumps(dossier), encoding="utf-8"
    )
    ledger_path = tmp_path / "ledger.csv"
    base_row = {field: "" for field in LEDGER_FIELDS}
    base_row.update(
        {
            "case_id": dossier["case_id"],
            "status": "reviewed",
            "model_sha256": dossier["model_sha256"],
            "experimental_sha256": dossier["experimental_sha256"],
            "media_sha256": dossier["media_sha256"],
            "target_fingerprint": dossier["target_fingerprint"],
            "chemistry_fingerprint": dossier["chemistry_fingerprint"],
            "reaction_ids": "R1",
        }
    )
    write_ledger([base_row], ledger_path)
    with pytest.raises(ValueError, match="Cannot accept"):
        record_human_decision(
            dossier["case_id"],
            "accept",
            ledger_path=ledger_path,
            evidence_dir=evidence_dir,
        )

    base_row["status"] = "awaiting_human"
    write_ledger([base_row], ledger_path)
    accepted = record_human_decision(
        dossier["case_id"],
        "accept",
        ledger_path=ledger_path,
        evidence_dir=evidence_dir,
        model_path=model_path,
        experimental_path=experimental_path,
        media_path=media_path,
    )
    assert accepted["status"] == "accepted"
    assert accepted["approved_by"] == "human_user"
    assert accepted["approved_at"]
    persisted = json.loads(
        (evidence_dir / f"{dossier['case_id']}.json").read_text(encoding="utf-8")
    )
    assert persisted["acceptance_live_provenance"]["model_sha256"] == sha256_file(
        model_path
    )


def test_illegal_state_jump_is_rejected() -> None:
    with pytest.raises(ValueError, match="Illegal"):
        assert_transition("queued", "accepted")


def test_generic_transition_cannot_create_acceptance(tmp_path) -> None:
    dossier = _valid_dossier()
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    (evidence_dir / f"{dossier['case_id']}.json").write_text(
        json.dumps(dossier), encoding="utf-8"
    )
    ledger_path = tmp_path / "ledger.csv"
    row = {field: "" for field in LEDGER_FIELDS}
    row.update(
        {
            "case_id": dossier["case_id"],
            "status": "awaiting_human",
            "model_sha256": dossier["model_sha256"],
            "target_fingerprint": dossier["target_fingerprint"],
            "chemistry_fingerprint": dossier["chemistry_fingerprint"],
        }
    )
    write_ledger([row], ledger_path)
    with pytest.raises(ValueError, match="record_human_decision"):
        transition_case_status(
            dossier["case_id"],
            "accepted",
            ledger_path=ledger_path,
            evidence_dir=evidence_dir,
        )


def test_changed_input_sha_requeues_case_without_reusing_simulation_state(tmp_path) -> None:
    context = _reaction_context()
    case = {
        "case_id": "EGC-abcdef123456",
        "category": "isozyme_redundancy",
        "gene_ids": ["g1"],
        "reaction_ids": ["R1"],
        "model_sha256": "model-a",
        "experimental_sha256": "experiment-a",
        "media_sha256": "medium-a",
        "target_fingerprint": target_fingerprint([context]),
        "chemistry_fingerprint": chemistry_fingerprint([context]),
        "claim_under_review": "test claim",
        "model_context": {"reactions": [context]},
    }
    ledger_path = tmp_path / "ledger.csv"
    evidence_dir = tmp_path / "evidence"
    merge_detected_cases(
        [case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    rows = read_ledger(ledger_path)
    rows[0]["status"] = "awaiting_human"
    write_ledger(rows, ledger_path)
    dossier_path = evidence_dir / f"{case['case_id']}.json"
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    dossier["primary_sources"] = [{"title": "reusable direct literature"}]
    dossier["adversarial_review"] = {
        "status": "complete",
        "verdict": "pass",
        "findings": [],
    }
    dossier["human_decision"] = {
        "decision": "accepted",
        "approved_by": "human_user",
        "approved_at": "earlier",
    }
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")

    changed = dict(case, model_sha256="model-b")
    merge_detected_cases(
        [changed], ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    refreshed = read_ledger(ledger_path)[0]
    assert refreshed["status"] == "queued"
    assert refreshed["previous_status"] == "awaiting_human"
    assert "model_sha256" in refreshed["stale_reason"]
    refreshed_dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    assert refreshed_dossier["model_sha256"] == "model-b"
    assert refreshed_dossier["primary_sources"] == [
        {"title": "reusable direct literature"}
    ]
    assert refreshed_dossier["adversarial_review"]["status"] == "not_run"
    assert refreshed_dossier["human_decision"]["decision"] == "pending"


def _make_same_sha_legacy_case(tmp_path):
    case = _case_packet("EGC-legacy000001", 17)
    ledger_path = tmp_path / "ledger.csv"
    evidence_dir = tmp_path / "evidence"
    merge_detected_cases(
        [case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    rows = read_ledger(ledger_path)
    rows[0].update(
        {
            "status": "needs_more_evidence",
            "chemistry_fingerprint": "",
            "human_decision": "deferred",
            "updated_at": "legacy-ledger-time",
        }
    )
    write_ledger(rows, ledger_path)

    dossier_path = evidence_dir / f"{case['case_id']}.json"
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    dossier.pop("chemistry_fingerprint")
    for context in dossier["model_context"]["reactions"]:
        context.pop("metabolite_chemistry")
    dossier.update(
        {
            "primary_sources": [{"title": "preserve this direct source"}],
            "identity_crosschecks": [{"database": "UniProt", "status": "match"}],
            "contradictions": [{"claim": "preserve this contradiction"}],
            "verdict": "needs_more_evidence",
            "confidence": "medium",
            "adversarial_review": {
                "status": "complete",
                "verdict": "needs_more_evidence",
                "findings": ["preserve skeptic review"],
            },
            "human_decision": {
                "decision": "deferred",
                "approved_by": "",
                "approved_at": "",
            },
            "workflow_status": "needs_more_evidence",
            "workflow_updated_at": "legacy-workflow-time",
            "chemistry_review": {
                "status": "not_run",
                "ready_for_activation": False,
            },
        }
    )
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")
    return case, ledger_path, evidence_dir, dossier_path


def test_same_sha_legacy_chemistry_backfill_preserves_evidence_and_state(
    tmp_path,
) -> None:
    case, ledger_path, evidence_dir, dossier_path = _make_same_sha_legacy_case(
        tmp_path
    )

    merge_detected_cases(
        [case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    row = read_ledger(ledger_path)[0]
    assert row["status"] == "needs_more_evidence"
    assert row["human_decision"] == "deferred"
    assert row["chemistry_fingerprint"] == case["chemistry_fingerprint"]
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    assert dossier["primary_sources"] == [
        {"title": "preserve this direct source"}
    ]
    assert dossier["adversarial_review"]["findings"] == [
        "preserve skeptic review"
    ]
    assert dossier["human_decision"]["decision"] == "deferred"
    assert dossier["workflow_status"] == "needs_more_evidence"
    assert dossier["workflow_updated_at"] == "legacy-workflow-time"
    assert dossier["chemistry_fingerprint"] == case["chemistry_fingerprint"]
    assert dossier["model_context"]["reactions"][0][
        "metabolite_chemistry"
    ] == case["model_context"]["reactions"][0]["metabolite_chemistry"]
    assert dossier["chemistry_fingerprint_migration"]["migration_type"] == (
        "same_sha_legacy_backfill"
    )
    assert dossier["chemistry_review"] == {
        "status": "not_run",
        "ready_for_activation": False,
    }
    assert "chemistry_fingerprint" not in dossier["chemistry_review"]
    assert not (evidence_dir / "archive").exists()


def test_same_sha_legacy_backfill_supports_empty_reaction_context(tmp_path) -> None:
    case = {
        "schema_version": "2.0",
        "case_id": "EGC-empty000001",
        "category": "inactive_reaction",
        "gene_ids": ["g-empty"],
        "reaction_ids": [],
        "model_sha256": "same-model-sha",
        "experimental_sha256": "same-experimental-sha",
        "media_sha256": "same-media-sha",
        "target_fingerprint": target_fingerprint([]),
        "chemistry_fingerprint": chemistry_fingerprint([]),
        "claim_under_review": "No model reaction is linked to this gene.",
        "model_context": {"reactions": []},
    }
    ledger_path = tmp_path / "ledger.csv"
    evidence_dir = tmp_path / "evidence"
    merge_detected_cases(
        [case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    rows = read_ledger(ledger_path)
    rows[0]["status"] = "reviewed"
    rows[0]["chemistry_fingerprint"] = ""
    write_ledger(rows, ledger_path)
    dossier_path = evidence_dir / f"{case['case_id']}.json"
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    dossier.pop("chemistry_fingerprint")
    dossier["primary_sources"] = [{"title": "preserved empty-context evidence"}]
    dossier["workflow_status"] = "reviewed"
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")

    merge_detected_cases(
        [case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    row = read_ledger(ledger_path)[0]
    assert row["status"] == "reviewed"
    assert row["reaction_ids"] == ""
    assert row["chemistry_fingerprint"] == chemistry_fingerprint([])
    migrated = json.loads(dossier_path.read_text(encoding="utf-8"))
    assert migrated["model_context"]["reactions"] == []
    assert migrated["primary_sources"] == [
        {"title": "preserved empty-context evidence"}
    ]
    assert migrated["workflow_status"] == "reviewed"
    assert migrated["chemistry_fingerprint"] == chemistry_fingerprint([])
    assert migrated["chemistry_fingerprint_migration"]["migration_type"] == (
        "same_sha_legacy_backfill"
    )
    assert not (evidence_dir / "archive").exists()


@pytest.mark.parametrize("conflict", ["partial_chemistry", "changed_bounds"])
def test_legacy_partial_or_conflicting_chemistry_is_archived_and_reset(
    tmp_path, conflict
) -> None:
    case, ledger_path, evidence_dir, dossier_path = _make_same_sha_legacy_case(
        tmp_path
    )
    dossier = json.loads(dossier_path.read_text(encoding="utf-8"))
    context = dossier["model_context"]["reactions"][0]
    if conflict == "partial_chemistry":
        context["metabolite_chemistry"] = {
            "a_c": case["model_context"]["reactions"][0][
                "metabolite_chemistry"
            ]["a_c"]
        }
    else:
        context["upper_bound"] = 999.0
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")

    merge_detected_cases(
        [case], ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    row = read_ledger(ledger_path)[0]
    assert row["status"] == "queued"
    assert row["previous_status"] == "needs_more_evidence"
    assert "chemistry_fingerprint" in row["stale_reason"]
    reset = json.loads(dossier_path.read_text(encoding="utf-8"))
    assert reset["primary_sources"] == []
    assert reset["human_decision"]["decision"] == "pending"
    archived = list((evidence_dir / "archive").glob("*.json"))
    assert len(archived) == 1
    archived_dossier = json.loads(archived[0].read_text(encoding="utf-8"))
    assert archived_dossier["primary_sources"] == [
        {"title": "preserve this direct source"}
    ]


def test_internally_inconsistent_fresh_fingerprint_writes_nothing(tmp_path) -> None:
    case, ledger_path, evidence_dir, _dossier_path = _make_same_sha_legacy_case(
        tmp_path
    )
    before = {
        path: path.read_bytes()
        for path in tmp_path.rglob("*")
        if path.is_file()
    }
    inconsistent = json.loads(json.dumps(case))
    inconsistent["target_fingerprint"] = "sha256:internally-wrong"

    with pytest.raises(ValueError, match="internally inconsistent"):
        merge_detected_cases(
            [inconsistent], ledger_path=ledger_path, evidence_dir=evidence_dir
        )

    after = {
        path: path.read_bytes()
        for path in tmp_path.rglob("*")
        if path.is_file()
    }
    assert after == before


def test_merge_detected_cases_rolls_back_entire_migration_bundle(
    tmp_path, monkeypatch
) -> None:
    case, ledger_path, evidence_dir, dossier_path = _make_same_sha_legacy_case(
        tmp_path
    )
    before_ledger = ledger_path.read_bytes()
    before_dossier = dossier_path.read_bytes()
    real_replace = evidence_module.os.replace
    failed = False

    def fail_once_on_ledger(source, destination):
        nonlocal failed
        if not failed and str(destination) == str(ledger_path):
            failed = True
            raise OSError("injected ledger replacement failure")
        return real_replace(source, destination)

    monkeypatch.setattr(evidence_module.os, "replace", fail_once_on_ledger)

    with pytest.raises(OSError, match="injected ledger replacement failure"):
        merge_detected_cases(
            [case], ledger_path=ledger_path, evidence_dir=evidence_dir
        )

    assert failed
    assert ledger_path.read_bytes() == before_ledger
    assert dossier_path.read_bytes() == before_dossier
    assert not (evidence_dir / "archive").exists()


def test_agent_cases_group_shared_reaction_and_rank_counterfactual_first() -> None:
    model = Model("mini")
    a = Metabolite("a_c", compartment="c", formula="C", charge=0)
    b = Metabolite("b_c", compartment="c", formula="C", charge=0)
    reaction = Reaction("R1")
    reaction.name = "shared enzyme"
    reaction.bounds = (0.0, 1000.0)
    reaction.add_metabolites({a: -1.0, b: 1.0})
    reaction.gene_reaction_rule = "g1 or g2"
    model.add_reactions([reaction])

    diagnostics = pd.DataFrame(
        [
            {
                "gene_id": gene_id,
                "category": "isozyme_redundancy",
                "reaction_ids": "R1",
                "ko_growth_ratio": 1.0,
                "all_linked_reactions_closed_growth_ratio": 0.0,
                "n_verified_single_bypasses": 0,
                "verified_single_bypasses": "",
                "wt_capacity": 1.0,
                "ko_capacity": 1.0,
            }
            for gene_id in ("g1", "g2")
        ]
    )
    per_gene = pd.DataFrame(
        [
            {
                "gene_id": gene_id,
                "essential_at_1pct": False,
                "essential_at_5pct": False,
                "essential_at_10pct": False,
                "essential_at_15pct": False,
            }
            for gene_id in ("g1", "g2")
        ]
    )
    summary = {
        "model": {"sha256": "model-sha"},
        "experimental": {"sha256": "experiment-sha"},
        "medium": {"sha256": "medium-sha"},
        "wt_growth": 1.0,
        "cutoff_curve": [
            {"cutoff_fraction_of_wt": value} for value in (0.01, 0.05, 0.10, 0.15)
        ],
    }
    cases = build_agent_cases(model, diagnostics, per_gene, summary, 0.10)

    assert len(cases) == 1
    assert cases[0]["gene_ids"] == ["g1", "g2"]
    assert cases[0]["reaction_ids"] == ["R1"]
    assert cases[0]["priority"] == 1
    assert cases[0]["target_fingerprint"].startswith("sha256:")
    assert cases[0]["chemistry_fingerprint"].startswith("sha256:")
    chemistry = cases[0]["model_context"]["reactions"][0][
        "metabolite_chemistry"
    ]
    assert chemistry["a_c"] == {
        "formula": "C",
        "charge": 0,
        "compartment": "c",
    }


def test_selected_batch_transition_changes_only_exact_three_cases_atomically(
    tmp_path,
) -> None:
    packets, extra, ledger_path, evidence_dir = _prepare_batch_state(tmp_path)

    updated = transition_selected_batch_to_researching(
        packets, ledger_path=ledger_path, evidence_dir=evidence_dir
    )

    assert [row["case_id"] for row in updated] == [
        packet["case_id"] for packet in packets
    ]
    rows = {row["case_id"]: row for row in read_ledger(ledger_path)}
    assert all(rows[packet["case_id"]]["status"] == "researching" for packet in packets)
    assert rows[extra["case_id"]]["status"] == "queued"
    for packet in packets:
        dossier = json.loads(
            (evidence_dir / f"{packet['case_id']}.json").read_text(encoding="utf-8")
        )
        assert dossier["workflow_status"] == "researching"


def test_metabolic_batch_filter_precedes_durable_merge_and_ranks_cases(
    tmp_path,
) -> None:
    output_dir = tmp_path / "results"
    output_dir.mkdir()
    ledger_path = tmp_path / "ledger.csv"
    evidence_dir = tmp_path / "evidence"
    unrelated = _case_packet(
        "EGC-unrelated0001", 9, category="isozyme_redundancy", priority=1
    )
    merge_detected_cases(
        [unrelated], ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    unrelated_row_before = {
        row["case_id"]: row for row in read_ledger(ledger_path)
    }[unrelated["case_id"]].copy()
    unrelated_dossier = evidence_dir / f"{unrelated['case_id']}.json"
    unrelated_dossier_before = unrelated_dossier.read_bytes()

    cases = [
        unrelated,
        _case_packet(
            "EGC-metabolic-low", 1,
            category="metabolic_bypass", priority=2, ko_growth_ratio=0.2,
        ),
        _case_packet(
            "EGC-metabolic-high", 2,
            category="metabolic_bypass", priority=2, ko_growth_ratio=0.8,
        ),
        _case_packet(
            "EGC-metabolic-unverified", 3,
            category="metabolic_bypass", priority=5, ko_growth_ratio=0.9,
        ),
    ]
    result = prepare_agent_case_files(
        cases,
        output_dir,
        3,
        ledger_path=ledger_path,
        evidence_dir=evidence_dir,
        case_category="metabolic_bypass",
        require_full_batch=True,
    )

    batch = json.loads((output_dir / "essentiality_agent_batch.json").read_text())
    assert [case["case_id"] for case in batch] == [
        "EGC-metabolic-high",
        "EGC-metabolic-low",
        "EGC-metabolic-unverified",
    ]
    assert {case["category"] for case in batch} == {"metabolic_bypass"}
    assert result["case_category"] == "metabolic_bypass"
    unrelated_row_after = {
        row["case_id"]: row for row in read_ledger(ledger_path)
    }[unrelated["case_id"]]
    assert unrelated_row_after == unrelated_row_before
    assert unrelated_dossier.read_bytes() == unrelated_dossier_before


def test_metabolic_batch_requires_three_before_any_durable_write(tmp_path) -> None:
    output_dir = tmp_path / "results"
    output_dir.mkdir()
    ledger_path = tmp_path / "ledger.csv"
    evidence_dir = tmp_path / "evidence"
    cases = [
        _case_packet(
            f"EGC-metabolic{i}", i,
            category="metabolic_bypass", priority=2,
        )
        for i in range(2)
    ]
    cases.append(
        _case_packet(
            "EGC-nutrient", 7, category="nutrient_bypass", priority=1
        )
    )

    with pytest.raises(ValueError, match="Only 2 queued metabolic_bypass"):
        prepare_agent_case_files(
            cases,
            output_dir,
            3,
            ledger_path=ledger_path,
            evidence_dir=evidence_dir,
            case_category="metabolic_bypass",
            require_full_batch=True,
        )

    assert not ledger_path.exists()
    assert not evidence_dir.exists()
    assert not (output_dir / "essentiality_agent_batch.json").exists()


def test_selected_batch_preflight_failure_writes_nothing(tmp_path) -> None:
    packets, _extra, ledger_path, evidence_dir = _prepare_batch_state(tmp_path)
    rows = read_ledger(ledger_path)
    rows_by_id = {row["case_id"]: row for row in rows}
    rows_by_id[packets[-1]["case_id"]]["status"] = "reviewed"
    write_ledger(rows, ledger_path)
    paths = [ledger_path] + [
        evidence_dir / f"{packet['case_id']}.json" for packet in packets
    ]
    before = {path: path.read_bytes() for path in paths}

    with pytest.raises(ValueError, match="must be queued"):
        transition_selected_batch_to_researching(
            packets, ledger_path=ledger_path, evidence_dir=evidence_dir
        )

    assert {path: path.read_bytes() for path in paths} == before


def test_complete_review_batch_moves_supported_pass_only_to_awaiting_human(
    tmp_path,
) -> None:
    packets, extra, ledger_path, evidence_dir = _prepare_batch_state(tmp_path)
    transition_selected_batch_to_researching(
        packets, ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    supported_dossier_path = evidence_dir / f"{packets[0]['case_id']}.json"
    supported_dossier = json.loads(
        supported_dossier_path.read_text(encoding="utf-8")
    )
    supported_dossier["chemistry_review"] = _verified_chemistry_review(
        packets[0]["model_sha256"],
        packets[0]["chemistry_fingerprint"],
        packets[0]["reaction_ids"],
    )
    supported_dossier["identity_review"] = _verified_identity_review(
        packets[0]["model_sha256"], packets[0]["gene_ids"]
    )
    supported_dossier_path.write_text(
        json.dumps(supported_dossier), encoding="utf-8"
    )
    reviewer_results = [
        _reviewer_result(packet, supported=index == 0)
        for index, packet in enumerate(packets)
    ]

    result = import_research_batch_results(
        packets,
        reviewer_results,
        _skeptic_batch(packets),
        ledger_path=ledger_path,
        evidence_dir=evidence_dir,
    )

    assert result["awaiting_human_case_ids"] == [packets[0]["case_id"]]
    assert result["reviewed_case_ids"] == [
        packets[1]["case_id"],
        packets[2]["case_id"],
    ]
    assert result["accepted_case_ids"] == []
    rows = {row["case_id"]: row for row in read_ledger(ledger_path)}
    assert rows[packets[0]["case_id"]]["status"] == "awaiting_human"
    assert rows[packets[1]["case_id"]]["status"] == "reviewed"
    assert rows[packets[2]["case_id"]]["status"] == "reviewed"
    assert rows[extra["case_id"]]["status"] == "queued"
    assert not any(row["status"] == "accepted" for row in rows.values())
    supported_dossier = json.loads(
        (evidence_dir / f"{packets[0]['case_id']}.json").read_text(
            encoding="utf-8"
        )
    )
    assert supported_dossier["adversarial_review"]["verdict"] == "pass"
    assert supported_dossier["workflow_status"] == "awaiting_human"


def test_incomplete_reviewer_batch_is_rejected_without_any_writes(tmp_path) -> None:
    packets, _extra, ledger_path, evidence_dir = _prepare_batch_state(tmp_path)
    transition_selected_batch_to_researching(
        packets, ledger_path=ledger_path, evidence_dir=evidence_dir
    )
    paths = [ledger_path] + [
        evidence_dir / f"{packet['case_id']}.json" for packet in packets
    ]
    before = {path: path.read_bytes() for path in paths}
    reviewer_results = [_reviewer_result(packet) for packet in packets[:-1]]

    with pytest.raises(ValueError, match="exactly cover"):
        import_research_batch_results(
            packets,
            reviewer_results,
            _skeptic_batch(packets),
            ledger_path=ledger_path,
            evidence_dir=evidence_dir,
        )

    assert {path: path.read_bytes() for path in paths} == before


@pytest.mark.parametrize("stale_input", ["model", "experimental", "media"])
def test_acceptance_recomputes_all_live_input_shas(tmp_path, stale_input) -> None:
    dossier = _valid_dossier()
    model_path, experimental_path, media_path = _write_live_acceptance_inputs(
        tmp_path, dossier
    )
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    dossier_path = evidence_dir / f"{dossier['case_id']}.json"
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")
    row = {field: "" for field in LEDGER_FIELDS}
    row.update(
        {
            "case_id": dossier["case_id"],
            "status": "awaiting_human",
            "reaction_ids": "R1",
            "model_sha256": dossier["model_sha256"],
            "experimental_sha256": dossier["experimental_sha256"],
            "media_sha256": dossier["media_sha256"],
            "target_fingerprint": dossier["target_fingerprint"],
            "chemistry_fingerprint": dossier["chemistry_fingerprint"],
        }
    )
    ledger_path = tmp_path / "ledger.csv"
    write_ledger([row], ledger_path)
    stale_path = {
        "model": model_path,
        "experimental": experimental_path,
        "media": media_path,
    }[stale_input]
    with stale_path.open("ab") as handle:
        handle.write(b"stale")
    before_ledger = ledger_path.read_bytes()
    before_dossier = dossier_path.read_bytes()

    with pytest.raises(ValueError, match="does not match the durable ledger"):
        record_human_decision(
            dossier["case_id"],
            "accept",
            ledger_path=ledger_path,
            evidence_dir=evidence_dir,
            model_path=model_path,
            experimental_path=experimental_path,
            media_path=media_path,
        )

    assert ledger_path.read_bytes() == before_ledger
    assert dossier_path.read_bytes() == before_dossier


def test_acceptance_recomputes_live_target_fingerprint(tmp_path) -> None:
    dossier = _valid_dossier()
    model_path, experimental_path, media_path = _write_live_acceptance_inputs(
        tmp_path, dossier
    )
    dossier["target_fingerprint"] = "sha256:stale-target"
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    dossier_path = evidence_dir / f"{dossier['case_id']}.json"
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")
    row = {field: "" for field in LEDGER_FIELDS}
    row.update(
        {
            "case_id": dossier["case_id"],
            "status": "awaiting_human",
            "reaction_ids": "R1",
            "model_sha256": dossier["model_sha256"],
            "experimental_sha256": dossier["experimental_sha256"],
            "media_sha256": dossier["media_sha256"],
            "target_fingerprint": dossier["target_fingerprint"],
            "chemistry_fingerprint": dossier["chemistry_fingerprint"],
        }
    )
    ledger_path = tmp_path / "ledger.csv"
    write_ledger([row], ledger_path)

    with pytest.raises(ValueError, match="Current target fingerprint"):
        record_human_decision(
            dossier["case_id"],
            "accept",
            ledger_path=ledger_path,
            evidence_dir=evidence_dir,
            model_path=model_path,
            experimental_path=experimental_path,
            media_path=media_path,
        )

    assert read_ledger(ledger_path)[0]["status"] == "awaiting_human"
    assert json.loads(dossier_path.read_text(encoding="utf-8"))["human_decision"][
        "decision"
    ] == "pending"


@pytest.mark.parametrize("field,value", [("formula", "C2"), ("charge", -1)])
def test_acceptance_recomputes_live_chemistry_fingerprint(
    tmp_path, field, value
) -> None:
    dossier = _valid_dossier()
    model_path, experimental_path, media_path = _write_live_acceptance_inputs(
        tmp_path, dossier
    )
    model = read_sbml_model(str(model_path))
    setattr(model.metabolites.get_by_id("a_c"), field, value)
    write_sbml_model(model, str(model_path))
    current_model_sha = sha256_file(model_path)
    dossier["model_sha256"] = current_model_sha
    dossier["chemistry_review"]["model_sha256"] = current_model_sha
    dossier["identity_review"]["model_sha256"] = current_model_sha

    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    dossier_path = evidence_dir / f"{dossier['case_id']}.json"
    dossier_path.write_text(json.dumps(dossier), encoding="utf-8")
    row = {field_name: "" for field_name in LEDGER_FIELDS}
    row.update(
        {
            "case_id": dossier["case_id"],
            "status": "awaiting_human",
            "reaction_ids": "R1",
            "model_sha256": current_model_sha,
            "experimental_sha256": dossier["experimental_sha256"],
            "media_sha256": dossier["media_sha256"],
            "target_fingerprint": dossier["target_fingerprint"],
            "chemistry_fingerprint": dossier["chemistry_fingerprint"],
        }
    )
    ledger_path = tmp_path / "ledger.csv"
    write_ledger([row], ledger_path)

    with pytest.raises(ValueError, match="Current chemistry fingerprint"):
        record_human_decision(
            dossier["case_id"],
            "accept",
            ledger_path=ledger_path,
            evidence_dir=evidence_dir,
            model_path=model_path,
            experimental_path=experimental_path,
            media_path=media_path,
        )
