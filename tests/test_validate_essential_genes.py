import json

import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model, write_sbml_model
from cobra.manipulation.delete import knock_out_model_genes

from scripts.gem_annotate.validate_essential_genes import (
    DEFAULT_EXPERIMENTAL,
    DEFAULT_MEDIA,
    IRON_EXCHANGE_ID,
    LEUCINE_EXCHANGE_ID,
    REPO_ROOT,
    apply_media,
    audit_translation_module,
    load_experimental,
    load_media,
    normalize_gene_id,
    run_single_gene_deletions,
    _verify_bypass_candidates,
)
from scripts.gem_annotate.essentiality_evidence import sha256_file
from scripts.gem_annotate import essentiality_evidence, patches as patches_module
from scripts.gem_annotate.patches import apply_curated_essentiality_patches


MODEL_PATH = REPO_ROOT / "model.xml"
SOURCE_XLSX = DEFAULT_EXPERIMENTAL.with_name("42003_2023_4996_MOESM10_ESM.xlsx")


def test_normalize_gene_id() -> None:
    assert normalize_gene_id("YALI1_A00309g") == "YALI1A00309g"
    assert normalize_gene_id("YALI1A00309g") == "YALI1A00309g"


def test_runtime_pseudo_gene_can_be_excluded_from_deletion_screen() -> None:
    model = Model("runtime-pseudo-gene")
    precursor = Metabolite("precursor_c", compartment="c")
    source = Reaction("SOURCE")
    source.bounds = (0.0, 1.0)
    source.add_metabolites({precursor: 1.0})
    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1.0)
    biomass.add_metabolites({precursor: -1.0})
    biomass.gene_reaction_rule = "real_gene or runtime_pseudo_gene"
    model.add_reactions([source, biomass])
    model.objective = biomass

    predictions, growth = run_single_gene_deletions(
        model,
        "glpk",
        excluded_gene_ids={"runtime_pseudo_gene"},
    )

    assert growth == pytest.approx(1.0)
    assert predictions["gene_id"].tolist() == ["real_gene"]


def test_positive_only_source_requires_explicit_flag() -> None:
    with pytest.raises(ValueError, match="positive-only"):
        load_experimental(SOURCE_XLSX, positive_only=False)


def test_reference_counts_and_model_intersection() -> None:
    experimental = load_experimental(DEFAULT_EXPERIMENTAL, positive_only=True)
    model = read_sbml_model(str(MODEL_PATH))
    experimental_ids = set(experimental["gene_id"])
    model_ids = {gene.id for gene in model.genes}

    assert len(experimental) == 1612
    assert experimental["gene_id"].nunique() == 1612
    assert len(experimental_ids & model_ids) == 322
    assert len(experimental_ids - model_ids) == 1290


def test_sd_leu_medium_and_wild_type_growth() -> None:
    model = read_sbml_model(str(MODEL_PATH))
    media = load_media(DEFAULT_MEDIA)
    active = apply_media(model, media)
    model.solver = "glpk"
    growth = model.slim_optimize()

    assert len(media) == 35
    assert LEUCINE_EXCHANGE_ID not in active
    assert active[IRON_EXCHANGE_ID] == pytest.approx(0.0000667)
    assert model.reactions.get_by_id(IRON_EXCHANGE_ID).lower_bound == pytest.approx(
        -0.0000667
    )
    assert 0.1 <= growth <= 2.0


def test_verified_bypass_excludes_globally_essential_dependency() -> None:
    model = Model("gene_specific_bypass")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    c = Metabolite("c_c", compartment="c")

    uptake = Reaction("UPTAKE")
    uptake.bounds = (0.0, 1.5)
    uptake.add_metabolites({a: 1.0})

    global_dependency = Reaction("GLOBAL")
    global_dependency.add_metabolites({a: -1.0, b: 1.0})

    target = Reaction("TARGET")
    target.gene_reaction_rule = "target_gene"
    target.add_metabolites({b: -1.0, c: 1.0})

    bypass = Reaction("ALT")
    bypass.add_metabolites({b: -2.0, c: 1.0})

    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1.0)
    biomass.add_metabolites({c: -1.0})

    model.add_reactions([uptake, global_dependency, target, bypass, biomass])
    model.objective = biomass
    wt_solution = model.optimize()

    verified = _verify_bypass_candidates(
        model,
        model.genes.get_by_id("target_gene"),
        wt_solution.fluxes,
        wt_growth=float(wt_solution.objective_value),
        lethal_growth=0.1,
        max_candidates=10,
    )

    assert verified == ["ALT"]


def test_translation_candidates_are_audited_but_not_connected() -> None:
    model = read_sbml_model(str(MODEL_PATH))
    apply_media(model, load_media(DEFAULT_MEDIA))
    rows, summary = audit_translation_module(model)

    assert len(rows) == 20
    assert rows["pair_complete"].all()
    assert rows["carrier_balanced_in_candidate"].all()
    assert summary["candidate_feasible"] is False
    assert summary["ready_to_connect"] is False


def test_empty_curated_patch_table_does_not_mutate_model(tmp_path) -> None:
    model = read_sbml_model(str(MODEL_PATH))
    before = (len(model.reactions), len(model.genes), model.slim_optimize())

    empty_table = tmp_path / "curated_model_patches.csv"
    empty_table.write_text(
        "patch_id,status,operation,target_id,value,evidence_url,rationale\n"
    )

    applied = apply_curated_essentiality_patches(model, str(empty_table))
    after = (len(model.reactions), len(model.genes), model.slim_optimize())

    assert applied == []
    assert after == before


def test_schema_v2_patch_requires_evidence_gate_fields(tmp_path) -> None:
    model = read_sbml_model(str(MODEL_PATH))
    patch_table = tmp_path / "curated_model_patches.csv"
    patch_table.write_text(
        "patch_id,status,operation,target_id,value,evidence_url,rationale,schema_version\n"
        "EGC-PATCH,accepted,set_gpr,R159,YALI1E11768g,https://example.org/source,test,2\n"
    )
    with pytest.raises(ValueError, match="schema-v2 gate fields"):
        apply_curated_essentiality_patches(model, str(patch_table))


def test_only_existing_patch_may_use_legacy_gate(tmp_path) -> None:
    model = read_sbml_model(str(MODEL_PATH))
    patch_table = tmp_path / "curated_model_patches.csv"
    patch_table.write_text(
        "patch_id,status,operation,target_id,value,evidence_url,rationale,schema_version\n"
        "NEW-PATCH,accepted,set_gpr,R159,YALI1E11768g,https://example.org/source,test,1\n"
    )
    with pytest.raises(ValueError, match="cannot use the legacy"):
        apply_curated_essentiality_patches(model, str(patch_table))


def _write_schema_v2_binding_fixture(
    tmp_path,
    monkeypatch,
    *,
    proposed_target_id: str,
    proposed_value: str,
    live_chemistry_fingerprint: str = "sha256:chemistry-test",
):
    repo_root = tmp_path / "repo"
    evidence_dir = repo_root / "data" / "essentiality" / "evidence"
    evidence_dir.mkdir(parents=True)
    case_id = "EGC-binding-test"
    fingerprint = "sha256:binding-test"
    expected_chemistry_fingerprint = "sha256:chemistry-test"
    approved_at = "2026-07-20T12:00:00-07:00"
    chemistry_dir = evidence_dir / "chemistry"
    chemistry_dir.mkdir()
    chemistry_audit_path = chemistry_dir / "audit.json"
    chemistry_audit_path.write_text('{"ready_for_activation": true}\n')
    chemistry_audit_sha256 = sha256_file(chemistry_audit_path)
    evidence_path = evidence_dir / f"{case_id}.json"
    evidence_path.write_text(
        json.dumps(
            {
                "case_id": case_id,
                "target_fingerprint": fingerprint,
                "chemistry_fingerprint": expected_chemistry_fingerprint,
                "chemistry_review": {
                    "audit_path": "data/essentiality/evidence/chemistry/audit.json",
                    "audit_sha256": chemistry_audit_sha256,
                    "audited_reaction_ids": ["R159"],
                },
                "human_decision": {
                    "approved_by": "human_user",
                    "approved_at": approved_at,
                },
                "proposed_operation": {
                    "operation": "set_gpr",
                    "target_id": proposed_target_id,
                    "value": proposed_value,
                },
                "model_context": {"reactions": [{"reaction_id": "R159"}]},
            }
        )
    )
    patch_table = tmp_path / "curated_model_patches.csv"
    patch_table.write_text(
        "patch_id,status,operation,target_id,value,evidence_url,rationale,"
        "schema_version,case_id,evidence_path,approved_by,approved_at,"
        "target_fingerprint\n"
        "EGC-PATCH,accepted,set_gpr,R159,YALI1E11768g,https://example.org/source,"
        f"test,2,{case_id},data/essentiality/evidence/{case_id}.json,human_user,"
        f"{approved_at},{fingerprint}\n"
    )

    monkeypatch.setattr(
        patches_module,
        "__file__",
        str(repo_root / "scripts" / "gem_annotate" / "patches.py"),
    )
    monkeypatch.setattr(
        essentiality_evidence,
        "require_valid_evidence_dossier",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        essentiality_evidence,
        "read_ledger",
        lambda _path: [
            {
                "case_id": case_id,
                "status": "accepted",
                "target_fingerprint": fingerprint,
                "chemistry_fingerprint": expected_chemistry_fingerprint,
                "approved_by": "human_user",
                "approved_at": approved_at,
            }
        ],
    )
    monkeypatch.setattr(
        essentiality_evidence,
        "target_fingerprint",
        lambda _contexts: fingerprint,
    )
    monkeypatch.setattr(
        essentiality_evidence,
        "chemistry_fingerprint",
        lambda _contexts: live_chemistry_fingerprint,
    )
    return patch_table


def test_schema_v2_patch_rejects_proposal_target_id_mismatch(
    tmp_path, monkeypatch
) -> None:
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R159")
    before = reaction.gene_reaction_rule
    patch_table = _write_schema_v2_binding_fixture(
        tmp_path,
        monkeypatch,
        proposed_target_id="R190",
        proposed_value="YALI1E11768g",
    )

    with pytest.raises(ValueError, match="target_id does not match evidence proposal"):
        apply_curated_essentiality_patches(model, str(patch_table))

    assert reaction.gene_reaction_rule == before


def test_schema_v2_patch_rejects_proposal_value_mismatch(
    tmp_path, monkeypatch
) -> None:
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R159")
    before = reaction.gene_reaction_rule
    patch_table = _write_schema_v2_binding_fixture(
        tmp_path,
        monkeypatch,
        proposed_target_id="R159",
        proposed_value="YALI1C33005g",
    )

    with pytest.raises(ValueError, match="value does not match evidence proposal"):
        apply_curated_essentiality_patches(model, str(patch_table))

    assert reaction.gene_reaction_rule == before


def test_schema_v2_patch_rejects_stale_live_chemistry_before_mutation(
    tmp_path, monkeypatch
) -> None:
    model = read_sbml_model(str(MODEL_PATH))
    reaction = model.reactions.get_by_id("R159")
    before = reaction.gene_reaction_rule
    patch_table = _write_schema_v2_binding_fixture(
        tmp_path,
        monkeypatch,
        proposed_target_id="R159",
        proposed_value="YALI1E11768g",
        live_chemistry_fingerprint="sha256:changed-chemistry",
    )

    with pytest.raises(ValueError, match="chemistry fingerprint changed"):
        apply_curated_essentiality_patches(model, str(patch_table))

    assert reaction.gene_reaction_rule == before


def test_schema_v2_post_patch_gate_rejects_mass_imbalance() -> None:
    model = Model("chemistry-gate")
    substrate = Metabolite("a_c", compartment="c", formula="C", charge=0)
    product = Metabolite("b_c", compartment="c", formula="CH", charge=0)
    reaction = Reaction("R1")
    reaction.add_metabolites({substrate: -1, product: 1})
    model.add_reactions([reaction])

    with pytest.raises(ValueError, match="post-patch mass/charge gate"):
        patches_module._assert_schema_v2_post_patch_balance(
            model, "R1", "EGC-PATCH"
        )


def test_current_input_sha256_values_are_stable() -> None:
    assert sha256_file(MODEL_PATH) == (
        "bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee"
    )
    assert sha256_file(DEFAULT_EXPERIMENTAL) == (
        "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b"
    )
    assert sha256_file(DEFAULT_MEDIA) == (
        "ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85"
    )


def test_cpa_ura2_partition_is_idempotent_and_growth_safe(tmp_path) -> None:
    model = read_sbml_model(str(MODEL_PATH))
    before_counts = (len(model.reactions), len(model.metabolites), len(model.genes))
    audited_reactions = ("R159", "R190", "R607")
    before_unbalanced = {
        reaction_id: bool(model.reactions.get_by_id(reaction_id).check_mass_balance())
        for reaction_id in audited_reactions
    }

    first = apply_curated_essentiality_patches(model)
    second = apply_curated_essentiality_patches(model)

    assert [row["patch_id"] for row in first] == ["EG-GPR-001"]
    assert [row["patch_id"] for row in second] == ["EG-GPR-001"]
    assert (len(model.reactions), len(model.metabolites), len(model.genes)) == before_counts
    assert {
        reaction_id: bool(model.reactions.get_by_id(reaction_id).check_mass_balance())
        for reaction_id in audited_reactions
    } == before_unbalanced

    r159 = model.reactions.get_by_id("R159")
    r190 = model.reactions.get_by_id("R190")
    r607 = model.reactions.get_by_id("R607")
    carbamoyl_phosphate = model.metabolites.get_by_id("m325[C_cy]")
    assert r159.gene_reaction_rule == "YALI1E11768g"
    assert r190.gene_reaction_rule == "YALI1C33005g and YALI1D09420g"
    assert r159.subsystem == "Pyrimidine metabolism"
    assert r190.subsystem == "Arginine and proline metabolism"
    assert carbamoyl_phosphate not in r159.metabolites
    assert carbamoyl_phosphate in r190.metabolites
    assert carbamoyl_phosphate in r607.metabolites

    # SBML stores subsystems through Groups; verify the partition survives a
    # complete pipeline write/read cycle instead of reverting to the old group.
    roundtrip_path = tmp_path / "cpa_ura2_roundtrip.xml"
    write_sbml_model(model, str(roundtrip_path))
    roundtrip = read_sbml_model(str(roundtrip_path))
    assert roundtrip.reactions.get_by_id("R159").subsystem == "Pyrimidine metabolism"
    assert (
        roundtrip.reactions.get_by_id("R190").subsystem
        == "Arginine and proline metabolism"
    )

    apply_media(model, load_media(DEFAULT_MEDIA))
    model.solver = "glpk"
    wt_growth = float(model.slim_optimize())
    assert wt_growth == pytest.approx(1.4492618988553403, rel=1e-9)

    expected_ratios = {
        "YALI1C33005g": 0.14210846955691958,
        "YALI1D09420g": 0.14210846955691955,
        "YALI1E11768g": 0.08072394846415512,
    }
    for gene_id, expected_ratio in expected_ratios.items():
        with model:
            knock_out_model_genes(model, [gene_id])
            ratio = float(model.slim_optimize()) / wt_growth
        assert ratio == pytest.approx(expected_ratio, rel=1e-8)
