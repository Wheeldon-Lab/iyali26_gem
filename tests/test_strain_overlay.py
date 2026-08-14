import copy
import json

import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import write_sbml_model

from scripts.gem_annotate.strain_overlay import (
    CONCENTRATION_RATIO_STATUS,
    LEU2_EXPECTED_GPR,
    LEU2_PLASMID_PSEUDO_GENE,
    RUNTIME_OVERRIDE_STATUS,
    StrainProfileError,
    apply_strain_overlay,
    derive_concentration_ratio_uptake,
    load_strain_profile,
)
from scripts.gem_annotate.validate_essential_genes import validate_essential_genes


def _profile(*, runtime_uracil_uptake: float = 0.2) -> dict:
    return {
        "schema_version": "1.0",
        "profile_id": "po1f_sd_leu_test",
        "strain": {
            "name": "PO1f",
            "reference_background": "W29/CLIB89",
            "genotype": "MatA, leu2-270, ura3-302, xpr2-322, axp-2",
            "assay_background": "NHEJ-intact acCRISPR fitness screen",
        },
        "operations": [
            {
                "type": "disable_reaction",
                "reaction_id": "R612",
                "expected_before": {
                    "lower_bound": 0.0,
                    "upper_bound": 1000.0,
                },
                "set_bounds": {
                    "lower_bound": 0.0,
                    "upper_bound": 0.0,
                },
                "locus": "URA3",
                "gene_id": "YALI1E31685g",
                "legacy_gene_id": "YALI0E26741g",
                "allele": "ura3-302",
                "protein_function": ("orotidine-5'-phosphate decarboxylase"),
                "evidence_status": "PO1f_ura3_genotype",
            },
            {
                "type": "plasmid_complement",
                "reaction_id": "R45",
                "expected_before": {
                    "gene_reaction_rule": LEU2_EXPECTED_GPR,
                },
                "set_gene_reaction_rule": LEU2_PLASMID_PSEUDO_GENE,
                "pseudo_gene": LEU2_PLASMID_PSEUDO_GENE,
                "locus": "LEU2",
                "gene_id": "YALI1C00464g",
                "legacy_gene_id": "YALI0C00407g",
                "allele": "leu2-270",
                "complement_source": "pCas9yl-GW/pLbCas12ayl-GW LEU2 marker",
                "protein_function": ("3-isopropylmalate dehydrogenase"),
                "evidence_status": "PO1f_LEU2_plasmid_complementation",
            },
        ],
        "medium": {
            "uptake_assertions": {
                "R1219": 0.0,
                "R1354": 0.01607,
            },
            "runtime_uptake_overrides": {
                "R1354": runtime_uracil_uptake,
            },
            "formulation": {
                "uracil_mg_per_l": 20.0,
                "uracil_millimolar": 0.17843,
                "uracil_molecular_weight_g_per_mol": 112.09,
                "glucose_g_per_l": 20.0,
                "glucose_molecular_weight_g_per_mol": 180.156,
                "glucose_uptake": 10.0,
                "supply_ratio_surrogate_uptake": 0.01607,
                "supply_ratio_surrogate_status": (CONCENTRATION_RATIO_STATUS),
                "runtime_override_status": RUNTIME_OVERRIDE_STATUS,
            },
        },
        "assay_confounded_loci": ["URA3", "LEU2"],
        "provenance_only_variants": ["PO1f genotype aliases"],
        "sources": ["https://example.test/po1f-profile-source"],
    }


def _toy_model() -> Model:
    model = Model("po1f-overlay-test")

    omp = Metabolite("omp_c", compartment="c")
    ump = Metabolite("ump_c", compartment="c")
    ura3 = Reaction("R612")
    ura3.bounds = (0.0, 1000.0)
    ura3.add_metabolites({omp: -1.0, ump: 1.0})
    ura3.gene_reaction_rule = "YALI1E31685g"

    precursor = Metabolite("leu_precursor_c", compartment="c")
    product = Metabolite("leu_product_c", compartment="c")
    leu2 = Reaction("R45")
    leu2.bounds = (0.0, 1000.0)
    leu2.add_metabolites({precursor: -1.0, product: 1.0})
    leu2.gene_reaction_rule = LEU2_EXPECTED_GPR

    leucine_external = Metabolite("leucine_e", compartment="e")
    leucine_exchange = Reaction("R1219")
    leucine_exchange.bounds = (0.0, 1000.0)
    leucine_exchange.add_metabolites({leucine_external: -1.0})

    uracil_external = Metabolite("uracil_e", compartment="e")
    uracil_exchange = Reaction("R1354")
    uracil_exchange.bounds = (-0.01607, 1000.0)
    uracil_exchange.add_metabolites({uracil_external: -1.0})

    glucose_external = Metabolite("glucose_e", compartment="e")
    glucose_exchange = Reaction("R1070")
    glucose_exchange.bounds = (-10.0, 1000.0)
    glucose_exchange.add_metabolites({glucose_external: -1.0})

    uracil_cytosolic = Metabolite("uracil_c", compartment="c")
    uracil_transport = Reaction("R935")
    uracil_transport.bounds = (0.0, 1000.0)
    uracil_transport.add_metabolites(
        {uracil_external: -1.0, uracil_cytosolic: 1.0}
    )
    uracil_transport.gene_reaction_rule = "transport_gene"

    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1.0)
    biomass.add_metabolites({uracil_cytosolic: -1.0})

    model.add_reactions(
        [
            ura3,
            leu2,
            leucine_exchange,
            uracil_exchange,
            glucose_exchange,
            uracil_transport,
            biomass,
        ]
    )
    model.objective = biomass
    return model


def test_load_profile_is_strict_and_rejects_duplicate_or_unknown_keys(
    tmp_path,
) -> None:
    profile_path = tmp_path / "profile.json"
    profile_path.write_text(
        json.dumps(_profile(), indent=2),
        encoding="utf-8",
    )
    loaded = load_strain_profile(profile_path)
    assert loaded == _profile()

    unknown = _profile()
    unknown["unexpected"] = True
    profile_path.write_text(json.dumps(unknown), encoding="utf-8")
    with pytest.raises(StrainProfileError, match="unknown keys"):
        load_strain_profile(profile_path)

    profile_path.write_text(
        '{"schema_version":"1.0","schema_version":"1.0"}',
        encoding="utf-8",
    )
    with pytest.raises(StrainProfileError, match="Duplicate JSON key"):
        load_strain_profile(profile_path)


def test_apply_overlay_changes_copy_and_preserves_baseline() -> None:
    baseline = _toy_model()
    baseline_signature = {
        "R612": baseline.reactions.get_by_id("R612").bounds,
        "R45": baseline.reactions.get_by_id("R45").gene_reaction_rule,
        "medium": dict(baseline.medium),
        "genes": sorted(gene.id for gene in baseline.genes),
    }
    overlay_model = baseline.copy()
    active_medium = dict(overlay_model.medium)

    audit = apply_strain_overlay(
        overlay_model,
        _profile(),
        active_medium=active_medium,
    )

    assert overlay_model.reactions.get_by_id("R612").bounds == (0.0, 0.0)
    assert (
        overlay_model.reactions.get_by_id("R45").gene_reaction_rule
        == LEU2_PLASMID_PSEUDO_GENE
    )
    pseudo_gene = overlay_model.genes.get_by_id(LEU2_PLASMID_PSEUDO_GENE)
    assert pseudo_gene.name == "LEU2 plasmid complementation"
    assert pseudo_gene.notes["complemented_locus"] == "LEU2"
    assert overlay_model.medium["R1354"] == pytest.approx(0.2)
    assert overlay_model.medium.get("R1219", 0.0) == 0.0

    assert baseline.reactions.get_by_id("R612").bounds == baseline_signature["R612"]
    assert (
        baseline.reactions.get_by_id("R45").gene_reaction_rule
        == baseline_signature["R45"]
    )
    assert dict(baseline.medium) == baseline_signature["medium"]
    assert sorted(gene.id for gene in baseline.genes) == baseline_signature["genes"]

    assert audit["medium"]["active_medium"]["R1354"] == pytest.approx(0.2)
    assert audit["medium"]["active_medium"]["R1219"] == 0.0
    assert audit["medium"]["formulation"]["uracil_mg_per_l"] == 20.0
    assert audit["model_layer_modified_in_memory"] is True
    assert audit["canonical_sbml_written"] is False
    json.dumps(audit, allow_nan=False)


def test_overlay_checks_expected_model_and_medium_state() -> None:
    wrong_bounds = _toy_model()
    wrong_bounds.reactions.get_by_id("R612").upper_bound = 1.0
    with pytest.raises(ValueError, match="expected.*before PO1f overlay"):
        apply_strain_overlay(wrong_bounds, _profile())

    wrong_medium = _toy_model()
    wrong_medium.reactions.get_by_id("R1219").lower_bound = -1.0
    with pytest.raises(ValueError, match="Medium assertion failed for R1219"):
        apply_strain_overlay(wrong_medium, _profile())

    mismatched_report = _toy_model()
    reported_medium = dict(mismatched_report.medium)
    reported_medium["R1354"] = 0.5
    with pytest.raises(ValueError, match="active_medium reports R1354"):
        apply_strain_overlay(
            mismatched_report,
            _profile(),
            active_medium=reported_medium,
        )


def test_concentration_ratio_formula_is_not_runtime_override() -> None:
    surrogate = derive_concentration_ratio_uptake(
        20.0,
        112.09,
        20.0,
        10.0,
        glucose_molecular_weight_g_per_mol=180.156,
    )
    assert surrogate == pytest.approx(0.01607, rel=5e-4)

    profile = _profile(runtime_uracil_uptake=1000.0)
    model = _toy_model()
    audit = apply_strain_overlay(model, profile)

    assert model.medium["R1354"] == pytest.approx(1000.0)
    assert audit["medium"]["runtime_uptake_overrides"]["R1354"] == 1000.0
    formulation = audit["medium"]["formulation"]
    assert formulation["supply_ratio_surrogate_uptake"] == 0.01607
    assert formulation["supply_ratio_surrogate_status"] == (
        "concentration_ratio_surrogate_not_measured_vmax"
    )
    assert formulation["runtime_override_status"] == (
        "static_fba_nonlimiting_not_measured_vmax"
    )


def test_profile_validation_does_not_mutate_input() -> None:
    profile = _profile()
    snapshot = copy.deepcopy(profile)

    apply_strain_overlay(_toy_model(), profile)

    assert profile == snapshot


def test_overlay_reapplication_is_idempotent() -> None:
    model = _toy_model()
    first = apply_strain_overlay(model, _profile())
    second = apply_strain_overlay(
        model,
        _profile(),
        active_medium=dict(model.medium),
    )

    assert first["overlay_already_applied"] is False
    assert second["overlay_already_applied"] is True
    assert second["model_layer_modified_in_memory"] is False
    assert all(not operation["changed"] for operation in second["operations"])
    assert model.reactions.get_by_id("R612").bounds == (0.0, 0.0)
    assert model.medium["R1354"] == pytest.approx(0.2)


def test_validator_records_profile_and_excludes_runtime_pseudo_gene(
    tmp_path,
) -> None:
    model_path = tmp_path / "model.xml"
    experimental_path = tmp_path / "experimental.csv"
    medium_path = tmp_path / "medium.csv"
    profile_path = tmp_path / "po1f.json"
    output_dir = tmp_path / "result"
    write_sbml_model(_toy_model(), str(model_path))
    experimental_path.write_text(
        "gene_id,function\n"
        "transport_gene,uracil transport fixture\n",
        encoding="utf-8",
    )
    medium_path.write_text(
        "exchange,uptake,comment\n"
        "R1354,0.01607,legacy uracil surrogate\n"
        "R1219,0,SD-Leu invariant\n",
        encoding="utf-8",
    )
    profile_path.write_text(
        json.dumps(_profile(runtime_uracil_uptake=1000.0), indent=2),
        encoding="utf-8",
    )

    summary = validate_essential_genes(
        experimental_path=experimental_path,
        model_path=model_path,
        media_path=medium_path,
        output_dir=output_dir,
        positive_only=True,
        solver="glpk",
        strain_profile_path=profile_path,
    )

    assert summary["wt_growth"] == pytest.approx(1.0)
    assert summary["primary"]["TP"] == 1
    assert summary["medium"]["uracil_uptake_bound"] == 1000.0
    assert summary["excluded_runtime_genes"] == [LEU2_PLASMID_PSEUDO_GENE]
    assert summary["strain_overlay"]["profile_sha256"]
    assert summary["simulation_context"]["fingerprint"]
    assert summary["simulation_context"]["simulation_context_fingerprint_version"] == "1"
    assert summary["simulation_context"]["strain_overlay_enabled"] is True
    assert summary["simulation_context"]["strain_overlay_effect_sha256"]
    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["inputs"]["strain_profile"]["path"] == str(
        profile_path.resolve()
    )
    assert manifest["strain_overlay"]["profile_id"] == "po1f_sd_leu_test"
    assert manifest["simulation_context"]["simulation_context_fingerprint"] == summary[
        "simulation_context"
    ]["simulation_context_fingerprint"]
