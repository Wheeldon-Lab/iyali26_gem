import csv
from pathlib import Path

import pytest
from cobra import Metabolite, Model, Reaction
from cobra.core import Group
from cobra.io import read_sbml_model, write_sbml_model
from cobra.manipulation.delete import knock_out_model_genes

from scripts.gem_annotate.essentiality_evidence import sha256_file, target_fingerprint
from scripts.gem_annotate.main import main as build_model
from scripts.gem_annotate.provisional_capacity import (
    MARKER_STATUS,
    apply_provisional_isozyme_capacities,
)
from scripts.gem_annotate.validate_essential_genes import REPO_ROOT


PROFILE_PATH = (
    REPO_ROOT / "data" / "essentiality" / "provisional_isozyme_capacities.csv"
)


def _toy_model() -> Model:
    model = Model("toy_provisional_capacity")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")

    uptake = Reaction("UPTAKE")
    uptake.bounds = (0.0, 10.0)
    uptake.add_metabolites({a: 1.0})

    isozyme = Reaction("R_ISO")
    isozyme.bounds = (0.0, 1000.0)
    isozyme.add_metabolites({a: -1.0, b: 1.0})
    isozyme.gene_reaction_rule = "g_backup or g_main"

    biomass = Reaction("BIOMASS")
    biomass.bounds = (0.0, 1000.0)
    biomass.add_metabolites({b: -1.0})

    model.add_reactions([uptake, isozyme, biomass])
    model.objective = biomass
    pathway = Group("toy_pathway", members=[isozyme])
    model.add_groups([pathway])
    return model


def _fingerprint(reaction: Reaction) -> str:
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


def _write_toy_profile(
    path: Path,
    model: Model,
    fingerprint: str | None = None,
    expected_gpr: str = "g_backup or g_main",
) -> None:
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
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(
            {
                "capacity_id": "PCAP-TOY",
                "status": "active_exploratory",
                "source_reaction_id": "R_ISO",
                "expected_gpr": expected_gpr,
                "primary_gpr": "g_main",
                "backup_gpr": "g_backup",
                "backup_reaction_id": "R_ISO__PCAP_BACKUP",
                "provisional_upper_bound": "0.5",
                "units": "mmol_gDW_h",
                "parameter_basis": "test_not_measured",
                "case_id": "EGC-toy",
                "validated_model_sha256": "toy-sha",
                "target_fingerprint": fingerprint
                or _fingerprint(model.reactions.get_by_id("R_ISO")),
                "requires_protein_abundance": "true",
                "requires_kcat": "true",
                "replacement_formula": "kcat*E*3600",
                "rationale": "test",
            }
        )


def test_global_split_preserves_total_capacity_and_is_idempotent(
    tmp_path: Path,
) -> None:
    model = _toy_model()
    profile = tmp_path / "profile.csv"
    _write_toy_profile(profile, model)
    baseline_growth = float(model.slim_optimize())

    first = apply_provisional_isozyme_capacities(
        model, profile, reference_model_sha256="toy-sha"
    )
    second = apply_provisional_isozyme_capacities(
        model, profile, reference_model_sha256="toy-sha"
    )

    primary = model.reactions.get_by_id("R_ISO")
    backup = model.reactions.get_by_id("R_ISO__PCAP_BACKUP")
    assert first[0]["outcome"] == "applied"
    assert second[0]["outcome"] == "already_applied"
    assert primary.gene_reaction_rule == "g_main"
    assert backup.gene_reaction_rule == "g_backup"
    assert primary.upper_bound + backup.upper_bound == pytest.approx(1000.0)
    assert len(model.reactions) == 4
    assert float(model.slim_optimize()) == pytest.approx(baseline_growth)

    with model:
        knock_out_model_genes(model, ["g_main"])
        assert float(model.slim_optimize()) / baseline_growth == pytest.approx(0.05)
    with model:
        knock_out_model_genes(model, ["g_backup"])
        assert float(model.slim_optimize()) / baseline_growth == pytest.approx(1.0)

    assert model.genes.get_by_id("g_main").notes["provisional_capacity_status"] == MARKER_STATUS
    assert model.genes.get_by_id("g_backup").notes["provisional_capacity_status"] == MARKER_STATUS
    assert (
        model.genes.get_by_id("g_backup").notes[
            "provisional_capacity_requires_protein_abundance"
        ]
        == "true"
    )
    assert (
        model.genes.get_by_id("g_backup").notes[
            "provisional_capacity_requires_kcat"
        ]
        == "true"
    )
    assert (
        model.genes.get_by_id("g_backup").notes[
            "provisional_capacity_reference_model_sha256"
        ]
        == "toy-sha"
    )
    assert (
        model.genes.get_by_id("g_backup").notes[
            "provisional_capacity_target_fingerprints"
        ]
        == _fingerprint(_toy_model().reactions.get_by_id("R_ISO"))
    )
    assert backup in model.groups.get_by_id("toy_pathway").members


def test_provisional_markers_survive_sbml_roundtrip(tmp_path: Path) -> None:
    model = _toy_model()
    profile = tmp_path / "profile.csv"
    _write_toy_profile(profile, model)
    apply_provisional_isozyme_capacities(
        model, profile, reference_model_sha256="toy-sha"
    )
    output = tmp_path / "capacity.xml"
    write_sbml_model(model, str(output))
    roundtrip = read_sbml_model(str(output))

    assert (
        roundtrip.genes.get_by_id("g_backup").notes["provisional_capacity_status"]
        == MARKER_STATUS
    )
    assert (
        roundtrip.reactions.get_by_id("R_ISO__PCAP_BACKUP").notes[
            "provisional_capacity_role"
        ]
        == "capacity_limited_backup_pool"
    )
    assert (
        roundtrip.reactions.get_by_id("R_ISO__PCAP_BACKUP")
        in roundtrip.groups.get_by_id("toy_pathway").members
    )
    assert roundtrip.notes["provisional_capacity_reference_model_sha256"] == "toy-sha"
    assert roundtrip.notes["provisional_capacity_profile_sha256"] == sha256_file(
        profile
    )


def test_stale_profile_fails_without_mutating_model(tmp_path: Path) -> None:
    model = _toy_model()
    profile = tmp_path / "stale.csv"
    _write_toy_profile(profile, model, fingerprint="sha256:stale")
    before = (
        len(model.reactions),
        model.reactions.get_by_id("R_ISO").gene_reaction_rule,
        model.reactions.get_by_id("R_ISO").bounds,
    )
    with pytest.raises(ValueError, match="fingerprint is stale"):
        apply_provisional_isozyme_capacities(
            model, profile, reference_model_sha256="toy-sha"
        )
    after = (
        len(model.reactions),
        model.reactions.get_by_id("R_ISO").gene_reaction_rule,
        model.reactions.get_by_id("R_ISO").bounds,
    )
    assert after == before


def test_reference_model_sha_is_a_gate(tmp_path: Path) -> None:
    model = _toy_model()
    profile = tmp_path / "profile.csv"
    _write_toy_profile(profile, model)
    before = (
        len(model.reactions),
        model.reactions.get_by_id("R_ISO").gene_reaction_rule,
        model.reactions.get_by_id("R_ISO").bounds,
    )

    with pytest.raises(ValueError, match="calibrated against model SHA"):
        apply_provisional_isozyme_capacities(
            model, profile, reference_model_sha256="different-sha"
        )

    assert before == (
        len(model.reactions),
        model.reactions.get_by_id("R_ISO").gene_reaction_rule,
        model.reactions.get_by_id("R_ISO").bounds,
    )


def test_already_split_model_rechecks_original_target_fingerprint(
    tmp_path: Path,
) -> None:
    model = _toy_model()
    profile = tmp_path / "profile.csv"
    _write_toy_profile(profile, model)
    apply_provisional_isozyme_capacities(
        model, profile, reference_model_sha256="toy-sha"
    )
    _write_toy_profile(profile, _toy_model(), fingerprint="sha256:tampered")

    with pytest.raises(ValueError, match="backup ID already exists"):
        apply_provisional_isozyme_capacities(
            model, profile, reference_model_sha256="toy-sha"
        )


def test_already_split_model_rechecks_original_groups(tmp_path: Path) -> None:
    model = _toy_model()
    profile = tmp_path / "profile.csv"
    _write_toy_profile(profile, model)
    apply_provisional_isozyme_capacities(
        model, profile, reference_model_sha256="toy-sha"
    )
    group = model.groups.get_by_id("toy_pathway")
    group.remove_members(
        [
            model.reactions.get_by_id("R_ISO"),
            model.reactions.get_by_id("R_ISO__PCAP_BACKUP"),
        ]
    )

    with pytest.raises(ValueError, match="backup ID already exists"):
        apply_provisional_isozyme_capacities(
            model, profile, reference_model_sha256="toy-sha"
        )


def test_nonzero_lower_bound_is_rejected_before_mutation(tmp_path: Path) -> None:
    model = _toy_model()
    source = model.reactions.get_by_id("R_ISO")
    source.lower_bound = 0.1
    profile = tmp_path / "profile.csv"
    _write_toy_profile(profile, model)

    with pytest.raises(ValueError, match="requires a zero lower bound"):
        apply_provisional_isozyme_capacities(
            model, profile, reference_model_sha256="toy-sha"
        )

    assert "R_ISO__PCAP_BACKUP" not in model.reactions
    assert source.bounds == (0.1, 1000.0)


def test_complex_gpr_is_rejected_by_parsed_logic(tmp_path: Path) -> None:
    model = _toy_model()
    source = model.reactions.get_by_id("R_ISO")
    source.gene_reaction_rule = "g_backup and(g_main)"
    profile = tmp_path / "profile.csv"
    _write_toy_profile(
        profile,
        model,
        expected_gpr=source.gene_reaction_rule,
    )

    with pytest.raises(ValueError, match="only simple OR isozyme partitions"):
        apply_provisional_isozyme_capacities(
            model, profile, reference_model_sha256="toy-sha"
        )

    assert "R_ISO__PCAP_BACKUP" not in model.reactions
    assert source.gene_reaction_rule == "g_backup and g_main"


def test_python_api_cannot_overwrite_canonical_model() -> None:
    with pytest.raises(ValueError, match="cannot overwrite canonical model.xml"):
        build_model(
            provisional_capacity_path=PROFILE_PATH,
            output_model_path=REPO_ROOT / "model.xml",
        )


def test_trna_biomass_overlay_cannot_overwrite_canonical_model() -> None:
    with pytest.raises(ValueError, match="cannot overwrite canonical model.xml"):
        build_model(
            trna_biomass_mode="split",
            output_model_path=REPO_ROOT / "model.xml",
        )


def test_current_profile_is_rejected_after_reference_model_changes() -> None:
    model = read_sbml_model(str(REPO_ROOT / "model.xml"))
    reaction_count = len(model.reactions)
    original_gprs = {
        reaction_id: model.reactions.get_by_id(reaction_id).gene_reaction_rule
        for reaction_id in ("R4", "R1846")
    }

    with pytest.raises(ValueError, match="calibrated against model SHA"):
        apply_provisional_isozyme_capacities(
            model,
            PROFILE_PATH,
            reference_model_sha256=sha256_file(REPO_ROOT / "model.xml"),
        )

    assert len(model.reactions) == reaction_count
    assert "R4__PCAP_BACKUP" not in model.reactions
    assert "R1846__PCAP_BACKUP" not in model.reactions
    assert {
        reaction_id: model.reactions.get_by_id(reaction_id).gene_reaction_rule
        for reaction_id in ("R4", "R1846")
    } == original_gprs
