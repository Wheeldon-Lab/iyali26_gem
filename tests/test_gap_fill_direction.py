from pathlib import Path

import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model
from memote.support import consistency

from scripts.gem_annotate import gaps as gaps_module
from scripts.gem_annotate.gap_fill_direction import (
    load_gap_fill_direction_curation,
)
from scripts.gem_annotate.gaps import add_gap_fill_reactions
from scripts.gem_annotate.locus_resolver import LocusCrosswalk


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def toy_gap_fill_resolver(monkeypatch, tmp_path: Path) -> LocusCrosswalk:
    iyli21 = tmp_path / "iyli21_genes_vs_S2.csv"
    iyli21.write_text(
        "model_gene,yali1_s2,yali0,category,n_reactions\n",
        encoding="utf-8",
    )
    metabolic = tmp_path / "s2_metabolic_genes.csv"
    metabolic.write_text(
        "yali1,yali0,geneid,metab_pathways,ec,verdict,in_model\n",
        encoding="utf-8",
    )
    resolver = LocusCrosswalk.from_csvs(iyli21, metabolic)
    monkeypatch.setattr(
        gaps_module, "load_default_locus_crosswalk", lambda: resolver
    )
    return resolver


def _toy_gap_fill_inputs(tmp_path: Path) -> tuple[Model, Path]:
    model = Model("gap-fill-direction-test")
    substrate = Metabolite("substrate", formula="C", compartment="C_cy")
    substrate.annotation = {"metanetx.chemical": "MNXM100"}
    product = Metabolite("product", formula="C", compartment="C_cy")
    product.annotation = {"metanetx.chemical": "MNXM101"}
    model.add_metabolites([substrate, product])

    candidates = tmp_path / "gap_fill_prioritized.csv"
    candidates.write_text(
        "priority,mnxr_id,bigg_reaction,gene_id,equation,ec_number,"
        "kegg_reaction\n"
        "P0,MNXR_TEST,TEST_c,,"
        '"1 MNXM100@MNXD1 = 1 MNXM101@MNXD1",1.1.1.1,R00001\n',
        encoding="utf-8",
    )
    return model, candidates


def _write_direction_table(tmp_path: Path, *, reaction: str = "TEST_c") -> Path:
    table = tmp_path / "gap_fill_direction_curation.csv"
    table.write_text(
        "schema_version,status,bigg_reaction,mnxr_id,stoichiometry_action,"
        "lower_bound,upper_bound,evidence_url,rationale\n"
        f"1,active,{reaction},MNXR_TEST,reverse,0,1000,"
        "https://example.org/reaction,Verified hydrolysis direction.\n",
        encoding="utf-8",
    )
    return table


def test_curated_direction_reverses_equation_and_sets_bounds(
    tmp_path: Path, toy_gap_fill_resolver: LocusCrosswalk
) -> None:
    model, candidates = _toy_gap_fill_inputs(tmp_path)
    curation = _write_direction_table(tmp_path)

    stats = add_gap_fill_reactions(
        model,
        candidates,
        cache_dir=tmp_path / "cache",
        direction_curation_path=curation,
    )

    reaction = model.reactions.get_by_id("TEST_c")
    assert reaction.bounds == (0.0, 1000.0)
    assert reaction.metabolites[model.metabolites.product] == -1.0
    assert reaction.metabolites[model.metabolites.substrate] == 1.0
    assert reaction.notes["gap_fill_direction_status"] == "active"
    assert stats["direction_curated"] == ["TEST_c"]
    assert stats["uncurated_direction"] == []


def test_uncurated_legacy_direction_is_reported(
    tmp_path: Path, toy_gap_fill_resolver: LocusCrosswalk
) -> None:
    model, candidates = _toy_gap_fill_inputs(tmp_path)

    stats = add_gap_fill_reactions(
        model,
        candidates,
        cache_dir=tmp_path / "cache",
    )

    reaction = model.reactions.get_by_id("TEST_c")
    assert reaction.bounds == (-1000.0, 1000.0)
    assert reaction.notes["gap_fill_direction_status"] == "legacy_unreviewed"
    assert stats["direction_curated"] == []
    assert stats["uncurated_direction"] == ["TEST_c"]


def test_curated_existing_reaction_prevents_sphpl_reinsertion(
    tmp_path: Path, toy_gap_fill_resolver: LocusCrosswalk
) -> None:
    model = Model("curated-existing-reaction-test")
    model.add_reactions([Reaction("R730")])
    candidates = tmp_path / "gap_fill_prioritized.csv"
    candidates.write_text(
        "priority,mnxr_id,bigg_reaction,gene_id,equation,ec_number,kegg_reaction\n"
        "P0,MNXR188844,SPHPL,YALI1E33285g,"
        '"1 MNXM1103529@MNXD1 = 1 MNXM187@MNXD1 + 1 MNXM528@MNXD1",4.1.2.27,R02464\n',
        encoding="utf-8",
    )

    stats = add_gap_fill_reactions(model, candidates, cache_dir=tmp_path / "cache")

    assert "SPHPL" not in model.reactions
    assert stats["skipped_curated_existing"] == ["SPHPL"]


def test_curated_existing_reaction_fails_closed_when_target_is_absent(
    tmp_path: Path,
    toy_gap_fill_resolver: LocusCrosswalk,
) -> None:
    candidates = tmp_path / "gap_fill_prioritized.csv"
    candidates.write_text(
        "priority,mnxr_id,bigg_reaction,gene_id,equation,ec_number,kegg_reaction\n"
        "P0,MNXR188844,SPHPL,YALI1E33285g,"
        '"1 MNXM1103529@MNXD1 = 1 MNXM187@MNXD1 + 1 MNXM528@MNXD1",4.1.2.27,R02464\n',
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="curated existing-reaction target is absent"):
        add_gap_fill_reactions(
            Model("missing-curated-target-test"),
            candidates,
            cache_dir=tmp_path / "cache",
        )


def test_stale_direction_row_is_rejected(
    tmp_path: Path, toy_gap_fill_resolver: LocusCrosswalk
) -> None:
    model, candidates = _toy_gap_fill_inputs(tmp_path)
    curation = _write_direction_table(tmp_path, reaction="OTHER_c")

    with pytest.raises(ValueError, match="does not match a P0"):
        add_gap_fill_reactions(
            model,
            candidates,
            cache_dir=tmp_path / "cache",
            direction_curation_path=curation,
        )


@pytest.mark.external_data
@pytest.mark.integration
def test_real_direction_table_contains_verified_egc_roots(
    external_data_file,
) -> None:
    rows = load_gap_fill_direction_curation(
        external_data_file("data/gap_fill_direction_curation.csv")
    )

    assert {
        reaction_id: (row.stoichiometry_action, row.lower_bound, row.upper_bound)
        for reaction_id, row in rows.items()
    } == {
        "R_NTP1": ("reverse", 0.0, 1000.0),
        "R_NDP1": ("reverse", 0.0, 1000.0),
        "R_NTP3pp": ("reverse", 0.0, 1000.0),
        "R_NTP7": ("keep", 0.0, 1000.0),
        "R_PGAM1_PhosHydro": ("reverse", 0.0, 1000.0),
        "R_CAT2p": ("keep", 0.0, 1000.0),
    }


@pytest.mark.integration
def test_built_model_applies_curated_directions_and_has_no_nucleotide_egc() -> None:
    model = read_sbml_model(str(REPO_ROOT / "model.xml"))

    for reaction_id in (
        "R_NTP1",
        "R_NDP1",
        "R_NTP3pp",
        "R_NTP7",
        "R_PGAM1_PhosHydro",
        "R_CAT2p",
    ):
        reaction = model.reactions.get_by_id(reaction_id)
        assert reaction.bounds == (0.0, 1000.0)
        assert reaction.notes["gap_fill_direction_status"] == "active"

    for energy_metabolite in ("MNXM3", "MNXM63", "MNXM51", "MNXM121", "MNXM423"):
        assert consistency.detect_energy_generating_cycles(
            model.copy(), energy_metabolite
        ) == []
