import json
from pathlib import Path

import pytest
from cobra import Metabolite, Model, Reaction

from scripts.gem_annotate import genes
from scripts.gem_annotate.gaps import _load_gene_cache, add_gap_fill_reactions
from scripts.gem_annotate.locus_resolver import (
    AmbiguousLocusMappingError,
    LocusCrosswalk,
    canonical_locus_key,
    load_default_locus_crosswalk,
    locus_spelling_variants,
)


def _write_crosswalks(
    tmp_path: Path, iyli21_rows: str, metabolic_rows: str
) -> LocusCrosswalk:
    iyli21 = tmp_path / "iyli21_genes_vs_S2.csv"
    iyli21.write_text(
        "model_gene,yali1_s2,yali0,category,n_reactions\n" + iyli21_rows,
        encoding="utf-8",
    )
    metabolic = tmp_path / "s2_metabolic_genes.csv"
    metabolic.write_text(
        "yali1,yali0,geneid,metab_pathways,ec,verdict,in_model\n"
        + metabolic_rows,
        encoding="utf-8",
    )
    return LocusCrosswalk.from_csvs(iyli21, metabolic)


def _model_with_gene(gene_id: str) -> Model:
    model = Model("locus-test")
    anchor = Reaction("ANCHOR")
    model.add_reactions([anchor])
    anchor.gene_reaction_rule = gene_id
    return model


def test_spelling_normalisation_preserves_assembly() -> None:
    assert canonical_locus_key("yali1_C32184G") == "yali1c32184g"
    assert locus_spelling_variants("YALI1C32184g") == {
        "YALI1C32184g",
        "YALI1_C32184g",
    }
    assert "YALI0C32184g" not in locus_spelling_variants("YALI1C32184g")
    assert genes._normalise_locus_tag("YALI0_C23364g") == {
        "YALI0C23364g",
        "YALI0_C23364g",
    }


def test_explicit_crosswalk_resolves_real_alias_but_rejects_same_suffix_alias() -> None:
    resolver = load_default_locus_crosswalk()
    lookup = resolver.build_lookup(["YALI1C32184g"])

    assert lookup["yali1c32184g"] == "YALI1C32184g"
    assert lookup["yali0c23364g"] == "YALI1C32184g"
    assert "yali0c32184g" not in lookup
    assert resolver.counterpart("YALI0C23364g") == "yali1c32184g"


def test_crosswalk_refuses_entire_ambiguous_component(tmp_path: Path) -> None:
    resolver = _write_crosswalks(
        tmp_path,
        "YALI1A00001g,YALI1_A00001g,YALI0A00009g,has_YALI0,1\n",
        "YALI1_A00002g,YALI0A00009g,1,,1.1.1.1,metabolic,no\n",
    )

    with pytest.raises(AmbiguousLocusMappingError, match="Ambiguous"):
        resolver.counterpart("YALI0A00009g")
    with pytest.raises(AmbiguousLocusMappingError, match="Ambiguous"):
        resolver.counterpart("YALI1A00001g")

    lookup = resolver.build_lookup(["YALI1A00001g", "YALI1A00002g"])
    assert "yali0a00009g" not in lookup


def test_r289_cross_assembly_mappings_are_explicit_and_one_to_one() -> None:
    resolver = load_default_locus_crosswalk()
    r289_pairs = {
        "YALI0C23364g": "YALI1C32184g",
        "YALI0E05929g": "YALI1E07121g",
        "YALI0E15081g": "YALI1E18121g",
        "YALI0D10549g": "YALI1D13201g",
    }
    lookup = resolver.build_lookup(r289_pairs.values())

    for yali0, yali1 in r289_pairs.items():
        assert lookup[canonical_locus_key(yali0)] == yali1
        assert resolver.counterpart(yali0) == canonical_locus_key(yali1)


def test_evidence_backed_competing_orf_crosswalk_is_excluded() -> None:
    resolver = load_default_locus_crosswalk()
    excluded_pair = (
        canonical_locus_key("YALI1_E18171g"),
        canonical_locus_key("YALI0E15125g"),
    )

    assert excluded_pair in resolver.excluded_pairs
    assert resolver.counterpart("YALI1_E18171g") is None
    assert resolver.counterpart("YALI0E15125g") is None
    lookup = resolver.build_lookup(["YALI1E18171g"])
    assert canonical_locus_key("YALI0E15125g") not in lookup


def test_tier_a_uses_crosswalk_but_not_fabricated_same_suffix(monkeypatch) -> None:
    entries = [
        {
            "primaryAccession": "GOOD",
            "genes": [
                {"orderedLocusNames": [{"value": "YALI0C23364g"}]}
            ],
        },
        {
            "primaryAccession": "BAD",
            "genes": [
                {"orderedLocusNames": [{"value": "YALI0C32184g"}]}
            ],
        },
    ]
    monkeypatch.setattr(genes, "_PROTEOME_IDS", ("toy",))
    monkeypatch.setattr(genes, "_fetch_proteome", lambda _proteome_id: entries)

    mapping = genes._tier_a(
        ["YALI1C32184g"], resolver=load_default_locus_crosswalk()
    )

    assert mapping == {"YALI1C32184g": {"uniprot": ["GOOD"]}}


def test_tier_b_queries_only_exact_safe_spellings(monkeypatch, tmp_path: Path) -> None:
    queries: list[str] = []

    class Response:
        def __init__(self, results: list[dict]):
            self._results = results

        def raise_for_status(self) -> None:
            return None

        def json(self) -> dict:
            return {"results": self._results}

    good_entry = {
        "primaryAccession": "GOOD",
        "genes": [{"orderedLocusNames": [{"value": "YALI0C23364g"}]}],
    }

    def request(_method: str, _url: str, *, params: dict, timeout: int) -> Response:
        del timeout
        query = params["query"]
        queries.append(query)
        return Response([good_entry] if "YALI0C23364g" in query else [])

    monkeypatch.setattr(genes, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(genes, "_request_with_retry", request)

    mapping = genes._tier_b(
        ["YALI1C32184g"], resolver=load_default_locus_crosswalk()
    )

    assert mapping == {"YALI1C32184g": {"uniprot": ["GOOD"]}}
    assert queries
    assert all(query.startswith('gene_exact:"') for query in queries)
    assert all('gene:"' not in query for query in queries)
    assert not any("YALI0C32184g" in query for query in queries)


def test_ncbi_tier_rejects_fabricated_same_suffix_alias(
    monkeypatch, tmp_path: Path
) -> None:
    class Response:
        def __init__(self, *, json_data: dict | None = None, content: bytes = b""):
            self._json_data = json_data or {}
            self.content = content

        def raise_for_status(self) -> None:
            return None

        def json(self) -> dict:
            return self._json_data

    xml = b"""<Entrezgene-Set><Entrezgene>
      <Gene-track><Gene-track_geneid>999</Gene-track_geneid></Gene-track>
      <Gene-ref><Gene-ref_locus-tag>YALI0C32184g</Gene-ref_locus-tag></Gene-ref>
    </Entrezgene></Entrezgene-Set>"""

    def request(_method: str, url: str, **_kwargs) -> Response:
        if "esearch" in url:
            return Response(json_data={"esearchresult": {"idlist": ["999"]}})
        return Response(content=xml)

    monkeypatch.setattr(genes, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(genes, "_request_with_retry", request)

    mapping = genes._tier_ncbi(
        ["YALI1C32184g"], resolver=load_default_locus_crosswalk()
    )

    assert mapping == {}


def test_gap_gene_cache_rejects_legacy_false_alias(tmp_path: Path) -> None:
    model = _model_with_gene("YALI1C32184g")
    cache_path = tmp_path / "gene_locus_tag_map.json"
    cache_path.write_text(
        json.dumps({"yali0c32184g": "YALI1C32184g"}), encoding="utf-8"
    )

    lookup = _load_gene_cache(cache_path, model)

    assert lookup["yali0c23364g"] == "YALI1C32184g"
    assert "yali0c32184g" not in lookup
    cached = json.loads(cache_path.read_text(encoding="utf-8"))
    assert cached["_meta"]["schema"] == "safe-locus-cache-v1"


def test_gap_fill_skips_unresolved_cross_assembly_gene_without_creating_it(
    tmp_path: Path,
) -> None:
    model = _model_with_gene("YALI1C32184g")
    substrate = Metabolite("substrate", formula="C", compartment="C_cy")
    substrate.annotation = {"metanetx.chemical": "MNXM100"}
    product = Metabolite("product", formula="C", compartment="C_cy")
    product.annotation = {"metanetx.chemical": "MNXM101"}
    model.add_metabolites([substrate, product])

    table = tmp_path / "gap_fill.csv"
    table.write_text(
        "priority,mnxr_id,bigg_reaction,gene_id,equation,ec_number,kegg_reaction\n"
        "P0,MNXR_SAFE,SAFE_c,YALI0C23364g,"
        '"1 MNXM100@MNXD1 = 1 MNXM101@MNXD1",1.1.1.1,R00001\n'
        "P0,MNXR_UNSAFE,UNSAFE_c,YALI0C32184g,"
        '"1 MNXM100@MNXD1 = 1 MNXM101@MNXD1",1.1.1.1,R00002\n',
        encoding="utf-8",
    )

    stats = add_gap_fill_reactions(model, table, cache_dir=tmp_path / "cache")

    assert stats["added"] == ["SAFE_c"]
    assert stats["skipped_unresolved_genes"] == ["UNSAFE_c"]
    assert model.reactions.get_by_id("SAFE_c").gene_reaction_rule == "YALI1C32184g"
    assert "UNSAFE_c" not in model.reactions
    assert "YALI0C32184g" not in model.genes
