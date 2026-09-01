import csv
import json
from pathlib import Path

import pytest
from cobra import Model, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.genes import apply_curated_gene_annotation_overrides


REPO_ROOT = Path(__file__).resolve().parents[1]


def _model_with_gene(gene_id: str = "YALI1E08382g") -> Model:
    model = Model("gene-override-test")
    reaction = Reaction("TEST_RXN")
    reaction.gene_reaction_rule = gene_id
    model.add_reactions([reaction])
    return model


def _write_override_table(path: Path, evidence_path: Path, duplicate: bool = False) -> None:
    fields = [
        "gene_id",
        "uniprot",
        "kegg.genes",
        "ncbigene",
        "refseq",
        "ec-code",
        "case_id",
        "evidence_path",
        "source_url",
        "reason",
    ]
    row = {
        "gene_id": "YALI1E08382g",
        "uniprot": "Q6C6R1",
        "kegg.genes": "yli:2912425",
        "ncbigene": "2912425",
        "refseq": "XP_503651.1",
        "ec-code": "2.4.2.19",
        "case_id": "EGC-test",
        "evidence_path": str(evidence_path),
        "source_url": "https://www.kegg.jp/entry/yli:2912425",
        "reason": "Correct an assembly identity collision.",
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(row)
        if duplicate:
            writer.writerow(row)


def test_override_replaces_stale_identity_and_preserves_other_annotations(tmp_path):
    model = _model_with_gene()
    gene = model.genes.get_by_id("YALI1E08382g")
    gene.annotation = {
        "sbo": "SBO:0000243",
        "uniprot": "Q6C6L8",
        "kegg.genes": "yli:2912487",
        "refseq": "XP_503694.1",
        "ec-code": "9.9.9.9",
    }
    evidence_path = tmp_path / "EGC-test.json"
    evidence_path.write_text(json.dumps({"case_id": "EGC-test"}), encoding="utf-8")
    table_path = tmp_path / "overrides.csv"
    _write_override_table(table_path, evidence_path)

    assert apply_curated_gene_annotation_overrides(model, table_path) == 1
    assert gene.annotation == {
        "sbo": "SBO:0000243",
        "uniprot": ["Q6C6R1"],
        "kegg.genes": ["yli:2912425"],
        "ncbigene": ["2912425"],
        "refseq": ["XP_503651.1"],
        "ec-code": ["2.4.2.19"],
    }
    assert apply_curated_gene_annotation_overrides(model, table_path) == 0


def test_override_rejects_duplicate_gene_rows(tmp_path):
    model = _model_with_gene()
    evidence_path = tmp_path / "EGC-test.json"
    evidence_path.write_text(json.dumps({"case_id": "EGC-test"}), encoding="utf-8")
    table_path = tmp_path / "overrides.csv"
    _write_override_table(table_path, evidence_path, duplicate=True)

    with pytest.raises(ValueError, match="Duplicate gene annotation override"):
        apply_curated_gene_annotation_overrides(model, table_path)


@pytest.mark.external_data
@pytest.mark.integration
def test_repository_override_corrects_base_model_identity(external_data_file):
    model = read_sbml_model(str(REPO_ROOT / "data" / "iyali26.xml"))

    table_path = external_data_file(
        "data/essentiality/curated_gene_annotation_overrides.csv"
    )
    assert apply_curated_gene_annotation_overrides(model, table_path) == 1
    annotation = model.genes.get_by_id("YALI1E08382g").annotation
    assert annotation["uniprot"] == ["Q6C6R1"]
    assert annotation["kegg.genes"] == ["yli:2912425"]
    assert annotation["ncbigene"] == ["2912425"]
    assert annotation["refseq"] == ["XP_503651.1"]
    assert annotation["ec-code"] == ["2.4.2.19"]
