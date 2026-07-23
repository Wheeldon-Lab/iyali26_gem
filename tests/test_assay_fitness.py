import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from scripts.gem_annotate.essentiality_evidence import sha256_file
from scripts.gem_annotate.prepare_essentiality_assays import (
    NORMALIZED_COLUMNS,
    SOURCE_COLUMNS,
    normalize_workbook,
    write_normalized_assays,
)
from scripts.gem_annotate.validate_essential_genes import (
    build_assay_fitness_table,
    build_run_manifest,
    load_assay_fitness,
    make_assay_fitness_summary,
    validate_essential_genes,
)


def _source_rows(gene_ids: list[str], calls: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Gene ID": gene_ids,
            "FS": [float(index - 2) for index in range(len(gene_ids))],
            "Raw p-value": [0.01 * (index + 1) for index in range(len(gene_ids))],
            "Corrected p-value": [0.02 * (index + 1) for index in range(len(gene_ids))],
            "Essentiality": calls,
        },
        columns=SOURCE_COLUMNS,
    )


def _write_source_workbook(
    path: Path,
    cas9: pd.DataFrame,
    cas12a: pd.DataFrame,
) -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        cas9.to_excel(writer, sheet_name="Cas9", index=False)
        cas12a.to_excel(writer, sheet_name="Cas12a", index=False)


def test_workbook_normalization_preserves_source_provenance(tmp_path: Path) -> None:
    workbook = tmp_path / "screens.xlsx"
    output = tmp_path / "normalized.csv"
    _write_source_workbook(
        workbook,
        _source_rows(
            ["YALI1_E10500g", "YALI1_C353g"],
            ["Essential", "Non-essential"],
        ),
        _source_rows(
            ["YALI1_A00309g", "YALI1_F34083g"],
            ["Non-essential", "Essential"],
        ),
    )

    normalized = write_normalized_assays(
        workbook,
        output,
        expected_rows={"Cas9": 2, "Cas12a": 2},
    )

    assert list(normalized.columns) == list(NORMALIZED_COLUMNS)
    assert normalized["gene_id"].tolist() == [
        "YALI1E10500g",
        "YALI1C353g",
        "YALI1A00309g",
        "YALI1F34083g",
    ]
    assert normalized["source_gene_id"].tolist()[0] == "YALI1_E10500g"
    assert normalized["experimental_call"].tolist() == [
        "essential",
        "nonessential",
        "nonessential",
        "essential",
    ]
    assert normalized.groupby("assay")["source_row"].apply(list).to_dict() == {
        "Cas12a": [2, 3],
        "Cas9": [2, 3],
    }
    assert normalized["source_sha256"].nunique() == 1
    assert normalized["source_sha256"].iloc[0] == sha256_file(workbook)
    roundtrip = load_assay_fitness(output, expected_rows=None)
    assert len(roundtrip) == 4
    assert not roundtrip.isna().any().any()


@pytest.mark.parametrize("defect", ["row_count", "missing", "duplicate", "columns"])
def test_workbook_normalization_fails_closed(tmp_path: Path, defect: str) -> None:
    workbook = tmp_path / "screens.xlsx"
    cas9 = _source_rows(
        ["YALI1_A00309g", "YALI1_B07538g"],
        ["Essential", "Non-essential"],
    )
    cas12a = _source_rows(
        ["YALI1_C24755g", "YALI1_D28221g"],
        ["Essential", "Non-essential"],
    )
    expected = {"Cas9": 2, "Cas12a": 2}
    if defect == "row_count":
        expected["Cas9"] = 3
    elif defect == "missing":
        cas9.loc[0, "FS"] = None
    elif defect == "duplicate":
        cas9.loc[1, "Gene ID"] = "YALI1_A00309g"
    elif defect == "columns":
        cas9 = cas9.rename(columns={"FS": "Fitness"})
    _write_source_workbook(workbook, cas9, cas12a)

    with pytest.raises(ValueError):
        normalize_workbook(workbook, expected_rows=expected)


def _normalized_fixture() -> pd.DataFrame:
    rows = []
    source_sha = "a" * 64
    calls = ["essential", "essential", "nonessential", "nonessential"]
    for assay in ("Cas9", "Cas12a"):
        for index, gene_id in enumerate(("g1", "g2", "g3", "g4"), start=2):
            rows.append(
                {
                    "gene_id": gene_id,
                    "source_gene_id": gene_id,
                    "assay": assay,
                    "fitness_score": float(index - 2),
                    "raw_p_value": 0.1,
                    "q_value": 0.2,
                    "experimental_call": calls[index - 2],
                    "source_sheet": assay,
                    "source_row": index,
                    "source_sha256": source_sha,
                }
            )
    return pd.DataFrame(rows, columns=NORMALIZED_COLUMNS)


def test_normalized_loader_rejects_duplicate_assay_gene_rows(tmp_path: Path) -> None:
    table = _normalized_fixture()
    table.loc[1, "gene_id"] = table.loc[0, "gene_id"]
    table.loc[1, "source_gene_id"] = table.loc[0, "source_gene_id"]
    path = tmp_path / "duplicate.csv"
    table.to_csv(path, index=False)

    with pytest.raises(ValueError, match="Duplicate assay/gene"):
        load_assay_fitness(path, expected_rows=None)


def test_assay_metrics_are_separate_and_directionally_correct() -> None:
    assays = _normalized_fixture()
    predictions = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "ko_status": ["optimal"] * 4,
            "ko_growth": [0.01, 0.34, 0.67, 1.0],
            "ko_growth_ratio": [0.01, 0.34, 0.67, 1.0],
        }
    )
    joined = build_assay_fitness_table(assays, predictions, primary_cutoff=0.1)
    summary = make_assay_fitness_summary(
        joined,
        predictions["gene_id"],
        primary_cutoff=0.1,
        growth_cutoffs=(0.1, 0.9),
    )

    assert summary["semantics"]["positive_only_recall_unchanged"] is True
    for assay in ("Cas9", "Cas12a"):
        metrics = summary["per_assay"][assay]
        assert metrics["model_overlap"] == 4
        assert metrics["spearman_fitness_vs_ko_growth_ratio"] == pytest.approx(1.0)
        assert metrics["linear_r_squared_fitness_vs_ko_growth_ratio"] == pytest.approx(
            1.0
        )
        assert metrics["assay_call_roc_auc"] == pytest.approx(1.0)

    proxy = summary["concordant_nonessential_proxy_safety"]
    assert proxy["source_proxy_genes"] == 2
    assert proxy["primary"] == {
        "cutoff_fraction_of_wt": 0.1,
        "model_proxy_genes": 2,
        "predicted_essential_count": 0,
        "safe_nonessential_count": 2,
        "safety_rate": 1.0,
    }
    assert proxy["cutoff_curve"][1]["safety_rate"] == pytest.approx(0.5)


def test_assay_metrics_treat_solver_noise_as_a_tie() -> None:
    assays = _normalized_fixture()
    first = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "ko_growth_ratio": [0.01, 0.5, 1.0, 1.0],
        }
    )
    second = first.copy()
    second.loc[2, "ko_growth_ratio"] += 2e-12
    second.loc[3, "ko_growth_ratio"] -= 3e-12

    summaries = []
    for predictions in (first, second):
        joined = build_assay_fitness_table(assays, predictions, primary_cutoff=0.1)
        summaries.append(
            make_assay_fitness_summary(
                joined,
                predictions["gene_id"],
                primary_cutoff=0.1,
                growth_cutoffs=(0.01, 0.05, 0.1, 0.15),
            )
        )

    assert summaries[0]["per_assay"] == summaries[1]["per_assay"]
    assert "solver-stable" in summaries[0]["semantics"][
        "ko_ratio_metric_quantization"
    ]


def test_run_manifest_records_inputs_solver_versions_git_and_cutoffs(
    tmp_path: Path,
) -> None:
    model_path = tmp_path / "model.xml"
    experimental_path = tmp_path / "experimental.csv"
    media_path = tmp_path / "medium.csv"
    assay_path = tmp_path / "assays.csv"
    for path, content in (
        (model_path, "model"),
        (experimental_path, "experimental"),
        (media_path, "medium"),
        (assay_path, "assays"),
    ):
        path.write_text(content)
    model = SimpleNamespace(
        solver=SimpleNamespace(
            interface=SimpleNamespace(__name__="optlang.gurobi_interface")
        )
    )

    manifest = build_run_manifest(
        model=model,
        model_path=model_path,
        experimental_path=experimental_path,
        media_path=media_path,
        assay_fitness_path=assay_path,
        assay_source_sha256s=["b" * 64],
        requested_solver="gurobi",
        primary_cutoff=0.1,
        growth_cutoffs=(0.01, 0.05, 0.1, 0.15),
        positive_only=True,
        repo_root=tmp_path,
        generated_at="2026-07-20T00:00:00+00:00",
    )

    assert manifest["inputs"]["model"]["sha256"] == sha256_file(model_path)
    assert manifest["inputs"]["assay_fitness"]["source_workbook_sha256"] == ["b" * 64]
    assert manifest["solver"]["requested"] == "gurobi"
    assert manifest["solver"]["actual"] == "gurobi"
    assert manifest["solver"]["backend_distribution"] == "gurobipy"
    assert manifest["software"]["cobra"]
    assert set(manifest["git"]) == {"commit", "branch", "dirty"}
    assert manifest["cutoffs"]["primary_fraction_of_wt"] == 0.1
    assert manifest["configuration"]["credentials_included"] is False
    serialized = json.dumps(manifest).casefold()
    assert "password" not in serialized
    assert "license key" not in serialized


def test_validator_defaults_to_gurobi() -> None:
    solver = inspect.signature(validate_essential_genes).parameters["solver"]
    assert solver.default == "gurobi"
