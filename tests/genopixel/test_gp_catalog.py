from __future__ import annotations

from pathlib import Path

from gp_catalog import GenoPixelCatalogStore, parse_listish


FIXTURE_XLSX = Path(__file__).resolve().parents[1] / "fixtures" / "catalog_fixture.xlsx"


def test_parse_listish_supports_brackets_scalars_and_semicolons() -> None:
    assert parse_listish("['kidney', 'blood']") == ["kidney", "blood"]
    assert parse_listish("Homo sapiens; Mus musculus") == ["Homo sapiens", "Mus musculus"]
    assert parse_listish("single value") == ["single value"]


def test_catalog_snapshot_normalizes_parent_rows_and_variants() -> None:
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)

    snapshot = store.get_snapshot()

    assert snapshot["source"]["all_rows"] == 3
    assert snapshot["source"]["multiple_rows"] == 2
    assert snapshot["source"]["multiple_parent_rows"] == 1
    assert snapshot["totals"]["multiple_datasets"] == 1

    datasets = snapshot["datasets"]
    assert [dataset["all_excel_row"] for dataset in datasets] == [2, 3, 4]

    single = datasets[0]
    assert single["primary_file"] == "atlas-1.h5ad"
    assert single["tissues"] == ["kidney", "blood"]
    assert single["diseases"] == ["normal"]
    assert single["variant_count"] == 0

    multiple = snapshot["dataset_details"][3]
    assert multiple["author"] == ""
    assert multiple["variant_count"] == 2
    assert multiple["publications"] == ["doi:multi"]
    assert [variant["multiple_excel_row"] for variant in multiple["variants"]] == [2, 3]
    assert multiple["variants"][0]["tissues"] == ["lung", "pleural effusion"]
    assert multiple["source_refs"]["all"]["excel_row"] == 3
    assert multiple["source_refs"]["multiple"][1]["excel_row"] == 3


def test_catalog_facets_include_normalized_values() -> None:
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)

    catalog = store.get_catalog_payload()

    assert catalog["facets"]["project"] == ["HCA", "cellxgene"]
    assert "Homo sapiens" in catalog["facets"]["organism"]
    assert "Mus musculus" in catalog["facets"]["organism"]
    assert catalog["facets"]["merged"] == ["multiple", "single"]
