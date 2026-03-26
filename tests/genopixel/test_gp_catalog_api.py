from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from fastapi.testclient import TestClient

import gp_catalog
import gp_catalog_api
import gp_runtime_state
from gp_catalog import GenoPixelCatalogStore
from gp_catalog_api import create_app
from gp_models import PlotResult


FIXTURE_XLSX = Path(__file__).resolve().parents[1] / "fixtures" / "catalog_fixture.xlsx"


def test_catalog_endpoint_returns_expected_shape() -> None:
    app = create_app(GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX))
    client = TestClient(app)

    response = client.get("/api/genopixel-catalog/catalog")

    assert response.status_code == 200
    payload = response.json()
    assert payload["source"]["all_rows"] == 3
    assert len(payload["datasets"]) == 3
    assert payload["datasets"][1]["variant_count"] == 2
    assert payload["datasets"][1]["author"] == ""


def test_catalog_routes_are_hidden_from_openapi_schema() -> None:
    app = create_app(GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX))
    client = TestClient(app)

    response = client.get("/openapi.json")

    assert response.status_code == 200
    paths = response.json()["paths"]
    assert "/generate_scanpy_plot" in paths
    assert "/api/genopixel-catalog/catalog" not in paths
    assert "/api/genopixel-catalog/datasets/{all_excel_row}/analyze" not in paths
    assert "/resolve_and_plot" not in paths


def test_dataset_detail_returns_variants_for_multiple_rows_only() -> None:
    app = create_app(GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX))
    client = TestClient(app)

    single = client.get("/api/genopixel-catalog/datasets/2")
    multiple = client.get("/api/genopixel-catalog/datasets/3")

    assert single.status_code == 200
    assert single.json()["dataset"]["variants"] == []

    assert multiple.status_code == 200
    assert len(multiple.json()["dataset"]["variants"]) == 2
    assert "resolved_h5ad_path" in multiple.json()["dataset"]["variants"][0]


def test_dataset_detail_returns_404_for_missing_row() -> None:
    app = create_app(GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX))
    client = TestClient(app)

    response = client.get("/api/genopixel-catalog/datasets/999")

    assert response.status_code == 404
    assert "999" in response.json()["detail"]


def test_catalog_endpoint_returns_503_when_workbook_is_missing(tmp_path) -> None:
    missing = tmp_path / "missing.xlsx"
    app = create_app(GenoPixelCatalogStore(xlsx_path=missing))
    client = TestClient(app)

    response = client.get("/api/genopixel-catalog/catalog")

    assert response.status_code == 503
    assert "Workbook not found" in response.json()["detail"]


def test_analyze_endpoint_loads_single_dataset(monkeypatch, tmp_path) -> None:
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(
        gp_runtime_state,
        "load_h5ad",
        lambda path, backed=False: SimpleNamespace(filename=str(path), backed=backed, n_obs=1000),
    )

    app = create_app(store)
    client = TestClient(app)

    response = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["message"] == "data is loaded, happy analysis"
    assert payload["all_excel_row"] == 2
    assert payload["h5ad_path"] == str(resolved_path.resolve())
    assert payload["active_dataset"]["title"]


def test_analyze_endpoint_force_reloads_same_dataset(monkeypatch, tmp_path) -> None:
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    load_calls = {"count": 0}

    def fake_load_h5ad(path, backed=False):
        load_calls["count"] += 1
        return SimpleNamespace(filename=str(path), backed=backed, n_obs=1111, file=None)

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(gp_runtime_state, "load_h5ad", fake_load_h5ad)

    app = create_app(store)
    client = TestClient(app)

    first = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )
    second = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )

    assert first.status_code == 200
    assert second.status_code == 200
    assert load_calls["count"] == 2
    assert second.json()["backed"] is False
    assert second.json()["active_dataset"]["backed"] is False


def test_analyze_endpoint_loads_selected_subdataset(monkeypatch, tmp_path) -> None:
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-variant.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(
        gp_runtime_state,
        "load_h5ad",
        lambda path, backed=False: SimpleNamespace(filename=str(path), backed=backed, n_obs=2000),
    )

    app = create_app(store)
    client = TestClient(app)
    detail = client.get("/api/genopixel-catalog/datasets/3").json()["dataset"]

    response = client.post(
        "/api/genopixel-catalog/datasets/3/analyze",
        json={
            "h5ad_path": detail["variants"][0]["resolved_h5ad_path"],
            "multiple_excel_row": detail["variants"][0]["multiple_excel_row"],
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["message"] == "data is loaded, happy analysis"
    assert payload["all_excel_row"] == 3
    assert payload["multiple_excel_row"] == detail["variants"][0]["multiple_excel_row"]


def test_active_dataset_endpoint_reflects_latest_loaded_dataset(monkeypatch, tmp_path) -> None:
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(
        gp_runtime_state,
        "load_h5ad",
        lambda path, backed=False: SimpleNamespace(filename=str(path), backed=backed, n_obs=3000),
    )

    app = create_app(store)
    client = TestClient(app)

    client.post("/api/genopixel-catalog/datasets/2/analyze", json={"h5ad_path": str(resolved_path.resolve())})
    response = client.get("/api/genopixel-runtime/active-dataset")

    assert response.status_code == 200
    payload = response.json()
    assert payload["loaded"] is True
    assert payload["all_excel_row"] == 2
    assert payload["h5ad_path"] == str(resolved_path.resolve())
    assert payload["total_cells"] is None or isinstance(payload["total_cells"], int)


def test_generate_scanpy_plot_requires_loaded_dataset(monkeypatch) -> None:
    monkeypatch.setenv("API_KEY", "test-key")
    app = create_app(GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX))
    client = TestClient(app)

    response = client.post(
        "/generate_scanpy_plot",
        headers={"Authorization": "Bearer test-key"},
        json={"plot_type": "umap", "color_json": "[]", "genes_json": "[]"},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["ok"] is False
    assert payload["status"] == "no_active_dataset"


def test_generate_scanpy_plot_uses_loaded_dataset_without_reloading(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("API_KEY", "test-key")
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    load_calls = {"count": 0}

    def fake_load_h5ad(path, backed=False):
        load_calls["count"] += 1
        return SimpleNamespace(filename=str(path), backed=backed, n_obs=4000)

    class FakePlotter:
        def run(self, adata, request):
            assert adata.filename == str(resolved_path.resolve())
            assert request.plot_type == "umap"
            output_file = tmp_path / "plot.png"
            output_file.write_text("plot", encoding="utf-8")
            return PlotResult(
                plot_type="umap",
                output_file=output_file,
                embedding_basis="X_umap",
                color_columns=["cell_type"],
                resolved_coloring_label="Cell types (cell_type)",
                display_plot_type="UMAP embedding",
            )

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(gp_runtime_state, "load_h5ad", fake_load_h5ad)
    monkeypatch.setattr(gp_catalog_api, "_PLOTTER", FakePlotter())

    app = create_app(store)
    client = TestClient(app)

    analyze_response = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )
    assert analyze_response.status_code == 200
    assert load_calls["count"] == 1

    plot_response = client.post(
        "/generate_scanpy_plot",
        headers={"Authorization": "Bearer test-key"},
        json={"plot_type": "umap", "color_json": "[\"cell_type\"]", "genes_json": "[]"},
    )

    assert plot_response.status_code == 200
    payload = plot_response.json()
    assert payload["ok"] is True
    assert payload["status"] == "success"
    assert load_calls["count"] == 1
    assert payload["plot"]["output_file"].endswith("plot.png")
    assert payload["canonical_response_markdown"] is not None


def test_generate_scanpy_plot_uses_public_url_in_markdown(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("API_KEY", "test-key")
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    def fake_load_h5ad(path, backed=False):
        return SimpleNamespace(filename=str(path), backed=backed, n_obs=2500)

    class FakePlotter:
        def run(self, adata, request):
            output_root = Path("/code/out/genopixel")
            output_root.mkdir(parents=True, exist_ok=True)
            output_file = output_root / "plot.png"
            output_file.write_text("plot", encoding="utf-8")
            return PlotResult(
                plot_type="umap",
                output_file=output_file,
                embedding_basis="X_umap",
                color_columns=[],
                display_plot_type="UMAP embedding",
            )

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(gp_runtime_state, "load_h5ad", fake_load_h5ad)
    monkeypatch.setattr(gp_catalog_api, "_PLOTTER", FakePlotter())

    app = create_app(store)
    client = TestClient(app)

    analyze_response = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )
    assert analyze_response.status_code == 200

    plot_response = client.post(
        "/generate_scanpy_plot",
        headers={"Authorization": "Bearer test-key"},
        json={"plot_type": "umap", "color_json": "[]", "genes_json": "[]"},
    )

    assert plot_response.status_code == 200
    payload = plot_response.json()
    assert payload["output_url"] == "http://localhost/assets/genopixel/plot.png"
    assert payload["output_markdown"] == "![Plot](http://localhost/assets/genopixel/plot.png)"
    assert payload["inline_markdown"] == payload["output_markdown"]
    assert "![Plot](http://localhost/assets/genopixel/plot.png)" in payload["canonical_response_markdown"]


def test_generate_scanpy_plot_returns_plain_umap_response_markdown(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("API_KEY", "test-key")
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    fake_adata = SimpleNamespace(filename=str(resolved_path.resolve()), backed=False, n_obs=12345)

    def fake_load_h5ad(path, backed=False):
        return fake_adata

    class FakePlotter:
        def run(self, adata, request):
            output_root = Path("/code/out/genopixel")
            output_root.mkdir(parents=True, exist_ok=True)
            output_file = output_root / "umap.png"
            output_file.write_text("plot", encoding="utf-8")
            return PlotResult(
                plot_type="umap",
                output_file=output_file,
                embedding_basis="X_umap",
                color_columns=["author_cell_type1"],
                resolved_coloring_label="Cell types (author_cell_type1)",
                display_plot_type="UMAP embedding",
            )

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(gp_runtime_state, "load_h5ad", fake_load_h5ad)
    monkeypatch.setattr(gp_catalog_api, "_PLOTTER", FakePlotter())

    app = create_app(store)
    client = TestClient(app)

    analyze_response = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )
    assert analyze_response.status_code == 200

    plot_response = client.post(
        "/generate_scanpy_plot",
        headers={"Authorization": "Bearer test-key"},
        json={"plot_type": "umap", "color_json": "[]", "genes_json": "[]"},
    )

    assert plot_response.status_code == 200
    payload = plot_response.json()
    assert payload["active_dataset"]["total_cells"] == 12345
    assert payload["canonical_response_markdown"] == payload["plain_umap_response_markdown"]
    assert payload["plain_umap_response_markdown"] is not None
    assert "**Dataset:**" in payload["plain_umap_response_markdown"]
    assert "- **Total cells:** 12,345 cells" in payload["plain_umap_response_markdown"]
    assert "- **Visualization type:** UMAP embedding" in payload["plain_umap_response_markdown"]
    assert "- **Coloring:** Cell types (author_cell_type1)" in payload["plain_umap_response_markdown"]
    assert "![Plot](http://localhost/assets/genopixel/umap.png)" in payload["plain_umap_response_markdown"]


def test_generate_scanpy_plot_returns_distribution_canonical_response(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("API_KEY", "test-key")
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    def fake_load_h5ad(path, backed=False):
        return SimpleNamespace(filename=str(path), backed=backed, n_obs=6789)

    class FakePlotter:
        def run(self, adata, request):
            output_root = Path("/code/out/genopixel")
            output_root.mkdir(parents=True, exist_ok=True)
            output_file = output_root / "violin.png"
            output_file.write_text("plot", encoding="utf-8")
            return PlotResult(
                plot_type="violin",
                output_file=output_file,
                resolved_genes=["OSMR"],
                resolved_groupby="author_cell_type1",
                display_plot_type="Violin plot",
            )

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(gp_runtime_state, "load_h5ad", fake_load_h5ad)
    monkeypatch.setattr(gp_catalog_api, "_PLOTTER", FakePlotter())

    app = create_app(store)
    client = TestClient(app)

    analyze_response = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )
    assert analyze_response.status_code == 200

    plot_response = client.post(
        "/generate_scanpy_plot",
        headers={"Authorization": "Bearer test-key"},
        json={"plot_type": "violin", "genes_json": "[\"OSMR\"]", "groupby": "author_cell_type1"},
    )

    assert plot_response.status_code == 200
    payload = plot_response.json()
    assert payload["canonical_response_markdown"] is not None
    assert "- **Plot type:** Violin plot" in payload["canonical_response_markdown"]
    assert "- **Genes:** OSMR" in payload["canonical_response_markdown"]
    assert "- **Grouping column:** author_cell_type1" in payload["canonical_response_markdown"]


def test_generate_scanpy_plot_returns_rank_genes_cache_status(monkeypatch, tmp_path) -> None:
    monkeypatch.setenv("API_KEY", "test-key")
    store = GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX)
    resolved_path = tmp_path / "atlas-1.h5ad"
    resolved_path.write_text("placeholder", encoding="utf-8")

    def fake_load_h5ad(path, backed=False):
        return SimpleNamespace(filename=str(path), backed=backed, n_obs=2222)

    class FakePlotter:
        def run(self, adata, request):
            output_root = Path("/code/out/genopixel")
            output_root.mkdir(parents=True, exist_ok=True)
            output_file = output_root / "rank_dotplot.png"
            output_file.write_text("plot", encoding="utf-8")
            return PlotResult(
                plot_type="rank_genes_groups_dotplot",
                output_file=output_file,
                resolved_groupby="author_cell_type1",
                display_plot_type="Rank-genes dot plot",
                rank_genes_groups_computed=False,
                rank_genes_groups_notice="Reusing cached rank_genes_groups computed with groupby 'author_cell_type1'.",
            )

    monkeypatch.setattr(store, "_resolve_h5ad_path", lambda file_value: resolved_path)
    monkeypatch.setattr(gp_runtime_state, "load_h5ad", fake_load_h5ad)
    monkeypatch.setattr(gp_catalog_api, "_PLOTTER", FakePlotter())

    app = create_app(store)
    client = TestClient(app)

    analyze_response = client.post(
        "/api/genopixel-catalog/datasets/2/analyze",
        json={"h5ad_path": str(resolved_path.resolve())},
    )
    assert analyze_response.status_code == 200

    plot_response = client.post(
        "/generate_scanpy_plot",
        headers={"Authorization": "Bearer test-key"},
        json={"plot_type": "rank_genes_groups_dotplot", "genes_json": "[]", "groupby": "author_cell_type1"},
    )

    assert plot_response.status_code == 200
    payload = plot_response.json()
    assert payload["plot"]["rank_genes_groups_computed"] is False
    assert "Reusing cached rank_genes_groups" in payload["plot"]["rank_genes_groups_notice"]
    assert payload["rank_genes_groups_computed"] is False
    assert "Rank-genes status" in payload["canonical_response_markdown"]


def test_generate_scanpy_plot_requires_auth(monkeypatch) -> None:
    monkeypatch.setenv("API_KEY", "test-key")
    app = create_app(GenoPixelCatalogStore(xlsx_path=FIXTURE_XLSX))
    client = TestClient(app)

    response = client.post("/generate_scanpy_plot", json={"plot_type": "umap"})

    assert response.status_code == 401
