from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from gp_scanpy_plotter import ScanpyPlotExecutor


def _active_payload(*, row: int, loaded_at: str, path: str = "/tmp/test.h5ad") -> dict[str, object]:
    return {
        "loaded": True,
        "all_excel_row": row,
        "multiple_excel_row": None,
        "h5ad_path": path,
        "backed": False,
        "loaded_at": loaded_at,
    }


def test_sync_active_dataset_defaults_to_author_cell_type_and_resets_on_dataset_change(tmp_path) -> None:
    executor = ScanpyPlotExecutor(output_dir=tmp_path)
    adata = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "author_cell_type": ["T cell", "B cell", "NK cell"],
                "manual_cell_type": ["lymphocyte", "lymphocyte", "lymphocyte"],
            }
        )
    )

    first_active = _active_payload(row=2, loaded_at="2026-03-11T10:00:00Z")
    executor.sync_active_dataset(first_active, adata)

    assert executor._resolve_cell_type_column(adata) == "author_cell_type"

    assert executor._resolve_cell_type_column(adata, "manual_cell_type") == "manual_cell_type"
    assert executor._resolve_cell_type_column(adata) == "manual_cell_type"

    # Same active dataset should keep the user-selected override.
    executor.sync_active_dataset(first_active, adata)
    assert executor._resolve_cell_type_column(adata) == "manual_cell_type"

    # New active dataset should reset to the default author cell-type column.
    second_active = _active_payload(row=3, loaded_at="2026-03-11T10:05:00Z")
    executor.sync_active_dataset(second_active, adata)
    assert executor._resolve_cell_type_column(adata) == "author_cell_type"


def test_sync_active_dataset_matches_default_author_cell_type_case_insensitively(tmp_path) -> None:
    executor = ScanpyPlotExecutor(output_dir=tmp_path)
    adata = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "Author_Cell_Type": ["stromal", "immune", "endothelial"],
                "cluster_id": ["0", "1", "2"],
            }
        )
    )

    executor.sync_active_dataset(_active_payload(row=2, loaded_at="2026-03-11T11:00:00Z"), adata)
    assert executor._resolve_cell_type_column(adata) == "Author_Cell_Type"
