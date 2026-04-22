from __future__ import annotations

import gc
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad

from gp_h5ad_loader import load_h5ad


class NoActiveDatasetError(RuntimeError):
    pass


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


class GenoPixelRuntimeState:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._active_adata: ad.AnnData | None = None
        self._active_h5ad_path: str | None = None
        self._active_all_excel_row: int | None = None
        self._active_multiple_excel_row: int | None = None
        self._active_title: str | None = None
        self._active_loaded_at: str | None = None
        self._active_backed: bool | None = None
        self._active_total_cells: int | None = None
        # Pending selection: recorded from DynamoDB even before the h5ad is loaded.
        self._pending_selection: dict[str, Any] | None = None

    def set_pending_selection(
        self,
        *,
        all_excel_row: int,
        multiple_excel_row: int | None,
        title: str,
        primary_file: str,
    ) -> None:
        with self._lock:
            self._pending_selection = {
                "all_excel_row": all_excel_row,
                "multiple_excel_row": multiple_excel_row,
                "title": title,
                "primary_file": primary_file,
            }

    def get_pending_selection(self) -> dict[str, Any] | None:
        with self._lock:
            if self._active_adata is not None:
                return None  # already loaded, not pending
            return self._pending_selection

    def load_active_dataset(
        self,
        *,
        h5ad_path: str,
        all_excel_row: int,
        multiple_excel_row: int | None,
        title: str,
        backed: bool,
        force_reload: bool = False,
    ) -> dict[str, Any]:
        normalized_path = str(Path(h5ad_path).expanduser().resolve())
        dataset_key = (all_excel_row, multiple_excel_row, normalized_path, backed)

        with self._lock:
            current_key = (
                self._active_all_excel_row,
                self._active_multiple_excel_row,
                self._active_h5ad_path,
                self._active_backed,
            )
            if force_reload or self._active_adata is None or current_key != dataset_key:
                self._clear_active_dataset_unlocked()
                self._active_adata = load_h5ad(Path(normalized_path), backed=backed)
                self._active_h5ad_path = normalized_path
                self._active_all_excel_row = all_excel_row
                self._active_multiple_excel_row = multiple_excel_row
                self._active_title = title
                self._active_backed = backed
                self._active_loaded_at = _utc_now_iso()
                self._active_total_cells = int(self._active_adata.n_obs)
                self._pending_selection = None

            return self.get_active_dataset_payload_unlocked()

    def get_active_dataset_payload(self) -> dict[str, Any]:
        with self._lock:
            return self.get_active_dataset_payload_unlocked()

    def get_active_dataset_payload_unlocked(self) -> dict[str, Any]:
        return {
            "loaded": self._active_adata is not None,
            "all_excel_row": self._active_all_excel_row,
            "multiple_excel_row": self._active_multiple_excel_row,
            "title": self._active_title,
            "h5ad_path": self._active_h5ad_path,
            "loaded_at": self._active_loaded_at,
            "backed": self._active_backed,
            "total_cells": self._active_total_cells,
        }

    def require_active_adata(self) -> tuple[ad.AnnData, dict[str, Any]]:
        with self._lock:
            if self._active_adata is None:
                raise NoActiveDatasetError(
                    "No dataset is currently loaded. Please load a dataset from the GenoPixel browser first."
                )
            return self._active_adata, self.get_active_dataset_payload_unlocked()

    def clear_active_dataset(self) -> None:
        with self._lock:
            self._clear_active_dataset_unlocked()

    def _clear_active_dataset_unlocked(self) -> None:
        self._close_active_adata_unlocked()
        self._active_adata = None
        self._active_h5ad_path = None
        self._active_all_excel_row = None
        self._active_multiple_excel_row = None
        self._active_title = None
        self._active_loaded_at = None
        self._active_backed = None
        self._active_total_cells = None
        gc.collect()

    def _close_active_adata_unlocked(self) -> None:
        if self._active_adata is None:
            return
        file_manager = getattr(self._active_adata, "file", None)
        if file_manager is not None:
            try:
                file_manager.close()
            except Exception:
                pass


RUNTIME_STATE = GenoPixelRuntimeState()
