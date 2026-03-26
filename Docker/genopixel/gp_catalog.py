from __future__ import annotations

import ast
import threading
from collections import Counter
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from gp_config import load_settings
from gp_runtime_state import RUNTIME_STATE

PARENT_MULTI_VALUE_COLUMNS = {
    "tissue": "tissues",
    "tissue_type": "tissue_types",
    "disease": "diseases",
    "organism": "organisms",
}

VARIANT_MULTI_VALUE_COLUMNS = {
    "tissue": "tissues",
    "tissue_type": "tissue_types",
    "disease": "diseases",
    "organism": "organisms",
}

FACET_FIELDS = (
    "project",
    "organisms",
    "tissues",
    "tissue_types",
    "diseases",
    "journal",
    "merged",
)


class CatalogLoadError(RuntimeError):
    pass


class DatasetLoadError(RuntimeError):
    pass


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _normalize_column_name(value: object) -> str:
    return str(value).strip().lower()


def _normalize_text(value: object) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _normalized_key(value: object) -> str:
    return _normalize_text(value).lower()


def _coerce_int(value: object) -> int | None:
    if value is None or pd.isna(value):
        return None
    try:
        return int(value)
    except Exception:
        return None


def parse_listish(value: object) -> list[str]:
    raw = _normalize_text(value)
    if not raw:
        return []

    if raw.startswith("["):
        try:
            parsed = ast.literal_eval(raw)
        except Exception:
            parsed = None
        if isinstance(parsed, list):
            return [item for item in (_normalize_text(entry) for entry in parsed) if item]
        if parsed is not None:
            scalar = _normalize_text(parsed)
            return [scalar] if scalar else []

    if ";" in raw:
        return [item for item in (_normalize_text(part) for part in raw.split(";")) if item]

    return [raw]


def _public_parent_record(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "all_excel_row": row["all_excel_row"],
        "project": row["project"],
        "doi": row["doi"],
        "cellxgene_doi": row["cellxgene_doi"],
        "title": row["title"],
        "author": row["author"],
        "year": row["year"],
        "journal": row["journal"],
        "tissues": row["tissues"],
        "tissue_types": row["tissue_types"],
        "diseases": row["diseases"],
        "organisms": row["organisms"],
        "cell_counts": row["cell_counts"],
        "merged": row["merged"],
        "primary_file": row["primary_file"],
        "variant_count": row["variant_count"],
    }


class GenoPixelCatalogStore:
    def __init__(self, xlsx_path: Path | None = None):
        settings = load_settings()
        self.xlsx_path = Path(xlsx_path or settings.metadata_xlsx).expanduser()
        self._lock = threading.Lock()
        self._snapshot: dict[str, Any] | None = None
        self._workbook_mtime: float | None = None
        self._last_reload_time: str | None = None
        self._cache_hits = 0
        self._settings = settings

    def get_health_payload(self) -> dict[str, Any]:
        workbook_exists = self.xlsx_path.exists()
        cache_loaded = self._snapshot is not None
        return {
            "ok": workbook_exists,
            "xlsx_path": str(self.xlsx_path),
            "workbook_exists": workbook_exists,
            "cache": {
                "loaded": cache_loaded,
                "cache_hits": self._cache_hits,
                "last_reload_time": self._last_reload_time,
                "workbook_mtime": self._snapshot["source"]["mtime"] if cache_loaded else None,
            },
            "active_dataset": RUNTIME_STATE.get_active_dataset_payload(),
        }

    def get_catalog_payload(self) -> dict[str, Any]:
        snapshot = self.get_snapshot()
        return {
            "generated_at": _utc_now_iso(),
            "source": snapshot["source"],
            "totals": snapshot["totals"],
            "facets": snapshot["facets"],
            "datasets": snapshot["datasets"],
        }

    def get_dataset_payload(self, all_excel_row: int) -> dict[str, Any]:
        snapshot = self.get_snapshot()
        details = snapshot["dataset_details"]
        if all_excel_row not in details:
            raise KeyError(all_excel_row)
        return {
            "generated_at": _utc_now_iso(),
            "source": snapshot["source"],
            "dataset": details[all_excel_row],
        }

    def analyze_dataset(
        self,
        all_excel_row: int,
        h5ad_path: str,
        multiple_excel_row: int | None = None,
    ) -> dict[str, Any]:
        snapshot = self.get_snapshot()
        details = snapshot["dataset_details"]
        if all_excel_row not in details:
            raise KeyError(all_excel_row)

        detail = details[all_excel_row]
        if detail["merged"] == "multiple" and multiple_excel_row is None:
            raise DatasetLoadError(
                f"Dataset row {all_excel_row} requires a sub-dataset selection before loading."
            )

        raw_h5ad_path = str(h5ad_path or "").strip()
        if not raw_h5ad_path:
            raise DatasetLoadError(f"Dataset row {all_excel_row} does not have a resolved h5ad path.")
        normalized_h5ad_path = str(Path(raw_h5ad_path).expanduser().resolve())

        if detail["merged"] == "single":
            expected_path = str(detail.get("resolved_h5ad_path") or "").strip()
            if expected_path and normalized_h5ad_path != str(Path(expected_path).expanduser().resolve()):
                raise DatasetLoadError(
                    f"Dataset row {all_excel_row} does not match the provided resolved h5ad path."
                )
        else:
            selected_variant = next(
                (variant for variant in detail.get("variants", []) if variant.get("multiple_excel_row") == multiple_excel_row),
                None,
            )
            if selected_variant is None:
                raise DatasetLoadError(
                    f"Dataset row {all_excel_row} does not include sub-dataset row {multiple_excel_row}."
                )
            expected_path = str(selected_variant.get("resolved_h5ad_path") or "").strip()
            if expected_path and normalized_h5ad_path != str(Path(expected_path).expanduser().resolve()):
                raise DatasetLoadError(
                    f"Sub-dataset row {multiple_excel_row} does not match the provided resolved h5ad path."
                )
        backed = bool(self._settings.default_backed)

        try:
            active_dataset = RUNTIME_STATE.load_active_dataset(
                h5ad_path=normalized_h5ad_path,
                all_excel_row=all_excel_row,
                multiple_excel_row=multiple_excel_row,
                title=detail.get("title") or f"Dataset {all_excel_row}",
                backed=backed,
                force_reload=True,
            )
        except Exception as exc:
            raise DatasetLoadError(f"Failed to load h5ad file '{normalized_h5ad_path}': {exc}") from exc

        return {
            "message": "data is loaded, happy analysis",
            "all_excel_row": all_excel_row,
            "multiple_excel_row": multiple_excel_row,
            "h5ad_path": active_dataset["h5ad_path"],
            "backed": backed,
            "active_dataset": {
                "title": active_dataset["title"],
                "loaded_at": active_dataset["loaded_at"],
                "backed": active_dataset["backed"],
            },
        }

    def _resolve_h5ad_path(self, file_value: str) -> Path:
        candidate = Path(str(file_value).strip()).expanduser()
        if candidate.is_absolute() and candidate.exists():
            return candidate.resolve()

        joined = (self._settings.h5ad_base_dir / candidate).expanduser()
        if joined.exists():
            return joined.resolve()

        files = sorted(self._settings.h5ad_base_dir.rglob("*.h5ad"))
        target_name = candidate.name.lower()
        target_stem = candidate.stem.lower()

        basename_matches = [path for path in files if path.name.lower() == target_name]
        if len(basename_matches) == 1:
            return basename_matches[0].resolve()
        if len(basename_matches) > 1:
            raise FileNotFoundError(f"Multiple h5ad files matched basename '{candidate.name}'.")

        stem_matches = [path for path in files if target_stem and target_stem in path.stem.lower()]
        if len(stem_matches) == 1:
            return stem_matches[0].resolve()
        if len(stem_matches) > 1:
            stem_counts = Counter(path.name for path in stem_matches)
            duplicates = ", ".join(sorted(stem_counts))
            raise FileNotFoundError(f"Multiple h5ad files matched stem '{candidate.stem}': {duplicates}")

        raise FileNotFoundError(
            f"Could not resolve h5ad file for value '{file_value}' under '{self._settings.h5ad_base_dir}'."
        )

    def _safe_resolved_h5ad_path(self, parent_row: dict[str, Any]) -> str | None:
        if parent_row["merged"] != "single" or not parent_row["primary_file"]:
            return None
        try:
            return str(self._resolve_h5ad_path(parent_row["primary_file"]))
        except Exception:
            return None

    def _safe_resolved_variant_h5ad_path(self, variant_row: dict[str, Any]) -> str | None:
        if not variant_row.get("file"):
            return None
        try:
            return str(self._resolve_h5ad_path(str(variant_row["file"])))
        except Exception:
            return None

    def get_snapshot(self, force_reload: bool = False) -> dict[str, Any]:
        with self._lock:
            if not self.xlsx_path.exists():
                raise CatalogLoadError(f"Workbook not found: {self.xlsx_path}")

            workbook_mtime = self.xlsx_path.stat().st_mtime
            if (
                not force_reload
                and self._snapshot is not None
                and self._workbook_mtime is not None
                and workbook_mtime == self._workbook_mtime
            ):
                self._cache_hits += 1
                return self._snapshot

            snapshot = self._load_snapshot(workbook_mtime)
            self._snapshot = snapshot
            self._workbook_mtime = workbook_mtime
            self._last_reload_time = snapshot["loaded_at"]
            return snapshot

    def _load_snapshot(self, workbook_mtime: float) -> dict[str, Any]:
        try:
            all_df = pd.read_excel(self.xlsx_path, sheet_name="all")
            multiple_df = pd.read_excel(self.xlsx_path, sheet_name="multiple")
        except Exception as exc:
            raise CatalogLoadError(f"Failed to load workbook '{self.xlsx_path}': {exc}") from exc

        all_df.columns = [_normalize_column_name(column) for column in all_df.columns]
        multiple_df.columns = [_normalize_column_name(column) for column in multiple_df.columns]

        variant_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for index, row in multiple_df.iterrows():
            variant_row = self._build_variant_row(index=index, row=row.to_dict())
            variant_row["resolved_h5ad_path"] = self._safe_resolved_variant_h5ad_path(variant_row)
            publication_key = _normalized_key(variant_row["publication"])
            if publication_key:
                variant_groups[publication_key].append(variant_row)

        parents: list[dict[str, Any]] = []
        details: dict[int, dict[str, Any]] = {}
        multiple_parent_rows = 0

        for index, row in all_df.iterrows():
            parent_row = self._build_parent_row(index=index, row=row.to_dict())
            publication_key = _normalized_key(parent_row["cellxgene_doi"])
            variants = [dict(variant) for variant in variant_groups.get(publication_key, [])]
            parent_row["variant_count"] = len(variants)
            parents.append(parent_row)

            if parent_row["merged"] == "multiple":
                multiple_parent_rows += 1

            details[parent_row["all_excel_row"]] = {
                **_public_parent_record(parent_row),
                "resolved_h5ad_path": self._safe_resolved_h5ad_path(parent_row),
                "publications": sorted(
                    {variant["publication"] for variant in variants if variant.get("publication")}
                ),
                "source_refs": {
                    "all": {
                        "sheet": "all",
                        "excel_row": parent_row["all_excel_row"],
                    },
                    "multiple": [
                        {
                            "sheet": "multiple",
                            "excel_row": variant["multiple_excel_row"],
                        }
                        for variant in variants
                    ],
                },
                "variants": variants,
            }

        facets = self._build_facets(parents)
        mtime_iso = datetime.fromtimestamp(workbook_mtime, tz=timezone.utc).replace(microsecond=0).isoformat().replace(
            "+00:00", "Z"
        )

        return {
            "loaded_at": _utc_now_iso(),
            "source": {
                "xlsx_path": str(self.xlsx_path),
                "mtime": mtime_iso,
                "all_rows": len(parents),
                "multiple_rows": len(multiple_df),
                "multiple_parent_rows": multiple_parent_rows,
            },
            "totals": {
                "datasets": len(parents),
                "variants": len(multiple_df),
                "single_datasets": sum(1 for row in parents if row["merged"] == "single"),
                "multiple_datasets": sum(1 for row in parents if row["merged"] == "multiple"),
            },
            "facets": facets,
            "datasets": [_public_parent_record(row) for row in parents],
            "dataset_details": details,
        }

    @staticmethod
    def _build_parent_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
        normalized = {key: row.get(key) for key in row}
        parent_row: dict[str, Any] = {
            "all_excel_row": int(index) + 2,
            "project": _normalize_text(normalized.get("project")),
            "doi": _normalize_text(normalized.get("doi")),
            "cellxgene_doi": _normalize_text(normalized.get("cellxgene_doi")),
            "title": _normalize_text(normalized.get("title")),
            "author": _normalize_text(normalized.get("author")),
            "year": _coerce_int(normalized.get("year")),
            "journal": _normalize_text(normalized.get("journal")),
            "cell_counts": _coerce_int(normalized.get("cell_counts")),
            "merged": _normalize_text(normalized.get("merged")).lower() or "single",
            "primary_file": _normalize_text(normalized.get("file")),
            "variant_count": 0,
        }
        for source, target in PARENT_MULTI_VALUE_COLUMNS.items():
            parent_row[target] = parse_listish(normalized.get(source))
        return parent_row

    @staticmethod
    def _build_variant_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
        normalized = {key: row.get(key) for key in row}
        variant_row: dict[str, Any] = {
            "multiple_excel_row": int(index) + 2,
            "publication": _normalize_text(normalized.get("publication")),
            "file": _normalize_text(normalized.get("file")),
            "cell_counts": _coerce_int(normalized.get("cell_counts")),
            "description": _normalize_text(normalized.get("description")),
        }
        for source, target in VARIANT_MULTI_VALUE_COLUMNS.items():
            variant_row[target] = parse_listish(normalized.get(source))
        return variant_row

    @staticmethod
    def _build_facets(parents: list[dict[str, Any]]) -> dict[str, list[str]]:
        values: dict[str, set[str]] = {field: set() for field in FACET_FIELDS}
        for row in parents:
            for field in FACET_FIELDS:
                current = row.get(field)
                if isinstance(current, list):
                    values[field].update(item for item in current if item)
                elif current:
                    values[field].add(str(current))

        return {
            "project": sorted(values["project"]),
            "organism": sorted(values["organisms"]),
            "tissue": sorted(values["tissues"]),
            "tissue_type": sorted(values["tissue_types"]),
            "disease": sorted(values["diseases"]),
            "journal": sorted(values["journal"]),
            "merged": sorted(values["merged"]),
        }
