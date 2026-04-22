"""
GenoPixel Catalog API — Lambda handler (Phase 3).

Routes:
  GET  /api/catalog                  → list datasets (parsed from Excel in S3)
  GET  /api/catalog/active-dataset   → current active dataset for user session
  GET  /api/catalog/{row}            → single dataset detail
  POST /api/catalog/{row}/analyze    → set active dataset for session (DynamoDB)

Configuration (env vars):
  ACTIVE_DATASET_TABLE   DynamoDB table name  (required)
  H5AD_S3_BUCKET         S3 bucket for h5ad files and metadata Excel  (required)
  METADATA_XLSX_S3_KEY   S3 key of the metadata Excel file
                         (default: metadata/metadata.xlsx)
  CORS_ALLOWED_ORIGINS   comma-separated list of allowed CORS origins
"""
from __future__ import annotations

import ast
import io
import json
import os
import time
from collections import defaultdict
from datetime import datetime, timezone
from typing import Any

import boto3
import pandas as pd

# ─── Configuration ─────────────────────────────────────────────────────────────

ACTIVE_DATASET_TABLE = os.environ["ACTIVE_DATASET_TABLE"]
H5AD_S3_BUCKET = os.environ.get("H5AD_S3_BUCKET", "")
METADATA_XLSX_S3_KEY = os.environ.get("METADATA_XLSX_S3_KEY", "metadata/metadata.xlsx")
CORS_ORIGINS = os.environ.get("CORS_ALLOWED_ORIGINS", "http://localhost:3000")

SESSION_TTL_SECONDS = 24 * 3600

# ─── AWS clients (module-level for Lambda container reuse) ─────────────────────

_s3 = boto3.client("s3")
_dynamodb = boto3.resource("dynamodb")
_table = _dynamodb.Table(ACTIVE_DATASET_TABLE)

# ─── In-memory catalog cache (reused across warm Lambda invocations) ───────────

_CATALOG_CACHE: dict[str, Any] | None = None
_CATALOG_CACHE_ETAG: str | None = None

# ─── Column-mapping constants (mirrors gp_catalog.py) ─────────────────────────

_PARENT_MULTI_VALUE_COLUMNS = {
    "tissue": "tissues",
    "tissue_type": "tissue_types",
    "disease": "diseases",
    "organism": "organisms",
}

_VARIANT_MULTI_VALUE_COLUMNS = {
    "tissue": "tissues",
    "tissue_type": "tissue_types",
    "disease": "diseases",
    "organism": "organisms",
}

_FACET_FIELDS = (
    "project",
    "organisms",
    "tissues",
    "tissue_types",
    "diseases",
    "journal",
    "merged",
)


# ─── Parsing helpers ───────────────────────────────────────────────────────────

def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _normalize_column_name(value: object) -> str:
    return str(value).strip().lower()


def _normalize_text(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()


def _normalized_key(value: object) -> str:
    return _normalize_text(value).lower()


def _coerce_int(value: object) -> int | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except Exception:
        pass
    try:
        return int(value)
    except Exception:
        return None


def _parse_listish(value: object) -> list[str]:
    """Parse a cell value that may be a Python list literal, semicolon-list, or plain string."""
    raw = _normalize_text(value)
    if not raw:
        return []
    if raw.startswith("["):
        try:
            parsed = ast.literal_eval(raw)
        except Exception:
            parsed = None
        if isinstance(parsed, list):
            return [item for item in (_normalize_text(e) for e in parsed) if item]
        if parsed is not None:
            scalar = _normalize_text(parsed)
            return [scalar] if scalar else []
    if ";" in raw:
        return [item for item in (_normalize_text(p) for p in raw.split(";")) if item]
    return [raw]


# ─── Row builders ──────────────────────────────────────────────────────────────

def _build_parent_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
    parent_row: dict[str, Any] = {
        "all_excel_row": index + 2,
        "project": _normalize_text(row.get("project")),
        "doi": _normalize_text(row.get("doi")),
        "cellxgene_doi": _normalize_text(row.get("cellxgene_doi")),
        "title": _normalize_text(row.get("title")),
        "author": _normalize_text(row.get("author")),
        "year": _coerce_int(row.get("year")),
        "journal": _normalize_text(row.get("journal")),
        "cell_counts": _coerce_int(row.get("cell_counts")),
        "merged": _normalize_text(row.get("merged")).lower() or "single",
        "primary_file": _normalize_text(row.get("file")),
        "variant_count": 0,
    }
    for source, target in _PARENT_MULTI_VALUE_COLUMNS.items():
        parent_row[target] = _parse_listish(row.get(source))
    return parent_row


def _build_variant_row(index: int, row: dict[str, Any]) -> dict[str, Any]:
    variant_row: dict[str, Any] = {
        "multiple_excel_row": index + 2,
        "publication": _normalize_text(row.get("publication")),
        "file": _normalize_text(row.get("file")),
        "cell_counts": _coerce_int(row.get("cell_counts")),
        "description": _normalize_text(row.get("description")),
    }
    for source, target in _VARIANT_MULTI_VALUE_COLUMNS.items():
        variant_row[target] = _parse_listish(row.get(source))
    return variant_row


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


def _build_facets(parents: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    """Return {facet_name: {value: count}} dicts for the frontend FacetPanel."""
    counts: dict[str, dict[str, int]] = {field: {} for field in _FACET_FIELDS}
    for row in parents:
        for field in _FACET_FIELDS:
            current = row.get(field)
            if isinstance(current, list):
                for item in current:
                    if item:
                        counts[field][item] = counts[field].get(item, 0) + 1
            elif current:
                val = str(current)
                counts[field][val] = counts[field].get(val, 0) + 1
    return {
        "project": counts["project"],
        "organism": counts["organisms"],
        "tissue": counts["tissues"],
        "tissue_type": counts["tissue_types"],
        "disease": counts["diseases"],
        "journal": counts["journal"],
        "merged": counts["merged"],
    }


# ─── Catalog loading: S3 download + ETag-based in-memory cache ────────────────

def _load_catalog() -> dict[str, Any]:
    """Download metadata Excel from S3, parse it, and cache by ETag."""
    global _CATALOG_CACHE, _CATALOG_CACHE_ETAG

    if not H5AD_S3_BUCKET:
        raise RuntimeError("H5AD_S3_BUCKET environment variable is not set")

    # Check ETag first — avoids re-download on warm Lambda invocations
    try:
        head = _s3.head_object(Bucket=H5AD_S3_BUCKET, Key=METADATA_XLSX_S3_KEY)
    except _s3.exceptions.ClientError as exc:  # type: ignore[attr-defined]
        code = exc.response.get("Error", {}).get("Code", "")
        if code in ("404", "NoSuchKey"):
            raise RuntimeError(
                f"Metadata Excel not found at s3://{H5AD_S3_BUCKET}/{METADATA_XLSX_S3_KEY}. "
                "Upload the workbook to S3 (Phase 6) to enable the catalog."
            ) from exc
        raise RuntimeError(f"S3 error checking metadata: {exc}") from exc

    etag = head.get("ETag", "")
    if _CATALOG_CACHE is not None and etag == _CATALOG_CACHE_ETAG:
        return _CATALOG_CACHE

    # Download into memory (avoids /tmp write latency)
    obj = _s3.get_object(Bucket=H5AD_S3_BUCKET, Key=METADATA_XLSX_S3_KEY)
    xlsx_bytes = obj["Body"].read()
    xlsx_io = io.BytesIO(xlsx_bytes)

    all_df = pd.read_excel(xlsx_io, sheet_name="all")
    xlsx_io.seek(0)
    multiple_df = pd.read_excel(xlsx_io, sheet_name="multiple")

    all_df.columns = [_normalize_column_name(c) for c in all_df.columns]
    multiple_df.columns = [_normalize_column_name(c) for c in multiple_df.columns]

    # Build variant lookup keyed by cellxgene_doi (publication key)
    variant_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for idx, row in multiple_df.iterrows():
        variant_row = _build_variant_row(index=int(idx), row=row.to_dict())
        pub_key = _normalized_key(variant_row["publication"])
        if pub_key:
            variant_groups[pub_key].append(variant_row)

    parents: list[dict[str, Any]] = []
    details: dict[int, dict[str, Any]] = {}

    for idx, row in all_df.iterrows():
        parent_row = _build_parent_row(index=int(idx), row=row.to_dict())
        pub_key = _normalized_key(parent_row["cellxgene_doi"])
        variants = list(variant_groups.get(pub_key, []))
        parent_row["variant_count"] = len(variants)
        parents.append(parent_row)

        details[parent_row["all_excel_row"]] = {
            **_public_parent_record(parent_row),
            "publications": sorted(
                {v["publication"] for v in variants if v.get("publication")}
            ),
            "source_refs": {
                "all": {"sheet": "all", "excel_row": parent_row["all_excel_row"]},
                "multiple": [
                    {"sheet": "multiple", "excel_row": v["multiple_excel_row"]}
                    for v in variants
                ],
            },
            "variants": variants,
        }

    last_modified = head.get("LastModified")
    mtime_iso = (
        last_modified.isoformat().replace("+00:00", "Z")
        if last_modified
        else _utc_now_iso()
    )

    snapshot: dict[str, Any] = {
        "loaded_at": _utc_now_iso(),
        "source": {
            "s3_bucket": H5AD_S3_BUCKET,
            "s3_key": METADATA_XLSX_S3_KEY,
            "etag": etag,
            "mtime": mtime_iso,
            "all_rows": len(parents),
            "multiple_rows": len(multiple_df),
        },
        "totals": {
            "datasets": len(parents),
            "variants": len(multiple_df),
            "single_datasets": sum(1 for r in parents if r["merged"] == "single"),
            "multiple_datasets": sum(1 for r in parents if r["merged"] == "multiple"),
        },
        "facets": _build_facets(parents),
        "datasets": [_public_parent_record(r) for r in parents],
        "dataset_details": details,
    }

    _CATALOG_CACHE = snapshot
    _CATALOG_CACHE_ETAG = etag
    return snapshot


# ─── DynamoDB helpers for active-dataset state ────────────────────────────────

def _get_active_dataset(session_id: str) -> dict[str, Any] | None:
    resp = _table.get_item(Key={"sessionId": session_id})
    item = resp.get("Item")
    if not item:
        return None
    return {k: v for k, v in item.items() if k not in ("sessionId", "ttl")}


def _set_active_dataset(session_id: str, payload: dict[str, Any]) -> None:
    _table.put_item(
        Item={
            "sessionId": session_id,
            "ttl": int(time.time()) + SESSION_TTL_SECONDS,
            **payload,
        }
    )


# ─── HTTP helpers ──────────────────────────────────────────────────────────────

def _response(status: int, body: dict, origin: str = "") -> dict:
    allowed = [o.strip() for o in CORS_ORIGINS.split(",")]
    acao = origin if origin in allowed else allowed[0]
    return {
        "statusCode": status,
        "headers": {
            "Content-Type": "application/json",
            "Access-Control-Allow-Origin": acao,
            "Access-Control-Allow-Headers": "Content-Type,Authorization",
            "Access-Control-Allow-Methods": "GET,POST,OPTIONS",
        },
        "body": json.dumps(body, default=str),
    }


def _session_id(event: dict) -> str:
    """Extract Cognito user sub from API Gateway requestContext."""
    claims = (
        (event.get("requestContext") or {})
        .get("authorizer", {})
        .get("claims", {})
    )
    return claims.get("sub") or claims.get("cognito:username") or "anonymous"


# ─── Lambda entry point ────────────────────────────────────────────────────────

def _filter_datasets(
    datasets: list[dict[str, Any]],
    qp: dict[str, str],
) -> list[dict[str, Any]]:
    """Apply search + facet filters from query-string params."""
    search = qp.get("search", "").strip().lower()
    organism = qp.get("organism", "").strip().lower()
    tissue = qp.get("tissue", "").strip().lower()
    tissue_type = qp.get("tissue_type", "").strip().lower()
    disease = qp.get("disease", "").strip().lower()
    project = qp.get("project", "").strip().lower()
    merged_filter = qp.get("merged", "").strip().lower()  # "true" | "false" | ""

    result = []
    for ds in datasets:
        if search:
            haystack = " ".join([
                ds.get("title", ""),
                ds.get("author", ""),
                ds.get("journal", ""),
                " ".join(ds.get("tissues", [])),
                " ".join(ds.get("diseases", [])),
                " ".join(ds.get("organisms", [])),
            ]).lower()
            if search not in haystack:
                continue
        if organism and not any(organism in v.lower() for v in ds.get("organisms", [])):
            continue
        if tissue and not any(tissue in v.lower() for v in ds.get("tissues", [])):
            continue
        if tissue_type and not any(tissue_type in v.lower() for v in ds.get("tissue_types", [])):
            continue
        if disease and not any(disease in v.lower() for v in ds.get("diseases", [])):
            continue
        if project and project != (ds.get("project") or "").lower():
            continue
        if merged_filter == "true" and ds.get("merged") != "multiple":
            continue
        if merged_filter == "false" and ds.get("merged") == "multiple":
            continue
        result.append(ds)
    return result


def handler(event: dict, context: Any) -> dict:  # noqa: ANN401
    method = event.get("httpMethod", "GET")
    path = event.get("path", "/")
    origin = (event.get("headers") or {}).get("origin", "")
    qp: dict[str, str] = event.get("queryStringParameters") or {}

    if method == "OPTIONS":
        return _response(200, {}, origin)

    session_id = _session_id(event)

    try:
        # GET /api/catalog/active-dataset
        if method == "GET" and path.rstrip("/").endswith("active-dataset"):
            active = _get_active_dataset(session_id)
            if active is None:
                return _response(200, {"loaded": False}, origin)
            return _response(200, {"loaded": True, **active}, origin)

        # Route on path segments (after stripping leading slash)
        path_parts = [p for p in path.split("/") if p]

        # GET /api/catalog  — last segment is "catalog" or path is empty
        if method == "GET" and (not path_parts or path_parts[-1] == "catalog"):
            catalog = _load_catalog()
            all_datasets = catalog["datasets"]

            filtered = _filter_datasets(all_datasets, qp)

            # Pagination
            try:
                page = max(1, int(qp.get("page", "1")))
            except ValueError:
                page = 1
            try:
                page_size = min(200, max(1, int(qp.get("page_size", "20"))))
            except ValueError:
                page_size = 20

            total = len(filtered)
            start = (page - 1) * page_size
            page_datasets = filtered[start: start + page_size]

            # Rebuild facets from the full filtered set (not just the current page)
            facets = catalog["facets"] if not any(qp.values()) else _build_facets(filtered)

            return _response(
                200,
                {
                    "datasets": page_datasets,
                    "total": total,
                    "page": page,
                    "page_size": page_size,
                    "facets": facets,
                },
                origin,
            )

        # POST /api/catalog/{row}/analyze
        if method == "POST" and path.rstrip("/").endswith("analyze"):
            row_str = path.rstrip("/").split("/")[-2]
            try:
                row = int(row_str)
            except ValueError:
                return _response(400, {"error": f"Invalid row number: {row_str!r}"}, origin)

            catalog = _load_catalog()
            details = catalog["dataset_details"]
            if row not in details:
                return _response(404, {"error": f"Dataset row {row} not found"}, origin)

            detail = details[row]

            body: dict[str, Any] = {}
            if event.get("body"):
                try:
                    body = json.loads(event["body"])
                except Exception:
                    pass

            multiple_excel_row: int | None = body.get("multiple_excel_row")
            if (
                detail.get("merged") == "multiple"
                and multiple_excel_row is None
                and detail.get("variants")
            ):
                return _response(
                    400,
                    {
                        "error": (
                            f"Dataset row {row} is a multi-variant dataset. "
                            "Provide multiple_excel_row in the request body."
                        )
                    },
                    origin,
                )

            payload: dict[str, Any] = {
                "all_excel_row": row,
                "title": detail.get("title") or f"Dataset {row}",
                "merged": detail.get("merged"),
                "primary_file": detail.get("primary_file", ""),
                "selected_at": _utc_now_iso(),
            }
            if multiple_excel_row is not None:
                payload["multiple_excel_row"] = multiple_excel_row
                for variant in detail.get("variants", []):
                    if variant.get("multiple_excel_row") == multiple_excel_row:
                        payload["primary_file"] = variant.get("file", payload["primary_file"])
                        break

            _set_active_dataset(session_id, payload)
            return _response(200, {"loaded": True, **payload}, origin)

        # GET /api/catalog/{row}
        if method == "GET":
            row_str = path_parts[-1] if path_parts else ""
            try:
                row = int(row_str)
            except ValueError:
                return _response(400, {"error": f"Invalid row number: {row_str!r}"}, origin)

            catalog = _load_catalog()
            details = catalog["dataset_details"]
            if row not in details:
                return _response(404, {"error": f"Dataset row {row} not found"}, origin)

            return _response(200, details[row], origin)

        return _response(404, {"error": "Not found"}, origin)

    except RuntimeError as exc:
        # Expected operational errors (file not in S3 yet, env not configured, etc.)
        return _response(503, {"error": str(exc)}, origin)
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR [{type(exc).__name__}]: {exc}")
        return _response(500, {"error": "Internal server error"}, origin)
