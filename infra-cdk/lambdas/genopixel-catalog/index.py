"""
GenoPixel Catalog API — stub handler (Phase 2).

Routes:
  GET  /api/catalog                  → list datasets
  GET  /api/catalog/active-dataset   → current active dataset for session
  GET  /api/catalog/{row}            → single dataset detail
  POST /api/catalog/{row}/analyze    → set active dataset for session

Phase 3 will replace these stubs with the full gp_catalog_api.py logic
(Excel metadata parsing, DynamoDB active-dataset persistence, h5ad validation).
"""
import json
import os


CORS_ORIGINS = os.environ.get("CORS_ALLOWED_ORIGINS", "http://localhost:3000")


def _response(status: int, body: dict, origin: str = "") -> dict:
    allowed = CORS_ORIGINS.split(",")
    acao = origin if origin in allowed else allowed[0]
    return {
        "statusCode": status,
        "headers": {
            "Content-Type": "application/json",
            "Access-Control-Allow-Origin": acao,
            "Access-Control-Allow-Headers": "Content-Type,Authorization",
            "Access-Control-Allow-Methods": "GET,POST,OPTIONS",
        },
        "body": json.dumps(body),
    }


def handler(event: dict, context) -> dict:  # type: ignore[type-arg]
    method = event.get("httpMethod", "GET")
    path = event.get("path", "/")
    origin = (event.get("headers") or {}).get("origin", "")

    # OPTIONS preflight
    if method == "OPTIONS":
        return _response(200, {}, origin)

    # GET /api/catalog/active-dataset
    if method == "GET" and path.rstrip("/").endswith("active-dataset"):
        return _response(200, {"active_dataset": None}, origin)

    # GET /api/catalog
    if method == "GET" and path.rstrip("/").endswith("catalog"):
        return _response(
            200,
            {
                "datasets": [],
                "total": 0,
                "message": "Stub — Phase 3 will port full catalog from Excel metadata",
            },
            origin,
        )

    # GET /api/catalog/{row}
    if method == "GET":
        row = path.rstrip("/").split("/")[-1]
        return _response(
            200,
            {"dataset": None, "row": row, "message": "Stub"},
            origin,
        )

    # POST /api/catalog/{row}/analyze
    if method == "POST" and path.rstrip("/").endswith("analyze"):
        return _response(
            200,
            {"success": True, "message": "Stub — Phase 3 will persist active dataset to DynamoDB"},
            origin,
        )

    return _response(404, {"error": "Not found"}, origin)
