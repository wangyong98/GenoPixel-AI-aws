# AI Agent Instructions – GenoPixel + Docker Compose

Purpose and architecture
- Native FastAPI GenoPixel runtime with Scanpy for catalog browsing and plotting.
- Docker Compose stack for Open WebUI and a Caddy reverse proxy.

Key files
- GenoPixel image: Docker/genopixel/Dockerfile
- GenoPixel catalog API: Docker/genopixel/gp_catalog_api.py
- Compose stack: docker-compose.yaml
- Reverse proxy: Caddyfile (serves /assets from repo out/)

Build images
- GenoPixel: built automatically by compose as `genopixel:latest`; rebuild with `docker compose build genopixel_tool_server`

Environment (.env at repo root)
- Required for Open WebUI: `WEBUI_ADMIN_EMAIL`, `WEBUI_ADMIN_PASS`, `WEBUI_ADMIN_NAME`, `WEBUI_HOST=localhost`, `WEBUI_BANNERS=[]`
- Optional: `CHAT_HEALTH_SERV_BASEURL=http://openwebui:8080/health`
- GenoPixel tool server: `GENOPIXEL_TOOL_HOST=genopixel_tool_server`, `GENOPIXEL_TOOL_PORT=18889`, `GENOPIXEL_TOOL_API_KEY=<key>`, `GENOPIXEL_TOOL_HEALTH_SERV_BASEURL=http://genopixel_tool_server:18889/health`

Compose launch
- Start web UI (will also start tool servers due to dependencies): `docker compose up -d caddy openwebui`
- Start everything: `docker compose up -d`
- Access Open WebUI via proxy at http://localhost (port 80) or directly at http://localhost:3000
- Static assets: Caddy serves repo `out/` at `http://localhost/assets/...`

Tool server details
- Health: `POST /health` with `Authorization: Bearer <API_KEY>` and JSON body `{}`.
- OpenAPI: `/openapi.json` (bearer may be required depending on config).
- In Open WebUI, the connection is named "genopixel-tool".

Quick checks (host)
```bash
# GenoPixel server health
curl -X POST -H "Authorization: Bearer ${GENOPIXEL_TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{}' http://localhost:${GENOPIXEL_TOOL_PORT}/health | head -c 200

# Sample GenoPixel tool
curl -X POST -H "Authorization: Bearer ${GENOPIXEL_TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{"plot_type":"umap","color_json":"[\"cell_type\"]","genes_json":"[]"}' \
  http://localhost:${GENOPIXEL_TOOL_PORT}/generate_scanpy_plot | jq
```

Troubleshooting
- WEBUI banners JSON error: If logs show a JSONDecodeError for `WEBUI_BANNERS`, set it to a valid JSON array in `.env` (e.g., `WEBUI_BANNERS=[]`) and restart `openwebui`.
- Tool server unhealthy:
  - Ensure `.env` API keys are set and healthchecks send an authorized POST with `{}`.
  - Interpreting responses:
    - 401 Unauthorized → API key missing/mismatch
    - 422 Unprocessable Content → missing JSON body; add `-H Content-Type` and `-d '{}'`
    - 405 Method Not Allowed → wrong HTTP method; use POST
- Logs
  - `docker compose logs --no-color --tail=100 genopixel_tool_server`
  - `docker compose logs --no-color --tail=100 openwebui`
- Rebuild after changes
  - GenoPixel server: `docker compose build genopixel_tool_server && docker compose up -d`
