# AI Agent Instructions – Agent‑Optimized Workflow (R + GenoPixel + Docker Compose)

Purpose and architecture
- Reproducible R runtime with preinstalled CRAN/Bioconductor packages.
- Example MCP stdio servers in R and GenoPixel.
- Docker Compose stack for Open WebUI and a Caddy reverse proxy.

Key files
- R image: Docker/r-mcp/Dockerfile
- R package bootstrap: Docker/r-mcp/initialize_r_packages.R
- R package list: Docker/r-mcp/packages.txt
- R MCP stdio server: Docker/r-mcp/r-mcp.R
- GenoPixel MCP stdio server: Docker/genopixel/server.py
- GenoPixel image: Docker/genopixel/Dockerfile
- Compose stack: docker-compose.yaml
- Reverse proxy: Caddyfile (serves /assets from repo out/)
- Launcher: /main.sh (starts the MCP OpenAPI Proxy over a stdio server)

Build images
- R: `docker build --pull --rm -t r-mcp:latest Docker/r-mcp` (rebuild after changing `packages.txt`; add `--no-cache` if system libs change)
- GenoPixel: built automatically by compose as `genopixel:latest`

Environment (.env at repo root)
- Required for Open WebUI: `WEBUI_ADMIN_EMAIL`, `WEBUI_ADMIN_PASS`, `WEBUI_ADMIN_NAME`, `WEBUI_HOST=localhost`, `WEBUI_BANNERS=[]`
- Optional: `CHAT_HEALTH_SERV_BASEURL=http://openwebui:8080/health`
- R tool server: `TOOL_HOST=r_tool_server`, `TOOL_PORT=18888`, `TOOL_API_KEY=<key>`, `TOOL_HEALTH_SERV_BASEURL=http://r_tool_server:18888/health`
- GenoPixel tool server: `GENOPIXEL_TOOL_HOST=genopixel_tool_server`, `GENOPIXEL_TOOL_PORT=18889`, `GENOPIXEL_TOOL_API_KEY=<key>`, `GENOPIXEL_TOOL_HEALTH_SERV_BASEURL=http://genopixel_tool_server:18889/health`

Compose launch
- Start web UI (will also start tool servers due to dependencies): `docker compose up -d caddy openwebui`
- Start everything: `docker compose up -d`
- Access Open WebUI via proxy at http://localhost (port 80) or directly at http://localhost:3000
- Static assets: Caddy serves repo `out/` at `http://localhost/assets/...`

Tool server details
- `/main.sh` launches the MCP OpenAPI Proxy via `uvx mcpo`, wrapping the stdio MCP server (`SERVER_CMD` defaults to R; Python service sets it to run server.py).
- Health: `POST /health` with `Authorization: Bearer <API_KEY>` and JSON body `{}`.
- OpenAPI: `/openapi.json` (bearer may be required depending on config).
- Example tool paths: `POST /dice-roll`, `POST /boxplot`.
- In Open WebUI, connections are named "r-tool" and "genopixel-tool".

Quick checks (host)
```bash
# R server health
curl -X POST -H "Authorization: Bearer ${TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{}' http://localhost:${TOOL_PORT}/health | head -c 200

# GenoPixel server health
curl -X POST -H "Authorization: Bearer ${GENOPIXEL_TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{}' http://localhost:${GENOPIXEL_TOOL_PORT}/health | head -c 200

# OpenAPI (R)
curl -H "Authorization: Bearer ${TOOL_API_KEY}" http://localhost:${TOOL_PORT}/openapi.json \
  | jq '.info, (.paths | keys)[:5]'

# Tools
curl -X POST -H "Authorization: Bearer ${TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{"sides":6,"seed":123,"include_quote":true}' http://localhost:${TOOL_PORT}/dice-roll | jq
curl -X POST -H "Authorization: Bearer ${TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{"n":100, "groups":4, "seed":42, "title":"Demo from R"}' http://localhost:${TOOL_PORT}/boxplot | jq
```

Run stdio servers directly (optional)
- R (host): `Rscript Docker/r-mcp/r-mcp.R`
- R (docker): `docker run --rm -i -v "$PWD":/work -w /work/Docker/r-mcp r-mcp:latest Rscript r-mcp.R`
- GenoPixel runs under compose; for direct stdio use, run `python Docker/genopixel/server.py` in an env with `mcp` installed.
- For Copilot Chat MCP: add a custom server pointing to the absolute path and appropriate command (e.g., `Rscript` for the R server).

Patterns for new R scripts
- CLI: `Rscript <script.R> <output_path> [args…]`; print `usage()` and exit non‑zero if required args are missing.
- Parse args via `commandArgs(trailingOnly = TRUE)`; set defaults; validate types; seed randomness.
- Suppress package messages with `suppressPackageStartupMessages({ library(pkg) })`.
- Ensure output directory exists (e.g., `dir.create(..., recursive = TRUE)`); write under `out/`.

Package management (R)
- Edit `Docker/r-mcp/packages.txt` (one per line; `#` comments OK).
- Installer uses `BiocManager` for CRAN + Bioconductor; install failures fail the image build.
- Override list path via `PACKAGES_FILE` during build if needed.

System dependencies (R image)
- Add Debian/Ubuntu libs in `Docker/r-mcp/Dockerfile` when R packages need headers (e.g., `libxml2-dev`, `libcurl4-openssl-dev`).
- Base image: `r-base:4.5.1`; keep apt clean‑ups to minimize size.

Conventions
- Avoid absolute host paths inside scripts; container working dir is `/work` when running ad‑hoc.
- Outputs are expected under `out/` on the host; tools save images to `out/boxplots/` and Caddy serves them under `/assets/boxplots/…`.

Troubleshooting
- WEBUI banners JSON error: If logs show a JSONDecodeError for `WEBUI_BANNERS`, set it to a valid JSON array in `.env` (e.g., `WEBUI_BANNERS=[]`) and restart `openwebui`.
- Tool server unhealthy:
  - Ensure `.env` API keys are set and healthchecks send an authorized POST with `{}`.
  - Interpreting responses:
    - 401 Unauthorized → API key missing/mismatch
    - 422 Unprocessable Content → missing JSON body; add `-H Content-Type` and `-d '{}'`
    - 405 Method Not Allowed → wrong HTTP method; use POST
- Logs
  - `docker compose logs --no-color --tail=100 r_tool_server`
  - `docker compose logs --no-color --tail=100 genopixel_tool_server`
  - `docker compose logs --no-color --tail=100 openwebui`
- Rebuild after changes
  - R server or R packages: `docker compose build r_tool_server && docker compose up -d`
  - GenoPixel server: `docker compose build genopixel_tool_server && docker compose up -d`
