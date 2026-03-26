# Agent‑Optimized Workflow

Containerized R and GenoPixel tooling with an optional web UI, orchestrated via Docker Compose. This repo includes:
- A reproducible R image with required CRAN/Bioconductor packages.
- An example MCP stdio server implemented in R.
- A native FastAPI GenoPixel runtime for catalog browsing and plotting.
- A docker-compose stack with Open WebUI and a Caddy reverse proxy.

## Repository Structure
- docker-compose.yaml — Compose stack for web UI and tools
- Docker/
  - r-mcp/
    - Dockerfile — Builds the `r-mcp:latest` image (R base + packages)
    - initialize_r_packages.R — Installs packages listed in `packages.txt`
    - packages.txt — CRAN/Bioconductor packages to include
    - r-mcp.R — MCP stdio server (consumed by the HTTP proxy)
    - btw-mcp-helper.R — Optional BTW MCP helper
  - genopixel/
    - Dockerfile — Builds `genopixel:latest` (GenoPixel runtime + Scanpy + catalog API stack)
    - gp_catalog.py — Excel normalization and cached catalog loader
    - gp_catalog_api.py — Unified FastAPI runtime for the dataset browser and plotting API
    - gp_runtime_state.py — Shared in-memory active dataset state
  - main.sh — Launcher for the R MCP OpenAPI Proxy
- web/genopixel-datasets/ — Static dataset browser UI served through Caddy
- openwebui/actions/ — Open WebUI action sources for manual import
- openwebui/functions/ — Open WebUI filter/action sources for runtime behavior
- openwebui/skills/ — Repo-managed GenoPixel skill artifacts mirrored into Open WebUI
- data/ — Persistent data mounted by services
- docs/ — Documentation assets
- .github/copilot-instructions.md — Maintainer and Copilot guidance
- out/ — Output artifacts written by scripts
- Caddyfile — Reverse proxy config that forwards port 80 to Open WebUI

## Prerequisites
- Docker Desktop (Compose v2) running on your machine
- macOS/Linux/Windows supported

## Quickstart
Create a minimal `.env` at the repo root and bring up the UI. The dataset browser API starts alongside the existing tool servers.

```env
# Open WebUI admin
WEBUI_ADMIN_EMAIL=admin@example.com
WEBUI_ADMIN_PASS=changeme
WEBUI_ADMIN_NAME=Admin
WEBUI_HOST=localhost
WEBUI_BANNERS=[]
CHAT_HEALTH_SERV_BASEURL=http://openwebui:8080/health

# R tool server
TOOL_HOST=r_tool_server
TOOL_PORT=18888
TOOL_API_KEY=local-dev-key
TOOL_HEALTH_SERV_BASEURL=http://r_tool_server:18888/health

# GenoPixel tool server
GENOPIXEL_TOOL_HOST=genopixel_tool_server
GENOPIXEL_TOOL_PORT=18889
GENOPIXEL_TOOL_API_KEY=local-dev-key
GENOPIXEL_TOOL_HEALTH_SERV_BASEURL=http://genopixel_tool_server:18889/health

# GenoPixel runtime settings
GENOPIXEL_METADATA_XLSX=/Volumes/cx10/Single_cell_data_0117_2026/Final_metadata/cellxgene_HCA_final_02182026.xlsx
GENOPIXEL_H5AD_BASE_DIR=/Volumes/cx10/Single_cell_data_0117_2026/cellxgene_final
GENOPIXEL_OUTPUT_DIR=/code/out/genopixel
GENOPIXEL_DEFAULT_BACKED=false
GENOPIXEL_TOOL_MEMORY_LIMIT=65536M
```

Start services and open http://localhost when healthy:

```bash
docker compose up -d caddy openwebui
```

## Build the R Image
Build once to create `r-mcp:latest` used for R-based tooling.

```bash
docker build --pull --rm -t r-mcp:latest Docker/r-mcp
```

Notes
- The image installs system libs and all R packages from Docker/r-mcp/packages.txt.
- Rebuild after editing packages.txt; use `--no-cache` if you changed system libs in the Dockerfile.

## Configure Environment
Copy the snippet below into a `.env` file at the repo root and adjust values. Only the first block is required to bring up Open WebUI. Ensure `WEBUI_BANNERS` is valid JSON (use `[]`).

```env
# Open WebUI admin
WEBUI_ADMIN_EMAIL=admin@example.com
WEBUI_ADMIN_PASS=changeme
WEBUI_ADMIN_NAME=Admin
WEBUI_HOST=localhost
WEBUI_BANNERS=[]

# Optional: Open WebUI health URL used by healthcheck
CHAT_HEALTH_SERV_BASEURL=http://openwebui:8080/health

# Tool server (HTTP proxy over R MCP)
TOOL_HOST=r_tool_server
TOOL_PORT=18888
TOOL_API_KEY=local-dev-key
TOOL_HEALTH_SERV_BASEURL=http://r_tool_server:18888/health

# GenoPixel runtime
GENOPIXEL_TOOL_HOST=genopixel_tool_server
GENOPIXEL_TOOL_PORT=18889
GENOPIXEL_TOOL_API_KEY=local-dev-key
GENOPIXEL_TOOL_HEALTH_SERV_BASEURL=http://genopixel_tool_server:18889/health

# GenoPixel runtime settings
GENOPIXEL_METADATA_XLSX=/Volumes/cx10/Single_cell_data_0117_2026/Final_metadata/cellxgene_HCA_final_02182026.xlsx
GENOPIXEL_H5AD_BASE_DIR=/Volumes/cx10/Single_cell_data_0117_2026/cellxgene_final
GENOPIXEL_OUTPUT_DIR=/code/out/genopixel
GENOPIXEL_DEFAULT_BACKED=false
GENOPIXEL_TOOL_MEMORY_LIMIT=65536M

# Optional integrations (only if you wire up corresponding services)
TOOLUNIVERSE_PORT=9999
OPENAI_API_KEY=
WEBUI_API_KEY=
FILE_USER=
FILE_PASS=
DB_PATH=
STORAGE_PATH=
FILE_SERV_BASEURL=
```

## Launch with Docker Compose

The compose stack defines four services:
- caddy — Reverse proxy on port 80. Also serves static files from `out/` at `/assets/*` and the dataset browser at `/apps/genopixel-datasets/`.
- openwebui — Chat UI (also mapped to host port 3000)
- r_tool_server — MCP OpenAPI Proxy over the R MCP stdio server
- genopixel_tool_server — Unified FastAPI GenoPixel runtime that serves both the dataset catalog and plotting API

Notes
- `r_tool_server` is launched by `main.sh`, which runs the MCP OpenAPI Proxy over the R stdio server (`Rscript Docker/r-mcp/r-mcp.R`).
- Images for external services use pull policy `if_not_present`. The `r_tool_server` image builds locally if missing.
- `GENOPIXEL_TOOL_MEMORY_LIMIT` controls the GenoPixel container limit. The default in this repo is `65536M`.

Start the web UI (will also start tool servers due to dependencies):

```bash
docker compose up -d caddy openwebui
```

Access
- Open WebUI via proxy: http://localhost
- Direct port (bypass proxy): http://localhost:3000

Stop the stack:

```bash
docker compose down
```

## Enable the Tool Server (HTTP)
Start all services. Compose will build the `r-mcp:latest` image if it’s not available locally.

```bash
docker compose up -d
```

Health and tools
- Health check uses `POST /health` with bearer auth and an empty JSON body `{}`.
- OpenAPI spec is served at `/openapi.json`.
- Example tool endpoint: `POST /dice-roll` (schema generated from the MCP tool definition).
- In Open WebUI, the tool connections appear as "r-tool" and "genopixel-tool".
- `genopixel-tool` exposes plotting and basic AnnData metadata inspection: `generate_scanpy_plot` and `print_adata_obs`.

Static artifacts
- Caddy serves the repository `out/` folder at `http://localhost/assets/…`. Tools save images under `out/boxplots/` and return public URLs like `http://localhost/assets/boxplots/<file>.png`.
- Caddy serves the dataset browser at `http://localhost/apps/genopixel-datasets/`.
- The browser reads catalog data from `http://localhost/api/genopixel-catalog/catalog`.
- Runtime readiness is available at `http://localhost/api/genopixel-runtime/active-dataset`.

Dataset browser
- The catalog UI uses the `all` sheet as the parent dataset list and shows `multiple` sheet rows inside the dataset detail drawer.
- Import [openwebui/actions/genopixel_dataset_browser.py](/Users/wangyong98/AI_test/ftowfic-agent-optimized-workflow-ff2a654b03e2/openwebui/actions/genopixel_dataset_browser.py) through Open WebUI Admin > Functions to add the `Browse Datasets` action.
- After importing the action, use it from chat to open the browser inline, or open the direct route above for a full-page view.

Runtime skills for GenoPixels
- Repo-managed GenoPixel skills live under [openwebui/skills/genopixel-tool-usage](/Users/wangyong98/AI_test/ftowfic-agent-optimized-workflow-ff2a654b03e2/openwebui/skills/genopixel-tool-usage) and [openwebui/skills/genopixel-plot-formatting](/Users/wangyong98/AI_test/ftowfic-agent-optimized-workflow-ff2a654b03e2/openwebui/skills/genopixel-plot-formatting).
- The runtime filter source lives at [openwebui/functions/genopixel_skill_injector.py](/Users/wangyong98/AI_test/ftowfic-agent-optimized-workflow-ff2a654b03e2/openwebui/functions/genopixel_skill_injector.py).
- Sync the mirrored skill records, the runtime filter, and the GenoPixels model linkage into Open WebUI with:

```bash
python3 openwebui/skills/scripts/sync_skills_to_openwebui.py
docker compose restart openwebui
```

- The sync script makes the mirrored skills publicly readable in Open WebUI, attaches the filter to the `GenoPixels` model, and clears the model-specific plotting system prompt so runtime behavior comes from the skill artifacts instead.

Quick verification from host
```bash
# Health (requires bearer)
curl -X POST -H "Authorization: Bearer ${TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{}' http://localhost:${TOOL_PORT}/health | head -c 200

# OpenAPI (often requires bearer)
curl -H "Authorization: Bearer ${TOOL_API_KEY}" http://localhost:${TOOL_PORT}/openapi.json | jq '.info, (.paths | keys)[:5]'

# Sample tool (dice-roll on R server)
curl -X POST -H "Authorization: Bearer ${TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{"sides":6,"seed":123,"include_quote":true}' http://localhost:${TOOL_PORT}/dice-roll | jq

# Sample tool (boxplot on R server)
curl -X POST -H "Authorization: Bearer ${TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{"n":100, "groups":4, "seed":42, "title":"Demo from R"}' \
  http://localhost:${TOOL_PORT}/boxplot | jq

# Sample GenoPixel tool
curl -X POST -H "Authorization: Bearer ${GENOPIXEL_TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{"plot_type":"umap","color_json":"[\"cell_type\"]","genes_json":"[]"}' \
  http://localhost:${GENOPIXEL_TOOL_PORT}/generate_scanpy_plot | jq

# Inspect observation columns from active dataset
curl -X POST -H "Authorization: Bearer ${GENOPIXEL_TOOL_API_KEY}" -H "Content-Type: application/json" \
  -d '{"command":"print(adata.obs)"}' \
  http://localhost:${GENOPIXEL_TOOL_PORT}/print_adata_obs | jq
```

## Optional: Run the R MCP Stdio Server
The provided R MCP server speaks stdio (not HTTP). You can run it directly for integration with MCP‑aware clients.

Host R (if you have R locally):

```bash
Rscript Docker/r-mcp/r-mcp.R
```

Dockerized (recommended):

```bash
docker run --rm -i \
  -v "$PWD":/work \
  -w /work/Docker/r-mcp \
  r-mcp:latest \
  Rscript r-mcp.R
```

Client tips
- MCP over stdio requires keeping STDIN open (use `-i` with Docker).
- For GitHub Copilot Chat (MCP), add a custom server with Command `Rscript` and Args pointing to the absolute path of `Docker/r-mcp/r-mcp.R`.

## Developing New R Scripts
- Follow CLI pattern: `Rscript <script.R> <output_path> [args…]`; print a friendly usage() and exit non‑zero on missing args.
- Parse args with `commandArgs(trailingOnly = TRUE)`; set defaults; validate types; seed randomness for reproducibility.
- Suppress package noise: `suppressPackageStartupMessages({ library(pkg) })`.
- Ensure output directories exist (e.g., `dir.create(..., recursive = TRUE)`) and prefer writing into `out/`.

## Package and System Dependencies
- R packages come from Docker/r-mcp/packages.txt and are installed using `BiocManager` (CRAN + Bioconductor). Any install failure fails the image build.
- Add required Debian/Ubuntu libs in Docker/r-mcp/Dockerfile when packages need headers (e.g., libxml2-dev, libcurl4-openssl-dev).

## Known Gaps / Next Steps
- If you plan to add `file_server` or `tool_universe`, define them in docker-compose.yaml and wire them into `openwebui.depends_on`.

## Troubleshooting
- Open WebUI banners JSON error: If logs show a JSONDecodeError for `WEBUI_BANNERS`, set it to a valid JSON array in `.env` (e.g., `WEBUI_BANNERS=[]`) and restart `openwebui`.

- Tool server unhealthy:
  - Ensure `.env` has `TOOL_API_KEY` set and compose healthcheck uses an authorized POST with an empty JSON body `{}`.
  - Quick local checks:
    ```bash
    # Health (authorized POST with empty JSON)
    curl -X POST -H "Authorization: Bearer ${TOOL_API_KEY}" -H "Content-Type: application/json" \
      -d '{}' http://localhost:${TOOL_PORT}/health | head -c 200

    # If 401 Unauthorized → API key missing/mismatch
    # If 422 Unprocessable Content → missing JSON body; add -H Content-Type and -d '{}'
    # If 405 Method Not Allowed → wrong HTTP method; use POST
    ```

- Open WebUI slow to become healthy: First start may run DB migrations. Tail logs and wait:
  ```bash
  docker compose logs --no-color --tail=100 openwebui
  ```

- Port conflicts (80/81/443): If those are busy on your host, edit the host port mappings under `proxymanager.ports` or stop the conflicting service, then restart compose.

- Rebuild after dependency changes: If you change packages or R server code, rebuild and restart:
  ```bash
  docker compose build r_tool_server
  docker compose up -d
  ```
