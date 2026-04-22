# GenoPixel agent — scientific single-cell genomics assistant.
# Base: Python 3.11 with uv (matches original GenoPixel Python version).
# Native libs: HDF5 + OpenBLAS required by anndata/numpy; nfs-common for EFS mounting.

FROM ghcr.io/astral-sh/uv:python3.11-bookworm-slim

# --- System dependencies ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    nfs-common \
    libhdf5-dev \
    libopenblas-dev \
    liblapack-dev \
    libxml2-dev \
    libglib2.0-0 \
    libgl1 \
    libfreetype6 \
    libpng16-16 \
    fonts-dejavu-core \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

ENV UV_SYSTEM_PYTHON=1 \
    UV_COMPILE_BYTECODE=1 \
    DOCKER_CONTAINER=1 \
    OTEL_PYTHON_LOG_CORRELATION=true \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg

# --- Python dependencies ---
COPY pyproject.toml .
COPY gateway/ gateway/
COPY tools/ agentcore_tools/

COPY patterns/strands-genopixel-agent/requirements.txt requirements.txt
RUN uv pip install --no-cache -r requirements.txt && \
    uv pip install --no-cache aws-opentelemetry-distro==0.16.0

RUN uv pip install --no-cache -e . --no-deps && \
    uv pip install --no-cache requests>=2.31.0

# --- GenoPixel source files from original Docker/genopixel/ ---
# gp_h5ad_loader, gp_models, gp_runtime_state, gp_plot_response_formatter are
# used directly by gp_tools.py.  gp_scanpy_plotter is NOT copied — the tools
# call scanpy directly to avoid the squidpy / spatial dependency.
COPY Docker/genopixel/gp_h5ad_loader.py .
COPY Docker/genopixel/gp_models.py .
COPY Docker/genopixel/gp_runtime_state.py .
COPY Docker/genopixel/gp_plot_response_formatter.py .

# --- Agent code ---
COPY patterns/strands-genopixel-agent/basic_agent.py .
COPY patterns/strands-genopixel-agent/tools/ tools/
COPY patterns/utils/ utils/

# --- EFS mount entrypoint ---
COPY patterns/strands-genopixel-agent/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

EXPOSE 8080

# Longer start-period: scientific Python imports + possible EFS mount take ~20s on cold start.
HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8080/ping', timeout=4)" || exit 1

ENTRYPOINT ["/entrypoint.sh"]
CMD ["opentelemetry-instrument", "python", "-m", "basic_agent"]
