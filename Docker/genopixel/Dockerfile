FROM python:3.11-slim

# Install curl and build essentials needed by the GenoPixel runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
  curl ca-certificates \
  build-essential \
  libopenblas-dev \
  liblapack-dev \
  libhdf5-dev \
  libxml2-dev \
  libxslt1-dev \
  zlib1g-dev \
  libgl1 \
  libglib2.0-0 \
    libfreetype6 \
    libpng16-16 \
    fonts-dejavu-core \
  && rm -rf /var/lib/apt/lists/*

# Ensure pip is up-to-date and install runtime dependencies
RUN python -m pip install --no-cache-dir -U pip \
  && python -m pip install --no-cache-dir \
    fastapi \
    uvicorn \
    matplotlib \
    numpy \
    pandas \
    openpyxl \
    requests \
    python-dotenv \
    anndata \
    scanpy \
    seaborn

WORKDIR /work
