"""
GenoPixel Strands tools — h5ad loading, dataset inspection, and plot generation.

These tools wrap the in-process AnnData/Scanpy stack.  The RUNTIME_STATE singleton
holds the currently-loaded AnnData object across requests within the same container
instance (warm start reuse).

h5ad resolution order:
  1. H5AD_BASE_DIR (EFS path, default /mnt/genopixel/h5ad)
  2. /tmp/genopixel/h5ad/<filename>  (cached S3 download)
  3. Download from S3 bucket H5AD_S3_BUCKET → /tmp/genopixel/h5ad/<filename>
"""
from __future__ import annotations

import base64
from datetime import datetime
import io
import json
import logging
import os
from pathlib import Path
from typing import Any

import boto3
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from strands import tool

matplotlib.use("Agg")
sc.settings.set_figure_params(dpi=150, facecolor="white", frameon=False)

# These modules are copied from Docker/genopixel/ into /app at build time.
from gp_h5ad_loader import load_h5ad  # type: ignore[import-untyped]
from gp_runtime_state import NoActiveDatasetError, RUNTIME_STATE  # type: ignore[import-untyped]

logger = logging.getLogger(__name__)

# ─── Gene name resolution (ported from GenoPixel-AI gp_scanpy_plotter.py) ─────

GENE_SYNONYM_MAP: dict[str, str] = {
    "TL1A": "TNFSF15",
    "TLA1": "TNFSF15",
    "DR3": "TNFRSF25",
    "DDR3": "TNFRSF25",
    "TNFALPHA": "TNF",
    "TNF-A": "TNF",
    "TNFA": "TNF",
    "IL-1B": "IL1B",
    "IL1-BETA": "IL1B",
    "IL-6": "IL6",
    "IL-8": "CXCL8",
    "CXCL-8": "CXCL8",
    "MCP1": "CCL2",
    "MCP-1": "CCL2",
    "RANTES": "CCL5",
    "MIP1A": "CCL3",
    "MIP-1A": "CCL3",
    "MIP1B": "CCL4",
    "MIP-1B": "CCL4",
    "IP10": "CXCL10",
    "IP-10": "CXCL10",
    "MIG": "CXCL9",
    "TARC": "CCL17",
    "MDC": "CCL22",
    "SCF": "KITLG",
    "FLT3L": "FLT3LG",
    "TRAIL": "TNFSF10",
    "BAFF": "TNFSF13B",
    "APRIL": "TNFSF13",
    "LIGHT": "TNFSF14",
    "4-1BB": "TNFRSF9",
    "4-1BBL": "TNFSF9",
    "CD137": "TNFRSF9",
    "CD137L": "TNFSF9",
    "OX40": "TNFRSF4",
    "OX40L": "TNFSF4",
    "ICOSL": "ICOSLG",
    "PD-L1": "CD274",
    "PDL1": "CD274",
    "PD-L2": "PDCD1LG2",
    "PDL2": "PDCD1LG2",
    "PD-1": "PDCD1",
    "TIM3": "HAVCR2",
    "TIM-3": "HAVCR2",
    "LAG3": "LAG3",
    "CTLA4": "CTLA4",
    "TGF-BETA1": "TGFB1",
    "TGF-β1": "TGFB1",
    "IFN-GAMMA": "IFNG",
    "IFN-γ": "IFNG",
    "GM-CSF": "CSF2",
    "G-CSF": "CSF3",
    "M-CSF": "CSF1",
}

_GENE_SYMBOL_COLUMNS = (
    "feature_name",
    "gene_name",
    "gene_names",
    "gene_symbols",
    "symbol",
    "hgnc_symbol",
)


def _find_gene_symbol_column(adata: Any) -> str | None:
    """Return the first adata.var column that stores human-readable gene symbols, or None."""
    for col in _GENE_SYMBOL_COLUMNS:
        if col in adata.var.columns:
            return col
    return None


def _has_feature_name_column(adata: Any) -> bool:
    return _find_gene_symbol_column(adata) is not None


def _gene_to_ens(df: pd.DataFrame, genes: list[str], column: str = "feature_name") -> list[str]:
    """Return var_names (Ensembl IDs) matching the given gene symbols via a var column."""
    if column not in df.columns:
        return []
    normalized_genes = [str(g).strip().lower() for g in genes if str(g).strip()]
    if not normalized_genes:
        return []
    normalized_set = set(normalized_genes)
    col_values = df[column].astype(str).str.strip().str.lower()
    return [str(item) for item in df.index[col_values.isin(normalized_set)].tolist()]


def _canonicalize_gene_tokens(genes: list[str]) -> list[str]:
    """Apply GENE_SYNONYM_MAP and deduplicate (case-insensitive)."""
    canonical: list[str] = []
    seen: set[str] = set()
    for gene in genes:
        token = str(gene).strip()
        if not token:
            continue
        mapped = GENE_SYNONYM_MAP.get(token.upper(), token)
        lowered = mapped.lower()
        if lowered not in seen:
            canonical.append(mapped)
            seen.add(lowered)
    return canonical


def _resolve_gene_names(adata: Any, genes: list[str]) -> list[str]:
    """
    Multi-tier gene name resolution matching GenoPixel-AI behaviour:
    1. Apply GENE_SYNONYM_MAP (e.g. "TNF-A" → "TNF", "PD-L1" → "CD274")
    2. If dataset has a gene-symbol column (feature_name / hgnc_symbol / …):
       translate user symbols → Ensembl var_names via case-insensitive lookup
    3. Fallback: direct case-insensitive var_names match
    4. Last resort: keep original token (lets scanpy produce a clear error)
    """
    if not genes:
        return []

    canonical = _canonicalize_gene_tokens(genes)

    if not _has_feature_name_column(adata):
        # No symbol column — deduplicate and return as-is
        seen: set[str] = set()
        out: list[str] = []
        for g in canonical:
            if g not in seen:
                out.append(g)
                seen.add(g)
        return out

    sym_col = _find_gene_symbol_column(adata) or "feature_name"
    resolved: list[str] = []
    seen_ids: set[str] = set()

    # Primary: translate symbols → Ensembl IDs via var column
    resolved_ids = _gene_to_ens(adata.var, canonical, column=sym_col)
    for gid in resolved_ids:
        if gid not in seen_ids:
            resolved.append(gid)
            seen_ids.add(gid)

    # Fallback: direct var_names match (case-insensitive)
    available_lower = {str(n).lower() for n in adata.var_names}
    for gene in canonical:
        if gene.lower() in available_lower and gene not in seen_ids:
            resolved.append(gene)
            seen_ids.add(gene)
        elif gene.lower() not in available_lower and gene not in seen_ids and not resolved_ids:
            # Last resort: pass original through (scanpy will give a clear error)
            resolved.append(gene)
            seen_ids.add(gene)

    return resolved


def _resolve_color_tokens(adata: Any, tokens: list[str]) -> list[str]:
    """
    Resolve a mixed list of obs-column names and gene symbols.
    Obs columns are kept verbatim; non-obs tokens are resolved via _resolve_gene_names.
    """
    if not tokens:
        return tokens
    obs_cols = set(adata.obs.columns)
    gene_tokens = [t for t in tokens if t not in obs_cols]
    resolved_genes = _resolve_gene_names(adata, gene_tokens) if gene_tokens else []

    # Rebuild in original order, substituting resolved genes
    gene_iter = iter(resolved_genes)
    result: list[str] = []
    gene_map: dict[str, str] = {}
    gi = 0
    for orig, resolved in zip(gene_tokens, resolved_genes):
        gene_map[orig] = resolved

    for token in tokens:
        if token in obs_cols:
            result.append(token)
        else:
            result.append(gene_map.get(token, token))
    return result


H5AD_BASE_DIR = Path(os.environ.get("H5AD_BASE_DIR", "/mnt/genopixel/h5ad"))
OUTPUT_DIR = Path(os.environ.get("OUTPUT_DIR", "/tmp/genopixel/out"))
H5AD_S3_BUCKET = os.environ.get("H5AD_S3_BUCKET", "")
_TMP_H5AD = Path("/tmp/genopixel/h5ad")


# ─── Internal helpers ─────────────────────────────────────────────────────────

def _is_truthy(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes", "on"}


def _is_under_base(path: Path, base: Path) -> bool:
    try:
        path.resolve().relative_to(base.resolve())
        return True
    except Exception:
        return False


def _should_use_backed(path: Path) -> bool:
    if _is_under_base(path, H5AD_BASE_DIR):
        return True
    # For S3 fallback files in /tmp, default to backed mode to avoid RAM spikes.
    return _is_truthy(os.environ.get("S3_FALLBACK_BACKED", "true"))


def _resolve_h5ad_path(filename: str) -> Path:
    """Resolve h5ad filename → local Path; download from S3 if needed."""
    candidate = Path(filename)

    # Absolute path that already exists
    if candidate.is_absolute() and candidate.exists():
        return candidate

    # Relative to EFS base dir
    for base in (H5AD_BASE_DIR, _TMP_H5AD):
        full = base / candidate.name
        if full.exists():
            return full
        full2 = base / filename
        if full2.exists():
            return full2

    # Fuzzy match in EFS base dir
    if H5AD_BASE_DIR.exists():
        name_l = candidate.name.lower()
        stem_l = candidate.stem.lower()
        hits = [p for p in H5AD_BASE_DIR.rglob("*.h5ad") if p.name.lower() == name_l]
        if len(hits) == 1:
            return hits[0]
        if not hits:
            hits = [p for p in H5AD_BASE_DIR.rglob("*.h5ad") if stem_l in p.stem.lower()]
            if len(hits) == 1:
                return hits[0]

    # S3 fallback
    if not H5AD_S3_BUCKET:
        raise FileNotFoundError(
            f"h5ad '{filename}' not found on EFS at {H5AD_BASE_DIR} "
            "and H5AD_S3_BUCKET is not configured. "
            "Upload the file to S3 or mount EFS first (Phase 6)."
        )

    fname = candidate.name
    tmp_path = _TMP_H5AD / fname
    if tmp_path.exists():
        logger.info("Using cached S3 download: %s", tmp_path)
        return tmp_path

    # Try both s3://<bucket>/h5ad/<file> and s3://<bucket>/<file>
    s3_client = boto3.client("s3")
    s3_key = f"h5ad/{fname}"
    try:
        s3_client.head_object(Bucket=H5AD_S3_BUCKET, Key=s3_key)
    except Exception:
        # Fall back to flat key (no prefix) in case files were uploaded without prefix
        s3_key = fname

    logger.info("Downloading s3://%s/%s → %s", H5AD_S3_BUCKET, s3_key, tmp_path)
    tmp_path.parent.mkdir(parents=True, exist_ok=True)
    s3_client.download_file(H5AD_S3_BUCKET, s3_key, str(tmp_path))
    return tmp_path


def _fig_to_markdown(title: str = "") -> str:
    """Capture current matplotlib figure as base64-encoded PNG markdown."""
    buf = io.BytesIO()
    plt.savefig(buf, dpi=150, bbox_inches="tight", format="png")
    plt.close("all")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    return f"![{title}](data:image/png;base64,{b64})"


def _output_dir() -> Path:
    out = OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)
    return out


# ─── Strands tools ─────────────────────────────────────────────────────────────

@tool
def get_active_dataset_info() -> str:
    """
    Return metadata about the dataset currently loaded in memory.

    Returns JSON with: loaded (bool), title, h5ad_path, total_cells, n_genes,
    obs_columns (list), embeddings (list of keys in adata.obsm).
    If no dataset is loaded, returns loaded=false with a helpful message.
    """
    payload = RUNTIME_STATE.get_active_dataset_payload()
    if not payload.get("loaded"):
        pending = RUNTIME_STATE.get_pending_selection()
        if pending and pending.get("primary_file"):
            return json.dumps({
                "loaded": False,
                "pending_selection": pending,
                "message": (
                    f"Dataset '{pending['title']}' is selected but not yet loaded into memory. "
                    f"Call load_dataset(h5ad_filename='{pending['primary_file']}', "
                    f"all_excel_row={pending['all_excel_row']}, "
                    f"title='{pending['title']}') to load it now."
                ),
            })
        return json.dumps({
            "loaded": False,
            "message": (
                "No dataset is loaded. Ask the user to select a dataset from the "
                "GenoPixel browser (Dataset tab), or call load_dataset with the h5ad filename."
            ),
        })

    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
        obs_cols = list(adata.obs.columns)
        embeddings = list(adata.obsm.keys())
        n_vars = int(adata.n_vars)
    except Exception:
        obs_cols = []
        embeddings = []
        n_vars = 0

    return json.dumps({
        "loaded": True,
        "title": payload.get("title"),
        "h5ad_path": payload.get("h5ad_path"),
        "total_cells": payload.get("total_cells"),
        "n_genes": n_vars,
        "obs_columns": obs_cols,
        "embeddings": embeddings,
        "all_excel_row": payload.get("all_excel_row"),
    })


@tool
def load_dataset(
    h5ad_filename: str,
    all_excel_row: int,
    title: str = "",
    multiple_excel_row: int = 0,
) -> str:
    """
    Load an h5ad dataset into memory for analysis.

    Resolves the file from EFS (/mnt/genopixel/h5ad) first, then falls back
    to downloading from the S3 data bucket. Uses backed=True for EFS files and
    (by default) for S3 fallback files in /tmp to reduce memory pressure.
    Set S3_FALLBACK_BACKED=false to force full in-memory load for S3 files.

    Args:
        h5ad_filename: Filename or relative path of the h5ad file (e.g. "dataset.h5ad")
        all_excel_row: Row number from the catalog Excel 'all' sheet (for tracking)
        title: Human-readable name shown in responses
        multiple_excel_row: Sub-dataset row from 'multiple' sheet; 0 = single dataset

    Returns:
        JSON with dataset summary or an error message.
    """
    try:
        path = _resolve_h5ad_path(h5ad_filename)
    except FileNotFoundError as exc:
        return json.dumps({"error": str(exc)})

    # EFS always uses backed mode; S3 fallback defaults to backed for stability.
    backed = _should_use_backed(path)

    try:
        payload = RUNTIME_STATE.load_active_dataset(
            h5ad_path=str(path),
            all_excel_row=all_excel_row,
            multiple_excel_row=multiple_excel_row if multiple_excel_row else None,
            title=title or h5ad_filename,
            backed=backed,
            force_reload=False,
        )
    except Exception as exc:
        return json.dumps({"error": f"Failed to load dataset: {exc}"})

    adata, _ = RUNTIME_STATE.require_active_adata()
    return json.dumps({
        "loaded": True,
        "title": payload["title"],
        "total_cells": payload["total_cells"],
        "n_genes": int(adata.n_vars),
        "obs_columns": list(adata.obs.columns),
        "embeddings": list(adata.obsm.keys()),
        "backed": payload["backed"],
    })


@tool
def get_obs_columns() -> str:
    """
    List all observation (cell metadata) columns in the loaded dataset.

    For categorical columns with ≤50 unique values, also returns the value list.
    For high-cardinality columns, returns the count of unique values only.

    Returns:
        JSON mapping column_name → {dtype, n_unique, values}
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    info: dict[str, Any] = {}
    for col in adata.obs.columns:
        dtype = str(adata.obs[col].dtype)
        n_unique = int(adata.obs[col].nunique())
        values: Any = f"{n_unique} unique values"
        if n_unique <= 50:
            values = sorted(adata.obs[col].dropna().astype(str).unique().tolist())
        info[col] = {"dtype": dtype, "n_unique": n_unique, "values": values}
    return json.dumps(info, default=str)


@tool
def get_obs_column_values(column: str) -> str:
    """
    Get all unique values and their cell counts for a specific obs column.

    Args:
        column: Name of the obs column (e.g. 'author_cell_type', 'disease', 'tissue')

    Returns:
        JSON mapping value → cell_count (top 200, sorted by count descending)
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    if column not in adata.obs.columns:
        avail = list(adata.obs.columns)
        return json.dumps({"error": f"Column '{column}' not found.", "available": avail})

    counts = adata.obs[column].value_counts().head(200)
    return json.dumps({str(k): int(v) for k, v in counts.items()})


@tool
def set_session_markers(genes: str) -> str:
    """
    Store a default list of marker genes for this session.
    These are used automatically in dot plots, heatmaps, and violin plots
    when no genes are specified explicitly.

    Args:
        genes: Comma-separated gene names (e.g. "CD3E,CD4,CD8A,FOXP3,IL2RA")

    Returns:
        Confirmation with the resolved gene list.
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    requested = [g.strip() for g in genes.split(",") if g.strip()]
    resolved = _resolve_gene_names(adata, requested)

    # Check which resolved names are actually in var_names
    var_names_set = set(adata.var_names)
    found = [g for g in resolved if g in var_names_set]
    not_found = [g for g in requested if not any(
        r in var_names_set for r in _resolve_gene_names(adata, [g])
    )]

    # Store on a module-level dict keyed by a sentinel (single-user per container)
    _SESSION_STATE["markers"] = found

    msg = f"Markers set: {', '.join(found)}"
    if not_found:
        msg += f"\nNot found in dataset: {', '.join(not_found)}"
    return msg


# Module-level session state (persists across warm invocations in same container)
_SESSION_STATE: dict[str, Any] = {"markers": []}


@tool
def generate_umap(
    color_by: str = "",
    genes: str = "",
    title: str = "",
) -> str:
    """
    Generate a UMAP embedding plot and return it as an embedded PNG image.

    Args:
        color_by: Obs column to color cells by (e.g. "author_cell_type", "disease").
                  Leave empty to use default (first available categorical column).
        genes: Comma-separated gene names to color by (e.g. "CD3E,CD4").
               When provided, overrides color_by — one panel per gene.
        title: Optional plot title.

    Returns:
        Markdown string with embedded PNG, e.g. ![title](data:image/png;base64,...)
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    if "X_umap" not in adata.obsm:
        return "No UMAP embedding found in this dataset (adata.obsm has no 'X_umap' key)."

    gene_list = _resolve_gene_names(adata, [g.strip() for g in genes.split(",") if g.strip()]) if genes else []
    color: list[str] = gene_list or ([color_by] if color_by else [])

    try:
        sc.pl.umap(adata, color=color or None, show=False, size=3,
                   title=title or None)
        return _fig_to_markdown(title or f"UMAP — {color[0] if color else 'unlabeled'}")
    except Exception as exc:
        plt.close("all")
        return f"UMAP plot failed: {exc}"


@tool
def generate_tsne(
    color_by: str = "",
    genes: str = "",
    title: str = "",
) -> str:
    """
    Generate a t-SNE embedding plot and return it as an embedded PNG image.

    Args:
        color_by: Obs column to color cells by.
        genes: Comma-separated gene names to color by.
        title: Optional plot title.

    Returns:
        Markdown string with embedded PNG.
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    if "X_tsne" not in adata.obsm:
        return "No t-SNE embedding found in this dataset (adata.obsm has no 'X_tsne' key)."

    gene_list = _resolve_gene_names(adata, [g.strip() for g in genes.split(",") if g.strip()]) if genes else []
    color: list[str] = gene_list or ([color_by] if color_by else [])

    try:
        sc.pl.tsne(adata, color=color or None, show=False, size=3,
                   title=title or None)
        return _fig_to_markdown(title or f"t-SNE — {color[0] if color else 'unlabeled'}")
    except Exception as exc:
        plt.close("all")
        return f"t-SNE plot failed: {exc}"


@tool
def generate_violin(
    genes: str,
    groupby: str = "",
    title: str = "",
) -> str:
    """
    Generate a violin plot showing gene expression distributions across cell groups.

    Args:
        genes: Comma-separated gene names (e.g. "CD3E,CD4,CD8A"). Required.
        groupby: Obs column to group cells by (default: author_cell_type if present).
        title: Optional plot title.

    Returns:
        Markdown string with embedded PNG.
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    raw_keys = [g.strip() for g in genes.split(",") if g.strip()]
    if not raw_keys:
        return "Please provide at least one gene name in the 'genes' argument."
    keys = _resolve_gene_names(adata, raw_keys)

    resolved_groupby: str | None = _resolve_groupby(adata, groupby)

    try:
        sc.pl.violin(adata, keys=keys, groupby=resolved_groupby,
                     rotation=45, show=False)
        return _fig_to_markdown(title or f"Violin — {', '.join(keys[:3])}")
    except Exception as exc:
        plt.close("all")
        return f"Violin plot failed: {exc}"


@tool
def generate_dotplot(
    genes: str = "",
    groupby: str = "",
    title: str = "",
) -> str:
    """
    Generate a dot plot showing mean expression and fraction of expressing cells.

    Args:
        genes: Comma-separated gene names. If empty, uses session markers set via set_session_markers.
        groupby: Obs column to group cells by (default: author_cell_type if present).
        title: Optional plot title.

    Returns:
        Markdown string with embedded PNG.
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    raw_keys = [g.strip() for g in genes.split(",") if g.strip()] or _SESSION_STATE["markers"]
    if not raw_keys:
        return "No genes specified and no session markers set. Use set_session_markers first or provide genes."
    keys = _resolve_gene_names(adata, raw_keys)

    resolved_groupby: str | None = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return "No suitable groupby column found. Pass a groupby argument (e.g. 'author_cell_type')."

    try:
        sc.pl.dotplot(adata, var_names=keys, groupby=resolved_groupby, show=False)
        return _fig_to_markdown(title or f"Dot plot — {resolved_groupby}")
    except Exception as exc:
        plt.close("all")
        return f"Dot plot failed: {exc}"


@tool
def generate_heatmap(
    genes: str = "",
    groupby: str = "",
    title: str = "",
) -> str:
    """
    Generate a heatmap of mean gene expression per cell group.

    Args:
        genes: Comma-separated gene names. If empty, uses session markers.
        groupby: Obs column to group cells by (default: author_cell_type if present).
        title: Optional plot title.

    Returns:
        Markdown string with embedded PNG.
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    raw_keys = [g.strip() for g in genes.split(",") if g.strip()] or _SESSION_STATE["markers"]
    if not raw_keys:
        return "No genes specified and no session markers set. Use set_session_markers first or provide genes."
    keys = _resolve_gene_names(adata, raw_keys)

    resolved_groupby: str | None = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return "No suitable groupby column found. Pass a groupby argument."

    try:
        sc.pl.heatmap(adata, var_names=keys, groupby=resolved_groupby, show=False)
        return _fig_to_markdown(title or f"Heatmap — {resolved_groupby}")
    except Exception as exc:
        plt.close("all")
        return f"Heatmap failed: {exc}"


@tool
def generate_cell_counts_barplot(
    groupby: str = "",
    title: str = "",
) -> str:
    """
    Generate a bar plot showing cell counts per category (e.g. cell type, disease).

    Args:
        groupby: Obs column to group by (default: author_cell_type if present).
        title: Optional plot title.

    Returns:
        Markdown string with embedded PNG.
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return "No suitable groupby column found. Pass a groupby argument."

    counts = adata.obs[resolved_groupby].value_counts().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.4), 5))
    counts.plot(kind="bar", ax=ax, color="steelblue", edgecolor="white")
    ax.set_xlabel(resolved_groupby)
    ax.set_ylabel("Cell count")
    ax.set_title(title or f"Cell counts by {resolved_groupby}")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.tight_layout()
    return _fig_to_markdown(title or f"Cell counts — {resolved_groupby}")


@tool
def generate_cell_type_proportion_barplot(
    groupby: str = "",
    sample_col: str = "",
    title: str = "",
) -> str:
    """
    Generate a stacked bar plot showing cell type proportions per sample or condition.

    Args:
        groupby: Cell-type obs column (default: author_cell_type if present).
        sample_col: Sample/donor obs column to stack by (auto-detected if empty).
        title: Optional plot title.

    Returns:
        Markdown string with embedded PNG.
    """
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return str(exc)

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return "No suitable groupby column found. Pass a groupby argument."

    # Auto-detect sample column
    resolved_sample = sample_col
    if not resolved_sample:
        for candidate in ("sample", "donor_id", "donor", "patient", "sample_id", "batch"):
            if candidate in adata.obs.columns:
                resolved_sample = candidate
                break
    if not resolved_sample:
        return (
            "Could not auto-detect a sample column. "
            "Pass sample_col (e.g. 'donor_id', 'sample', 'batch')."
        )

    proportions = (
        adata.obs.groupby([resolved_sample, resolved_groupby])
        .size()
        .unstack(fill_value=0)
        .div(adata.obs.groupby(resolved_sample).size(), axis=0)
    )

    fig, ax = plt.subplots(figsize=(max(8, len(proportions) * 0.5), 5))
    proportions.plot(kind="bar", stacked=True, ax=ax, legend=True, edgecolor="none")
    ax.set_xlabel(resolved_sample)
    ax.set_ylabel("Proportion")
    ax.set_title(title or f"Cell type proportions by {resolved_sample}")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7, frameon=False)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.tight_layout()
    return _fig_to_markdown(title or f"Proportions — {resolved_sample}")


# ─── Internal helper ──────────────────────────────────────────────────────────

def _resolve_groupby(adata: Any, requested: str) -> str | None:
    """Return requested groupby if valid, otherwise fall back to author_cell_type."""
    if requested and requested in adata.obs.columns:
        return requested
    for fallback in ("author_cell_type", "cell_type", "celltype", "cluster", "louvain", "leiden"):
        if fallback in adata.obs.columns:
            return fallback
    # Return first categorical column as last resort
    for col in adata.obs.columns:
        if str(adata.obs[col].dtype) in ("category", "object"):
            return col
    return None


# ─── GenoPixel compatibility tool layer ──────────────────────────────────────
# These wrappers expose the same function-tool names used in GenoPixel-AI
# OpenWebUI (gp_catalog_api operation_id values) while keeping this AgentCore
# runtime independent from OpenWebUI.


def _json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, default=str)


def _parse_string_list(value: str | list[str] | None) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v).strip() for v in value if str(v).strip()]
    raw = str(value).strip()
    if not raw:
        return []
    try:
        parsed = json.loads(raw)
    except Exception:
        parsed = None
    if isinstance(parsed, list):
        return [str(v).strip() for v in parsed if str(v).strip()]
    if isinstance(parsed, str):
        raw = parsed
    return [part.strip() for part in raw.split(",") if part.strip()]


def _has_data_image(text: str) -> bool:
    normalized = str(text or "").replace("\\n", "\n").replace('\\"', '"')
    return "data:image/" in normalized and "base64," in normalized


def _no_active_dataset_response(message: str) -> str:
    return _json({"ok": False, "status": "no_active_dataset", "message": message})


def _plot_error_response(message: str, active: dict[str, Any] | None = None) -> str:
    payload: dict[str, Any] = {"ok": False, "status": "plot_error", "message": message}
    if active is not None:
        payload["active_dataset"] = active
    return _json(payload)


def _obs_filter_error_response(message: str, active: dict[str, Any] | None = None) -> str:
    payload: dict[str, Any] = {"ok": False, "status": "obs_filter_error", "message": message}
    if active is not None:
        payload["active_dataset"] = active
    return _json(payload)


def _success_plot_response(
    *,
    active: dict[str, Any],
    inline_markdown: str,
    plot: dict[str, Any] | None = None,
    extra: dict[str, Any] | None = None,
) -> str:
    payload: dict[str, Any] = {
        "ok": True,
        "status": "success",
        "active_dataset": active,
        "inline_markdown": inline_markdown,
    }
    if plot is not None:
        payload["plot"] = plot
    if extra:
        payload.update(extra)
    return _json(payload)


def _apply_obs_filter_compat(adata: Any, obs_filter_json: str) -> Any:
    raw = str(obs_filter_json or "").strip()
    if not raw or raw in {"{}", "null", "None"}:
        return adata
    try:
        parsed = json.loads(raw)
    except Exception as exc:
        raise ValueError(
            "obs_filter_json must be a JSON object like {\"author_cell_type\": [\"T cell\"]}."
        ) from exc
    if not isinstance(parsed, dict):
        raise ValueError("obs_filter_json must be a JSON object mapping obs columns to value lists.")

    filtered = adata
    for col, values in parsed.items():
        col_name = str(col)
        if col_name not in filtered.obs.columns:
            raise ValueError(f"obs filter column '{col_name}' was not found in adata.obs.")
        allowed = _parse_string_list(values if isinstance(values, list) else str(values))
        if not allowed:
            continue
        allowed_norm = {str(v).casefold() for v in allowed}
        mask = filtered.obs[col_name].astype(str).str.casefold().isin(allowed_norm)
        filtered = filtered[mask].copy()
    return filtered


def _active_adata_with_filter(
    obs_filter_json: str = "{}",
) -> tuple[Any | None, dict[str, Any] | None, str | None]:
    try:
        adata, active = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return None, None, _no_active_dataset_response(str(exc))

    try:
        filtered = _apply_obs_filter_compat(adata, obs_filter_json)
    except ValueError as exc:
        return None, active, _obs_filter_error_response(str(exc), active)
    return filtered, active, None


def _resolve_embedding_basis(adata: Any, basis: str) -> tuple[str | None, list[str]]:
    available = [str(k) for k in adata.obsm.keys()]
    requested = str(basis or "").strip()
    if not requested:
        return None, available

    if requested in available:
        return requested[2:] if requested.startswith("X_") else requested, available

    prefixed = f"X_{requested}"
    if prefixed in available:
        return requested, available

    req_l = requested.lower()
    exact = [k for k in available if k.lower() == req_l or k.lower().replace("x_", "", 1) == req_l]
    if len(exact) == 1:
        k = exact[0]
        return k[2:] if k.startswith("X_") else k, available

    partial = [k for k in available if req_l in k.lower().replace("x_", "", 1)]
    if len(partial) == 1:
        k = partial[0]
        return k[2:] if k.startswith("X_") else k, available

    return None, available


def _safe_rank_genes_groups_key(adata: Any, key: str) -> str | None:
    resolved_key = str(key or "rank_genes_groups").strip() or "rank_genes_groups"
    payload = adata.uns.get(resolved_key)
    if isinstance(payload, dict) and "names" in payload:
        return resolved_key
    return None


@tool
def set_markers(markers_json: str) -> str:
    """
    Store a session-level marker list. Compatible with GenoPixel set_markers.
    """
    raw_markers = _parse_string_list(markers_json)
    if not raw_markers:
        return _json({"ok": False, "status": "error", "message": "markers_json must contain at least one gene."})
    try:
        adata, _ = RUNTIME_STATE.require_active_adata()
        markers = _resolve_gene_names(adata, raw_markers)
    except Exception:
        markers = raw_markers  # no dataset loaded yet — store as-is
    _SESSION_STATE["markers"] = markers
    return _json(
        {
            "ok": True,
            "status": "success",
            "markers": markers,
            "count": len(markers),
            "message": (
                f"Stored {len(markers)} marker gene(s). They will be used as the default gene list "
                "for violin, dotplot, matrixplot, and heatmap plots."
            ),
        }
    )


@tool
def get_markers() -> str:
    """
    Return current session-level markers. Compatible with GenoPixel get_markers.
    """
    markers = [str(v) for v in _SESSION_STATE.get("markers", []) if str(v).strip()]
    return _json({"ok": True, "markers": markers, "count": len(markers)})


@tool
def log_unmet_request(user_request: str, active_dataset: str = "") -> str:
    """
    Log unmet requests for future feature planning.
    """
    record = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "user_request": str(user_request or "").strip(),
        "active_dataset": str(active_dataset or "").strip(),
    }
    try:
        log_path = _output_dir().parent / "unmet_requests.jsonl"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(record, ensure_ascii=False) + "\n")
    except Exception as exc:
        return _json({"ok": False, "status": "log_error", "message": str(exc)})
    return _json(
        {
            "ok": True,
            "status": "logged",
            "message": "Your request has been noted and sent to the developers for consideration.",
        }
    )


@tool
def generate_umap_plot(
    color_json: str = "[]",
    obs_filter_json: str = "{}",
    title: str = "",
    size: float | None = None,
    color_map: str = "",
    palette: str = "",
    legend_loc: str = "right margin",
    add_outline: bool = False,
    edges: bool = False,
    vmin: str | None = None,
    vmax: str | None = None,
    groups_json: str = "",
    ncols: int = 4,
) -> str:
    """
    UMAP plot with GenoPixel-compatible name and response schema.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    if "X_umap" not in adata.obsm:
        try:
            if "X_pca" not in adata.obsm:
                sc.tl.pca(adata, svd_solver="arpack")
            if "neighbors" not in adata.uns:
                sc.pp.neighbors(adata)
            sc.tl.umap(adata)
        except Exception as exc:
            return _plot_error_response(f"Could not compute UMAP embedding: {exc}", active)

    colors = _resolve_color_tokens(adata, _parse_string_list(color_json))
    groups = _parse_string_list(groups_json) if groups_json else []
    kwargs: dict[str, Any] = {"show": False}
    if title:
        kwargs["title"] = title
    if size is not None:
        kwargs["size"] = float(size)
    if color_map:
        kwargs["color_map"] = color_map
    if palette:
        kwargs["palette"] = palette
    if legend_loc:
        kwargs["legend_loc"] = legend_loc
    if add_outline:
        kwargs["add_outline"] = bool(add_outline)
    if edges:
        kwargs["edges"] = bool(edges)
    if vmin not in (None, ""):
        kwargs["vmin"] = vmin
    if vmax not in (None, ""):
        kwargs["vmax"] = vmax
    if groups:
        kwargs["groups"] = groups
    if ncols:
        kwargs["ncols"] = int(ncols)

    try:
        sc.pl.umap(adata, color=colors or None, **kwargs)
        inline_markdown = _fig_to_markdown(title or "UMAP")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    plot = {
        "plot_type": "umap",
        "resolved_coloring_label": colors[0] if colors else None,
        "resolved_groupby": None,
        "resolved_genes": None,
    }
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot=plot)


@tool
def generate_tsne_plot(
    color_json: str = "[]",
    obs_filter_json: str = "{}",
    title: str = "",
    size: float | None = None,
    color_map: str = "",
    palette: str = "",
    legend_loc: str = "right margin",
    add_outline: bool = False,
    edges: bool = False,
    vmin: str | None = None,
    vmax: str | None = None,
    groups_json: str = "",
    ncols: int = 4,
) -> str:
    """
    tSNE plot with GenoPixel-compatible name and response schema.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    if "X_tsne" not in adata.obsm:
        try:
            if "X_pca" not in adata.obsm:
                sc.tl.pca(adata, svd_solver="arpack")
            sc.tl.tsne(adata)
        except Exception as exc:
            return _plot_error_response(f"Could not compute tSNE embedding: {exc}", active)

    colors = _resolve_color_tokens(adata, _parse_string_list(color_json))
    groups = _parse_string_list(groups_json) if groups_json else []
    kwargs: dict[str, Any] = {"show": False}
    if title:
        kwargs["title"] = title
    if size is not None:
        kwargs["size"] = float(size)
    if color_map:
        kwargs["color_map"] = color_map
    if palette:
        kwargs["palette"] = palette
    if legend_loc:
        kwargs["legend_loc"] = legend_loc
    if add_outline:
        kwargs["add_outline"] = bool(add_outline)
    if edges:
        kwargs["edges"] = bool(edges)
    if vmin not in (None, ""):
        kwargs["vmin"] = vmin
    if vmax not in (None, ""):
        kwargs["vmax"] = vmax
    if groups:
        kwargs["groups"] = groups
    if ncols:
        kwargs["ncols"] = int(ncols)

    try:
        sc.pl.tsne(adata, color=colors or None, **kwargs)
        inline_markdown = _fig_to_markdown(title or "tSNE")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    plot = {
        "plot_type": "tsne",
        "resolved_coloring_label": colors[0] if colors else None,
        "resolved_groupby": None,
        "resolved_genes": None,
    }
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot=plot)


@tool
def generate_embedding_plot(
    basis: str,
    color_json: str = "[]",
    obs_filter_json: str = "{}",
    components: str = "",
    title: str = "",
    legend_loc: str = "right margin",
    size: float | None = None,
    color_map: str = "",
    palette: str = "",
    vmin: str | None = None,
    vmax: str | None = None,
    add_outline: bool = False,
    edges: bool = False,
    ncols: int = 4,
) -> str:
    """
    Generic embedding plot (UMAP/PCA/etc.) by basis key.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_basis, available = _resolve_embedding_basis(adata, basis)
    if not resolved_basis:
        return _json(
            {
                "ok": False,
                "status": "no_embedding",
                "message": (
                    f"Embedding '{basis}' was not found in adata.obsm. "
                    f"Available embeddings: {', '.join(available) if available else '(none)'}"
                ),
                "active_dataset": active,
            }
        )

    colors = _resolve_color_tokens(adata, _parse_string_list(color_json))
    kwargs: dict[str, Any] = {"show": False}
    if title:
        kwargs["title"] = title
    if legend_loc:
        kwargs["legend_loc"] = legend_loc
    if size is not None:
        kwargs["size"] = float(size)
    if color_map:
        kwargs["color_map"] = color_map
    if palette:
        kwargs["palette"] = palette
    if vmin not in (None, ""):
        kwargs["vmin"] = vmin
    if vmax not in (None, ""):
        kwargs["vmax"] = vmax
    if add_outline:
        kwargs["add_outline"] = bool(add_outline)
    if edges:
        kwargs["edges"] = bool(edges)
    if ncols:
        kwargs["ncols"] = int(ncols)
    if components:
        kwargs["components"] = components

    try:
        sc.pl.embedding(adata, basis=resolved_basis, color=colors or None, **kwargs)
        inline_markdown = _fig_to_markdown(title or f"Embedding — {resolved_basis}")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    plot = {
        "plot_type": "embedding",
        "embedding_basis": resolved_basis,
        "resolved_coloring_label": colors[0] if colors else None,
    }
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot=plot)


@tool
def generate_diffmap_plot(
    color_json: str = "[]",
    obs_filter_json: str = "{}",
    components: str = "",
    title: str = "",
    legend_loc: str = "right margin",
    size: float | None = None,
    color_map: str = "",
    palette: str = "",
) -> str:
    """
    Diffusion-map embedding plot. Requires X_diffmap in adata.obsm.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    if "X_diffmap" not in adata.obsm:
        return _json(
            {
                "ok": False,
                "status": "no_diffmap_embedding",
                "message": "Diffusion map embedding was not found in adata.obsm (missing 'X_diffmap').",
                "active_dataset": active,
            }
        )

    colors = _resolve_color_tokens(adata, _parse_string_list(color_json))
    kwargs: dict[str, Any] = {"show": False}
    if title:
        kwargs["title"] = title
    if legend_loc:
        kwargs["legend_loc"] = legend_loc
    if size is not None:
        kwargs["size"] = float(size)
    if color_map:
        kwargs["color_map"] = color_map
    if palette:
        kwargs["palette"] = palette
    if components:
        kwargs["components"] = components

    try:
        sc.pl.diffmap(adata, color=colors or None, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Diffmap")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    plot = {
        "plot_type": "diffmap",
        "resolved_coloring_label": colors[0] if colors else None,
    }
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot=plot)


@tool
def generate_violin_plot(
    keys_json: str = "[]",
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
    log: bool = False,
) -> str:
    """
    Violin plot wrapper compatible with GenoPixel generate_violin_plot.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    raw_keys = _parse_string_list(keys_json) or [str(v) for v in _SESSION_STATE.get("markers", []) if str(v).strip()]
    if not raw_keys:
        return _json(
            {
                "ok": False,
                "status": "error",
                "message": "No genes specified and no session markers are set. Provide keys_json or call set_markers first.",
            }
        )
    keys = _resolve_gene_names(adata, raw_keys)

    resolved_groupby = _resolve_groupby(adata, groupby)
    try:
        sc.pl.violin(adata, keys=keys, groupby=resolved_groupby, show=False, log=bool(log), rotation=45)
        inline_markdown = _fig_to_markdown(title or "Violin")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "violin", "resolved_genes": keys, "resolved_groupby": resolved_groupby},
        extra={"resolved_genes": keys, "resolved_groupby": resolved_groupby},
    )


@tool
def generate_dotplot_plot(
    markers_json: str = "[]",
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = False,
    standard_scale: str = "",
    cmap: str = "Reds",
    expression_cutoff: float = 0.0,
) -> str:
    """
    Dot plot wrapper compatible with GenoPixel generate_dotplot_plot.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    raw_var_names = _parse_string_list(markers_json) or [str(v) for v in _SESSION_STATE.get("markers", []) if str(v).strip()]
    if not raw_var_names:
        return _json(
            {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }
        )
    var_names = _resolve_gene_names(adata, raw_var_names)

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found for dotplot.", active)

    kwargs: dict[str, Any] = {
        "show": False,
        "swap_axes": bool(swap_axes),
        "expression_cutoff": float(expression_cutoff),
        "cmap": cmap or "Reds",
    }
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    try:
        sc.pl.dotplot(adata, var_names=var_names, groupby=resolved_groupby, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Dot plot")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "dotplot", "resolved_genes": var_names, "resolved_groupby": resolved_groupby},
        extra={"var_names_used": var_names, "groupby_used": resolved_groupby},
    )


@tool
def generate_heatmap_plot(
    markers_json: str = "[]",
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = True,
    standard_scale: str = "",
    log: bool = False,
) -> str:
    """
    Heatmap wrapper compatible with GenoPixel generate_heatmap_plot.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    raw_var_names = _parse_string_list(markers_json) or [str(v) for v in _SESSION_STATE.get("markers", []) if str(v).strip()]
    if not raw_var_names:
        return _json(
            {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }
        )
    var_names = _resolve_gene_names(adata, raw_var_names)

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found for heatmap.", active)

    kwargs: dict[str, Any] = {
        "show": False,
        "swap_axes": bool(swap_axes),
        "log": bool(log),
    }
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    try:
        sc.pl.heatmap(adata, var_names=var_names, groupby=resolved_groupby, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Heatmap")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "heatmap", "resolved_genes": var_names, "resolved_groupby": resolved_groupby},
        extra={"var_names_used": var_names, "groupby_used": resolved_groupby},
    )


@tool
def cell_counts_barplot(
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
) -> str:
    """
    Cell counts bar plot wrapper compatible with GenoPixel cell_counts_barplot.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found.", active)

    try:
        counts = adata.obs[resolved_groupby].value_counts().sort_values(ascending=False)
        fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.4), 5))
        counts.plot(kind="bar", ax=ax, color="steelblue", edgecolor="white")
        ax.set_xlabel(resolved_groupby)
        ax.set_ylabel("Cell count")
        ax.set_title(title or f"Cell counts by {resolved_groupby}")
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.tight_layout()
        inline_markdown = _fig_to_markdown(title or "Cell counts")
    except Exception as exc:
        plt.close("all")
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "cell_counts_barplot", "resolved_groupby": resolved_groupby},
        extra={"groupby": resolved_groupby},
    )


@tool
def cell_type_proportion_barplot(
    groupby: str = "author_cell_type",
    sample_col: str = "",
    obs_filter_json: str = "{}",
    title: str = "",
) -> str:
    """
    Stacked proportion bar plot wrapper compatible with GenoPixel.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found.", active)

    resolved_sample = str(sample_col or "").strip()
    if not resolved_sample:
        for candidate in ("sample_id", "donor_id", "sample", "donor", "patient", "batch"):
            if candidate in adata.obs.columns:
                resolved_sample = candidate
                break
    if not resolved_sample or resolved_sample not in adata.obs.columns:
        return _plot_error_response(
            "Could not auto-detect a sample column. Please pass sample_col explicitly.",
            active,
        )

    try:
        proportions = (
            adata.obs.groupby([resolved_sample, resolved_groupby])
            .size()
            .unstack(fill_value=0)
            .div(adata.obs.groupby(resolved_sample).size(), axis=0)
        )

        fig, ax = plt.subplots(figsize=(max(8, len(proportions) * 0.5), 5))
        proportions.plot(kind="bar", stacked=True, ax=ax, legend=True, edgecolor="none")
        ax.set_xlabel(resolved_sample)
        ax.set_ylabel("Proportion")
        ax.set_title(title or f"Cell type proportions by {resolved_sample}")
        ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7, frameon=False)
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.tight_layout()
        inline_markdown = _fig_to_markdown(title or "Cell type proportions")
    except Exception as exc:
        plt.close("all")
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "cell_type_proportion_barplot", "resolved_groupby": resolved_groupby},
        extra={"groupby": resolved_groupby},
    )


@tool
def generate_matrixplot(
    markers_json: str = "[]",
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = False,
    standard_scale: str = "",
    cmap: str = "",
) -> str:
    """
    Matrix plot wrapper compatible with GenoPixel generate_matrixplot.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    raw_var_names = _parse_string_list(markers_json) or [str(v) for v in _SESSION_STATE.get("markers", []) if str(v).strip()]
    if not raw_var_names:
        return _json(
            {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }
        )
    var_names = _resolve_gene_names(adata, raw_var_names)

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found for matrixplot.", active)

    kwargs: dict[str, Any] = {"show": False, "swap_axes": bool(swap_axes)}
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    if cmap:
        kwargs["cmap"] = cmap
    try:
        sc.pl.matrixplot(adata, var_names=var_names, groupby=resolved_groupby, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Matrix plot")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "matrixplot", "resolved_genes": var_names, "resolved_groupby": resolved_groupby},
        extra={"var_names_used": var_names, "groupby_used": resolved_groupby},
    )


@tool
def generate_stacked_violin(
    markers_json: str = "[]",
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = False,
    standard_scale: str = "",
    cmap: str = "",
) -> str:
    """
    Stacked violin wrapper compatible with GenoPixel generate_stacked_violin.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    raw_var_names = _parse_string_list(markers_json) or [str(v) for v in _SESSION_STATE.get("markers", []) if str(v).strip()]
    if not raw_var_names:
        return _json(
            {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }
        )
    var_names = _resolve_gene_names(adata, raw_var_names)

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found for stacked violin.", active)

    kwargs: dict[str, Any] = {"show": False, "swap_axes": bool(swap_axes)}
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    if cmap:
        kwargs["cmap"] = cmap
    try:
        sc.pl.stacked_violin(adata, var_names=var_names, groupby=resolved_groupby, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Stacked violin")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "stacked_violin", "resolved_genes": var_names, "resolved_groupby": resolved_groupby},
        extra={"var_names_used": var_names, "groupby_used": resolved_groupby},
    )


@tool
def generate_tracksplot(
    markers_json: str = "[]",
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
    log: bool = False,
) -> str:
    """
    Tracksplot wrapper compatible with GenoPixel generate_tracksplot.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    raw_var_names = _parse_string_list(markers_json) or [str(v) for v in _SESSION_STATE.get("markers", []) if str(v).strip()]
    if not raw_var_names:
        return _json(
            {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }
        )
    var_names = _resolve_gene_names(adata, raw_var_names)

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found for tracksplot.", active)

    try:
        sc.pl.tracksplot(adata, var_names=var_names, groupby=resolved_groupby, show=False, log=bool(log))
        inline_markdown = _fig_to_markdown(title or "Tracksplot")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "tracksplot", "resolved_genes": var_names, "resolved_groupby": resolved_groupby},
        extra={"var_names_used": var_names, "groupby_used": resolved_groupby},
    )


@tool
def generate_dendrogram(
    groupby: str = "author_cell_type",
    obs_filter_json: str = "{}",
    title: str = "",
    orientation: str = "top",
) -> str:
    """
    Dendrogram wrapper compatible with GenoPixel generate_dendrogram.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found for dendrogram.", active)

    try:
        sc.tl.dendrogram(adata, groupby=resolved_groupby)
        sc.pl.dendrogram(adata, groupby=resolved_groupby, orientation=orientation or "top", show=False)
        inline_markdown = _fig_to_markdown(title or "Dendrogram")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "dendrogram", "resolved_groupby": resolved_groupby},
        extra={"groupby_used": resolved_groupby},
    )


@tool
def generate_clustermap() -> str:
    """
    Clustermap endpoint is intentionally disabled, matching GenoPixel behavior.
    """
    return _json({"ok": False, "status": "disabled", "message": "Clustermap plotting is temporarily disabled."})


@tool
def generate_correlation_matrix_plot(
    groupby: str = "",
    obs_filter_json: str = "{}",
    title: str = "",
    show_correlation_numbers: bool = False,
    dendrogram: bool | None = None,
    cmap: str = "",
) -> str:
    """
    Correlation matrix wrapper compatible with GenoPixel generate_correlation_matrix_plot.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_groupby = _resolve_groupby(adata, groupby)
    if not resolved_groupby:
        return _plot_error_response("No suitable groupby column found for correlation matrix.", active)

    kwargs: dict[str, Any] = {"show": False, "groupby": resolved_groupby, "show_correlation_numbers": bool(show_correlation_numbers)}
    if dendrogram is not None:
        kwargs["dendrogram"] = bool(dendrogram)
    if cmap:
        kwargs["cmap"] = cmap

    try:
        sc.pl.correlation_matrix(adata, **kwargs)
        if title:
            plt.title(title)
        inline_markdown = _fig_to_markdown(title or "Correlation matrix")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "correlation_matrix", "resolved_groupby": resolved_groupby},
    )


@tool
def check_rank_genes_groups() -> str:
    """
    Check whether rank_genes_groups results exist and plot a quick panel when available.
    """
    try:
        adata, active = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return _no_active_dataset_response(str(exc))

    key = _safe_rank_genes_groups_key(adata, "rank_genes_groups")
    if not key:
        return _json(
            {
                "ok": False,
                "status": "not_available",
                "message": "rank_genes_groups results are not available in the active dataset.",
                "active_dataset": active,
            }
        )

    try:
        sc.pl.rank_genes_groups(adata, key=key, n_genes=5, show=False)
        inline_markdown = _fig_to_markdown("Rank genes groups")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _json(
        {
            "ok": True,
            "status": "available",
            "active_dataset": active,
            "groupby_used": _resolve_groupby(adata, "author_cell_type"),
            "inline_markdown": inline_markdown,
        }
    )


@tool
def generate_rank_genes_groups_plot(
    groups_json: str = "",
    n_genes: int = 5,
    key: str = "rank_genes_groups",
    obs_filter_json: str = "{}",
    title: str = "",
    ncols: int = 3,
    sharey: bool = True,
) -> str:
    """
    Rank genes groups score-panel plot. Does not auto-compute rank_genes_groups.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_key = _safe_rank_genes_groups_key(adata, key)
    if not resolved_key:
        return _json(
            {
                "ok": False,
                "status": "no_rank_genes_groups",
                "message": f"rank_genes_groups key '{key}' not found in adata.uns.",
                "active_dataset": active,
            }
        )
    groups = _parse_string_list(groups_json) if groups_json else None
    try:
        sc.pl.rank_genes_groups(
            adata,
            key=resolved_key,
            groups=groups,
            n_genes=int(n_genes),
            ncols=int(ncols),
            sharey=bool(sharey),
            show=False,
        )
        inline_markdown = _fig_to_markdown(title or "Rank genes groups")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "rank_genes_groups"})


@tool
def generate_rank_genes_groups_violin(
    groups_json: str = "",
    n_genes: int = 5,
    gene_names_json: str = "",
    key: str = "rank_genes_groups",
    obs_filter_json: str = "{}",
    title: str = "",
    split: bool = True,
) -> str:
    """
    Rank genes violin plot. Does not auto-compute rank_genes_groups.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_key = _safe_rank_genes_groups_key(adata, key)
    if not resolved_key:
        return _json(
            {
                "ok": False,
                "status": "no_rank_genes_groups",
                "message": f"rank_genes_groups key '{key}' not found in adata.uns.",
                "active_dataset": active,
            }
        )
    groups = _parse_string_list(groups_json) if groups_json else None
    raw_gene_names = _parse_string_list(gene_names_json) if gene_names_json else None
    gene_names = _resolve_gene_names(adata, raw_gene_names) if raw_gene_names else None
    try:
        sc.pl.rank_genes_groups_violin(
            adata,
            key=resolved_key,
            groups=groups,
            n_genes=int(n_genes),
            gene_names=gene_names,
            split=bool(split),
            show=False,
        )
        inline_markdown = _fig_to_markdown(title or "Rank genes violin")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "rank_genes_groups_violin"})


@tool
def generate_rank_genes_groups_dotplot_plot(
    groups_json: str = "",
    n_genes: int = 5,
    groupby: str = "",
    key: str = "rank_genes_groups",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = False,
    dendrogram: bool = False,
    standard_scale: str = "",
    values_to_plot: str = "",
) -> str:
    """
    Rank genes dotplot. Does not auto-compute rank_genes_groups.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_key = _safe_rank_genes_groups_key(adata, key)
    if not resolved_key:
        return _json(
            {
                "ok": False,
                "status": "no_rank_genes_groups",
                "message": f"rank_genes_groups key '{key}' not found in adata.uns.",
                "active_dataset": active,
            }
        )

    groups = _parse_string_list(groups_json) if groups_json else None
    kwargs: dict[str, Any] = {
        "key": resolved_key,
        "groups": groups,
        "n_genes": int(n_genes),
        "show": False,
        "swap_axes": bool(swap_axes),
        "dendrogram": bool(dendrogram),
    }
    resolved_groupby = _resolve_groupby(adata, groupby)
    if resolved_groupby:
        kwargs["groupby"] = resolved_groupby
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    if values_to_plot:
        kwargs["values_to_plot"] = values_to_plot
    try:
        sc.pl.rank_genes_groups_dotplot(adata, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Rank genes dotplot")
    except Exception as exc:
        return _plot_error_response(str(exc), active)
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "rank_genes_groups_dotplot"})


@tool
def generate_rank_genes_groups_matrixplot_plot(
    groups_json: str = "",
    n_genes: int = 5,
    groupby: str = "",
    key: str = "rank_genes_groups",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = False,
    dendrogram: bool = False,
    standard_scale: str = "",
    values_to_plot: str = "",
) -> str:
    """
    Rank genes matrix plot. Does not auto-compute rank_genes_groups.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_key = _safe_rank_genes_groups_key(adata, key)
    if not resolved_key:
        return _json(
            {
                "ok": False,
                "status": "no_rank_genes_groups",
                "message": f"rank_genes_groups key '{key}' not found in adata.uns.",
                "active_dataset": active,
            }
        )

    groups = _parse_string_list(groups_json) if groups_json else None
    kwargs: dict[str, Any] = {
        "key": resolved_key,
        "groups": groups,
        "n_genes": int(n_genes),
        "show": False,
        "swap_axes": bool(swap_axes),
        "dendrogram": bool(dendrogram),
    }
    resolved_groupby = _resolve_groupby(adata, groupby)
    if resolved_groupby:
        kwargs["groupby"] = resolved_groupby
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    if values_to_plot:
        kwargs["values_to_plot"] = values_to_plot
    try:
        sc.pl.rank_genes_groups_matrixplot(adata, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Rank genes matrixplot")
    except Exception as exc:
        return _plot_error_response(str(exc), active)
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "rank_genes_groups_matrixplot"})


@tool
def generate_rank_genes_groups_heatmap_plot(
    groups_json: str = "",
    n_genes: int = 5,
    groupby: str = "",
    key: str = "rank_genes_groups",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = False,
    standard_scale: str = "",
) -> str:
    """
    Rank genes heatmap. Does not auto-compute rank_genes_groups.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_key = _safe_rank_genes_groups_key(adata, key)
    if not resolved_key:
        return _json(
            {
                "ok": False,
                "status": "no_rank_genes_groups",
                "message": f"rank_genes_groups key '{key}' not found in adata.uns.",
                "active_dataset": active,
            }
        )

    groups = _parse_string_list(groups_json) if groups_json else None
    kwargs: dict[str, Any] = {
        "key": resolved_key,
        "groups": groups,
        "n_genes": int(n_genes),
        "show": False,
        "swap_axes": bool(swap_axes),
    }
    resolved_groupby = _resolve_groupby(adata, groupby)
    if resolved_groupby:
        kwargs["groupby"] = resolved_groupby
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    try:
        sc.pl.rank_genes_groups_heatmap(adata, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Rank genes heatmap")
    except Exception as exc:
        return _plot_error_response(str(exc), active)
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "rank_genes_groups_heatmap"})


@tool
def generate_rank_genes_groups_tracksplot_plot(
    groups_json: str = "",
    n_genes: int = 5,
    groupby: str = "",
    key: str = "rank_genes_groups",
    obs_filter_json: str = "{}",
    title: str = "",
    dendrogram: bool = False,
) -> str:
    """
    Rank genes tracks plot. Does not auto-compute rank_genes_groups.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_key = _safe_rank_genes_groups_key(adata, key)
    if not resolved_key:
        return _json(
            {
                "ok": False,
                "status": "no_rank_genes_groups",
                "message": f"rank_genes_groups key '{key}' not found in adata.uns.",
                "active_dataset": active,
            }
        )

    groups = _parse_string_list(groups_json) if groups_json else None
    kwargs: dict[str, Any] = {
        "key": resolved_key,
        "groups": groups,
        "n_genes": int(n_genes),
        "show": False,
        "dendrogram": bool(dendrogram),
    }
    resolved_groupby = _resolve_groupby(adata, groupby)
    if resolved_groupby:
        kwargs["groupby"] = resolved_groupby
    try:
        sc.pl.rank_genes_groups_tracksplot(adata, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Rank genes tracksplot")
    except Exception as exc:
        return _plot_error_response(str(exc), active)
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "rank_genes_groups_tracksplot"})


@tool
def generate_rank_genes_groups_stacked_violin(
    groups_json: str = "",
    n_genes: int = 5,
    key: str = "rank_genes_groups",
    obs_filter_json: str = "{}",
    title: str = "",
    swap_axes: bool = False,
    standard_scale: str = "",
    cmap: str = "Blues",
) -> str:
    """
    Rank genes stacked violin. Does not auto-compute rank_genes_groups.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_key = _safe_rank_genes_groups_key(adata, key)
    if not resolved_key:
        return _json(
            {
                "ok": False,
                "status": "no_rank_genes_groups",
                "message": f"rank_genes_groups key '{key}' not found in adata.uns.",
                "active_dataset": active,
            }
        )

    plot_fn = getattr(sc.pl, "rank_genes_groups_stacked_violin", None)
    if not callable(plot_fn):
        return _plot_error_response("This Scanpy version does not provide rank_genes_groups_stacked_violin.", active)

    groups = _parse_string_list(groups_json) if groups_json else None
    kwargs: dict[str, Any] = {
        "key": resolved_key,
        "groups": groups,
        "n_genes": int(n_genes),
        "show": False,
        "swap_axes": bool(swap_axes),
        "cmap": cmap or "Blues",
    }
    if standard_scale:
        kwargs["standard_scale"] = standard_scale
    try:
        plot_fn(adata, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Rank genes stacked violin")
    except Exception as exc:
        return _plot_error_response(str(exc), active)
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "rank_genes_groups_stacked_violin"})


@tool
def print_adata() -> str:
    """
    Return active AnnData summary payload.
    """
    try:
        adata, active = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return _no_active_dataset_response(str(exc))

    obs_cols = [str(c) for c in adata.obs.columns]
    var_cols = [str(c) for c in adata.var.columns]
    obsm_keys = [str(k) for k in adata.obsm.keys()]
    obsp_keys = [str(k) for k in adata.obsp.keys()]
    uns_keys = [str(k) for k in adata.uns.keys()]
    layers_keys = [str(k) for k in adata.layers.keys()]

    cell_type_keywords = (
        "cell_type",
        "celltype",
        "cell_label",
        "annotation",
        "cluster",
        "leiden",
        "louvain",
        "subtype",
    )
    annotation_col: str | None = None
    for col in obs_cols:
        lc = col.lower().replace(" ", "_").replace("-", "_")
        if any(token in lc for token in cell_type_keywords):
            annotation_col = col
            break
    ann = annotation_col or (obs_cols[0] if obs_cols else "cell_type")
    has_rgg = any(isinstance(adata.uns.get(k), dict) and "names" in adata.uns[k] for k in adata.uns)

    suggested_next_steps = [
        {"label": "UMAP", "prompt": "show me a UMAP"},
        {"label": "Cell counts", "prompt": f"plot cell counts by {ann}"},
        {"label": "Violin plot", "prompt": f"use violin plot to show the expression of any gene across {ann}"},
    ]
    if has_rgg:
        suggested_next_steps.append({"label": "Markers/top genes", "prompt": "show top marker genes"})

    return _json(
        {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "summary": str(adata),
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "obs_columns": obs_cols,
            "var_columns": var_cols,
            "obsm_keys": obsm_keys,
            "obsp_keys": obsp_keys,
            "uns_keys": uns_keys,
            "layers": layers_keys,
            "suggested_next_steps": suggested_next_steps,
        }
    )


@tool
def print_adata_obs(command: str = "print(adata.obs)") -> str:
    """
    List adata.obs columns with categorical previews.
    """
    try:
        adata, active = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return _no_active_dataset_response(str(exc))

    normalized_command = "".join(str(command or "").lower().split())
    if normalized_command != "print(adata.obs)":
        return _json(
            {
                "ok": False,
                "status": "unsupported_command",
                "message": "Only print(adata.obs) is supported in this endpoint.",
                "command": command,
                "supported_commands": ["print(adata.obs)"],
            }
        )

    categorical_cols: list[dict[str, Any]] = []
    numerical_cols: list[str] = []
    other_cols: list[str] = []

    for col in adata.obs.columns:
        col_str = str(col)
        series = adata.obs[col]
        dtype = str(series.dtype)
        if dtype == "category" or dtype == "object":
            try:
                unique_vals = sorted(series.dropna().astype(str).unique().tolist())
            except Exception:
                unique_vals = []
            max_values = 50
            categorical_cols.append(
                {
                    "column": col_str,
                    "n_unique": len(unique_vals),
                    "values": unique_vals[:max_values],
                    "values_truncated": len(unique_vals) > max_values,
                }
            )
        elif dtype.startswith(("int", "float", "uint")):
            numerical_cols.append(col_str)
        else:
            other_cols.append(col_str)

    n_cells, n_cols = int(adata.n_obs), int(len(adata.obs.columns))
    md_lines = [f"**Dataset OBS — {n_cells:,} cells x {n_cols} columns**", ""]
    if categorical_cols:
        md_lines.append("**Categorical columns:**")
        for c in sorted(categorical_cols, key=lambda row: str(row["column"])):
            vals = [str(v) for v in c["values"]]
            preview = ", ".join(vals[:8])
            if len(vals) > 8:
                preview += f", ... (+{len(vals) - 8} more)"
            md_lines.append(f"- **{c['column']}** ({c['n_unique']} unique): {preview}")
        md_lines.append("")
    if numerical_cols:
        md_lines.append("**Numerical columns:**")
        md_lines.append("- " + ", ".join(numerical_cols))
        md_lines.append("")
    if other_cols:
        md_lines.append("**Other columns:**")
        md_lines.append("- " + ", ".join(other_cols))

    return _json(
        {
            "ok": True,
            "status": "success",
            "command": "print(adata.obs)",
            "active_dataset": active,
            "obs_shape": [n_cells, n_cols],
            "categorical_columns": sorted(categorical_cols, key=lambda row: str(row["column"])),
            "numerical_columns": numerical_cols,
            "other_columns": other_cols,
            "inline_markdown": "\n".join(md_lines).strip(),
        }
    )


@tool
def get_obs_unique_values(column: str = "") -> str:
    """
    Return unique values for one obs column or all categorical columns.
    """
    try:
        adata, active = RUNTIME_STATE.require_active_adata()
    except NoActiveDatasetError as exc:
        return _no_active_dataset_response(str(exc))

    requested = str(column).strip()

    def col_unique(series: Any) -> list[str]:
        try:
            vals = sorted(series.dropna().unique().tolist(), key=lambda v: str(v))
            return [str(v) for v in vals]
        except Exception:
            return [str(v) for v in series.dropna().unique().tolist()]

    if requested:
        lower_to_col = {str(c).lower(): str(c) for c in adata.obs.columns}
        resolved = requested if requested in adata.obs.columns else lower_to_col.get(requested.lower())
        if not resolved:
            return _json(
                {
                    "ok": False,
                    "status": "column_not_found",
                    "message": f"Column '{requested}' not found in adata.obs.",
                    "available_columns": [str(c) for c in adata.obs.columns],
                    "active_dataset": active,
                }
            )
        vals = col_unique(adata.obs[resolved])
        md = f"**{resolved}** — {len(vals)} unique values:\n\n" + "\n".join(f"- {v}" for v in vals)
        return _json(
            {
                "ok": True,
                "status": "success",
                "active_dataset": active,
                "column": resolved,
                "n_unique": len(vals),
                "values": vals,
                "inline_markdown": md,
            }
        )

    results = []
    for c in adata.obs.columns:
        series = adata.obs[c]
        if str(series.dtype) in {"category", "object"}:
            vals = col_unique(series)
            results.append({"column": str(c), "n_unique": len(vals), "values": vals})

    md_lines = [f"**All categorical obs columns ({len(results)} total):**", ""]
    for item in results:
        md_lines.append(f"**{item['column']}** ({item['n_unique']} unique): " + ", ".join(item["values"]))

    return _json(
        {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "columns": results,
            "inline_markdown": "\n".join(md_lines),
        }
    )


@tool
def generate_highest_expr_genes(
    n_top: int = 30,
    obs_filter_json: str = "{}",
    title: str = "",
    log: bool = False,
) -> str:
    """
    Plot highest expressed genes in the active dataset.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    try:
        sc.pl.highest_expr_genes(adata, n_top=int(n_top), show=False, log=bool(log))
        if title:
            plt.title(title)
        inline_markdown = _fig_to_markdown(title or "Highest expressed genes")
    except Exception as exc:
        return _plot_error_response(str(exc), active)
    return _success_plot_response(active=active, inline_markdown=inline_markdown, plot={"plot_type": "highest_expr_genes"})


@tool
def obs_count_table(
    row_col: str,
    col_col: str,
    obs_filter_json: str = "{}",
) -> str:
    """
    Return a contingency table from two obs columns.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    row_name = str(row_col or "").strip()
    col_name = str(col_col or "").strip()
    if not row_name or not col_name:
        return _json(
            {
                "ok": False,
                "status": "error",
                "message": "row_col and col_col are required.",
                "active_dataset": active,
            }
        )
    if row_name not in adata.obs.columns or col_name not in adata.obs.columns:
        return _json(
            {
                "ok": False,
                "status": "column_not_found",
                "message": "row_col or col_col was not found in adata.obs.",
                "available_columns": [str(c) for c in adata.obs.columns],
                "active_dataset": active,
            }
        )

    table = adata.obs.groupby([row_name, col_name], dropna=False).size().unstack(fill_value=0)
    preview_rows = min(20, len(table))
    md_lines = [
        f"**Obs count table: {row_name} x {col_name}**",
        "",
        f"Rows: {len(table):,}, Columns: {len(table.columns):,}",
        "",
    ]
    for idx, row in table.head(preview_rows).iterrows():
        row_bits = ", ".join(f"{str(k)}={int(v)}" for k, v in row.items())
        md_lines.append(f"- **{idx}**: {row_bits}")
    if len(table) > preview_rows:
        md_lines.append(f"- ... ({len(table) - preview_rows} more rows)")

    return _json(
        {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "row_col": row_name,
            "col_col": col_name,
            "shape": [int(table.shape[0]), int(table.shape[1])],
            "table": {
                str(idx): {str(k): int(v) for k, v in row.items()}
                for idx, row in table.iterrows()
            },
            "inline_markdown": "\n".join(md_lines),
        }
    )


@tool
def generate_spatial_scatter(
    color_json: str = "[]",
    obs_filter_json: str = "{}",
    title: str = "",
    ncols: int | None = None,
    size: float | None = None,
    cmap: str | None = None,
    palette: str | None = None,
    legend_loc: str | None = None,
    groups_json: str = "[]",
) -> str:
    """
    Spatial scatter fallback using an embedding basis named 'spatial' or 'X_spatial'.
    """
    adata, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    assert adata is not None and active is not None

    resolved_basis, available = _resolve_embedding_basis(adata, "spatial")
    if not resolved_basis:
        return _json(
            {
                "ok": False,
                "status": "no_spatial_embedding",
                "message": (
                    "No spatial embedding found. Expected 'spatial' or 'X_spatial' in adata.obsm. "
                    f"Available embeddings: {', '.join(available) if available else '(none)'}"
                ),
                "active_dataset": active,
            }
        )

    colors = _resolve_color_tokens(adata, _parse_string_list(color_json))
    groups = _parse_string_list(groups_json)
    kwargs: dict[str, Any] = {"show": False}
    if title:
        kwargs["title"] = title
    if ncols is not None:
        kwargs["ncols"] = int(ncols)
    if size is not None:
        kwargs["size"] = float(size)
    if cmap:
        kwargs["color_map"] = str(cmap)
    if palette:
        kwargs["palette"] = str(palette)
    if legend_loc:
        kwargs["legend_loc"] = str(legend_loc)
    if groups:
        kwargs["groups"] = groups

    try:
        sc.pl.embedding(adata, basis=resolved_basis, color=colors or None, **kwargs)
        inline_markdown = _fig_to_markdown(title or "Spatial scatter")
    except Exception as exc:
        return _plot_error_response(str(exc), active)

    return _success_plot_response(
        active=active,
        inline_markdown=inline_markdown,
        plot={"plot_type": "spatial_scatter", "embedding_basis": resolved_basis},
    )


@tool
def generate_nhood_enrichment(
    cluster_key: str = "",
    mode: str = "zscore",
    obs_filter_json: str = "{}",
    title: str = "",
) -> str:
    """
    Neighborhood enrichment is not supported in this lightweight AgentCore runtime.
    """
    _, active, err = _active_adata_with_filter(obs_filter_json)
    if err:
        return err
    payload: dict[str, Any] = {
        "ok": False,
        "status": "not_available",
        "message": (
            "Neighborhood enrichment requires Squidpy and spatial graph preprocessing, "
            "which are not enabled in this AgentCore mirror runtime."
        ),
        "requested_cluster_key": str(cluster_key or ""),
        "requested_mode": str(mode or ""),
        "requested_title": str(title or ""),
    }
    if active is not None:
        payload["active_dataset"] = active
    return _json(payload)


@tool
def generate_scanpy_plot(
    plot_type: str = "umap",
    color_json: str = "[]",
    genes_json: str = "[]",
    groupby: str = "",
    gene_symbols_column: str = "",
    title: str = "",
    obs_filter_json: str = "{}",
) -> str:
    """
    Generic Scanpy plot endpoint compatible with GenoPixel generate_scanpy_plot.
    """
    kind = str(plot_type or "umap").strip().lower()
    if kind == "umap":
        return generate_umap_plot(color_json=color_json or genes_json, title=title, obs_filter_json=obs_filter_json)
    if kind == "tsne":
        return generate_tsne_plot(color_json=color_json or genes_json, title=title, obs_filter_json=obs_filter_json)
    if kind == "violin":
        keys_json = genes_json or color_json
        return generate_violin_plot(keys_json=keys_json, groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind == "dotplot":
        return generate_dotplot_plot(markers_json=genes_json, groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind == "heatmap":
        return generate_heatmap_plot(markers_json=genes_json, groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind == "matrixplot":
        return generate_matrixplot(markers_json=genes_json, groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind == "stacked_violin":
        return generate_stacked_violin(markers_json=genes_json, groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind == "tracksplot":
        return generate_tracksplot(markers_json=genes_json, groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind == "cell_counts_barplot":
        return cell_counts_barplot(groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind == "cell_type_proportion_barplot":
        return cell_type_proportion_barplot(groupby=groupby, title=title, obs_filter_json=obs_filter_json)
    if kind in {"gene_cell_embedding", "embedding"}:
        basis = gene_symbols_column or "umap"
        return generate_embedding_plot(basis=basis, color_json=color_json or genes_json, title=title, obs_filter_json=obs_filter_json)

    return _json(
        {
            "ok": False,
            "status": "plot_error",
            "message": (
                f"Unsupported plot type '{plot_type}'. "
                "Supported in this runtime: umap, tsne, violin, dotplot, heatmap, "
                "matrixplot, stacked_violin, tracksplot, cell_counts_barplot, "
                "cell_type_proportion_barplot, gene_cell_embedding."
            ),
        }
    )


# Keep original lightweight tools and expose GenoPixel-compatible tool names.
ALL_GENOPIXEL_TOOLS = [
    get_active_dataset_info,
    load_dataset,
    get_obs_columns,
    get_obs_column_values,
    set_session_markers,
    generate_umap,
    generate_tsne,
    generate_violin,
    generate_dotplot,
    generate_heatmap,
    generate_cell_counts_barplot,
    generate_cell_type_proportion_barplot,
    generate_scanpy_plot,
    generate_heatmap_plot,
    set_markers,
    get_markers,
    log_unmet_request,
    generate_violin_plot,
    cell_counts_barplot,
    cell_type_proportion_barplot,
    generate_dotplot_plot,
    check_rank_genes_groups,
    generate_rank_genes_groups_violin,
    generate_rank_genes_groups_plot,
    generate_rank_genes_groups_dotplot_plot,
    generate_rank_genes_groups_tracksplot_plot,
    generate_correlation_matrix_plot,
    generate_rank_genes_groups_matrixplot_plot,
    generate_rank_genes_groups_heatmap_plot,
    generate_embedding_plot,
    generate_diffmap_plot,
    generate_umap_plot,
    generate_tsne_plot,
    generate_dendrogram,
    generate_clustermap,
    generate_matrixplot,
    generate_stacked_violin,
    generate_tracksplot,
    print_adata,
    print_adata_obs,
    get_obs_unique_values,
    generate_rank_genes_groups_stacked_violin,
    generate_highest_expr_genes,
    obs_count_table,
    generate_spatial_scatter,
    generate_nhood_enrichment,
]
