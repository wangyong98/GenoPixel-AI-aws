from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Any

from fastapi import Depends, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from pydantic import BaseModel, Field

from gp_catalog import CatalogLoadError, DatasetLoadError, GenoPixelCatalogStore
from gp_models import PlotRequest
from gp_plot_response_formatter import build_canonical_response_markdown
from gp_runtime_state import NoActiveDatasetError, RUNTIME_STATE
from gp_scanpy_plotter import ScanpyPlotExecutor


CATALOG_PREFIX = "/api/genopixel-catalog"
RUNTIME_PREFIX = "/api/genopixel-runtime"
START_TIME = time.time()
_PLOTTER: ScanpyPlotExecutor | None = None
_bearer = HTTPBearer(auto_error=False)


class AnalyzeDatasetRequest(BaseModel):
    h5ad_path: str
    multiple_excel_row: int | None = None


class GeneratePlotRequest(BaseModel):
    plot_type: str = Field(
        default="umap",
        description=(
            "Type of Scanpy plot to generate for the currently loaded dataset. "
            "Use 'umap' when the user asks for a UMAP plot."
        ),
    )
    color_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of observation columns or genes to color by. "
            "Use [] when the user does not specify a coloring field."
        ),
    )
    genes_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of genes used by gene-based plot types. "
            "Use [] for a plain UMAP unless the user explicitly asks for genes."
        ),
    )
    groupby: str = Field(
        default="",
        description="Optional observation column used by grouped plot types such as dotplot or violin.",
    )
    gene_symbols_column: str = Field(
        default="",
        description="Optional column name that contains gene symbols if the dataset uses an alternate gene label column.",
    )
    title: str = Field(
        default="",
        description="Optional human-readable plot title.",
    )


class HeatmapPlotRequest(BaseModel):
    markers_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of genes to display in the heatmap. "
            "If empty, the session markers set by set_markers are used. "
            "Example: '[\"C1QA\",\"PSAP\",\"CD79A\",\"CD79B\",\"CST3\",\"LYZ\"]'."
        ),
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Observation column to group cells by on one axis. Default: 'author_cell_type'.",
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata.")
    log: bool = Field(default=False, description="Plot on a logarithmic scale.")
    num_categories: int = Field(
        default=7,
        description="Number of categories when groupby is continuous (non-categorical).",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols (leave empty to use var_names index).",
    )
    layer: str = Field(default="", description="AnnData layer to plot. Leave empty to use X or raw.")
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'var' (gene) or 'obs' (cell). Leave empty to skip.",
    )
    swap_axes: bool = Field(
        default=True,
        description="Swap axes: cell types on x-axis and genes on y-axis. Default True — gives category labels more room to spread out.",
    )
    show_gene_labels: bool | None = Field(
        default=None,
        description="Show gene labels on the plot. Auto-detected (hidden when >50 genes) when None.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed from gene/category count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed from gene/category count when None.",
    )
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None, description="Center of the color scale (useful for diverging colormaps)."
    )
    title: str = Field(default="", description="Optional plot title.")


class SetMarkersRequest(BaseModel):
    markers_json: str | list[str] = Field(
        description=(
            "JSON array or comma-separated list of gene markers to store as the session default. "
            "Examples: '[\"C1QA\",\"PSAP\",\"CD79A\"]' or 'C1QA, PSAP, CD79A'. "
            "These will be used automatically by any tool that accepts genes when no genes are specified."
        ),
    )


class ViolinPlotRequest(BaseModel):
    keys_json: str | list[str] = Field(
        description=(
            "JSON array or comma-separated list of genes or obs fields to plot "
            "(e.g. '[\"OSMR\", \"TNF\"]'). Required."
        ),
    )
    groupby: str = Field(
        default="author_cell_type",
        description=(
            "Observation column to group by on the x-axis. "
            "Defaults to 'author_cell_type'."
        ),
    )
    rotation: float = Field(
        default=45.0,
        description="Rotation angle in degrees for x-axis tick labels. Default: 45.",
    )
    log: bool = Field(default=False, description="Plot on a logarithmic axis.")
    use_raw: bool | None = Field(
        default=None,
        description="Use the raw attribute of adata. Defaults to False when groupby is set.",
    )
    stripplot: bool = Field(
        default=True,
        description="Overlay a strip plot (individual data points) on the violin.",
    )
    jitter: float | bool = Field(
        default=True,
        description="Add jitter to the strip plot. Pass a float to control the amount.",
    )
    size: int = Field(default=1, description="Size of the jitter points.")
    layer: str = Field(
        default="",
        description="AnnData layer to plot. Leave empty to use X or raw.",
    )
    density_norm: str = Field(
        default="width",
        description="How to scale violin width: 'width' (same width), 'area' (same area), or 'count' (proportional to N).",
    )
    order_json: str = Field(
        default="",
        description="JSON array specifying the order of categories on the x-axis (e.g. '[\"B cell\", \"T cell\"]').",
    )
    multi_panel: bool | None = Field(
        default=None,
        description="Display each key in a separate panel. Auto-detected when None.",
    )
    xlabel: str = Field(default="", description="X-axis label.")
    ylabel: str = Field(default="", description="Y-axis label.")
    title: str = Field(default="", description="Optional plot title.")


class CellCountsBarplotRequest(BaseModel):
    groupby: str = Field(
        default="author_cell_type",
        description=(
            "Observation column from adata.obs to count cells by. "
            "Examples: 'author_cell_type', 'cell_type', 'disease', 'tissue', 'donor_id'. "
            "Use print_adata_obs to discover available columns."
        ),
    )
    title: str = Field(
        default="",
        description="Optional plot title. Defaults to 'Cell Counts by <groupby>'.",
    )


class DotplotPlotRequest(BaseModel):
    markers_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of genes to display in the dot plot. "
            "If empty, the session markers set by set_markers are used. "
            "Example: '[\"C1QA\",\"PSAP\",\"CD79A\"]'."
        ),
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Observation column to group cells by. Default: 'author_cell_type'.",
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata.")
    log: bool = Field(default=False, description="Plot on a logarithmic scale.")
    num_categories: int = Field(
        default=7,
        description="Number of categories when groupby is continuous.",
    )
    categories_order_json: str = Field(
        default="",
        description="JSON array specifying display order of groupby categories (e.g. '[\"B cell\",\"T cell\"]').",
    )
    expression_cutoff: float = Field(
        default=0.0,
        description="Expression value cutoff. Genes expressed below this in a group are not shown as expressed.",
    )
    mean_only_expressed: bool = Field(
        default=False,
        description="If True, compute mean expression only over cells that express the gene (> expression_cutoff).",
    )
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'var' (gene) or 'obs' (cell). Leave empty to skip.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    layer: str = Field(default="", description="AnnData layer to plot. Leave empty to use X or raw.")
    swap_axes: bool = Field(
        default=False,
        description="Swap axes: groups on x-axis, genes on y-axis. Useful when there are many groups.",
    )
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None, description="Center of the color scale (useful for diverging colormaps)."
    )
    cmap: str = Field(default="Reds", description="Colormap for dot color. Default: 'Reds'.")
    dot_max: float | None = Field(
        default=None,
        description="Maximum dot size. Values above this are clipped to this size.",
    )
    dot_min: float | None = Field(
        default=None,
        description="Minimum dot size. Values below this are not shown.",
    )
    smallest_dot: float = Field(
        default=0.0,
        description="Smallest dot size in points. Useful to make dots with very low expression still visible.",
    )
    colorbar_title: str = Field(
        default="",
        description="Title for the color bar. Defaults to 'Mean expression in group'.",
    )
    size_title: str = Field(
        default="",
        description="Title for the size legend. Defaults to 'Fraction of cells in group (%)'.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed from gene/category count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed from gene/category count when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class RankGenesGroupsViolinRequest(BaseModel):
    groups_json: str = Field(
        default="",
        description=(
            "JSON array of specific group names to display (e.g. '[\"B cell\",\"T cell\"]'). "
            "Leave empty to display all groups."
        ),
    )
    n_genes: int = Field(
        default=20,
        description="Number of top-ranked genes to show per group. Default: 20.",
    )
    gene_names_json: str = Field(
        default="",
        description=(
            "JSON array of specific gene names to plot instead of the top-ranked genes "
            "(e.g. '[\"CD3E\",\"CD79A\"]'). Leave empty to use the ranked results."
        ),
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata.")
    key: str = Field(
        default="rank_genes_groups",
        description="Key in adata.uns where ranking results are stored. Default: 'rank_genes_groups'.",
    )
    split: bool = Field(
        default=True,
        description="Split violin by group (True) or show a single violin per gene (False). Default: True.",
    )
    density_norm: str = Field(
        default="width",
        description="Violin density normalization: 'width', 'area', or 'count'. Default: 'width'.",
    )
    strip: bool = Field(default=True, description="Overlay individual data points as a strip plot.")
    jitter: float | bool = Field(default=True, description="Add jitter to strip plot points.")
    size: int = Field(default=1, description="Point size for strip plot.")
    title: str = Field(default="", description="Optional plot title.")


class RankGenesGroupsStackedViolinRequest(BaseModel):
    groups_json: str = Field(
        default="",
        description=(
            "JSON array of specific group names to display (e.g. '[\"B cell\",\"T cell\"]'). "
            "Leave empty to display all groups."
        ),
    )
    n_genes: int = Field(
        default=10,
        description="Number of top-ranked genes to show per group. Default: 10.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    var_names_json: str = Field(
        default="",
        description=(
            "JSON array of specific gene names to override the ranked genes "
            "(e.g. '[\"CD3E\",\"CD79A\"]'). Leave empty to use the ranked results."
        ),
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Minimum log fold change threshold. Genes below this are excluded. Default: None.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="Key in adata.uns where ranking results are stored. Default: 'rank_genes_groups'.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap genes and groups axes. Useful when there are many cell types. Default: False.",
    )
    cmap: str = Field(
        default="Blues",
        description="Colormap for the violin fill. Default: 'Blues'.",
    )
    stripplot: bool = Field(
        default=False,
        description="Overlay individual data points as a strip plot. Default: False.",
    )
    jitter: bool = Field(
        default=False,
        description="Add jitter to strip plot points. Default: False.",
    )
    size: int = Field(default=1, description="Point size for strip plot. Default: 1.")
    row_palette: str = Field(
        default="",
        description="Palette for the row colors (gene rows). Leave empty for default.",
    )
    yticklabels: bool = Field(
        default=False,
        description="Show y-axis tick labels. Default: False.",
    )
    standard_scale: str = Field(
        default="",
        description="Standardize data: 'var' (per gene) or 'obs' (per cell). Leave empty for none.",
    )
    vmin: float | None = Field(default=None, description="Minimum color scale value.")
    vmax: float | None = Field(default=None, description="Maximum color scale value.")
    vcenter: float | None = Field(default=None, description="Center of the color scale (for diverging colormaps).")
    colorbar_title: str = Field(default="", description="Title for the colorbar.")
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed from gene/group count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed from gene/group count when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class RankGenesGroupsPlotRequest(BaseModel):
    groups_json: str = Field(
        default="",
        description=(
            "JSON array of specific group names to display (e.g. '[\"B cell\",\"T cell\"]'). "
            "Leave empty to display all groups."
        ),
    )
    n_genes: int = Field(
        default=20,
        description="Number of top-ranked genes to show per group panel. Default: 20.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description=(
            "Key in adata.uns where the ranking results are stored. "
            "Default: 'rank_genes_groups'. Change only if results were stored under a different key."
        ),
    )
    fontsize: int = Field(default=8, description="Font size for gene labels in each panel. Default: 8.")
    ncols: int = Field(default=4, description="Number of group panels per row. Default: 4.")
    sharey: bool = Field(
        default=True,
        description="Share the y-axis range across all panels for easy comparison. Default: True.",
    )
    title: str = Field(default="", description="Optional plot title.")


class RankGenesGroupsHeatmapRequest(BaseModel):
    groups_json: str = Field(
        default="",
        description=(
            "JSON array of specific group names to include (e.g. '[\"B cell\",\"T cell\"]'). "
            "Leave empty to display all groups."
        ),
    )
    n_genes: int | None = Field(
        default=None,
        description=(
            "Number of top genes to show per group. Use a negative value to show down-regulated genes. "
            "Defaults to 10 when not provided."
        ),
    )
    groupby: str = Field(
        default="",
        description="Observation column to group cells by. Inferred from rank_genes_groups params when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to auto-detect or use var_names index.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Minimum log fold-change threshold — genes below this value are excluded.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="Key in adata.uns where rank_genes_groups results are stored.",
    )
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'var' (gene) or 'obs' (cell). Leave empty to skip.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes so genes appear on the x-axis and cell groups on the y-axis.",
    )
    show_gene_labels: bool | None = Field(
        default=None,
        description="Show gene labels on the heatmap. Auto-detected when None.",
    )
    cmap: str = Field(
        default="",
        description="Matplotlib colormap name (e.g. 'bwr', 'viridis', 'RdBu_r'). Leave empty for scanpy default.",
    )
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None,
        description="Center of the color scale — useful for diverging colormaps like 'bwr'.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class RankGenesGroupsDotplotRequest(BaseModel):
    groups_json: str = Field(
        default="",
        description=(
            "JSON array of specific group names to include (e.g. '[\"B cell\",\"T cell\"]'). "
            "Leave empty to display all groups."
        ),
    )
    n_genes: int | None = Field(
        default=None,
        description=(
            "Number of top genes to show per group. Use a negative value to show down-regulated genes. "
            "Defaults to 5 when not provided."
        ),
    )
    groupby: str = Field(
        default="",
        description="Observation column to group cells by. Inferred from rank_genes_groups params when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to auto-detect or use var_names index.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Minimum log fold-change threshold — genes below this value are excluded.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="Key in adata.uns where rank_genes_groups results are stored.",
    )
    values_to_plot: str = Field(
        default="",
        description=(
            "Metric to encode as dot color instead of mean expression. "
            "One of: 'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'. "
            "Leave empty to use mean expression."
        ),
    )
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'var' (gene) or 'obs' (cell). Leave empty to skip.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes so genes appear on the x-axis and cell groups on the y-axis.",
    )
    cmap: str = Field(
        default="",
        description="Matplotlib colormap name (e.g. 'Reds', 'Blues', 'viridis'). Leave empty for scanpy default.",
    )
    dot_max: float | None = Field(
        default=None,
        description="Maximum dot size (fraction of largest dot). Between 0 and 1.",
    )
    dot_min: float | None = Field(
        default=None,
        description="Minimum dot size (fraction of largest dot). Between 0 and 1.",
    )
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None,
        description="Center of the color scale — useful for diverging colormaps.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class RankGenesGroupsMatrixplotRequest(BaseModel):
    groups_json: str = Field(
        default="",
        description=(
            "JSON array of specific group names to include (e.g. '[\"B cell\",\"T cell\"]'). "
            "Leave empty to display all groups."
        ),
    )
    n_genes: int | None = Field(
        default=None,
        description=(
            "Number of top genes to show per group. Use a negative value to show down-regulated genes. "
            "Defaults to 3 when not provided."
        ),
    )
    groupby: str = Field(
        default="",
        description="Observation column to group cells by. Inferred from rank_genes_groups params when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to auto-detect or use var_names index.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Minimum log fold-change threshold — genes below this value are excluded.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="Key in adata.uns where rank_genes_groups results are stored.",
    )
    values_to_plot: str = Field(
        default="",
        description=(
            "Metric to display instead of mean expression. "
            "One of: 'scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'log10_pvals', 'log10_pvals_adj'. "
            "Leave empty to use mean expression."
        ),
    )
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'var' (gene) or 'obs' (cell). Leave empty to skip.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes so genes appear on the x-axis and cell groups on the y-axis.",
    )
    cmap: str = Field(
        default="",
        description="Matplotlib colormap name (e.g. 'bwr', 'viridis', 'RdBu_r'). Leave empty for scanpy default.",
    )
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None,
        description="Center of the color scale — useful for diverging colormaps like 'bwr'.",
    )
    colorbar_title: str = Field(
        default="",
        description="Title for the color bar (e.g. 'Mean expression\\nin group').",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class RankGenesGroupsTracksplotRequest(BaseModel):
    groups_json: str = Field(
        default="",
        description=(
            "JSON array of specific group names to include (e.g. '[\"B cell\",\"T cell\"]'). "
            "Leave empty to display all groups."
        ),
    )
    n_genes: int | None = Field(
        default=None,
        description=(
            "Number of top genes to show per group. Use a negative value to show down-regulated genes. "
            "Defaults to 10 when not provided."
        ),
    )
    groupby: str = Field(
        default="",
        description="Observation column to group cells by. Inferred from rank_genes_groups params when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to auto-detect or use var_names index.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Minimum log fold-change threshold — genes below this value are excluded.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="Key in adata.uns where rank_genes_groups results are stored.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Use the raw attribute of adata if present. Auto-detected when None.",
    )
    log: bool = Field(
        default=False,
        description="Plot values on a logarithmic scale.",
    )
    layer: str = Field(
        default="",
        description="AnnData layer to use instead of X or raw. Leave empty to use default.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class CorrelationMatrixRequest(BaseModel):
    groupby: str = Field(
        default="",
        description="Observation column to group cells by for the correlation matrix. Inferred when empty.",
    )
    show_correlation_numbers: bool = Field(
        default=False,
        description="Overlay the correlation value on each cell of the matrix.",
    )
    dendrogram: bool | None = Field(
        default=None,
        description="Add a hierarchical clustering dendrogram. Auto-detected when None.",
    )
    cmap: str = Field(
        default="",
        description="Matplotlib colormap name (e.g. 'RdBu_r', 'coolwarm'). Leave empty for scanpy default.",
    )
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None,
        description="Center of the color scale — useful for diverging colormaps.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed as a square based on category count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed as a square based on category count when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class EmbeddingPlotRequest(BaseModel):
    basis: str = Field(
        description=(
            "Name of the embedding to plot. Examples: 'umap', 'tsne', 'pca', 'diffmap', 'draw_graph_fa'. "
            "Can be specified with or without the 'X_' prefix. "
            "Call print_adata_obs or check active dataset info to see which embeddings are available."
        ),
    )
    color_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of observation columns or gene names to color by. "
            "Examples: '[\"author_cell_type\"]' or '[\"CD3E\",\"CD79A\"]'. "
            "Leave empty to auto-color by the inferred cell type column."
        ),
    )
    components: str = Field(
        default="",
        description=(
            "Components to plot, e.g. '1,2' or '2,3'. "
            "Leave empty to use the default (first two components)."
        ),
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata for gene expression.")
    layer: str = Field(default="", description="AnnData layer to use for coloring. Leave empty for X or raw.")
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position: 'right margin', 'on data', 'best', 'upper right', etc. Default: 'right margin'.",
    )
    legend_fontsize: float | None = Field(default=None, description="Legend font size in points. Auto-scaled when None.")
    legend_fontweight: str = Field(default="bold", description="Legend font weight. Default: 'bold'.")
    colorbar_loc: str = Field(default="right", description="Colorbar position for continuous color keys. Default: 'right'.")
    color_map: str = Field(
        default="",
        description="Colormap for continuous variables (e.g. 'viridis', 'RdBu_r'). Leave empty for default.",
    )
    palette: str = Field(
        default="",
        description="Color palette for categorical variables (e.g. 'tab20'). Leave empty for default.",
    )
    na_color: str = Field(default="lightgray", description="Color for cells with missing values. Default: 'lightgray'.")
    na_in_legend: bool = Field(default=True, description="Include NA category in the legend.")
    size: float | None = Field(
        default=None,
        description="Point size. Auto-calculated from cell count when None (120000 / n_cells).",
    )
    frameon: bool | None = Field(default=None, description="Show plot frame. Uses global settings when None.")
    vmin: str | None = Field(
        default=None,
        description="Lower color scale limit. Supports percentile syntax e.g. 'p1.5'. Leave empty for auto.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color scale limit. Supports percentile syntax e.g. 'p98'. Leave empty for auto.",
    )
    vcenter: float | None = Field(default=None, description="Center value for diverging colormaps.")
    add_outline: bool = Field(default=False, description="Add a border outline around each group of dots.")
    sort_order: bool = Field(default=True, description="Plot higher-value points on top of lower-value ones.")
    edges: bool = Field(default=False, description="Overlay neighborhood graph edges on the embedding.")
    edges_width: float = Field(default=0.1, description="Width of graph edges when edges=True.")
    edges_color: str = Field(default="grey", description="Color of graph edges when edges=True.")
    groups_json: str = Field(
        default="",
        description="JSON array of category values to highlight; all others are greyed out (e.g. '[\"B cell\"]').",
    )
    projection: str = Field(default="2d", description="Projection type: '2d' or '3d'. Default: '2d'.")
    ncols: int = Field(default=4, description="Number of panels per row when multiple color keys are given.")
    title: str = Field(default="", description="Optional plot title.")


class DiffmapPlotRequest(BaseModel):
    color_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of observation columns or gene names to color by. "
            "Examples: '[\"author_cell_type\"]' or '[\"CD3E\",\"CD79A\"]'. "
            "Leave empty to auto-color by the inferred cell type column."
        ),
    )
    components: str = Field(
        default="",
        description=(
            "Diffusion components to plot, e.g. '1,2' (DC1 vs DC2) or '2,3' (DC2 vs DC3). "
            "Leave empty to use the default (first two components)."
        ),
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata for gene expression.")
    layer: str = Field(default="", description="AnnData layer to use for coloring. Leave empty for X or raw.")
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position: 'right margin', 'on data', 'best', 'upper right', etc. Default: 'right margin'.",
    )
    legend_fontsize: float | None = Field(default=None, description="Legend font size in points. Auto-scaled when None.")
    legend_fontweight: str = Field(default="bold", description="Legend font weight. Default: 'bold'.")
    colorbar_loc: str = Field(default="right", description="Colorbar position for continuous color keys. Default: 'right'.")
    color_map: str = Field(
        default="",
        description="Colormap for continuous variables (e.g. 'viridis', 'RdBu_r'). Leave empty for default.",
    )
    palette: str = Field(
        default="",
        description="Color palette for categorical variables (e.g. 'tab20'). Leave empty for default.",
    )
    na_color: str = Field(default="lightgray", description="Color for cells with missing values. Default: 'lightgray'.")
    na_in_legend: bool = Field(default=True, description="Include NA category in the legend.")
    size: float | None = Field(
        default=None,
        description="Point size. Auto-calculated from cell count when None (120000 / n_cells).",
    )
    frameon: bool | None = Field(default=None, description="Show plot frame. Uses global settings when None.")
    vmin: str | None = Field(
        default=None,
        description="Lower color scale limit. Supports percentile syntax e.g. 'p1.5'. Leave empty for auto.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color scale limit. Supports percentile syntax e.g. 'p98'. Leave empty for auto.",
    )
    vcenter: float | None = Field(default=None, description="Center value for diverging colormaps.")
    add_outline: bool = Field(default=False, description="Add a border outline around each group of dots.")
    sort_order: bool = Field(default=True, description="Plot higher-value points on top of lower-value ones.")
    edges: bool = Field(default=False, description="Overlay neighborhood graph edges on the embedding.")
    edges_width: float = Field(default=0.1, description="Width of graph edges when edges=True.")
    edges_color: str = Field(default="grey", description="Color of graph edges when edges=True.")
    groups_json: str = Field(
        default="",
        description="JSON array of category values to highlight; all others are greyed out (e.g. '[\"B cell\"]').",
    )
    ncols: int = Field(default=4, description="Number of panels per row when multiple color keys are given.")
    title: str = Field(default="", description="Optional plot title.")


class UmapPlotRequest(BaseModel):
    color_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of observation columns or gene names to color by. "
            "Examples: '[\"author_cell_type\"]' or '[\"CD3E\",\"CD79A\"]'. "
            "Leave empty to auto-color by the inferred cell type column."
        ),
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata for gene expression.")
    layer: str = Field(default="", description="AnnData layer to use for coloring. Leave empty for X or raw.")
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position: 'right margin', 'on data', 'best', 'upper right', etc. Default: 'right margin'.",
    )
    legend_fontsize: float | None = Field(default=None, description="Legend font size in points. Auto-scaled when None.")
    legend_fontweight: str = Field(default="bold", description="Legend font weight. Default: 'bold'.")
    colorbar_loc: str = Field(default="right", description="Colorbar position for continuous color keys. Default: 'right'.")
    color_map: str = Field(
        default="",
        description="Colormap for continuous variables (e.g. 'viridis', 'RdBu_r'). Leave empty for default.",
    )
    palette: str = Field(
        default="",
        description="Color palette for categorical variables (e.g. 'tab20'). Leave empty for default.",
    )
    na_color: str = Field(default="lightgray", description="Color for cells with missing values. Default: 'lightgray'.")
    na_in_legend: bool = Field(default=True, description="Include NA category in the legend.")
    size: float | None = Field(
        default=None,
        description="Point size. Auto-calculated from cell count when None (120000 / n_cells).",
    )
    frameon: bool | None = Field(default=None, description="Show plot frame. Uses global settings when None.")
    vmin: str | None = Field(
        default=None,
        description="Lower color scale limit. Supports percentile syntax e.g. 'p1.5'. Leave empty for auto.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color scale limit. Supports percentile syntax e.g. 'p98'. Leave empty for auto.",
    )
    vcenter: float | None = Field(default=None, description="Center value for diverging colormaps.")
    add_outline: bool = Field(default=False, description="Add a border outline around each group of dots.")
    sort_order: bool = Field(default=True, description="Plot higher-value points on top of lower-value ones.")
    edges: bool = Field(default=False, description="Overlay neighborhood graph edges on the embedding.")
    edges_width: float = Field(default=0.1, description="Width of graph edges when edges=True.")
    edges_color: str = Field(default="grey", description="Color of graph edges when edges=True.")
    groups_json: str = Field(
        default="",
        description="JSON array of category values to highlight; all others are greyed out (e.g. '[\"B cell\"]').",
    )
    ncols: int = Field(default=4, description="Number of panels per row when multiple color keys are given.")
    title: str = Field(default="", description="Optional plot title.")


class TsnePlotRequest(BaseModel):
    color_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of observation columns or gene names to color by. "
            "Examples: '[\"author_cell_type\"]' or '[\"CD3E\",\"CD79A\"]'. "
            "Leave empty to auto-color by the inferred cell type column."
        ),
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata for gene expression.")
    layer: str = Field(default="", description="AnnData layer to use for coloring. Leave empty for X or raw.")
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position: 'right margin', 'on data', 'best', 'upper right', etc. Default: 'right margin'.",
    )
    legend_fontsize: float | None = Field(default=None, description="Legend font size in points. Auto-scaled when None.")
    legend_fontweight: str = Field(default="bold", description="Legend font weight. Default: 'bold'.")
    colorbar_loc: str = Field(default="right", description="Colorbar position. Default: 'right'.")
    color_map: str = Field(
        default="",
        description="Colormap for continuous variables (e.g. 'viridis', 'RdBu_r'). Leave empty for default.",
    )
    palette: str = Field(
        default="",
        description="Color palette for categorical variables (e.g. 'tab20'). Leave empty for default.",
    )
    na_color: str = Field(default="lightgray", description="Color for cells with missing values. Default: 'lightgray'.")
    na_in_legend: bool = Field(default=True, description="Show NA category in the legend.")
    size: float | None = Field(
        default=None,
        description="Point size. Auto-calculated from cell count when None (120000 / n_cells).",
    )
    frameon: bool | None = Field(default=None, description="Show plot frame. Uses global settings when None.")
    vmin: str | None = Field(
        default=None,
        description="Lower color scale limit. Supports percentile syntax e.g. 'p1.5'. Leave empty for auto.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color scale limit. Supports percentile syntax e.g. 'p98'. Leave empty for auto.",
    )
    vcenter: float | None = Field(default=None, description="Center value for diverging colormaps.")
    add_outline: bool = Field(default=False, description="Add a border outline around each group of dots.")
    sort_order: bool = Field(default=True, description="Plot higher-value points on top of lower-value ones.")
    edges: bool = Field(default=False, description="Overlay neighborhood graph edges on the embedding.")
    edges_width: float = Field(default=0.1, description="Width of graph edges when edges=True.")
    edges_color: str = Field(default="grey", description="Color of graph edges when edges=True.")
    groups_json: str = Field(
        default="",
        description="JSON array of category values to highlight; all others are greyed out (e.g. '[\"B cell\"]').",
    )
    ncols: int = Field(default=4, description="Number of panels per row when multiple color keys are given.")
    title: str = Field(default="", description="Optional plot title.")


class DendrogramRequest(BaseModel):
    groupby: str = Field(
        default="author_cell_type",
        description="Observation column to build the dendrogram from. Default: 'author_cell_type'.",
    )
    dendrogram_key: str = Field(
        default="",
        description=(
            "Key in adata.uns where the precomputed dendrogram is stored. "
            "Defaults to 'dendrogram_{groupby}'. Leave empty to use the default."
        ),
    )
    orientation: str = Field(
        default="top",
        description="Direction the tree grows from: 'top', 'bottom', 'left', or 'right'. Default: 'top'.",
    )
    remove_labels: bool = Field(
        default=False,
        description="Hide category labels on the dendrogram leaves.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed from category count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed from category count when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class ClustermapRequest(BaseModel):
    markers_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of genes to subset before clustering. "
            "Strongly recommended — running on all genes is very slow for large datasets. "
            "If empty, session markers set by set_markers are used. "
            "If neither is provided, all genes in the dataset are used (may be slow)."
        ),
    )
    obs_keys: str = Field(
        default="",
        description=(
            "Optional observation column to color-code rows by (e.g. 'author_cell_type', 'disease'). "
            "Only one key is supported. Leave empty for no row coloring."
        ),
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata.")
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'row' or 'col' before clustering. Leave empty to skip.",
    )
    z_score: int | None = Field(
        default=None,
        description="Compute z-score along rows (0) or columns (1) before clustering. Leave None to skip.",
    )
    method: str = Field(
        default="average",
        description="Linkage method for hierarchical clustering (e.g. 'average', 'complete', 'ward').",
    )
    metric: str = Field(
        default="euclidean",
        description="Distance metric for clustering (e.g. 'euclidean', 'correlation', 'cosine').",
    )
    cmap: str = Field(default="viridis", description="Colormap for the heatmap cells. Default: 'viridis'.")
    figsize_width: float | None = Field(default=None, description="Figure width in inches.")
    figsize_height: float | None = Field(default=None, description="Figure height in inches.")
    title: str = Field(default="", description="Optional plot title.")


class MatrixplotRequest(BaseModel):
    markers_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of genes to display in the matrix plot. "
            "If empty, the session markers set by set_markers are used. "
            "Example: '[\"C1QA\",\"PSAP\",\"CD79A\"]'."
        ),
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Observation column to group cells by. Default: 'author_cell_type'.",
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata.")
    log: bool = Field(default=False, description="Plot on a logarithmic scale.")
    num_categories: int = Field(
        default=7,
        description="Number of categories when groupby is continuous.",
    )
    categories_order_json: str = Field(
        default="",
        description="JSON array specifying display order of groupby categories (e.g. '[\"B cell\",\"T cell\"]').",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    layer: str = Field(default="", description="AnnData layer to plot. Leave empty to use X or raw.")
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'var' (gene) or 'obs' (cell). Leave empty to skip.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes: cell types on x-axis, genes on y-axis. Useful when there are many cell types.",
    )
    cmap: str = Field(default="viridis", description="Colormap for cell fill color. Default: 'viridis'.")
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None, description="Center of the color scale (useful for diverging colormaps)."
    )
    colorbar_title: str = Field(
        default="",
        description="Title for the color bar. Defaults to 'Mean expression in group'.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed from gene/category count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed from gene/category count when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class StackedViolinRequest(BaseModel):
    markers_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of genes to display in the stacked violin plot. "
            "If empty, the session markers set by set_markers are used. "
            "Example: '[\"C1QA\",\"PSAP\",\"CD79A\"]'."
        ),
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Observation column to group cells by. Default: 'author_cell_type'.",
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata.")
    log: bool = Field(default=False, description="Plot on a logarithmic scale.")
    num_categories: int = Field(
        default=7,
        description="Number of categories when groupby is continuous.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    layer: str = Field(default="", description="AnnData layer to plot. Leave empty to use X or raw.")
    standard_scale: str = Field(
        default="",
        description="Standardize values to [0,1] per 'var' (gene) or 'obs' (cell). Leave empty to skip.",
    )
    categories_order_json: str = Field(
        default="",
        description="JSON array specifying display order of groupby categories (e.g. '[\"B cell\",\"T cell\"]').",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes: genes on x-axis, groups on y-axis. Useful when there are many groups.",
    )
    vmin: float | None = Field(default=None, description="Lower color scale limit.")
    vmax: float | None = Field(default=None, description="Upper color scale limit.")
    vcenter: float | None = Field(
        default=None, description="Center of the color scale (useful for diverging colormaps)."
    )
    cmap: str = Field(default="Blues", description="Colormap for violin fill color. Default: 'Blues'.")
    stripplot: bool = Field(default=False, description="Overlay a strip plot (individual data points) on each violin.")
    jitter: float | bool = Field(default=False, description="Add jitter to the strip plot.")
    size: int = Field(default=1, description="Size of the jitter points.")
    row_palette: str = Field(
        default="",
        description="Color palette for violin rows (e.g. 'tab20'). Leave empty to use cmap.",
    )
    yticklabels: bool = Field(default=False, description="Show y-axis tick labels on each violin row.")
    colorbar_title: str = Field(
        default="",
        description="Title for the color bar. Defaults to 'Median expression in group'.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed from gene/category count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed from gene/category count when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class TracksplotRequest(BaseModel):
    markers_json: str | list[str] = Field(
        default="[]",
        description=(
            "JSON array or comma-separated list of genes to display in the tracks plot. "
            "If empty, the session markers set by set_markers are used. "
            "Example: '[\"C1QA\",\"PSAP\",\"CD79A\"]'."
        ),
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Observation column to group cells by. Default: 'author_cell_type'.",
    )
    use_raw: bool | None = Field(default=None, description="Use the raw attribute of adata.")
    log: bool = Field(default=False, description="Plot on a logarithmic scale.")
    dendrogram: bool = Field(
        default=False,
        description="Add a dendrogram based on hierarchical clustering of the groupby categories.",
    )
    gene_symbols: str = Field(
        default="",
        description="Column in .var that contains gene symbols. Leave empty to use var_names index.",
    )
    layer: str = Field(default="", description="AnnData layer to plot. Leave empty to use X or raw.")
    figsize_width: float | None = Field(
        default=None,
        description="Figure width in inches. Auto-computed from gene/category count when None.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Figure height in inches. Auto-computed from gene/category count when None.",
    )
    title: str = Field(default="", description="Optional plot title.")


class ExecuteDatasetCommandRequest(BaseModel):
    command: str = Field(
        default="print(adata.obs)",
        description="Observation inspection command. Supported value: `print(adata.obs)`.",
    )


def _settings():
    from gp_config import load_settings

    return load_settings()


def _tool_api_key() -> str:
    return str(os.getenv("API_KEY", "")).strip()


def _require_api_key(credentials: HTTPAuthorizationCredentials | None = Depends(_bearer)) -> None:
    expected = _tool_api_key()
    if not expected:
        return
    provided = credentials.credentials if credentials else ""
    if not credentials or credentials.scheme.lower() != "bearer" or provided != expected:
        raise HTTPException(status_code=401, detail="Unauthorized")


def _public_assets_url(file_path: str) -> str | None:
    candidate = Path(str(file_path)).expanduser()
    for base in (Path("/code/out"), Path("/work/out")):
        try:
            rel = candidate.resolve().relative_to(base)
            return f"http://localhost/assets/{rel.as_posix()}"
        except Exception:
            continue
    return None


def _output_markdown(output_file: str, output_url: str | None) -> str:
    target = str(output_url or output_file).strip()
    if not target:
        return ""
    return f"![Plot]({target})"


def _is_plain_umap_request(payload: GeneratePlotRequest, plot_request: PlotRequest) -> bool:
    return (
        plot_request.plot_type == "umap"
        and not plot_request.color
        and not plot_request.genes
        and not (plot_request.groupby or "").strip()
        and not (plot_request.gene_symbols_column or "").strip()
        and not (payload.title or "").strip()
    )


def _parse_string_list(value: str | list[str]) -> list[str]:
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]

    raw = str(value or "").strip()
    if not raw:
        return []
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        parsed = None
    if isinstance(parsed, list):
        return [str(item).strip() for item in parsed if str(item).strip()]
    return [item.strip() for item in raw.split(",") if item.strip()]


def _build_plot_request(payload: GeneratePlotRequest) -> PlotRequest:
    return PlotRequest(
        plot_type=str(payload.plot_type or "umap").strip() or "umap",
        color=_parse_string_list(payload.color_json),
        genes=_parse_string_list(payload.genes_json),
        groupby=str(payload.groupby or "").strip() or None,
        gene_symbols_column=str(payload.gene_symbols_column or "").strip() or None,
        title=str(payload.title or "").strip() or None,
    )


def _normalize_dataset_command(command: str) -> str:
    return "".join(str(command or "").lower().split())


def _plotter() -> ScanpyPlotExecutor:
    global _PLOTTER
    if _PLOTTER is None:
        _PLOTTER = ScanpyPlotExecutor(_settings().output_dir)
    return _PLOTTER


def create_app(store: GenoPixelCatalogStore | None = None) -> FastAPI:
    catalog_store = store or GenoPixelCatalogStore()
    app = FastAPI(title="GenoPixel Runtime API", version="0.2.0")

    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=False,
        allow_methods=["GET", "POST"],
        allow_headers=["*"],
    )

    @app.exception_handler(CatalogLoadError)
    async def handle_catalog_load_error(_, exc: CatalogLoadError) -> JSONResponse:
        return JSONResponse(status_code=503, content={"detail": str(exc)})

    @app.exception_handler(DatasetLoadError)
    async def handle_dataset_load_error(_, exc: DatasetLoadError) -> JSONResponse:
        return JSONResponse(status_code=500, content={"detail": str(exc)})

    @app.get(f"{CATALOG_PREFIX}/health", include_in_schema=False)
    async def catalog_health() -> dict:
        return catalog_store.get_health_payload()

    @app.get(f"{CATALOG_PREFIX}/catalog", include_in_schema=False)
    async def catalog() -> dict:
        return catalog_store.get_catalog_payload()

    @app.get(f"{CATALOG_PREFIX}/datasets/{{all_excel_row}}", include_in_schema=False)
    async def dataset_detail(all_excel_row: int) -> dict:
        try:
            return catalog_store.get_dataset_payload(all_excel_row)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=f"Dataset row {all_excel_row} was not found.") from exc

    @app.post(f"{CATALOG_PREFIX}/datasets/{{all_excel_row}}/analyze", include_in_schema=False)
    async def analyze_dataset(all_excel_row: int, payload: AnalyzeDatasetRequest) -> dict:
        try:
            return catalog_store.analyze_dataset(
                all_excel_row,
                h5ad_path=payload.h5ad_path,
                multiple_excel_row=payload.multiple_excel_row,
            )
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=f"Dataset row {all_excel_row} was not found.") from exc

    @app.get(f"{RUNTIME_PREFIX}/active-dataset", include_in_schema=False)
    async def active_dataset() -> dict:
        return RUNTIME_STATE.get_active_dataset_payload()

    @app.post("/health", dependencies=[Depends(_require_api_key)], include_in_schema=False)
    async def tool_health() -> dict:
        now = time.time()
        return {
            "ok": True,
            "name": "genopixel-runtime",
            "time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(now)),
            "uptime_seconds": int(now - START_TIME),
            "versions": {
                "python": sys.version.split()[0],
            },
            "active_dataset": RUNTIME_STATE.get_active_dataset_payload(),
        }

    @app.post(
        "/generate_scanpy_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_scanpy_plot",
        summary="Generate a Scanpy plot for the currently loaded GenoPixel dataset",
        description=(
            "Use this tool whenever the user asks for a plot or visualization of the dataset that was already "
            "loaded from the GenoPixel browser. Do not ask for file paths. For a plain UMAP request, call this "
            "tool with plot_type='umap', color_json='[]', and genes_json='[]' unless the user asked for a "
            "specific coloring field or gene."
        ),
    )
    async def generate_scanpy_plot(payload: GeneratePlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {
                "ok": False,
                "status": "no_active_dataset",
                "message": str(exc),
            }

        plot_request = _build_plot_request(payload)
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        try:
            result = plotter.run(adata, plot_request)
        except Exception as exc:
            return {
                "ok": False,
                "status": "plot_error",
                "message": str(exc),
                "active_dataset": active,
            }
        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        plot_payload: dict[str, Any] = {
            "plot_type": result.plot_type,
            "embedding_basis": result.embedding_basis,
            "color_columns": result.color_columns,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
            "resolved_coloring_label": result.resolved_coloring_label,
            "display_plot_type": result.display_plot_type,
            "rank_genes_groups_computed": result.rank_genes_groups_computed,
            "rank_genes_groups_notice": result.rank_genes_groups_notice,
            "output_file": output_file,
            "output_markdown": output_markdown,
        }
        if output_url:
            plot_payload["output_url"] = output_url

        canonical_response_markdown = build_canonical_response_markdown(active, plot_payload, output_markdown)
        plain_umap_markdown = None
        if _is_plain_umap_request(payload, plot_request):
            plain_umap_markdown = canonical_response_markdown

        return {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "plot": plot_payload,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "output_url": output_url,
            "canonical_response_markdown": canonical_response_markdown,
            "plain_umap_response_markdown": plain_umap_markdown,
            "rank_genes_groups_computed": result.rank_genes_groups_computed,
            "rank_genes_groups_notice": result.rank_genes_groups_notice,
        }

    @app.post(
        "/generate_heatmap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_heatmap_plot",
        summary="Generate a customizable heatmap for the active dataset",
        description=(
            "Plot gene expression as a heatmap grouped by any observation column. "
            "Genes are taken from markers_json; if empty, session markers set by set_markers are used. "
            "If neither is provided, the tool returns an error asking for genes. "
            "groupby defaults to 'author_cell_type'. "
            "All scanpy.pl.heatmap parameters are exposed."
        ),
    )
    async def generate_heatmap_plot(payload: HeatmapPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        var_names = _parse_string_list(payload.markers_json)
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        # Fall back to session markers if no genes provided
        if not var_names:
            var_names = plotter.get_markers()
        if not var_names:
            return {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_heatmap(
                adata,
                var_names,
                groupby=str(payload.groupby or "").strip() or None,
                use_raw=payload.use_raw,
                log=bool(payload.log),
                num_categories=int(payload.num_categories),
                dendrogram=bool(payload.dendrogram),
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                layer=str(payload.layer or "").strip() or None,
                standard_scale=str(payload.standard_scale or "").strip() or None,
                swap_axes=bool(payload.swap_axes),
                show_gene_labels=payload.show_gene_labels,
                figsize=figsize,
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "var_names_used": result.resolved_genes,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/set_markers",
        dependencies=[Depends(_require_api_key)],
        operation_id="set_markers",
        summary="Set a session-level gene marker list",
        description=(
            "Store a list of genes as the session default 'markers'. "
            "Once set, any tool that accepts genes (violin, dotplot, matrixplot, heatmap) "
            "will automatically use these markers when no genes are explicitly specified. "
            "Example: markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']. "
            "Call this once at the start of a session to avoid repeating the gene list."
        ),
    )
    async def set_markers(payload: SetMarkersRequest) -> dict[str, Any]:
        markers = _parse_string_list(payload.markers_json)
        if not markers:
            return {"ok": False, "status": "error", "message": "markers_json must contain at least one gene."}
        plotter = _plotter()
        plotter.set_markers(markers)
        return {
            "ok": True,
            "status": "success",
            "markers": markers,
            "count": len(markers),
            "message": f"Stored {len(markers)} marker gene(s). They will be used as the default gene list for violin, dotplot, matrixplot, and heatmap plots.",
        }

    @app.post(
        "/get_markers",
        dependencies=[Depends(_require_api_key)],
        operation_id="get_markers",
        summary="Get the current session-level gene marker list",
        description="Return the gene markers currently stored in this session.",
    )
    async def get_markers() -> dict[str, Any]:
        plotter = _plotter()
        markers = plotter.get_markers()
        return {
            "ok": True,
            "markers": markers,
            "count": len(markers),
        }

    @app.post(
        "/generate_violin_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_violin_plot",
        summary="Generate a customizable violin plot for the active dataset",
        description=(
            "Plot gene expression or obs values as violins grouped by any observation column. "
            "Exposes all scanpy.pl.violin parameters. "
            "groupby defaults to 'author_cell_type'. rotation defaults to 45 degrees. "
            "Use print_adata_obs to discover available obs columns."
        ),
    )
    async def generate_violin_plot(payload: ViolinPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        keys = _parse_string_list(payload.keys_json)
        # Fall back to session markers if no keys provided
        if not keys:
            keys = _plotter().get_markers()
        if not keys:
            return {
                "ok": False,
                "status": "error",
                "message": "No genes specified and no session markers are set. Provide keys_json or call set_markers first.",
            }

        order = _parse_string_list(payload.order_json) if payload.order_json else None
        layer = str(payload.layer or "").strip() or None
        ylabel = str(payload.ylabel or "").strip() or None
        xlabel = str(payload.xlabel or "").strip()

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        try:
            result = plotter.run_violin(
                adata,
                keys,
                groupby=str(payload.groupby or "").strip() or None,
                rotation=float(payload.rotation),
                log=bool(payload.log),
                use_raw=payload.use_raw,
                stripplot=bool(payload.stripplot),
                jitter=payload.jitter,
                size=int(payload.size),
                layer=layer,
                density_norm=str(payload.density_norm or "width"),
                order=order,
                multi_panel=payload.multi_panel,
                xlabel=xlabel,
                ylabel=ylabel,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/cell_counts_barplot",
        dependencies=[Depends(_require_api_key)],
        operation_id="cell_counts_barplot",
        summary="Plot cell counts for a categorical observation column",
        description=(
            "Count cells by any adata.obs column (e.g. 'author_cell_type', 'disease', 'tissue', 'donor_id') "
            "and generate a bar plot sorted from most to fewest cells. "
            "Figure size adapts automatically to the number of categories and label length. "
            "Use print_adata_obs first to see which columns are available."
        ),
    )
    async def cell_counts_barplot(payload: CellCountsBarplotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {
                "ok": False,
                "status": "no_active_dataset",
                "message": str(exc),
            }

        from gp_models import PlotRequest as _PlotRequest
        plot_request = _PlotRequest(
            plot_type="cell_counts_barplot",
            groupby=str(payload.groupby or "").strip() or None,
            title=str(payload.title or "").strip() or None,
        )
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        try:
            result = plotter.run(adata, plot_request)
        except Exception as exc:
            return {
                "ok": False,
                "status": "plot_error",
                "message": str(exc),
                "active_dataset": active,
            }

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_dotplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_dotplot_plot",
        summary="Generate a customizable dot plot for the active dataset",
        description=(
            "Plot gene expression as a dot plot grouped by any observation column. "
            "Dot size encodes the fraction of cells expressing the gene; "
            "dot color encodes mean expression. "
            "Genes are taken from markers_json; if empty, session markers set by set_markers are used. "
            "If neither is provided, the tool returns an error asking for genes. "
            "groupby defaults to 'author_cell_type'. "
            "All scanpy.pl.dotplot parameters are exposed."
        ),
    )
    async def generate_dotplot_plot(payload: DotplotPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        var_names = _parse_string_list(payload.markers_json)
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        if not var_names:
            var_names = plotter.get_markers()
        if not var_names:
            return {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        categories_order = _parse_string_list(payload.categories_order_json) if payload.categories_order_json else None

        try:
            result = plotter.run_dotplot(
                adata,
                var_names,
                groupby=str(payload.groupby or "").strip() or None,
                use_raw=payload.use_raw,
                log=bool(payload.log),
                num_categories=int(payload.num_categories),
                categories_order=categories_order,
                expression_cutoff=float(payload.expression_cutoff),
                mean_only_expressed=bool(payload.mean_only_expressed),
                standard_scale=str(payload.standard_scale or "").strip() or None,
                dendrogram=bool(payload.dendrogram),
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                layer=str(payload.layer or "").strip() or None,
                swap_axes=bool(payload.swap_axes),
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                cmap=str(payload.cmap or "Reds").strip() or "Reds",
                dot_max=payload.dot_max,
                dot_min=payload.dot_min,
                smallest_dot=float(payload.smallest_dot),
                colorbar_title=str(payload.colorbar_title or "").strip() or None,
                size_title=str(payload.size_title or "").strip() or None,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "var_names_used": result.resolved_genes,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_rank_genes_groups_violin",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_violin",
        summary="Plot ranked marker genes as violins from precomputed sc.tl.rank_genes_groups results",
        description=(
            "Plot expression distributions of top-ranked genes per group as violin plots. "
            "With split=True (default), each violin is split to compare expression in the target group vs all others. "
            "This tool does NOT auto-compute rank_genes_groups — if results are missing it returns "
            "'no_rank_genes_groups' with instructions. "
            "Use gene_names_json to override ranked genes with a specific list. "
            "All scanpy.pl.rank_genes_groups_violin parameters are exposed."
        ),
    )
    async def generate_rank_genes_groups_violin(payload: RankGenesGroupsViolinRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None
        gene_names = _parse_string_list(payload.gene_names_json) if payload.gene_names_json else None

        try:
            result = plotter.run_rank_genes_groups_violin(
                adata,
                groups=groups,
                n_genes=int(payload.n_genes),
                gene_names=gene_names,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                use_raw=payload.use_raw,
                key=str(payload.key or "rank_genes_groups").strip() or "rank_genes_groups",
                split=bool(payload.split),
                density_norm=str(payload.density_norm or "width").strip() or "width",
                strip=bool(payload.strip),
                jitter=payload.jitter,
                size=int(payload.size),
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_rank_genes_groups" if "not found" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_rank_genes_groups_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_plot",
        summary="Plot ranked marker genes per group from precomputed sc.tl.rank_genes_groups results",
        description=(
            "Plot the top-ranked genes for each group as computed by sc.tl.rank_genes_groups. "
            "Each panel shows score vs. gene for one group. "
            "This tool does NOT auto-compute rank_genes_groups — if results are not present in adata.uns "
            "it returns a 'no_rank_genes_groups' error with instructions on how to compute them. "
            "Use groups_json to restrict to specific groups. "
            "All scanpy.pl.rank_genes_groups parameters are exposed."
        ),
    )
    async def generate_rank_genes_groups_plot(payload: RankGenesGroupsPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None

        try:
            result = plotter.run_rank_genes_groups(
                adata,
                groups=groups,
                n_genes=int(payload.n_genes),
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                key=str(payload.key or "rank_genes_groups").strip() or "rank_genes_groups",
                fontsize=int(payload.fontsize),
                ncols=int(payload.ncols),
                sharey=bool(payload.sharey),
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_rank_genes_groups" if "not found" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_rank_genes_groups_dotplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_dotplot_plot",
        summary="Plot ranked marker genes as a customizable dot plot",
        description=(
            "Render sc.pl.rank_genes_groups_dotplot for the active dataset. "
            "Requires precomputed sc.tl.rank_genes_groups results in adata.uns. "
            "Exposes all scanpy.pl.rank_genes_groups_dotplot parameters including "
            "n_genes, groups, values_to_plot, standard_scale, dendrogram, swap_axes, "
            "cmap, dot_max/dot_min, vmin/vmax/vcenter, and figsize."
        ),
    )
    async def generate_rank_genes_groups_dotplot_plot(
        payload: RankGenesGroupsDotplotRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None
        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_rank_genes_groups_dotplot(
                adata,
                groups=groups or None,
                n_genes=payload.n_genes,
                groupby=str(payload.groupby or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                min_logfoldchange=payload.min_logfoldchange,
                key=str(payload.key or "rank_genes_groups").strip() or "rank_genes_groups",
                values_to_plot=str(payload.values_to_plot or "").strip() or None,
                standard_scale=str(payload.standard_scale or "").strip() or None,
                dendrogram=bool(payload.dendrogram),
                swap_axes=bool(payload.swap_axes),
                cmap=str(payload.cmap or "").strip() or None,
                dot_max=payload.dot_max,
                dot_min=payload.dot_min,
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_rank_genes_groups" if "not found" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_rank_genes_groups_tracksplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_tracksplot_plot",
        summary="Plot ranked marker genes as a customizable tracks plot",
        description=(
            "Render sc.pl.rank_genes_groups_tracksplot for the active dataset. "
            "Requires precomputed sc.tl.rank_genes_groups results in adata.uns. "
            "Exposes all scanpy.pl.rank_genes_groups_tracksplot parameters including "
            "n_genes, groups, dendrogram, use_raw, log, layer, and figsize."
        ),
    )
    async def generate_rank_genes_groups_tracksplot_plot(
        payload: RankGenesGroupsTracksplotRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None
        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_rank_genes_groups_tracksplot(
                adata,
                groups=groups or None,
                n_genes=payload.n_genes,
                groupby=str(payload.groupby or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                min_logfoldchange=payload.min_logfoldchange,
                key=str(payload.key or "rank_genes_groups").strip() or "rank_genes_groups",
                dendrogram=bool(payload.dendrogram),
                use_raw=payload.use_raw,
                log=bool(payload.log),
                layer=str(payload.layer or "").strip() or None,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_rank_genes_groups" if "not found" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_correlation_matrix_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_correlation_matrix_plot",
        summary="Plot a correlation matrix between cell groups for the active dataset",
        description=(
            "Render sc.pl.correlation_matrix for the active dataset. "
            "Computes and displays pairwise correlations between groups defined by groupby. "
            "Exposes all scanpy.pl.correlation_matrix parameters including "
            "show_correlation_numbers, dendrogram, cmap, vmin/vmax/vcenter, and figsize."
        ),
    )
    async def generate_correlation_matrix_plot(
        payload: CorrelationMatrixRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_correlation_matrix(
                adata,
                groupby=str(payload.groupby or "").strip() or None,
                show_correlation_numbers=bool(payload.show_correlation_numbers),
                dendrogram=payload.dendrogram,
                cmap=str(payload.cmap or "").strip() or None,
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_rank_genes_groups_matrixplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_matrixplot_plot",
        summary="Plot ranked marker genes as a customizable matrix plot",
        description=(
            "Render sc.pl.rank_genes_groups_matrixplot for the active dataset. "
            "Requires precomputed sc.tl.rank_genes_groups results in adata.uns. "
            "Exposes all scanpy.pl.rank_genes_groups_matrixplot parameters including "
            "n_genes, groups, values_to_plot, standard_scale, dendrogram, swap_axes, "
            "cmap, colorbar_title, vmin/vmax/vcenter, and figsize."
        ),
    )
    async def generate_rank_genes_groups_matrixplot_plot(
        payload: RankGenesGroupsMatrixplotRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None
        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_rank_genes_groups_matrixplot(
                adata,
                groups=groups or None,
                n_genes=payload.n_genes,
                groupby=str(payload.groupby or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                min_logfoldchange=payload.min_logfoldchange,
                key=str(payload.key or "rank_genes_groups").strip() or "rank_genes_groups",
                values_to_plot=str(payload.values_to_plot or "").strip() or None,
                standard_scale=str(payload.standard_scale or "").strip() or None,
                dendrogram=bool(payload.dendrogram),
                swap_axes=bool(payload.swap_axes),
                cmap=str(payload.cmap or "").strip() or None,
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                colorbar_title=str(payload.colorbar_title or "").strip() or None,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_rank_genes_groups" if "not found" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_rank_genes_groups_heatmap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_heatmap_plot",
        summary="Plot ranked marker genes as a customizable heatmap",
        description=(
            "Render sc.pl.rank_genes_groups_heatmap for the active dataset. "
            "Requires precomputed sc.tl.rank_genes_groups results in adata.uns. "
            "Exposes all scanpy.pl.rank_genes_groups_heatmap parameters including "
            "n_genes, groups, standard_scale, swap_axes, show_gene_labels, cmap, vmin/vmax/vcenter, and figsize."
        ),
    )
    async def generate_rank_genes_groups_heatmap_plot(
        payload: RankGenesGroupsHeatmapRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None
        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_rank_genes_groups_heatmap(
                adata,
                groups=groups or None,
                n_genes=payload.n_genes,
                groupby=str(payload.groupby or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                min_logfoldchange=payload.min_logfoldchange,
                key=str(payload.key or "rank_genes_groups").strip() or "rank_genes_groups",
                standard_scale=str(payload.standard_scale or "").strip() or None,
                swap_axes=bool(payload.swap_axes),
                show_gene_labels=payload.show_gene_labels,
                cmap=str(payload.cmap or "").strip() or None,
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_rank_genes_groups" if "not found" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_embedding_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_embedding_plot",
        summary="Generate a plot for any named embedding in the active dataset",
        description=(
            "Generic embedding plotter — the user specifies the basis name (e.g. 'umap', 'tsne', 'pca', 'diffmap', 'draw_graph_fa'). "
            "Use this when the user asks for an embedding that is not UMAP, tSNE, or diffmap, "
            "or when they want to specify the exact basis by name. "
            "If the requested basis does not exist in adata.obsm, the tool returns an error listing available embeddings. "
            "All scanpy.pl.embedding parameters are exposed."
        ),
    )
    async def generate_embedding_plot(payload: EmbeddingPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        color = _parse_string_list(payload.color_json)
        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None

        try:
            result = plotter.run_embedding(
                adata,
                basis=str(payload.basis).strip(),
                color=color or None,
                components=str(payload.components or "").strip() or None,
                use_raw=payload.use_raw,
                layer=str(payload.layer or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                legend_loc=str(payload.legend_loc or "right margin").strip() or "right margin",
                legend_fontsize=payload.legend_fontsize,
                legend_fontweight=str(payload.legend_fontweight or "bold").strip() or "bold",
                colorbar_loc=str(payload.colorbar_loc or "right").strip() or "right",
                color_map=str(payload.color_map or "").strip() or None,
                palette=str(payload.palette or "").strip() or None,
                na_color=str(payload.na_color or "lightgray").strip() or "lightgray",
                na_in_legend=bool(payload.na_in_legend),
                size=payload.size,
                frameon=payload.frameon,
                vmin=str(payload.vmin).strip() if payload.vmin is not None else None,
                vmax=str(payload.vmax).strip() if payload.vmax is not None else None,
                vcenter=payload.vcenter,
                add_outline=bool(payload.add_outline),
                sort_order=bool(payload.sort_order),
                edges=bool(payload.edges),
                edges_width=float(payload.edges_width),
                edges_color=str(payload.edges_color or "grey").strip() or "grey",
                groups=groups,
                projection=str(payload.projection or "2d").strip() or "2d",
                ncols=int(payload.ncols),
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_embedding" if "not found in adata.obsm" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_coloring_label": result.resolved_coloring_label,
            "resolved_genes": result.color_columns,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "basis_used": result.embedding_basis,
            "color_used": result.color_columns,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_diffmap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_diffmap_plot",
        summary="Generate a customizable diffusion map embedding plot for the active dataset",
        description=(
            "Plot a diffusion map (diffusion pseudo-time embedding) scatter plot. "
            "Requires sc.tl.diffmap to have been run on the dataset first (needs precomputed neighbors). "
            "If the embedding is missing, the tool returns a clear error. "
            "Color cells by any observation column or gene expression. "
            "Use the components parameter to select which diffusion components to plot (e.g. '1,2' or '2,3'). "
            "All scanpy.pl.diffmap parameters are exposed."
        ),
    )
    async def generate_diffmap_plot(payload: DiffmapPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        color = _parse_string_list(payload.color_json)
        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None

        try:
            result = plotter.run_diffmap(
                adata,
                color=color or None,
                components=str(payload.components or "").strip() or None,
                use_raw=payload.use_raw,
                layer=str(payload.layer or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                legend_loc=str(payload.legend_loc or "right margin").strip() or "right margin",
                legend_fontsize=payload.legend_fontsize,
                legend_fontweight=str(payload.legend_fontweight or "bold").strip() or "bold",
                colorbar_loc=str(payload.colorbar_loc or "right").strip() or "right",
                color_map=str(payload.color_map or "").strip() or None,
                palette=str(payload.palette or "").strip() or None,
                na_color=str(payload.na_color or "lightgray").strip() or "lightgray",
                na_in_legend=bool(payload.na_in_legend),
                size=payload.size,
                frameon=payload.frameon,
                vmin=str(payload.vmin).strip() if payload.vmin is not None else None,
                vmax=str(payload.vmax).strip() if payload.vmax is not None else None,
                vcenter=payload.vcenter,
                add_outline=bool(payload.add_outline),
                sort_order=bool(payload.sort_order),
                edges=bool(payload.edges),
                edges_width=float(payload.edges_width),
                edges_color=str(payload.edges_color or "grey").strip() or "grey",
                groups=groups,
                ncols=int(payload.ncols),
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            # Catch the missing embedding error specifically for a clean status
            msg = str(exc)
            if "X_diffmap" in msg:
                return {"ok": False, "status": "no_diffmap_embedding", "message": msg, "active_dataset": active}
            return {"ok": False, "status": "plot_error", "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_coloring_label": result.resolved_coloring_label,
            "resolved_genes": result.color_columns,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "color_used": result.color_columns,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_umap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_umap_plot",
        summary="Generate a customizable UMAP embedding plot for the active dataset",
        description=(
            "Plot a UMAP (Uniform Manifold Approximation and Projection) scatter plot. "
            "If no UMAP embedding has been computed, it is calculated automatically (falls back to tSNE if UMAP fails). "
            "Color cells by any observation column (e.g. 'author_cell_type') or gene expression. "
            "Legend is placed on the right margin by default with auto-scaled font size. "
            "All scanpy.pl.umap parameters are exposed."
        ),
    )
    async def generate_umap_plot(payload: UmapPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        color = _parse_string_list(payload.color_json)
        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None

        try:
            result = plotter.run_umap(
                adata,
                color=color or None,
                use_raw=payload.use_raw,
                layer=str(payload.layer or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                legend_loc=str(payload.legend_loc or "right margin").strip() or "right margin",
                legend_fontsize=payload.legend_fontsize,
                legend_fontweight=str(payload.legend_fontweight or "bold").strip() or "bold",
                colorbar_loc=str(payload.colorbar_loc or "right").strip() or "right",
                color_map=str(payload.color_map or "").strip() or None,
                palette=str(payload.palette or "").strip() or None,
                na_color=str(payload.na_color or "lightgray").strip() or "lightgray",
                na_in_legend=bool(payload.na_in_legend),
                size=payload.size,
                frameon=payload.frameon,
                vmin=str(payload.vmin).strip() if payload.vmin is not None else None,
                vmax=str(payload.vmax).strip() if payload.vmax is not None else None,
                vcenter=payload.vcenter,
                add_outline=bool(payload.add_outline),
                sort_order=bool(payload.sort_order),
                edges=bool(payload.edges),
                edges_width=float(payload.edges_width),
                edges_color=str(payload.edges_color or "grey").strip() or "grey",
                groups=groups,
                ncols=int(payload.ncols),
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_coloring_label": result.resolved_coloring_label,
            "resolved_genes": result.color_columns,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "color_used": result.color_columns,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_tsne_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_tsne_plot",
        summary="Generate a customizable tSNE embedding plot for the active dataset",
        description=(
            "Plot a tSNE (t-distributed stochastic neighbor embedding) scatter plot. "
            "If no tSNE embedding has been computed, it is calculated automatically. "
            "Color cells by any observation column (e.g. 'author_cell_type') or gene expression. "
            "Legend is placed on the right margin by default with auto-scaled font size. "
            "All scanpy.pl.tsne parameters are exposed."
        ),
    )
    async def generate_tsne_plot(payload: TsnePlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        color = _parse_string_list(payload.color_json)
        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None

        try:
            result = plotter.run_tsne(
                adata,
                color=color or None,
                use_raw=payload.use_raw,
                layer=str(payload.layer or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                legend_loc=str(payload.legend_loc or "right margin").strip() or "right margin",
                legend_fontsize=payload.legend_fontsize,
                legend_fontweight=str(payload.legend_fontweight or "bold").strip() or "bold",
                colorbar_loc=str(payload.colorbar_loc or "right").strip() or "right",
                color_map=str(payload.color_map or "").strip() or None,
                palette=str(payload.palette or "").strip() or None,
                na_color=str(payload.na_color or "lightgray").strip() or "lightgray",
                na_in_legend=bool(payload.na_in_legend),
                size=payload.size,
                frameon=payload.frameon,
                vmin=str(payload.vmin).strip() if payload.vmin is not None else None,
                vmax=str(payload.vmax).strip() if payload.vmax is not None else None,
                vcenter=payload.vcenter,
                add_outline=bool(payload.add_outline),
                sort_order=bool(payload.sort_order),
                edges=bool(payload.edges),
                edges_width=float(payload.edges_width),
                edges_color=str(payload.edges_color or "grey").strip() or "grey",
                groups=groups,
                ncols=int(payload.ncols),
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_coloring_label": result.resolved_coloring_label,
            "resolved_genes": result.color_columns,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "color_used": result.color_columns,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_dendrogram",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_dendrogram",
        summary="Generate a dendrogram showing hierarchical relationships between cell groups",
        description=(
            "Plot a dendrogram of groupby categories based on their expression profiles. "
            "If sc.tl.dendrogram has not been run yet for the given groupby, it is computed automatically. "
            "groupby defaults to 'author_cell_type'. "
            "Use orientation to control which direction the tree grows ('top', 'bottom', 'left', 'right')."
        ),
    )
    async def generate_dendrogram(payload: DendrogramRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_dendrogram(
                adata,
                groupby=str(payload.groupby or "").strip() or None,
                dendrogram_key=str(payload.dendrogram_key or "").strip() or None,
                orientation=str(payload.orientation or "top").strip() or "top",
                remove_labels=bool(payload.remove_labels),
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_clustermap",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_clustermap",
        summary="Generate a clustermap (hierarchically clustered heatmap) for the active dataset",
        description=(
            "Plot a seaborn clustermap — cells and genes are both hierarchically clustered and reordered. "
            "Unlike heatmap/dotplot, there is no fixed groupby axis; clustering arranges rows and columns automatically. "
            "Providing a gene list via markers_json is strongly recommended — plotting all genes is very slow. "
            "If markers_json is empty, session markers set by set_markers are used. "
            "obs_keys can optionally color-code rows by a categorical observation column."
        ),
    )
    async def generate_clustermap(payload: ClustermapRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        var_names = _parse_string_list(payload.markers_json)
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        if not var_names:
            var_names = plotter.get_markers()
        # Note: unlike other tools, clustermap does NOT error if no genes — it falls back to all genes.
        # But we warn in the response.

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_clustermap(
                adata,
                var_names,
                obs_keys=str(payload.obs_keys or "").strip() or None,
                use_raw=payload.use_raw,
                standard_scale=str(payload.standard_scale or "").strip() or None,
                z_score=payload.z_score,
                method=str(payload.method or "average").strip() or "average",
                metric=str(payload.metric or "euclidean").strip() or "euclidean",
                cmap=str(payload.cmap or "viridis").strip() or "viridis",
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "var_names_used": result.resolved_genes,
            "obs_keys_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if not var_names:
            response["warning"] = "No genes specified — clustered all genes. This may be slow for large datasets."
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_matrixplot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_matrixplot",
        summary="Generate a customizable matrix plot for the active dataset",
        description=(
            "Plot mean gene expression as a color-coded matrix — rows are genes, columns are cell type groups "
            "(or swapped with swap_axes). Each cell shows the mean expression of that gene in that group. "
            "Genes are taken from markers_json; if empty, session markers set by set_markers are used. "
            "If neither is provided, the tool returns an error asking for genes. "
            "groupby defaults to 'author_cell_type'. "
            "All scanpy.pl.matrixplot parameters are exposed."
        ),
    )
    async def generate_matrixplot(payload: MatrixplotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        var_names = _parse_string_list(payload.markers_json)
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        if not var_names:
            var_names = plotter.get_markers()
        if not var_names:
            return {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        categories_order = _parse_string_list(payload.categories_order_json) if payload.categories_order_json else None

        try:
            result = plotter.run_matrixplot(
                adata,
                var_names,
                groupby=str(payload.groupby or "").strip() or None,
                use_raw=payload.use_raw,
                log=bool(payload.log),
                num_categories=int(payload.num_categories),
                categories_order=categories_order,
                dendrogram=bool(payload.dendrogram),
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                layer=str(payload.layer or "").strip() or None,
                standard_scale=str(payload.standard_scale or "").strip() or None,
                swap_axes=bool(payload.swap_axes),
                cmap=str(payload.cmap or "viridis").strip() or "viridis",
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                colorbar_title=str(payload.colorbar_title or "").strip() or None,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "var_names_used": result.resolved_genes,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_stacked_violin",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_stacked_violin",
        summary="Generate a customizable stacked violin plot for the active dataset",
        description=(
            "Plot gene expression as stacked violins — one violin per gene per group — allowing "
            "comparison of expression distributions across many cell types at once. "
            "Genes are taken from markers_json; if empty, session markers set by set_markers are used. "
            "If neither is provided, the tool returns an error asking for genes. "
            "groupby defaults to 'author_cell_type'. "
            "All scanpy.pl.stacked_violin parameters are exposed."
        ),
    )
    async def generate_stacked_violin(payload: StackedViolinRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        var_names = _parse_string_list(payload.markers_json)
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        if not var_names:
            var_names = plotter.get_markers()
        if not var_names:
            return {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        categories_order = _parse_string_list(payload.categories_order_json) if payload.categories_order_json else None

        try:
            result = plotter.run_stacked_violin(
                adata,
                var_names,
                groupby=str(payload.groupby or "").strip() or None,
                use_raw=payload.use_raw,
                log=bool(payload.log),
                num_categories=int(payload.num_categories),
                dendrogram=bool(payload.dendrogram),
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                layer=str(payload.layer or "").strip() or None,
                standard_scale=str(payload.standard_scale or "").strip() or None,
                categories_order=categories_order,
                swap_axes=bool(payload.swap_axes),
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                cmap=str(payload.cmap or "Blues").strip() or "Blues",
                stripplot=bool(payload.stripplot),
                jitter=payload.jitter,
                size=int(payload.size),
                row_palette=str(payload.row_palette or "").strip() or None,
                yticklabels=bool(payload.yticklabels),
                colorbar_title=str(payload.colorbar_title or "").strip() or None,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "var_names_used": result.resolved_genes,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/generate_tracksplot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_tracksplot",
        summary="Generate a customizable tracks plot for the active dataset",
        description=(
            "Plot gene expression as a tracks plot (one track per group) grouped by any observation column. "
            "Each row represents a groupby category; each column is a gene. "
            "Genes are taken from markers_json; if empty, session markers set by set_markers are used. "
            "If neither is provided, the tool returns an error asking for genes. "
            "groupby defaults to 'author_cell_type'. "
            "All scanpy.pl.tracksplot parameters are exposed."
        ),
    )
    async def generate_tracksplot(payload: TracksplotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        var_names = _parse_string_list(payload.markers_json)
        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        if not var_names:
            var_names = plotter.get_markers()
        if not var_names:
            return {
                "ok": False,
                "status": "no_genes",
                "message": (
                    "No genes specified and no session markers are set. "
                    "Please provide markers_json or call set_markers first."
                ),
            }

        figsize: tuple[float, float] | None = None
        if payload.figsize_width is not None and payload.figsize_height is not None:
            figsize = (float(payload.figsize_width), float(payload.figsize_height))

        try:
            result = plotter.run_tracksplot(
                adata,
                var_names,
                groupby=str(payload.groupby or "").strip() or None,
                use_raw=payload.use_raw,
                log=bool(payload.log),
                dendrogram=bool(payload.dendrogram),
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                layer=str(payload.layer or "").strip() or None,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "var_names_used": result.resolved_genes,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    @app.post(
        "/print_adata_obs",
        dependencies=[Depends(_require_api_key)],
        operation_id="print_adata_obs",
        summary="Inspect the active dataset observation table",
        description=(
            "List AnnData observation metadata using the command `print(adata.obs)`. "
            "The tool returns the available `.obs` column names for the active dataset."
        ),
    )
    async def print_adata_obs(payload: ExecuteDatasetCommandRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {
                "ok": False,
                "status": "no_active_dataset",
                "message": str(exc),
            }

        normalized_command = _normalize_dataset_command(payload.command)
        if normalized_command != "print(adata.obs)":
            return {
                "ok": False,
                "status": "unsupported_command",
                "message": "Only print(adata.obs) is supported in this endpoint.",
                "command": payload.command,
                "supported_commands": ["print(adata.obs)"],
            }

        obs_columns = [str(column) for column in adata.obs.columns]
        return {
            "ok": True,
            "status": "success",
            "command": "print(adata.obs)",
            "active_dataset": active,
            "obs_columns": obs_columns,
            "obs_shape": [int(adata.n_obs), len(obs_columns)],
        }

    @app.post(
        "/generate_rank_genes_groups_stacked_violin",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_stacked_violin",
        summary="Plot top-ranked genes per group as a stacked violin from precomputed sc.tl.rank_genes_groups results",
        description=(
            "Plot expression distributions of top-ranked genes per group as stacked violins. "
            "Each row is a gene; each violin shows the expression distribution within a cell type group. "
            "This tool does NOT auto-compute rank_genes_groups — if results are missing it returns "
            "'no_rank_genes_groups' with instructions. "
            "Use var_names_json to override ranked genes with a specific list. "
            "Use swap_axes=True when there are many cell types to avoid label overlap. "
            "All scanpy.pl.rank_genes_groups_stacked_violin parameters are exposed."
        ),
    )
    async def generate_rank_genes_groups_stacked_violin(
        payload: RankGenesGroupsStackedViolinRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        groups = _parse_string_list(payload.groups_json) if payload.groups_json else None
        var_names = _parse_string_list(payload.var_names_json) if payload.var_names_json else None
        figsize = (
            (float(payload.figsize_width), float(payload.figsize_height))
            if payload.figsize_width and payload.figsize_height
            else None
        )

        try:
            result = plotter.run_rank_genes_groups_stacked_violin(
                adata,
                groups=groups,
                n_genes=int(payload.n_genes),
                groupby=None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                var_names=var_names,
                min_logfoldchange=payload.min_logfoldchange,
                key=str(payload.key or "rank_genes_groups").strip() or "rank_genes_groups",
                swap_axes=bool(payload.swap_axes),
                cmap=str(payload.cmap or "Blues").strip() or "Blues",
                stripplot=bool(payload.stripplot),
                jitter=bool(payload.jitter),
                size=int(payload.size),
                row_palette=str(payload.row_palette or "").strip() or None,
                yticklabels=bool(payload.yticklabels),
                standard_scale=str(payload.standard_scale or "").strip() or None,
                vmin=payload.vmin,
                vmax=payload.vmax,
                vcenter=payload.vcenter,
                colorbar_title=str(payload.colorbar_title or "").strip() or None,
                figsize=figsize,
                title=str(payload.title or "").strip() or None,
            )
        except ValueError as exc:
            msg = str(exc)
            status = "no_rank_genes_groups" if "not found" in msg else "plot_error"
            return {"ok": False, "status": status, "message": msg, "active_dataset": active}
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
            "resolved_groupby": result.resolved_groupby,
            "rank_genes_groups_notice": result.rank_genes_groups_notice,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "output_file": output_file,
            "output_markdown": output_markdown,
            "inline_markdown": output_markdown,
            "canonical_response_markdown": canonical_response_markdown,
        }
        if output_url:
            response["output_url"] = output_url
        return response

    return app


app = create_app()
