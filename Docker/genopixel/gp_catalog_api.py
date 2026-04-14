from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd
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


class _ObsFilterMixin(BaseModel):
    obs_filter_json: str = Field(
        default="{}",
        description="Subset cells: JSON {obs_col:[values]}, AND logic, case-insensitive.",
    )


class AnalyzeDatasetRequest(BaseModel):
    h5ad_path: str
    multiple_excel_row: int | None = None


class GeneratePlotRequest(_ObsFilterMixin):
    plot_type: str = Field(
        default="umap",
        description="Plot type e.g. 'umap'.",
    )
    color_json: str | list[str] = Field(
        default="[]",
        description="JSON array of obs columns or genes to color by; [] for none.",
    )
    genes_json: str | list[str] = Field(
        default="[]",
        description="JSON array of genes; [] unless user specifies genes.",
    )
    groupby: str = Field(
        default="",
        description="Groupby column.",
    )
    gene_symbols_column: str = Field(
        default="",
        description="var column with gene symbols; empty for auto.",
    )
    title: str = Field(default="")


class HeatmapPlotRequest(_ObsFilterMixin):
    markers_json: str | list[str] = Field(
        default="[]",
        description="JSON array or comma-separated genes. Empty=session markers.",
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    use_raw: bool | None = Field(default=None)
    log: bool = Field(default=False, description="Log scale.")
    num_categories: int = Field(
        default=7,
        description="Bins when groupby is continuous.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    layer: str = Field(default="", description="Layer; empty=X.")
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    swap_axes: bool = Field(
        default=True,
        description="Swap axes (genes/groups).",
    )
    show_gene_labels: bool | None = Field(
        default=None,
        description="Show gene labels; auto-detected when None.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    title: str = Field(default="")


class LogUnmetRequestRequest(BaseModel):
    user_request: str = Field(
        description="The user's unfulfilled request."
    )
    active_dataset: str = Field(default="")


class SetMarkersRequest(BaseModel):
    markers_json: str | list[str] = Field(
        description="JSON array or comma-separated gene list to store as session defaults.",
    )


class ViolinPlotRequest(_ObsFilterMixin):
    keys_json: str | list[str] = Field(
        description="JSON array or comma-separated genes/obs fields. Required.",
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    rotation: float = Field(default=45.0, include_in_schema=False)
    log: bool = Field(default=False, description="Log scale.")
    use_raw: bool | None = Field(
        default=None,
        description="Use adata.raw; default False.",
    )
    stripplot: bool = Field(default=True, include_in_schema=False)
    jitter: float | bool = Field(default=True, include_in_schema=False)
    size: int = Field(default=1, include_in_schema=False)
    layer: str = Field(
        default="",
        description="Layer; empty=X.",
    )
    density_norm: str = Field(default="width", include_in_schema=False)
    order_json: str = Field(
        default="",
        description="JSON array: x-axis category order.",
    )
    multi_panel: bool | None = Field(default=None, include_in_schema=False)
    xlabel: str = Field(default="", include_in_schema=False)
    ylabel: str = Field(default="", include_in_schema=False)
    title: str = Field(default="")


class CellCountsBarplotRequest(_ObsFilterMixin):
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column (e.g. 'author_cell_type', 'disease').",
    )
    title: str = Field(default="")


class CellTypeProportionBarplotRequest(_ObsFilterMixin):
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    sample_col: str = Field(
        default="",
        description="Sample/donor column; empty=auto-detect.",
    )
    title: str = Field(default="")


class DotplotPlotRequest(_ObsFilterMixin):
    markers_json: str | list[str] = Field(
        default="[]",
        description="JSON array or comma-separated genes. Empty=session markers.",
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    use_raw: bool | None = Field(default=None)
    log: bool = Field(default=False, description="Log scale.")
    num_categories: int = Field(
        default=7,
        description="Bins when groupby is continuous.",
    )
    categories_order_json: str = Field(
        default="",
        description="JSON array: display order of categories.",
    )
    expression_cutoff: float = Field(
        default=0.0,
        description="Min expression cutoff.",
    )
    mean_only_expressed: bool = Field(default=False, include_in_schema=False)
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    layer: str = Field(default="", description="Layer; empty=X.")
    swap_axes: bool = Field(
        default=False,
        description="Swap axes (genes/groups).",
    )
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    cmap: str = Field(default="Reds", description="Colormap name.")
    dot_max: float | None = Field(default=None, include_in_schema=False)
    dot_min: float | None = Field(default=None, include_in_schema=False)
    smallest_dot: float = Field(default=0.0, include_in_schema=False)
    colorbar_title: str = Field(default="", include_in_schema=False)
    size_title: str = Field(default="", include_in_schema=False)
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class RankGenesGroupsViolinRequest(_ObsFilterMixin):
    groups_json: str = Field(
        default="",
        description="JSON array of groups to show; empty=all.",
    )
    n_genes: int = Field(
        default=5,
        description="Top genes per group.",
    )
    gene_names_json: str = Field(
        default="",
        description="JSON gene list override; empty=use ranked.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    use_raw: bool | None = Field(default=None)
    key: str = Field(
        default="rank_genes_groups",
        description="adata.uns key for results.",
    )
    split: bool = Field(default=True, include_in_schema=False)
    density_norm: str = Field(default="width", include_in_schema=False)
    strip: bool = Field(default=True, include_in_schema=False)
    jitter: float | bool = Field(default=True, include_in_schema=False)
    size: int = Field(default=1, include_in_schema=False)
    title: str = Field(default="")


class RankGenesGroupsStackedViolinRequest(_ObsFilterMixin):
    groups_json: str = Field(
        default="",
        description="JSON array of groups to show; empty=all.",
    )
    n_genes: int = Field(
        default=5,
        description="Top genes per group.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    var_names_json: str = Field(
        default="",
        description="JSON gene list override; empty=use ranked.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Min log fold-change filter.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="adata.uns key for results.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes (genes/groups).",
    )
    cmap: str = Field(
        default="Blues",
        description="Colormap name.",
    )
    stripplot: bool = Field(default=False, include_in_schema=False)
    jitter: bool = Field(default=False, include_in_schema=False)
    size: int = Field(default=1, include_in_schema=False)
    row_palette: str = Field(default="", include_in_schema=False)
    yticklabels: bool = Field(default=False, include_in_schema=False)
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    colorbar_title: str = Field(default="", include_in_schema=False)
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class RankGenesGroupsPlotRequest(_ObsFilterMixin):
    groups_json: str = Field(
        default="",
        description="JSON array of groups to show; empty=all.",
    )
    n_genes: int = Field(
        default=5,
        description="Top genes per group.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="adata.uns key for results.",
    )
    fontsize: int = Field(default=11, include_in_schema=False)
    ncols: int = Field(default=3, include_in_schema=False)
    sharey: bool = Field(default=True, include_in_schema=False)
    title: str = Field(default="")


class RankGenesGroupsHeatmapRequest(_ObsFilterMixin):
    groups_json: str = Field(
        default="",
        description="JSON array of groups to show; empty=all.",
    )
    n_genes: int = Field(
        default=5,
        description="Top genes per group; negative for down-regulated.",
    )
    groupby: str = Field(
        default="",
        description="Groupby column; inferred when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Min log fold-change filter.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="adata.uns key for results.",
    )
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes (genes/groups).",
    )
    show_gene_labels: bool | None = Field(
        default=None,
        description="Show gene labels; auto-detected when None.",
    )
    cmap: str = Field(
        default="",
        description="Colormap name; empty=scanpy default.",
    )
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class RankGenesGroupsDotplotRequest(_ObsFilterMixin):
    groups_json: str = Field(
        default="",
        description="JSON array of groups to show; empty=all.",
    )
    n_genes: int = Field(
        default=5,
        description="Top genes per group; negative for down-regulated.",
    )
    groupby: str = Field(
        default="",
        description="Groupby column; inferred when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Min log fold-change filter.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="adata.uns key for results.",
    )
    values_to_plot: str = Field(
        default="",
        description="Color metric: 'scores','logfoldchanges','pvals','pvals_adj'. Empty=mean expr.",
    )
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes (genes/groups).",
    )
    cmap: str = Field(
        default="",
        description="Colormap name; empty=scanpy default.",
    )
    dot_max: float | None = Field(default=None, include_in_schema=False)
    dot_min: float | None = Field(default=None, include_in_schema=False)
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class RankGenesGroupsMatrixplotRequest(_ObsFilterMixin):
    groups_json: str = Field(
        default="",
        description="JSON array of groups to show; empty=all.",
    )
    n_genes: int = Field(
        default=5,
        description="Top genes per group; negative for down-regulated.",
    )
    groupby: str = Field(
        default="",
        description="Groupby column; inferred when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Min log fold-change filter.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="adata.uns key for results.",
    )
    values_to_plot: str = Field(
        default="",
        description="Color metric: 'scores','logfoldchanges','pvals','pvals_adj'. Empty=mean expr.",
    )
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes (genes/groups).",
    )
    cmap: str = Field(
        default="",
        description="Colormap name; empty=scanpy default.",
    )
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    colorbar_title: str = Field(default="", include_in_schema=False)
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class RankGenesGroupsTracksplotRequest(_ObsFilterMixin):
    groups_json: str = Field(
        default="",
        description="JSON array of groups to show; empty=all.",
    )
    n_genes: int = Field(
        default=5,
        description="Top genes per group; negative for down-regulated.",
    )
    groupby: str = Field(
        default="",
        description="Groupby column; inferred when empty.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    min_logfoldchange: float | None = Field(
        default=None,
        description="Min log fold-change filter.",
    )
    key: str = Field(
        default="rank_genes_groups",
        description="adata.uns key for results.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Use adata.raw; default False.",
    )
    log: bool = Field(default=False)
    layer: str = Field(
        default="",
        description="Layer; empty=X.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class CorrelationMatrixRequest(_ObsFilterMixin):
    groupby: str = Field(
        default="",
        description="Groupby column; inferred when empty.",
    )
    show_correlation_numbers: bool = Field(default=False, include_in_schema=False)
    dendrogram: bool | None = Field(
        default=None,
        description="Add dendrogram; auto when None.",
    )
    cmap: str = Field(
        default="",
        description="Colormap name; empty=scanpy default.",
    )
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class EmbeddingPlotRequest(_ObsFilterMixin):
    basis: str = Field(
        description="Embedding name e.g. 'umap','pca','draw_graph_fa' (X_ prefix optional).",
    )
    color_json: str | list[str] = Field(
        default="[]",
        description="JSON array of obs columns or genes to color by; empty=auto.",
    )
    components: str = Field(
        default="",
        description="Components e.g. '1,2'; empty=default.",
    )
    use_raw: bool | None = Field(default=None)
    layer: str = Field(default="", description="Layer; empty=X.")
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position.",
    )
    legend_fontsize: float | None = Field(default=None, include_in_schema=False)
    legend_fontweight: str = Field(default="bold", include_in_schema=False)
    colorbar_loc: str = Field(default="right", include_in_schema=False)
    color_map: str = Field(
        default="",
        description="Colormap; empty=default.",
    )
    palette: str = Field(
        default="",
        description="Palette; empty=default.",
    )
    na_color: str = Field(default="lightgray", include_in_schema=False)
    na_in_legend: bool = Field(default=True, include_in_schema=False)
    size: float | None = Field(default=None, include_in_schema=False)
    frameon: bool | None = Field(default=None, include_in_schema=False)
    vmin: str | None = Field(
        default=None,
        description="Lower color limit; percentile syntax e.g. 'p1.5'.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color limit; percentile syntax e.g. 'p98'.",
    )
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    add_outline: bool = Field(default=False, include_in_schema=False)
    sort_order: bool = Field(default=True, include_in_schema=False)
    edges: bool = Field(default=False, include_in_schema=False)
    edges_width: float = Field(default=0.1, include_in_schema=False)
    edges_color: str = Field(default="grey", include_in_schema=False)
    groups_json: str = Field(
        default="",
        description="JSON categories to highlight; others greyed out.",
    )
    projection: str = Field(default="2d", include_in_schema=False)
    ncols: int = Field(default=4, include_in_schema=False)
    title: str = Field(default="")


class DiffmapPlotRequest(_ObsFilterMixin):
    color_json: str | list[str] = Field(
        default="[]",
        description="JSON array of obs columns or genes to color by; empty=auto.",
    )
    components: str = Field(
        default="",
        description="Components e.g. '1,2'; empty=default.",
    )
    use_raw: bool | None = Field(default=None)
    layer: str = Field(default="", description="Layer; empty=X.")
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position.",
    )
    legend_fontsize: float | None = Field(default=None, include_in_schema=False)
    legend_fontweight: str = Field(default="bold", include_in_schema=False)
    colorbar_loc: str = Field(default="right", include_in_schema=False)
    color_map: str = Field(
        default="",
        description="Colormap; empty=default.",
    )
    palette: str = Field(
        default="",
        description="Palette; empty=default.",
    )
    na_color: str = Field(default="lightgray", include_in_schema=False)
    na_in_legend: bool = Field(default=True, include_in_schema=False)
    size: float | None = Field(default=None, include_in_schema=False)
    frameon: bool | None = Field(default=None, include_in_schema=False)
    vmin: str | None = Field(
        default=None,
        description="Lower color limit; percentile syntax e.g. 'p1.5'.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color limit; percentile syntax e.g. 'p98'.",
    )
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    add_outline: bool = Field(default=False, include_in_schema=False)
    sort_order: bool = Field(default=True, include_in_schema=False)
    edges: bool = Field(default=False, include_in_schema=False)
    edges_width: float = Field(default=0.1, include_in_schema=False)
    edges_color: str = Field(default="grey", include_in_schema=False)
    groups_json: str = Field(
        default="",
        description="JSON categories to highlight; others greyed out.",
    )
    ncols: int = Field(default=4, include_in_schema=False)
    title: str = Field(default="")


class UmapPlotRequest(_ObsFilterMixin):
    color_json: str | list[str] = Field(
        default="[]",
        description="JSON array of obs columns or genes to color by; empty=auto.",
    )
    use_raw: bool | None = Field(default=None)
    layer: str = Field(default="", description="Layer; empty=X.")
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position.",
    )
    legend_fontsize: float | None = Field(default=None, include_in_schema=False)
    legend_fontweight: str = Field(default="bold", include_in_schema=False)
    colorbar_loc: str = Field(default="right", include_in_schema=False)
    color_map: str = Field(
        default="",
        description="Colormap; empty=default.",
    )
    palette: str = Field(
        default="",
        description="Palette; empty=default.",
    )
    na_color: str = Field(default="lightgray", include_in_schema=False)
    na_in_legend: bool = Field(default=True, include_in_schema=False)
    size: float | None = Field(default=None, include_in_schema=False)
    frameon: bool | None = Field(default=None, include_in_schema=False)
    vmin: str | None = Field(
        default=None,
        description="Lower color limit; percentile syntax e.g. 'p1.5'.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color limit; percentile syntax e.g. 'p98'.",
    )
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    add_outline: bool = Field(default=False, include_in_schema=False)
    sort_order: bool = Field(default=True, include_in_schema=False)
    edges: bool = Field(default=False, include_in_schema=False)
    edges_width: float = Field(default=0.1, include_in_schema=False)
    edges_color: str = Field(default="grey", include_in_schema=False)
    groups_json: str = Field(
        default="",
        description="JSON categories to highlight; others greyed out.",
    )
    ncols: int = Field(default=4, include_in_schema=False)
    title: str = Field(default="")


class TsnePlotRequest(_ObsFilterMixin):
    color_json: str | list[str] = Field(
        default="[]",
        description="JSON array of obs columns or genes to color by; empty=auto.",
    )
    use_raw: bool | None = Field(default=None)
    layer: str = Field(default="", description="Layer; empty=X.")
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    legend_loc: str = Field(
        default="right margin",
        description="Legend position.",
    )
    legend_fontsize: float | None = Field(default=None, include_in_schema=False)
    legend_fontweight: str = Field(default="bold", include_in_schema=False)
    colorbar_loc: str = Field(default="right", include_in_schema=False)
    color_map: str = Field(
        default="",
        description="Colormap; empty=default.",
    )
    palette: str = Field(
        default="",
        description="Palette; empty=default.",
    )
    na_color: str = Field(default="lightgray", include_in_schema=False)
    na_in_legend: bool = Field(default=True, include_in_schema=False)
    size: float | None = Field(default=None, include_in_schema=False)
    frameon: bool | None = Field(default=None, include_in_schema=False)
    vmin: str | None = Field(
        default=None,
        description="Lower color limit; percentile syntax e.g. 'p1.5'.",
    )
    vmax: str | None = Field(
        default=None,
        description="Upper color limit; percentile syntax e.g. 'p98'.",
    )
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    add_outline: bool = Field(default=False, include_in_schema=False)
    sort_order: bool = Field(default=True, include_in_schema=False)
    edges: bool = Field(default=False, include_in_schema=False)
    edges_width: float = Field(default=0.1, include_in_schema=False)
    edges_color: str = Field(default="grey", include_in_schema=False)
    groups_json: str = Field(
        default="",
        description="JSON categories to highlight; others greyed out.",
    )
    ncols: int = Field(default=4, include_in_schema=False)
    title: str = Field(default="")


class DendrogramRequest(_ObsFilterMixin):
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    dendrogram_key: str = Field(
        default="",
        description="adata.uns dendrogram key; empty=default.",
    )
    orientation: str = Field(default="top", include_in_schema=False)
    remove_labels: bool = Field(default=False, include_in_schema=False)
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class ClustermapRequest(_ObsFilterMixin):
    markers_json: str | list[str] = Field(
        default="[]",
        description="JSON array or comma-separated genes. Empty=session markers.",
    )
    obs_keys: str = Field(
        default="",
        description="obs column for row colors; empty=none.",
    )
    use_raw: bool | None = Field(default=None)
    standard_scale: str = Field(
        default="",
        description="Normalize: 'row' or 'col'. Leave empty to skip.",
    )
    z_score: int | None = Field(default=None, include_in_schema=False)
    method: str = Field(default="average", include_in_schema=False)
    metric: str = Field(default="euclidean", include_in_schema=False)
    cmap: str = Field(default="viridis", description="Colormap name.")
    figsize_width: float | None = Field(default=None, description="Width inches.")
    figsize_height: float | None = Field(default=None, description="Height inches.")
    title: str = Field(default="")


class MatrixplotRequest(_ObsFilterMixin):
    markers_json: str | list[str] = Field(
        default="[]",
        description="JSON array or comma-separated genes. Empty=session markers.",
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    use_raw: bool | None = Field(default=None)
    log: bool = Field(default=False, description="Log scale.")
    num_categories: int = Field(
        default=7,
        description="Bins when groupby is continuous.",
    )
    categories_order_json: str = Field(
        default="",
        description="JSON array: display order of categories.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    layer: str = Field(default="", description="Layer; empty=X.")
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes (genes/groups).",
    )
    cmap: str = Field(default="viridis", description="Colormap name.")
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    colorbar_title: str = Field(default="", include_in_schema=False)
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class StackedViolinRequest(_ObsFilterMixin):
    markers_json: str | list[str] = Field(
        default="[]",
        description="JSON array or comma-separated genes. Empty=session markers.",
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    use_raw: bool | None = Field(default=None)
    log: bool = Field(default=False, description="Log scale.")
    num_categories: int = Field(
        default=7,
        description="Bins when groupby is continuous.",
    )
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    layer: str = Field(default="", description="Layer; empty=X.")
    standard_scale: str = Field(
        default="",
        description="Normalize: 'var' or 'obs'. Leave empty to skip.",
    )
    categories_order_json: str = Field(
        default="",
        description="JSON array: display order of categories.",
    )
    swap_axes: bool = Field(
        default=False,
        description="Swap axes (genes/groups).",
    )
    vmin: float | None = Field(default=None)
    vmax: float | None = Field(default=None)
    vcenter: float | None = Field(default=None, description="Diverging colormap center.")
    cmap: str = Field(default="Blues", description="Colormap name.")
    stripplot: bool = Field(default=False, include_in_schema=False)
    jitter: float | bool = Field(default=False, include_in_schema=False)
    size: int = Field(default=1, include_in_schema=False)
    row_palette: str = Field(default="", include_in_schema=False)
    yticklabels: bool = Field(default=False, include_in_schema=False)
    colorbar_title: str = Field(default="", include_in_schema=False)
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class TracksplotRequest(_ObsFilterMixin):
    markers_json: str | list[str] = Field(
        default="[]",
        description="JSON array or comma-separated genes. Empty=session markers.",
    )
    groupby: str = Field(
        default="author_cell_type",
        description="Groupby column.",
    )
    use_raw: bool | None = Field(default=None)
    log: bool = Field(default=False, description="Log scale.")
    dendrogram: bool = Field(
        default=False,
        description="Add hierarchical dendrogram.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    layer: str = Field(default="", description="Layer; empty=X.")
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )
    title: str = Field(default="")


class ScatterPlotRequest(_ObsFilterMixin):
    x: str = Field(
        default="",
        description="x-axis: obs column, gene, or empty if basis set.",
    )
    y: str = Field(
        default="",
        description="y-axis: obs column, gene, or empty if basis set.",
    )
    color_json: str = Field(
        default="[]",
        description="JSON array of obs columns or genes to color by.",
    )
    basis: str = Field(
        default="",
        description="Embedding name; when set, x/y from embedding coords.",
    )
    use_raw: bool | None = Field(
        default=None,
        description="Use adata.raw; default False.",
    )
    layer: str = Field(
        default="",
        description="Layer; empty=X.",
    )
    groups_json: str = Field(
        default="[]",
        description="JSON categories to restrict coloring to.",
    )
    components: str = Field(
        default="",
        description="Components e.g. '1,2'; only if basis is set.",
    )
    sort_order: bool = Field(default=True, include_in_schema=False)
    legend_loc: str = Field(
        default="right margin",
        description="Legend position.",
    )
    size: float | None = Field(default=None, include_in_schema=False)
    color_map: str = Field(
        default="",
        description="Colormap; empty=default.",
    )
    frameon: bool | None = Field(default=None, include_in_schema=False)
    title: str = Field(default="")
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed.",
    )


class HighestExprGenesRequest(_ObsFilterMixin):
    n_top: int = Field(
        default=30,
        description="Number of top genes.",
    )
    layer: str = Field(
        default="",
        description="Layer; empty=X.",
    )
    gene_symbols: str = Field(
        default="",
        description="var column with gene symbols; leave empty for auto.",
    )
    log: bool = Field(default=False)
    title: str = Field(default="")


class ExecuteDatasetCommandRequest(BaseModel):
    command: str = Field(
        default="print(adata.obs)",
        description="Command; supported: 'print(adata.obs)'.",
    )


class ObsCountTableRequest(_ObsFilterMixin):
    row_col: str = Field(
        description="obs column to use as rows (e.g. 'author_cell_type').",
    )
    col_col: str = Field(
        description="obs column to use as columns / second groupby key (e.g. 'sample').",
    )


class SpatialScatterRequest(_ObsFilterMixin):
    color_json: str | list[str] = Field(
        default="[]",
        description="JSON array of obs columns or gene names/ENS IDs to color by; [] for auto cell type.",
    )
    figsize_width: float | None = Field(default=None, description="Width inches; auto-computed.")
    figsize_height: float | None = Field(default=None, description="Height inches; auto-computed.")
    title: str = Field(default="")
    # Layout
    wspace: float | None = Field(default=None, description="Width space between panels (e.g. 0.1).")
    hspace: float | None = Field(default=None, description="Height space between panels.")
    ncols: int | None = Field(default=None, description="Number of panels per row.")
    # Shape / size
    shape: str | None = Field(default=None, description="Point shape: 'circle', 'square', or 'hex'.")
    size: float | None = Field(default=None, description="Size of the scatter point/shape.")
    alpha: float | None = Field(default=None, description="Alpha (opacity) for scatter points.")
    # Color / style
    cmap: str | None = Field(default=None, description="Colormap for continuous annotations (e.g. 'viridis').")
    palette: str | None = Field(default=None, description="Palette for discrete annotations.")
    na_color: str | None = Field(default=None, description="Color for NA values.")
    # Image
    img: bool | None = Field(default=None, description="Whether to plot the tissue image overlay.")
    img_alpha: float | None = Field(default=None, description="Alpha for the underlying tissue image.")
    # Crop
    crop_coord: list[int] | None = Field(default=None, description="[left, right, top, bottom] crop in pixel space.")
    # Legend / colorbar
    legend_loc: str | None = Field(default=None, description="Legend location (e.g. 'right margin', 'on data').")
    legend_fontsize: float | None = Field(default=None, description="Legend font size.")
    colorbar: bool | None = Field(default=None, description="Whether to show the colorbar.")
    frameon: bool | None = Field(default=None, description="Whether to draw a frame around panels.")
    # Outline
    outline: bool | None = Field(default=None, description="Whether to draw a thin border around points.")
    # Groups
    groups_json: str = Field(default="[]", description="JSON array of group values to show (others shown as NA).")
    # Layer
    layer: str | None = Field(default=None, description="adata.layers key to use instead of X.")
    use_raw: bool | None = Field(default=None, description="Whether to use adata.raw.")


class ObsUniqueValuesRequest(BaseModel):
    column: str = Field(
        default="",
        description="adata.obs column name to get unique values for. Leave empty to list all categorical columns.",
    )


class NhoodEnrichmentRequest(_ObsFilterMixin):
    cluster_key: str = Field(
        default="",
        description="adata.obs column for cluster labels; empty=auto cell type column.",
    )
    mode: str = Field(
        default="zscore",
        description="Enrichment mode: 'zscore' or 'count'.",
    )
    figsize_width: float | None = Field(
        default=None,
        description="Width inches; auto-computed from number of clusters.",
    )
    figsize_height: float | None = Field(
        default=None,
        description="Height inches; auto-computed from number of clusters.",
    )
    title: str = Field(default="")


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


def _apply_obs_filter(adata: ad.AnnData, obs_filter_json: str) -> ad.AnnData:
    """Return a view of adata restricted to cells matching obs_filter for a single plot.

    Returns the original adata unchanged when obs_filter_json is empty/null.
    Raises ValueError with a descriptive message when a column or value is not found.
    Column names and values are matched case-insensitively.
    The view is temporary — the active dataset in RUNTIME_STATE is never modified.
    """
    raw = str(obs_filter_json or "").strip()
    if not raw or raw in ("{}", "null", "none", ""):
        return adata
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        raise ValueError(f"obs_filter_json is not valid JSON: {raw!r}")
    if not isinstance(parsed, dict):
        raise ValueError(
            "obs_filter_json must be a JSON object, e.g. "
            '{\"author_cell_type\": [\"B cells\", \"T cells\"]}'
        )
    if not parsed:
        return adata

    lower_to_col = {str(c).lower(): str(c) for c in adata.obs.columns}
    mask = pd.Series([True] * adata.n_obs, index=adata.obs.index)

    for requested_col, values in parsed.items():
        col = lower_to_col.get(str(requested_col).strip().lower())
        if col is None:
            available = ", ".join(str(c) for c in list(adata.obs.columns)[:30])
            raise ValueError(
                f"obs_filter column '{requested_col}' not found in adata.obs. "
                f"Available columns: {available}"
            )
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, list):
            raise ValueError(
                f"obs_filter values for '{requested_col}' must be a list of strings, "
                f"e.g. [\"B cells\", \"T cells\"]."
            )
        str_values = {str(v).strip().lower() for v in values if str(v).strip()}
        series_lower = adata.obs[col].astype(str).str.strip().str.lower()
        col_mask = series_lower.isin(str_values)
        if not col_mask.any():
            available_vals = adata.obs[col].astype(str).unique().tolist()[:20]
            raise ValueError(
                f"obs_filter: no cells match '{col}' in {list(values)}. "
                f"Available values (up to 20): {available_vals}"
            )
        mask = mask & col_mask

    if not mask.any():
        raise ValueError(
            "obs_filter resulted in 0 cells after applying all filters. "
            "Check that your filter combinations are not mutually exclusive."
        )

    return adata[mask.values]


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
        summary="Generate a Scanpy plot",
        description="Plot the loaded dataset. For plain UMAP: plot_type='umap', color_json='[]', genes_json='[]'.",
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
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "resolved_genes": result.resolved_genes,
            "resolved_groupby": result.resolved_groupby,
            "resolved_coloring_label": result.resolved_coloring_label,
            "rank_genes_groups_computed": result.rank_genes_groups_computed,
            "rank_genes_groups_notice": result.rank_genes_groups_notice,
        }
        if output_url:
            plot_payload["output_url"] = output_url

        canonical_response_markdown = build_canonical_response_markdown(active, plot_payload, output_markdown)

        return {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "plot": plot_payload,
            "inline_markdown": canonical_response_markdown or output_markdown,
        }

    @app.post(
        "/generate_heatmap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_heatmap_plot",
        summary="Heatmap plot",
        description="Heatmap grouped by obs column. Genes from markers_json or session markers.",
    )
    async def generate_heatmap_plot(payload: HeatmapPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/set_markers",
        dependencies=[Depends(_require_api_key)],
        operation_id="set_markers",
        summary="Set a session-level gene marker list",
        description="Store genes as session default markers for all gene-based tools.",
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
        "/log_unmet_request",
        dependencies=[Depends(_require_api_key)],
        operation_id="log_unmet_request",
        summary="Log an unmet user request",
        description=(
            "Call this when the user asks for an analysis or plot that no available GenoPixel tool can fulfill. "
            "Appends the request to a log file for developer review and future feature planning."
        ),
    )
    async def log_unmet_request(payload: LogUnmetRequestRequest) -> dict[str, Any]:
        import datetime

        record = {
            "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
            "user_request": str(payload.user_request or "").strip(),
            "active_dataset": str(payload.active_dataset or "").strip(),
        }
        try:
            log_path = Path(_settings().output_dir).parent / "unmet_requests.jsonl"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(record, ensure_ascii=False) + "\n")
            return {
                "ok": True,
                "status": "logged",
                "message": "Your request has been noted and sent to the developers for consideration.",
            }
        except Exception as exc:
            return {
                "ok": False,
                "status": "log_error",
                "message": str(exc),
            }

    @app.post(
        "/generate_violin_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_violin_plot",
        summary="Violin plot",
        description="Violin plot grouped by obs column.",
    )
    async def generate_violin_plot(payload: ViolinPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/cell_counts_barplot",
        dependencies=[Depends(_require_api_key)],
        operation_id="cell_counts_barplot",
        summary="Cell counts bar plot",
        description="Bar plot of cell counts by obs column, sorted descending.",
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
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/cell_type_proportion_barplot",
        dependencies=[Depends(_require_api_key)],
        operation_id="cell_type_proportion_barplot",
        summary="Cell type proportion bar plot",
        description="Stacked bar plot of cell type proportions per sample.",
    )
    async def cell_type_proportion_barplot(payload: CellTypeProportionBarplotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {
                "ok": False,
                "status": "no_active_dataset",
                "message": str(exc),
            }
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

        from gp_models import PlotRequest as _PlotRequest
        sample_col = str(payload.sample_col or "").strip()
        plot_request = _PlotRequest(
            plot_type="cell_type_proportion_barplot",
            groupby=str(payload.groupby or "").strip() or None,
            color=[sample_col] if sample_col else [],
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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_dotplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_dotplot_plot",
        summary="Dot plot",
        description="Dot plot: dot size=fraction expressing, color=mean expression. Genes from markers_json or session markers.",
    )
    async def generate_dotplot_plot(payload: DotplotPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/check_rank_genes_groups",
        dependencies=[Depends(_require_api_key)],
        operation_id="check_rank_genes_groups",
        summary="Rank genes groups plot",
        description="Check if rank_genes_groups results exist and plot them if so.",
    )
    async def check_rank_genes_groups() -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        rgg = adata.uns.get("rank_genes_groups")
        if not isinstance(rgg, dict) or "names" not in rgg:
            return {
                "ok": False,
                "status": "not_available",
                "message": "rank_genes_groups results are not available in the active dataset.",
                "active_dataset": active,
            }

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        try:
            result = plotter.run_rank_genes_groups(adata, n_genes=5)
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
            "status": "available",
            "active_dataset": active,
            "groupby_used": result.resolved_groupby,
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_rank_genes_groups_violin",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_violin",
        summary="Rank genes violin",
        description="Violins of top-ranked genes per group. Requires precomputed rank_genes_groups.",
    )
    async def generate_rank_genes_groups_violin(payload: RankGenesGroupsViolinRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_rank_genes_groups_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_plot",
        summary="Rank genes groups",
        description="Score panels per group for top-ranked genes. Requires precomputed rank_genes_groups.",
    )
    async def generate_rank_genes_groups_plot(payload: RankGenesGroupsPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_rank_genes_groups_dotplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_dotplot_plot",
        summary="Rank genes dot plot",
        description="Dot plot of top-ranked genes. Requires precomputed rank_genes_groups.",
    )
    async def generate_rank_genes_groups_dotplot_plot(
        payload: RankGenesGroupsDotplotRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_rank_genes_groups_tracksplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_tracksplot_plot",
        summary="Rank genes tracks plot",
        description="Tracks plot of top-ranked genes. Requires precomputed rank_genes_groups.",
    )
    async def generate_rank_genes_groups_tracksplot_plot(
        payload: RankGenesGroupsTracksplotRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_correlation_matrix_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_correlation_matrix_plot",
        summary="Correlation matrix",
        description="Pairwise correlation matrix between groupby categories.",
    )
    async def generate_correlation_matrix_plot(
        payload: CorrelationMatrixRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_rank_genes_groups_matrixplot_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_matrixplot_plot",
        summary="Rank genes matrix plot",
        description="Matrix plot of top-ranked genes. Requires precomputed rank_genes_groups.",
    )
    async def generate_rank_genes_groups_matrixplot_plot(
        payload: RankGenesGroupsMatrixplotRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_rank_genes_groups_heatmap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_heatmap_plot",
        summary="Rank genes heatmap",
        description="Heatmap of top-ranked genes. Requires precomputed rank_genes_groups.",
    )
    async def generate_rank_genes_groups_heatmap_plot(
        payload: RankGenesGroupsHeatmapRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_embedding_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_embedding_plot",
        summary="Embedding plot",
        description="Generic embedding scatter plot; specify basis name (e.g. 'pca', 'draw_graph_fa').",
    )
    async def generate_embedding_plot(payload: EmbeddingPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_diffmap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_diffmap_plot",
        summary="Diffusion map plot",
        description="Diffusion map scatter plot. Requires precomputed sc.tl.diffmap.",
    )
    async def generate_diffmap_plot(payload: DiffmapPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_umap_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_umap_plot",
        summary="UMAP plot",
        description="UMAP scatter plot; auto-computes if missing.",
    )
    async def generate_umap_plot(payload: UmapPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_tsne_plot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_tsne_plot",
        summary="tSNE plot",
        description="tSNE scatter plot; auto-computes if missing.",
    )
    async def generate_tsne_plot(payload: TsnePlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_dendrogram",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_dendrogram",
        summary="Dendrogram",
        description="Dendrogram of groupby categories; auto-computes if missing.",
    )
    async def generate_dendrogram(payload: DendrogramRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_clustermap",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_clustermap",
        include_in_schema=False,
    )
    async def generate_clustermap(payload: ClustermapRequest) -> dict[str, Any]:
        return {"ok": False, "status": "disabled", "message": "Clustermap plotting is temporarily disabled."}

    @app.post(
        "/generate_matrixplot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_matrixplot",
        summary="Matrix plot",
        description="Matrix of mean gene expression per group. Genes from markers_json or session markers.",
    )
    async def generate_matrixplot(payload: MatrixplotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_stacked_violin",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_stacked_violin",
        summary="Stacked violin plot",
        description="Stacked violins per gene per group. Genes from markers_json or session markers.",
    )
    async def generate_stacked_violin(payload: StackedViolinRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_tracksplot",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_tracksplot",
        summary="Tracks plot",
        description="Tracks plot: each row=group, each col=gene. Genes from markers_json or session markers.",
    )
    async def generate_tracksplot(payload: TracksplotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/print_adata",
        dependencies=[Depends(_require_api_key)],
        operation_id="print_adata",
        summary="Dataset summary",
        description="Return AnnData summary: dims, obs/var columns, embeddings, uns keys, layers.",
    )
    async def print_adata() -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        obs_cols = [str(c) for c in adata.obs.columns]
        var_cols = [str(c) for c in adata.var.columns]
        obsm_keys = [str(k) for k in adata.obsm.keys()]
        obsp_keys = [str(k) for k in adata.obsp.keys()]
        uns_keys = [str(k) for k in adata.uns.keys()]
        layers_keys = [str(k) for k in adata.layers.keys()]

        summary_text = str(adata)

        # Detect best annotation column for contextual suggestions
        _CELL_TYPE_KEYWORDS = (
            "cell_type", "celltype", "cell_label", "celllabel",
            "annotation", "annot", "cluster", "leiden", "louvain",
            "subtype", "lineage", "phenotype", "identity", "ident",
            "label", "class", "state", "population",
        )
        annotation_col: str | None = None
        for col in obs_cols:
            lc = col.lower().replace(" ", "_").replace("-", "_")
            if any(kw in lc for kw in _CELL_TYPE_KEYWORDS):
                annotation_col = col
                break

        ann = annotation_col or (obs_cols[0] if obs_cols else "cell_type")
        has_rgg = any(
            isinstance(adata.uns.get(k), dict) and "names" in adata.uns[k]
            for k in adata.uns
        )
        suggested_next_steps = [
            {"label": "UMAP", "prompt": "show me a UMAP"},
            {"label": "Cell counts", "prompt": f"plot cell counts by {ann}"},
            {"label": "Violin plot", "prompt": f"use violin plot to show the expression of any gene across {ann}"},
        ]
        if has_rgg:
            suggested_next_steps.append({"label": "Markers/top genes", "prompt": "show top marker genes"})

        return {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "summary": summary_text,
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

    @app.post(
        "/print_adata_obs",
        dependencies=[Depends(_require_api_key)],
        operation_id="print_adata_obs",
        summary="Inspect obs table",
        description="List adata.obs column names.",
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

        _CELL_TYPE_KEYWORDS = (
            "cell_type", "celltype", "cell_label", "celllabel",
            "annotation", "annot", "cluster", "leiden", "louvain",
            "subtype", "lineage", "phenotype", "identity", "ident",
            "label", "class", "state", "population",
        )

        def _is_cell_type_col(col: str) -> bool:
            lc = col.lower().replace(" ", "_").replace("-", "_")
            return any(kw in lc for kw in _CELL_TYPE_KEYWORDS)

        categorical_cols: list[dict] = []
        numerical_cols: list[str] = []
        other_cols: list[str] = []

        for col in adata.obs.columns:
            col_str = str(col)
            series = adata.obs[col]
            if hasattr(series, "cat") or str(series.dtype) == "category" or series.dtype == object:
                try:
                    unique_vals = sorted(series.dropna().unique().tolist(), key=lambda v: str(v))
                    unique_vals_str = [str(v) for v in unique_vals]
                except Exception:
                    unique_vals_str = []
                _MAX_VALUES = 50
                categorical_cols.append({
                    "column": col_str,
                    "n_unique": len(unique_vals_str),
                    "values": unique_vals_str[:_MAX_VALUES],
                    "values_truncated": len(unique_vals_str) > _MAX_VALUES,
                    "is_annotation": _is_cell_type_col(col_str),
                })
            elif hasattr(series, "dtype") and str(series.dtype).startswith(("int", "float", "uint")):
                numerical_cols.append(col_str)
            else:
                other_cols.append(col_str)

        # Sort categorical: annotation/cell-type columns first, then the rest alphabetically
        annotation_cols = [c for c in categorical_cols if c["is_annotation"]]
        non_annotation_cols = [c for c in categorical_cols if not c["is_annotation"]]
        annotation_cols.sort(key=lambda c: c["column"])
        non_annotation_cols.sort(key=lambda c: c["column"])
        ordered_categorical = annotation_cols + non_annotation_cols

        # Strip internal flag before returning
        for c in ordered_categorical:
            del c["is_annotation"]

        # Build a human-readable markdown summary for inline display
        n_cells, n_cols = int(adata.n_obs), int(len(adata.obs.columns))
        md_lines: list[str] = [
            f"**Dataset OBS — {n_cells:,} cells × {n_cols} columns**",
            "",
        ]

        if ordered_categorical:
            md_lines.append("**Categorical columns:**")
            for c in ordered_categorical:
                vals = c["values"]
                preview = ", ".join(str(v) for v in vals[:8])
                if len(vals) > 8:
                    preview += f", … (+{len(vals) - 8} more)"
                md_lines.append(f"- **{c['column']}** ({c['n_unique']} unique): {preview}")
            md_lines.append("")

        if numerical_cols:
            md_lines.append("**Numerical columns:**")
            md_lines.append("- " + ", ".join(numerical_cols))
            md_lines.append("")

        if other_cols:
            md_lines.append("**Other columns:**")
            md_lines.append("- " + ", ".join(other_cols))

        inline_markdown = "\n".join(md_lines).strip()

        return {
            "ok": True,
            "status": "success",
            "command": "print(adata.obs)",
            "active_dataset": active,
            "obs_shape": [n_cells, n_cols],
            "categorical_columns": ordered_categorical,
            "numerical_columns": numerical_cols,
            "other_columns": other_cols,
            "inline_markdown": inline_markdown,
        }

    @app.post(
        "/get_obs_unique_values",
        dependencies=[Depends(_require_api_key)],
        operation_id="get_obs_unique_values",
        summary="Unique values for an obs column",
        description=(
            "Return ALL unique values for a given adata.obs column, sorted alphabetically. "
            "Use this whenever the user asks for 'list of unique values', 'what values does X have', "
            "'full list of diseases', 'all cell types', etc. "
            "Pass column='' to get all categorical columns with their full value lists."
        ),
    )
    async def get_obs_unique_values(payload: ObsUniqueValuesRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        col = str(payload.column).strip()

        def _col_unique(series: "pd.Series") -> list[str]:
            try:
                vals = sorted(series.dropna().unique().tolist(), key=lambda v: str(v))
                return [str(v) for v in vals]
            except Exception:
                return [str(v) for v in series.dropna().unique().tolist()]

        # Single column
        if col:
            # Case-insensitive match
            col_map = {str(c).lower(): str(c) for c in adata.obs.columns}
            resolved = col if col in adata.obs.columns else col_map.get(col.lower())
            if not resolved:
                return {
                    "ok": False,
                    "status": "column_not_found",
                    "message": f"Column '{col}' not found in adata.obs.",
                    "available_columns": [str(c) for c in adata.obs.columns],
                    "active_dataset": active,
                }
            vals = _col_unique(adata.obs[resolved])
            md = f"**{resolved}** — {len(vals)} unique values:\n\n" + "\n".join(f"- {v}" for v in vals)
            return {
                "ok": True,
                "status": "success",
                "active_dataset": active,
                "column": resolved,
                "n_unique": len(vals),
                "values": vals,
                "inline_markdown": md,
            }

        # All categorical columns
        results = []
        for c in adata.obs.columns:
            series = adata.obs[c]
            if hasattr(series, "cat") or str(series.dtype) in ("category", "object"):
                vals = _col_unique(series)
                results.append({"column": str(c), "n_unique": len(vals), "values": vals})

        md_lines = [f"**All categorical obs columns ({len(results)} total):**", ""]
        for r in results:
            md_lines.append(f"**{r['column']}** ({r['n_unique']} unique): " + ", ".join(r["values"]))
        return {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "columns": results,
            "inline_markdown": "\n".join(md_lines),
        }

    @app.post(
        "/generate_rank_genes_groups_stacked_violin",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_rank_genes_groups_stacked_violin",
        summary="Rank genes stacked violin",
        description="Stacked violins of top-ranked genes. Requires precomputed rank_genes_groups.",
    )
    async def generate_rank_genes_groups_stacked_violin(
        payload: RankGenesGroupsStackedViolinRequest,
    ) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

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
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_highest_expr_genes",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_highest_expr_genes",
        summary="Highest expressed genes",
        description=(
            "Plots sc.pl.highest_expr_genes — for each gene, computes the fraction of counts "
            "assigned to that gene within each cell, then shows the n_top genes with the highest "
            "mean fraction as horizontal boxplots. "
            "Useful for quality control: expect mitochondrial genes, actin, ribosomal proteins, and MALAT1. "
            "Does not require any prior computation — works directly on the raw count matrix."
        ),
    )
    async def generate_highest_expr_genes(payload: HighestExprGenesRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        try:
            result = plotter.run_highest_expr_genes(
                adata,
                n_top=int(payload.n_top),
                layer=str(payload.layer or "").strip() or None,
                gene_symbols=str(payload.gene_symbols or "").strip() or None,
                log=bool(payload.log),
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
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/generate_scatter",
        summary="Scatter plot",
        description=(
            "Generates a scatter plot using sc.pl.scatter for the active dataset. "
            "Supports QC scatter plots (e.g. n_counts vs pct_counts_mt), embedding-based scatter "
            "(basis='umap'), and any combination of obs columns or genes on x/y axes. "
            "Color by observation columns or gene expression."
        ),
    )
    async def generate_scatter(payload: ScatterPlotRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        try:
            color_list: list[str] = []
            raw_color = str(payload.color_json or "").strip()
            if raw_color and raw_color != "[]":
                try:
                    parsed = json.loads(raw_color)
                    color_list = [str(c) for c in parsed] if isinstance(parsed, list) else [str(parsed)]
                except json.JSONDecodeError:
                    color_list = [c.strip() for c in raw_color.split(",") if c.strip()]

            groups_list: list[str] | None = None
            raw_groups = str(payload.groups_json or "").strip()
            if raw_groups and raw_groups != "[]":
                try:
                    parsed_g = json.loads(raw_groups)
                    groups_list = [str(g) for g in parsed_g] if isinstance(parsed_g, list) else [str(parsed_g)]
                except json.JSONDecodeError:
                    groups_list = [g.strip() for g in raw_groups.split(",") if g.strip()]

            result = plotter.run_scatter(
                adata,
                x=str(payload.x or "").strip() or None,
                y=str(payload.y or "").strip() or None,
                color=color_list or None,
                basis=str(payload.basis or "").strip() or None,
                use_raw=payload.use_raw,
                layer=str(payload.layer or "").strip() or None,
                groups=groups_list,
                components=str(payload.components or "").strip() or None,
                sort_order=bool(payload.sort_order),
                legend_loc=str(payload.legend_loc or "right margin").strip() or "right margin",
                size=payload.size,
                color_map=str(payload.color_map or "").strip() or None,
                frameon=payload.frameon,
                title=str(payload.title or "").strip() or None,
                figsize_width=payload.figsize_width,
                figsize_height=payload.figsize_height,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
        }, output_markdown)

        response: dict[str, Any] = {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "inline_markdown": canonical_response_markdown or output_markdown,
        }
        return response

    @app.post(
        "/obs_count_table",
        dependencies=[Depends(_require_api_key)],
        operation_id="obs_count_table",
        summary="Cell count table for two OBS categories",
        description=(
            "Groups the active dataset's obs by two categorical columns and returns n_cells per group. "
            "If the result has more than 10 rows, only the top-5 and bottom-5 rows are returned."
        ),
    )
    async def obs_count_table(payload: ObsCountTableRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}

        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

        obs_cols = list(adata.obs.columns)
        for col in (payload.row_col, payload.col_col):
            if col not in obs_cols:
                return {
                    "ok": False,
                    "status": "column_not_found",
                    "message": f"Column '{col}' not found in adata.obs.",
                    "available_columns": obs_cols,
                    "active_dataset": active,
                }

        counts: pd.DataFrame = (
            adata.obs
            .groupby([payload.row_col, payload.col_col], observed=True)
            .size()
            .reset_index(name="n_cells")
        )
        counts = counts.sort_values("n_cells", ascending=False).reset_index(drop=True)

        total_rows = len(counts)
        truncated = total_rows > 10
        if truncated:
            display_df = pd.concat([counts.head(5), counts.tail(5)], ignore_index=True)
        else:
            display_df = counts

        # Save full CSV for download
        output_dir = Path(os.environ.get("OUTPUT_DIR", "/code/out/genopixel"))
        csv_dir = output_dir / "count_tables"
        csv_dir.mkdir(parents=True, exist_ok=True)
        csv_filename = f"count_table_{payload.row_col}_x_{payload.col_col}_{int(time.time())}.csv"
        csv_path = csv_dir / csv_filename
        counts.to_csv(csv_path, index=False)
        csv_url = _public_assets_url(str(csv_path))

        # Build markdown table (truncated for display)
        header = (
            f"**Cell count table — `{payload.row_col}` × `{payload.col_col}`** "
            f"({total_rows} groups, {int(counts['n_cells'].sum()):,} cells total"
            + (", showing top-5 and bottom-5)" if truncated else ")")
        )
        if csv_url:
            header += f"  \n[Download full CSV ↓]({csv_url})"

        md_lines: list[str] = [
            header,
            "",
            f"| {payload.row_col} | {payload.col_col} | n_cells |",
            f"|{'---'}|{'---'}|{'---'}|",
        ]
        if truncated:
            for _, row in display_df.head(5).iterrows():
                md_lines.append(f"| {row[payload.row_col]} | {row[payload.col_col]} | {int(row['n_cells']):,} |")
            md_lines.append(f"| … | … | … ({total_rows - 10} rows omitted) |")
            for _, row in display_df.tail(5).iterrows():
                md_lines.append(f"| {row[payload.row_col]} | {row[payload.col_col]} | {int(row['n_cells']):,} |")
        else:
            for _, row in display_df.iterrows():
                md_lines.append(f"| {row[payload.row_col]} | {row[payload.col_col]} | {int(row['n_cells']):,} |")

        return {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "row_col": payload.row_col,
            "col_col": payload.col_col,
            "total_groups": total_rows,
            "total_cells": int(counts["n_cells"].sum()),
            "truncated": truncated,
            "csv_url": csv_url,
            "rows": display_df.to_dict(orient="records"),
            "inline_markdown": "\n".join(md_lines),
        }

    @app.post(
        "/generate_spatial_scatter",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_spatial_scatter",
        summary="Spatial scatter plot",
        description=(
            "Generates a spatial scatter plot for spatial transcriptomics datasets "
            "using squidpy sq.pl.spatial_scatter. Only valid when adata.uns has 'spatial' "
            "and adata.obsm has 'spatial'. Color by adata.obs columns or gene names/ENS IDs. "
            "Automatically selects the cell type column when color is not specified. "
            "If adata.uns['spatial'] has library entries (image + spot info), the tissue image "
            "is overlaid; otherwise spots are plotted without an image background. "
            "Supports layout (wspace, hspace, ncols), shape, size, alpha, cmap, palette, "
            "legend_loc, colorbar, frameon, outline, groups, layer, and crop_coord."
        ),
    )
    async def generate_spatial_scatter(payload: SpatialScatterRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

        if "spatial" not in adata.uns or "spatial" not in adata.obsm:
            return {
                "ok": False,
                "status": "not_spatial",
                "message": (
                    "Active dataset is not a spatial transcriptomics dataset. "
                    "Requires adata.uns['spatial'] and adata.obsm['spatial']."
                ),
                "active_dataset": active,
            }

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        color = _parse_string_list(payload.color_json)

        groups = _parse_string_list(payload.groups_json)
        crop = tuple(payload.crop_coord) if payload.crop_coord and len(payload.crop_coord) == 4 else None  # type: ignore[arg-type]

        try:
            result = plotter.run_spatial_scatter(
                adata,
                color=color or None,
                figsize_width=payload.figsize_width,
                figsize_height=payload.figsize_height,
                title=str(payload.title or "").strip() or None,
                wspace=payload.wspace,
                hspace=payload.hspace,
                ncols=payload.ncols,
                shape=payload.shape,
                size=payload.size,
                alpha=payload.alpha,
                cmap=payload.cmap,
                palette=payload.palette,
                na_color=payload.na_color,
                img=payload.img,
                img_alpha=payload.img_alpha,
                crop_coord=crop,
                legend_loc=payload.legend_loc,
                legend_fontsize=payload.legend_fontsize,
                colorbar=payload.colorbar,
                frameon=payload.frameon,
                outline=payload.outline,
                groups=groups or None,
                layer=payload.layer,
                use_raw=payload.use_raw,
            )
        except Exception as exc:
            return {"ok": False, "status": "plot_error", "message": str(exc), "active_dataset": active}

        output_file = str(result.output_file.resolve())
        output_url = _public_assets_url(output_file)
        output_markdown = _output_markdown(output_file, output_url)
        canonical_response_markdown = build_canonical_response_markdown(active, {
            "plot_type": result.plot_type,
            "display_plot_type": result.display_plot_type,
        }, output_markdown)

        return {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "color_used": result.color_columns,
            "inline_markdown": canonical_response_markdown or output_markdown,
        }

    @app.post(
        "/generate_nhood_enrichment",
        dependencies=[Depends(_require_api_key)],
        operation_id="generate_nhood_enrichment",
        summary="Neighborhood enrichment analysis and plot",
        description=(
            "Computes and plots neighborhood enrichment for spatial transcriptomics datasets "
            "using squidpy. Runs sq.gr.spatial_neighbors, sq.gr.nhood_enrichment, and "
            "sq.pl.nhood_enrichment in sequence. Only valid when adata.uns has 'spatial' "
            "and adata.obsm has 'spatial'. cluster_key defaults to the auto-detected cell "
            "type column. mode can be 'zscore' (default) or 'count'."
        ),
    )
    async def generate_nhood_enrichment(payload: NhoodEnrichmentRequest) -> dict[str, Any]:
        try:
            adata, active = RUNTIME_STATE.require_active_adata()
        except NoActiveDatasetError as exc:
            return {"ok": False, "status": "no_active_dataset", "message": str(exc)}
        try:
            adata = _apply_obs_filter(adata, payload.obs_filter_json)
        except ValueError as exc:
            return {"ok": False, "status": "obs_filter_error", "message": str(exc), "active_dataset": active}

        if "spatial" not in adata.uns or "spatial" not in adata.obsm:
            return {
                "ok": False,
                "status": "not_spatial",
                "message": (
                    "Active dataset is not a spatial transcriptomics dataset. "
                    "Requires adata.uns['spatial'] and adata.obsm['spatial']."
                ),
                "active_dataset": active,
            }

        plotter = _plotter()
        sync_active_dataset = getattr(plotter, "sync_active_dataset", None)
        if callable(sync_active_dataset):
            sync_active_dataset(active, adata)

        try:
            result = plotter.run_nhood_enrichment(
                adata,
                cluster_key=str(payload.cluster_key or "").strip() or None,
                mode=str(payload.mode or "zscore").strip() or "zscore",
                figsize_width=payload.figsize_width,
                figsize_height=payload.figsize_height,
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
        }, output_markdown)

        return {
            "ok": True,
            "status": "success",
            "active_dataset": active,
            "cluster_key_used": result.resolved_groupby,
            "inline_markdown": canonical_response_markdown or output_markdown,
        }

    return app


app = create_app()
