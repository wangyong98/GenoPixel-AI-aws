from __future__ import annotations

from datetime import datetime
import math
from pathlib import Path
import re

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
from matplotlib import rc_context

from gp_models import PlotRequest, PlotResult


ALLOWED_PLOTS = {
    "umap",
    "tsne",
    "gene_cell_embedding",
    "violin",
    "dotplot",
    "matrixplot",
    "heatmap",
    "rank_genes_groups_dotplot",
    "rank_genes_groups_matrixplot",
    "rank_genes_groups_stacked_violin",
    "rank_genes_groups_heatmap",
    "correlation_matrix",
    "cell_counts_barplot",
    "cell_type_proportion_barplot",
    "spatial_scatter",
}

GENE_SYNONYM_MAP = {
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


class ScanpyPlotExecutor:
    EMBEDDING_FIGSIZE = (10, 8)
    EMBEDDING_FIGSIZE_WITH_RIGHT_LEGEND = (16, 8)
    EMBEDDING_LEGEND_FONTSIZE = 7
    EMBEDDING_LEGEND_MIN_FONTSIZE = 5
    EMBEDDING_LEGEND_ONDATA_THRESHOLD = 20
    UMAP_TSNE_POINT_SIZE = 3.0
    DOTPLOT_HEIGHT = 8
    DOTPLOT_MIN_WIDTH = 10
    DOTPLOT_MAX_WIDTH = 24
    DOTPLOT_BASE_WIDTH = 6
    DOTPLOT_WIDTH_PER_GENE = 0.7
    DOTPLOT_WIDTH_PER_GROUP = 0.2
    _RANK_GENES_CACHE_KEY = "_genopixel_rank_genes_groups_cache"
    _RANK_GENES_PREP_WARNING = (
        "rank_genes_groups cache not found for active dataset. Starting sc.tl.rank_genes_groups "
        "(this can take a long time)."
    )
    _PLACEHOLDER_MARKER_TOKENS = {
        "obsorgene",
        "obs",
        "gene",
        "genes",
        "marker",
        "markers",
        "obscolumn",
        "groupby",
        "column",
        "none",
        "null",
        "na",
    }
    DEFAULT_ACTIVE_CELL_TYPE_COLUMN = "author_cell_type"

    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._active_dataset_key: tuple[object, ...] | None = None
        self._session_cell_type_column: str | None = None
        self._session_markers: list[str] = []
        self._rank_genes_groups_computed_this_run: bool = False
        self._rank_genes_groups_notice: str | None = None
        self._rank_genes_groups_used_groupby: str | None = None

    def set_markers(self, markers: list[str]) -> None:
        self._session_markers = [str(m).strip() for m in markers if str(m).strip()]

    def get_markers(self) -> list[str]:
        return list(self._session_markers)

    def sync_active_dataset(self, active_dataset: dict[str, object] | None, adata: ad.AnnData | None = None) -> None:
        key = self._build_active_dataset_key(active_dataset)
        if key == self._active_dataset_key:
            return

        self._active_dataset_key = key
        self._session_cell_type_column = self._resolve_default_active_cell_type_column(adata)

    @staticmethod
    def _build_active_dataset_key(active_dataset: dict[str, object] | None) -> tuple[object, ...] | None:
        if not isinstance(active_dataset, dict):
            return None
        if not active_dataset.get("loaded"):
            return None
        return (
            active_dataset.get("all_excel_row"),
            active_dataset.get("multiple_excel_row"),
            active_dataset.get("h5ad_path"),
            active_dataset.get("backed"),
            active_dataset.get("loaded_at"),
        )

    def _resolve_default_active_cell_type_column(self, adata: ad.AnnData | None) -> str:
        default_column = self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        if adata is None:
            return default_column
        if default_column in adata.obs.columns:
            return default_column

        lower_to_original = {str(column).lower(): str(column) for column in adata.obs.columns}
        return lower_to_original.get(default_column.lower(), default_column)

    def run(self, adata: ad.AnnData, request: PlotRequest) -> PlotResult:
        self._rank_genes_groups_computed_this_run = False
        self._rank_genes_groups_notice = None
        self._rank_genes_groups_used_groupby = None

        requested_plot_type = request.plot_type.lower().strip()
        if requested_plot_type not in ALLOWED_PLOTS:
            raise ValueError(f"Unsupported plot type '{requested_plot_type}'. Allowed: {sorted(ALLOWED_PLOTS)}")

        actual_plot_type = requested_plot_type
        embedding_basis: str | None = None
        color_columns: list[str] | None = None

        plt.close("all")
        if requested_plot_type == "umap":
            actual_plot_type, embedding_basis, color_columns = self._plot_umap(adata, request)
        elif requested_plot_type == "tsne":
            actual_plot_type, embedding_basis, color_columns = self._plot_tsne(adata, request)
        elif requested_plot_type == "gene_cell_embedding":
            actual_plot_type, embedding_basis, color_columns = self._plot_gene_cell_embedding(adata, request)
        elif requested_plot_type == "violin":
            self._plot_violin(adata, request)
        elif requested_plot_type == "dotplot":
            self._plot_dotplot(adata, request)
        elif requested_plot_type == "matrixplot":
            self._plot_matrixplot(adata, request)
        elif requested_plot_type == "heatmap":
            self._plot_heatmap(adata, request)
        elif requested_plot_type == "rank_genes_groups_dotplot":
            self._plot_rank_genes_groups_dotplot(adata, request)
        elif requested_plot_type == "rank_genes_groups_matrixplot":
            self._plot_rank_genes_groups_matrixplot(adata, request)
        elif requested_plot_type == "rank_genes_groups_stacked_violin":
            self._plot_rank_genes_groups_stacked_violin(adata, request)
        elif requested_plot_type == "rank_genes_groups_heatmap":
            self._plot_rank_genes_groups_heatmap(adata, request)
        elif requested_plot_type == "correlation_matrix":
            self._plot_correlation_matrix(adata, request)
        elif requested_plot_type == "cell_counts_barplot":
            self._plot_cell_counts_barplot(adata, request)
        elif requested_plot_type == "cell_type_proportion_barplot":
            self._plot_cell_type_proportion_barplot(adata, request)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"{actual_plot_type}_{timestamp}.png"

        fig = plt.gcf()
        if actual_plot_type not in {"violin", "cell_counts_barplot", "cell_type_proportion_barplot"}:
            fig.tight_layout()
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")
        metadata = self._derive_plot_metadata(
            adata=adata,
            request=request,
            requested_plot_type=requested_plot_type,
            actual_plot_type=actual_plot_type,
            embedding_basis=embedding_basis,
            color_columns=color_columns,
        )
        return PlotResult(
            plot_type=actual_plot_type,
            output_file=out_file,
            embedding_basis=embedding_basis,
            color_columns=metadata["color_columns"],
            resolved_genes=metadata["resolved_genes"],
            resolved_groupby=metadata["resolved_groupby"],
            resolved_coloring_label=metadata["resolved_coloring_label"],
            display_plot_type=metadata["display_plot_type"],
            rank_genes_groups_computed=self._rank_genes_groups_computed_this_run if requested_plot_type.startswith("rank_genes_groups_") else None,
            rank_genes_groups_notice=self._rank_genes_groups_notice if requested_plot_type.startswith("rank_genes_groups_") else None,
        )

    def run_violin(
        self,
        adata: ad.AnnData,
        keys: list[str],
        *,
        groupby: str | None = None,
        rotation: float = 45.0,
        log: bool = False,
        use_raw: bool | None = None,
        stripplot: bool = True,
        jitter: float | bool = True,
        size: int = 1,
        layer: str | None = None,
        density_norm: str = "width",
        order: list[str] | None = None,
        multi_panel: bool | None = None,
        xlabel: str = "",
        ylabel: str | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized violin plot with dynamic sizing and label alignment."""
        plt.close("all")

        if not keys:
            keys = list(self._session_markers)
        plot_keys = self._resolve_gene_names(adata, keys)
        label_map = self._resolve_feature_labels_for_var_keys(adata, plot_keys)

        resolved_groupby = self._resolve_cell_type_column(
            adata, groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        )

        n_categories = 0
        max_label_len = 0
        if resolved_groupby and resolved_groupby in adata.obs.columns:
            n_categories = int(adata.obs[resolved_groupby].nunique(dropna=True))
            non_null = adata.obs[resolved_groupby].dropna().astype(str)
            max_label_len = int(non_null.str.len().max()) if not non_null.empty else 0

        fig_width = max(10.0, min(n_categories * 0.6, 36.0)) if n_categories else 12.0

        if n_categories > 40 or max_label_len > 30:
            xlabel_fontsize = 6
        elif n_categories > 20 or max_label_len > 20:
            xlabel_fontsize = 8
        else:
            xlabel_fontsize = 10

        sc_kwargs: dict[str, object] = {
            "keys": plot_keys,
            "show": False,
            "log": log,
            "stripplot": stripplot,
            "jitter": jitter,
            "size": size,
            "density_norm": density_norm,
            "rotation": rotation,
        }
        if resolved_groupby:
            sc_kwargs["groupby"] = resolved_groupby
        if use_raw is not None:
            sc_kwargs["use_raw"] = use_raw
        else:
            sc_kwargs["use_raw"] = False
        if layer:
            sc_kwargs["layer"] = layer
        if order:
            sc_kwargs["order"] = order
        if multi_panel is not None:
            sc_kwargs["multi_panel"] = multi_panel
        if xlabel:
            sc_kwargs["xlabel"] = xlabel
        if ylabel:
            sc_kwargs["ylabel"] = ylabel

        with rc_context({"figure.figsize": (fig_width, 6.0)}):
            sc.pl.violin(adata, **sc_kwargs)

        self._relabel_violin_axes(label_map)

        fig = plt.gcf()
        for ax in fig.axes:
            plt.setp(ax.get_xticklabels(), rotation=rotation, ha="right", fontsize=xlabel_fontsize)
        if title:
            fig.axes[0].set_title(str(title))

        # Remove legend — x-axis labels already identify each group.
        for ax in fig.axes:
            leg = ax.get_legend()
            if leg is not None:
                leg.remove()

        bottom_margin = min(0.08 + max_label_len * 0.008, 0.45)
        fig.subplots_adjust(bottom=bottom_margin)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"violin_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="violin",
            output_file=out_file,
            resolved_genes=plot_keys or None,
            resolved_groupby=resolved_groupby,
            display_plot_type="Violin plot",
        )

    def run_heatmap(
        self,
        adata: ad.AnnData,
        var_names: list[str],
        *,
        groupby: str | None = None,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        dendrogram: bool = False,
        gene_symbols: str | None = None,
        layer: str | None = None,
        standard_scale: str | None = None,
        swap_axes: bool = True,
        show_gene_labels: bool | None = None,
        figsize: tuple[float, float] | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized heatmap with automatic figsize and gene resolution."""
        plt.close("all")

        if not var_names:
            var_names = list(self._session_markers)
        if not var_names:
            raise ValueError(
                "No genes provided and no session markers are set. "
                "Call set_markers first or pass genes explicitly."
            )

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = (
            self._resolve_gene_names(adata, var_names)
        )

        resolved_groupby = self._resolve_cell_type_column(
            adata, groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        )
        if not resolved_groupby:
            raise ValueError("heatmap requires a groupby column; none found.")

        # Auto-compute figsize if not given: width ~ gene count, height ~ category count
        if figsize is None:
            gene_count = len(plot_genes)
            cat_count = int(adata.obs[resolved_groupby].nunique(dropna=True)) if resolved_groupby in adata.obs.columns else 10
            if swap_axes:
                width = max(8.0, min(cat_count * 0.45 + 4, 28.0))
                height = max(4.0, min(gene_count * 0.35 + 2, 24.0))
            else:
                width = max(8.0, min(gene_count * 0.45 + 4, 28.0))
                height = max(4.0, min(cat_count * 0.35 + 2, 24.0))
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": resolved_groupby,
            "show": False,
            "log": log,
            "num_categories": num_categories,
            "dendrogram": dendrogram,
            "swap_axes": swap_axes,
            "figsize": figsize,
        }
        if use_gene_symbols:
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if show_gene_labels is not None:
            sc_kwargs["show_gene_labels"] = show_gene_labels
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter

        sc.pl.heatmap(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()

        # Scale label font size so category labels don't overlap.
        # The groupby labels appear on y-axis (default) or x-axis (swap_axes).
        # Pick font size based on category count and max label length.
        cat_count = (
            int(adata.obs[resolved_groupby].nunique(dropna=True))
            if resolved_groupby in adata.obs.columns
            else 0
        )
        max_cat_len = 0
        if resolved_groupby in adata.obs.columns:
            max_cat_len = int(adata.obs[resolved_groupby].dropna().astype(str).str.len().max())

        if cat_count > 30 or max_cat_len > 25:
            label_fontsize = 6
        elif cat_count > 15 or max_cat_len > 15:
            label_fontsize = 8
        else:
            label_fontsize = 9

        for ax in fig.axes:
            if swap_axes:
                # groupby on x-axis — rotate so labels don't overlap
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=label_fontsize)
            else:
                # groupby on y-axis
                plt.setp(ax.get_yticklabels(), fontsize=label_fontsize)

        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"heatmap_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="heatmap",
            output_file=out_file,
            resolved_genes=plot_genes or None,
            resolved_groupby=resolved_groupby,
            display_plot_type="Heatmap",
        )

    def run_dotplot(
        self,
        adata: ad.AnnData,
        var_names: list[str],
        *,
        groupby: str | None = None,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: list[str] | None = None,
        expression_cutoff: float = 0.0,
        mean_only_expressed: bool = False,
        standard_scale: str | None = None,
        dendrogram: bool = False,
        gene_symbols: str | None = None,
        layer: str | None = None,
        swap_axes: bool = False,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        cmap: str = "Reds",
        dot_max: float | None = None,
        dot_min: float | None = None,
        smallest_dot: float = 0.0,
        colorbar_title: str | None = None,
        size_title: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized dot plot with automatic figsize and gene resolution."""
        plt.close("all")

        if not var_names:
            var_names = list(self._session_markers)
        if not var_names:
            raise ValueError(
                "No genes provided and no session markers are set. "
                "Call set_markers first or pass genes explicitly."
            )

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = (
            self._resolve_gene_names(adata, var_names)
        )

        resolved_groupby = self._resolve_cell_type_column(
            adata, groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        )
        if not resolved_groupby:
            raise ValueError("dotplot requires a groupby column; none found.")

        if figsize is None:
            figsize = self._compute_dotplot_figsize(adata, groupby=resolved_groupby, gene_count=len(plot_genes))

        sc_kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": resolved_groupby,
            "show": False,
            "log": log,
            "num_categories": num_categories,
            "expression_cutoff": expression_cutoff,
            "mean_only_expressed": mean_only_expressed,
            "dendrogram": dendrogram,
            "swap_axes": swap_axes,
            "cmap": cmap,
            "smallest_dot": smallest_dot,
            "figsize": figsize,
        }
        if use_gene_symbols:
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if categories_order:
            sc_kwargs["categories_order"] = categories_order
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if dot_max is not None:
            sc_kwargs["dot_max"] = dot_max
        if dot_min is not None:
            sc_kwargs["dot_min"] = dot_min
        if colorbar_title:
            sc_kwargs["colorbar_title"] = colorbar_title
        if size_title:
            sc_kwargs["size_title"] = size_title

        sc.pl.dotplot(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()

        # Rotate labels when swap_axes puts groupby on x-axis
        if swap_axes:
            cat_count = (
                int(adata.obs[resolved_groupby].nunique(dropna=True))
                if resolved_groupby in adata.obs.columns
                else 0
            )
            max_cat_len = 0
            if resolved_groupby in adata.obs.columns:
                max_cat_len = int(adata.obs[resolved_groupby].dropna().astype(str).str.len().max())
            if cat_count > 30 or max_cat_len > 25:
                label_fontsize = 6
            elif cat_count > 15 or max_cat_len > 15:
                label_fontsize = 8
            else:
                label_fontsize = 9
            for ax in fig.axes:
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=label_fontsize)

        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"dotplot_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="dotplot",
            output_file=out_file,
            resolved_genes=plot_genes or None,
            resolved_groupby=resolved_groupby,
            display_plot_type="Dot plot",
        )

    def run_tracksplot(
        self,
        adata: ad.AnnData,
        var_names: list[str],
        *,
        groupby: str | None = None,
        use_raw: bool | None = None,
        log: bool = False,
        dendrogram: bool = False,
        gene_symbols: str | None = None,
        layer: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized tracksplot with automatic figsize and gene resolution."""
        plt.close("all")

        if not var_names:
            var_names = list(self._session_markers)
        if not var_names:
            raise ValueError(
                "No genes provided and no session markers are set. "
                "Call set_markers first or pass genes explicitly."
            )

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = (
            self._resolve_gene_names(adata, var_names)
        )

        resolved_groupby = self._resolve_cell_type_column(
            adata, groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        )
        if not resolved_groupby:
            raise ValueError("tracksplot requires a groupby column; none found.")

        if figsize is None:
            gene_count = len(plot_genes)
            cat_count = (
                int(adata.obs[resolved_groupby].nunique(dropna=True))
                if resolved_groupby in adata.obs.columns
                else 10
            )
            width = max(10.0, min(gene_count * 0.5 + 4, 28.0))
            height = max(4.0, min(cat_count * 0.4 + 2, 24.0))
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": resolved_groupby,
            "show": False,
            "log": log,
            "dendrogram": dendrogram,
            "figsize": figsize,
        }
        if use_gene_symbols:
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer

        sc.pl.tracksplot(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()

        # Scale category label font size on y-axis (groupby labels appear there)
        cat_count = (
            int(adata.obs[resolved_groupby].nunique(dropna=True))
            if resolved_groupby in adata.obs.columns
            else 0
        )
        max_cat_len = 0
        if resolved_groupby in adata.obs.columns:
            max_cat_len = int(adata.obs[resolved_groupby].dropna().astype(str).str.len().max())
        if cat_count > 30 or max_cat_len > 25:
            label_fontsize = 6
        elif cat_count > 15 or max_cat_len > 15:
            label_fontsize = 8
        else:
            label_fontsize = 9
        for ax in fig.axes:
            plt.setp(ax.get_yticklabels(), fontsize=label_fontsize)

        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"tracksplot_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="tracksplot",
            output_file=out_file,
            resolved_genes=plot_genes or None,
            resolved_groupby=resolved_groupby,
            display_plot_type="Tracks plot",
        )

    def run_stacked_violin(
        self,
        adata: ad.AnnData,
        var_names: list[str],
        *,
        groupby: str | None = None,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        dendrogram: bool = False,
        gene_symbols: str | None = None,
        layer: str | None = None,
        standard_scale: str | None = None,
        categories_order: list[str] | None = None,
        swap_axes: bool = False,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        cmap: str = "Blues",
        stripplot: bool = False,
        jitter: bool | float = False,
        size: int = 1,
        row_palette: str | None = None,
        yticklabels: bool = False,
        colorbar_title: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized stacked violin plot with automatic figsize and gene resolution."""
        plt.close("all")

        if not var_names:
            var_names = list(self._session_markers)
        if not var_names:
            raise ValueError(
                "No genes provided and no session markers are set. "
                "Call set_markers first or pass genes explicitly."
            )

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = (
            self._resolve_gene_names(adata, var_names)
        )

        resolved_groupby = self._resolve_cell_type_column(
            adata, groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        )
        if not resolved_groupby:
            raise ValueError("stacked_violin requires a groupby column; none found.")

        if figsize is None:
            gene_count = len(plot_genes)
            cat_count = (
                int(adata.obs[resolved_groupby].nunique(dropna=True))
                if resolved_groupby in adata.obs.columns
                else 10
            )
            if swap_axes:
                width = max(8.0, min(cat_count * 0.45 + 4, 28.0))
                height = max(4.0, min(gene_count * 0.5 + 2, 28.0))
            else:
                width = max(8.0, min(gene_count * 0.5 + 4, 28.0))
                height = max(4.0, min(cat_count * 0.45 + 2, 28.0))
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": resolved_groupby,
            "show": False,
            "log": log,
            "num_categories": num_categories,
            "dendrogram": dendrogram,
            "swap_axes": swap_axes,
            "cmap": cmap,
            "stripplot": stripplot,
            "jitter": jitter,
            "size": size,
            "yticklabels": yticklabels,
            "figsize": figsize,
        }
        if use_gene_symbols:
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if categories_order:
            sc_kwargs["categories_order"] = categories_order
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if row_palette:
            sc_kwargs["row_palette"] = row_palette
        if colorbar_title:
            sc_kwargs["colorbar_title"] = colorbar_title
        if title:
            sc_kwargs["title"] = title

        sc.pl.stacked_violin(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()

        # Rotate x-axis labels when swap_axes puts groupby on x-axis
        if swap_axes:
            cat_count = (
                int(adata.obs[resolved_groupby].nunique(dropna=True))
                if resolved_groupby in adata.obs.columns
                else 0
            )
            max_cat_len = 0
            if resolved_groupby in adata.obs.columns:
                max_cat_len = int(adata.obs[resolved_groupby].dropna().astype(str).str.len().max())
            if cat_count > 30 or max_cat_len > 25:
                label_fontsize = 6
            elif cat_count > 15 or max_cat_len > 15:
                label_fontsize = 8
            else:
                label_fontsize = 9
            for ax in fig.axes:
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=label_fontsize)
        else:
            # groupby on y-axis — apply font sizing
            cat_count = (
                int(adata.obs[resolved_groupby].nunique(dropna=True))
                if resolved_groupby in adata.obs.columns
                else 0
            )
            max_cat_len = 0
            if resolved_groupby in adata.obs.columns:
                max_cat_len = int(adata.obs[resolved_groupby].dropna().astype(str).str.len().max())
            if cat_count > 30 or max_cat_len > 25:
                label_fontsize = 6
            elif cat_count > 15 or max_cat_len > 15:
                label_fontsize = 8
            else:
                label_fontsize = 9
            for ax in fig.axes:
                plt.setp(ax.get_yticklabels(), fontsize=label_fontsize)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"stacked_violin_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="stacked_violin",
            output_file=out_file,
            resolved_genes=plot_genes or None,
            resolved_groupby=resolved_groupby,
            display_plot_type="Stacked violin plot",
        )

    def run_matrixplot(
        self,
        adata: ad.AnnData,
        var_names: list[str],
        *,
        groupby: str | None = None,
        use_raw: bool | None = None,
        log: bool = False,
        num_categories: int = 7,
        categories_order: list[str] | None = None,
        dendrogram: bool = False,
        gene_symbols: str | None = None,
        layer: str | None = None,
        standard_scale: str | None = None,
        swap_axes: bool = False,
        cmap: str = "viridis",
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        colorbar_title: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized matrix plot with automatic figsize and gene resolution."""
        plt.close("all")

        if not var_names:
            var_names = list(self._session_markers)
        if not var_names:
            raise ValueError(
                "No genes provided and no session markers are set. "
                "Call set_markers first or pass genes explicitly."
            )

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = (
            self._resolve_gene_names(adata, var_names)
        )

        resolved_groupby = self._resolve_cell_type_column(
            adata, groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        )
        if not resolved_groupby:
            raise ValueError("matrixplot requires a groupby column; none found.")

        if figsize is None:
            gene_count = len(plot_genes)
            cat_count = (
                int(adata.obs[resolved_groupby].nunique(dropna=True))
                if resolved_groupby in adata.obs.columns
                else 10
            )
            if swap_axes:
                width = max(8.0, min(cat_count * 0.45 + 4, 28.0))
                height = max(4.0, min(gene_count * 0.35 + 2, 24.0))
            else:
                width = max(8.0, min(gene_count * 0.45 + 4, 28.0))
                height = max(4.0, min(cat_count * 0.35 + 2, 24.0))
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": resolved_groupby,
            "show": False,
            "log": log,
            "num_categories": num_categories,
            "dendrogram": dendrogram,
            "swap_axes": swap_axes,
            "cmap": cmap,
            "figsize": figsize,
        }
        if use_gene_symbols:
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if categories_order:
            sc_kwargs["categories_order"] = categories_order
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if colorbar_title:
            sc_kwargs["colorbar_title"] = colorbar_title

        sc.pl.matrixplot(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()

        # Scale label font size based on category count and max label length
        cat_count = (
            int(adata.obs[resolved_groupby].nunique(dropna=True))
            if resolved_groupby in adata.obs.columns
            else 0
        )
        max_cat_len = 0
        if resolved_groupby in adata.obs.columns:
            max_cat_len = int(adata.obs[resolved_groupby].dropna().astype(str).str.len().max())
        if cat_count > 30 or max_cat_len > 25:
            label_fontsize = 6
        elif cat_count > 15 or max_cat_len > 15:
            label_fontsize = 8
        else:
            label_fontsize = 9

        for ax in fig.axes:
            if swap_axes:
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=label_fontsize)
            else:
                plt.setp(ax.get_yticklabels(), fontsize=label_fontsize)

        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"matrixplot_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="matrixplot",
            output_file=out_file,
            resolved_genes=plot_genes or None,
            resolved_groupby=resolved_groupby,
            display_plot_type="Matrix plot",
        )

    def run_clustermap(
        self,
        adata: ad.AnnData,
        var_names: list[str],
        *,
        obs_keys: str | None = None,
        use_raw: bool | None = None,
        standard_scale: str | None = None,
        z_score: int | None = None,
        method: str = "average",
        metric: str = "euclidean",
        cmap: str = "viridis",
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Clustermap (seaborn-backed) with optional gene subsetting and session marker fallback."""
        plt.close("all")

        if not var_names:
            var_names = list(self._session_markers)

        # Subset adata to requested genes if provided
        if var_names:
            use_gene_symbols = self._has_feature_name_column(adata)
            plot_genes = (
                self._resolve_gene_names(adata, var_names)
            )
            # Intersect with actual var_names present in adata
            valid_genes = [g for g in plot_genes if g in adata.var_names]
            if valid_genes:
                adata = adata[:, valid_genes]
            else:
                plot_genes = []
        else:
            plot_genes = list(adata.var_names)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "method": method,
            "metric": metric,
            "cmap": cmap,
        }
        if obs_keys:
            sc_kwargs["obs_keys"] = obs_keys
        sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if standard_scale is not None:
            sc_kwargs["standard_scale"] = standard_scale
        if z_score is not None:
            sc_kwargs["z_score"] = z_score
        if figsize is not None:
            sc_kwargs["figsize"] = figsize

        sc.pl.clustermap(adata, **sc_kwargs)

        fig = plt.gcf()
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"clustermap_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="clustermap",
            output_file=out_file,
            resolved_genes=plot_genes or None,
            resolved_groupby=obs_keys or None,
            display_plot_type="Cluster map",
        )

    def run_dendrogram(
        self,
        adata: ad.AnnData,
        *,
        groupby: str | None = None,
        dendrogram_key: str | None = None,
        orientation: str = "top",
        remove_labels: bool = False,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Dendrogram plot — auto-computes sc.tl.dendrogram if not already cached."""
        plt.close("all")

        resolved_groupby = self._resolve_cell_type_column(
            adata, groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN
        )
        if not resolved_groupby:
            raise ValueError("dendrogram requires a groupby column; none found.")

        # Auto-compute dendrogram if not already stored in uns
        key = dendrogram_key or f"dendrogram_{resolved_groupby}"
        if key not in adata.uns:
            sc.tl.dendrogram(adata, resolved_groupby)

        if figsize is None:
            cat_count = (
                int(adata.obs[resolved_groupby].nunique(dropna=True))
                if resolved_groupby in adata.obs.columns
                else 10
            )
            if orientation in {"top", "bottom"}:
                width = max(6.0, min(cat_count * 0.4 + 2, 24.0))
                figsize = (width, 3.0)
            else:
                figsize = (4.0, max(4.0, min(cat_count * 0.4 + 2, 20.0)))

        sc_kwargs: dict[str, object] = {
            "groupby": resolved_groupby,
            "orientation": orientation,
            "remove_labels": remove_labels,
            "show": False,
        }
        if dendrogram_key:
            sc_kwargs["dendrogram_key"] = dendrogram_key

        fig, ax = plt.subplots(figsize=figsize)
        sc_kwargs["ax"] = ax
        sc.pl.dendrogram(adata, **sc_kwargs)

        if title:
            ax.set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"dendrogram_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="dendrogram",
            output_file=out_file,
            resolved_groupby=resolved_groupby,
            display_plot_type="Dendrogram",
        )

    def run_tsne(
        self,
        adata: ad.AnnData,
        *,
        color: list[str] | None = None,
        use_raw: bool | None = None,
        layer: str | None = None,
        gene_symbols: str | None = None,
        legend_loc: str = "right margin",
        legend_fontsize: float | None = None,
        legend_fontweight: str = "bold",
        colorbar_loc: str = "right",
        color_map: str | None = None,
        palette: str | None = None,
        na_color: str = "lightgray",
        na_in_legend: bool = True,
        size: float | None = None,
        frameon: bool | None = None,
        vmin: str | float | None = None,
        vmax: str | float | None = None,
        vcenter: float | None = None,
        add_outline: bool = False,
        sort_order: bool = True,
        edges: bool = False,
        edges_width: float = 0.1,
        edges_color: str = "grey",
        groups: list[str] | None = None,
        ncols: int = 4,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized tSNE plot — auto-computes embedding if missing."""
        plt.close("all")

        # Ensure tSNE embedding exists — use partial matching (handles X_tsne_perplexity30 etc.)
        basis = self._find_embedding_basis(adata, "tsne")
        if basis is None:
            self._ensure_embedding(adata, "tsne")
            basis = self._find_embedding_basis(adata, "tsne")
        if basis is None:
            available = self._available_embeddings(adata)
            raise ValueError(
                f"No tSNE embedding found or could be computed. "
                f"Available embeddings: {available or ['none']}."
            )

        # Resolve color columns
        resolved_color, gene_symbols_column = self._resolve_embedding_color(
            adata, color or [], preferred_gene_symbols_column=gene_symbols
        )
        if not resolved_color:
            inferred = self._resolve_cell_type_column(adata)
            if inferred:
                resolved_color = [inferred]

        # Legend and figsize — reuse existing smart sizing
        override_legend = legend_loc != "right margin"
        if override_legend:
            legend_options: dict[str, object] = {
                "legend_loc": legend_loc,
                "legend_fontsize": legend_fontsize or self.EMBEDDING_LEGEND_FONTSIZE,
                "legend_fontweight": legend_fontweight,
            }
            figure_size = self.EMBEDDING_FIGSIZE
        else:
            legend_options = self._embedding_legend_options(adata, resolved_color)
            legend_options["legend_fontweight"] = legend_fontweight
            if legend_fontsize is not None:
                legend_options["legend_fontsize"] = legend_fontsize
            figure_size = self._compute_embedding_figsize(adata, resolved_color, legend_options)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "sort_order": sort_order,
            "na_color": na_color,
            "na_in_legend": na_in_legend,
            "add_outline": add_outline,
            "edges": edges,
            "edges_width": edges_width,
            "edges_color": edges_color,
            "ncols": ncols,
            "colorbar_loc": colorbar_loc,
            **legend_options,
        }
        if resolved_color:
            sc_kwargs["color"] = resolved_color
        if gene_symbols_column:
            sc_kwargs["gene_symbols"] = gene_symbols_column
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if size is not None:
            sc_kwargs["size"] = float(size)
        if frameon is not None:
            sc_kwargs["frameon"] = bool(frameon)
        if color_map:
            sc_kwargs["color_map"] = color_map
        if palette:
            sc_kwargs["palette"] = palette
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if groups:
            sc_kwargs["groups"] = groups
        if title:
            sc_kwargs["title"] = title

        # Use sc.pl.tsne only when the matched basis is the canonical X_tsne;
        # otherwise fall back to sc.pl.embedding so the correct basis key is used.
        bare_basis = basis[2:] if basis.startswith("X_") else basis
        is_canonical = bare_basis == "tsne"
        with rc_context({"figure.figsize": figure_size}):
            if is_canonical:
                sc.pl.tsne(adata, **sc_kwargs)
            else:
                sc.pl.embedding(adata, basis=bare_basis, **sc_kwargs)

        fig = plt.gcf()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"tsne_{timestamp}.png"
        fig.tight_layout()
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        coloring_label = self._build_coloring_label(
            adata, resolved_color, self._resolve_cell_type_column(adata)
        )
        return PlotResult(
            plot_type="tsne",
            output_file=out_file,
            embedding_basis=basis,
            color_columns=resolved_color or None,
            resolved_coloring_label=coloring_label,
            display_plot_type="tSNE embedding",
        )

    def run_umap(
        self,
        adata: ad.AnnData,
        *,
        color: list[str] | None = None,
        use_raw: bool | None = None,
        layer: str | None = None,
        gene_symbols: str | None = None,
        legend_loc: str = "right margin",
        legend_fontsize: float | None = None,
        legend_fontweight: str = "bold",
        colorbar_loc: str = "right",
        color_map: str | None = None,
        palette: str | None = None,
        na_color: str = "lightgray",
        na_in_legend: bool = True,
        size: float | None = None,
        frameon: bool | None = None,
        vmin: str | float | None = None,
        vmax: str | float | None = None,
        vcenter: float | None = None,
        add_outline: bool = False,
        sort_order: bool = True,
        edges: bool = False,
        edges_width: float = 0.1,
        edges_color: str = "grey",
        groups: list[str] | None = None,
        ncols: int = 4,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized UMAP plot — uses pre-computed embedding only, never computes on the fly."""
        plt.close("all")

        # Find a pre-computed UMAP or tSNE embedding — never compute one
        basis = self._find_first_matching_embedding_basis(adata, ["umap", "tsne"])
        if basis is None:
            obsm_keys = list(adata.obsm.keys())
            keys_str = "\n".join(f"  - {k}" for k in obsm_keys) if obsm_keys else "  (none)"
            raise ValueError(
                f"No UMAP or tSNE embedding found in this dataset.\n\n"
                f"Available embeddings in adata.obsm:\n{keys_str}\n\n"
                f"Please choose one of the above keys to plot, or ask for a different plot type."
            )

        # Resolve color columns
        resolved_color, gene_symbols_column = self._resolve_embedding_color(
            adata, color or [], preferred_gene_symbols_column=gene_symbols
        )
        if not resolved_color:
            inferred = self._resolve_cell_type_column(adata)
            if inferred:
                resolved_color = [inferred]

        # Legend and figsize — reuse existing smart sizing
        override_legend = legend_loc != "right margin"
        if override_legend:
            legend_options: dict[str, object] = {
                "legend_loc": legend_loc,
                "legend_fontsize": legend_fontsize or self.EMBEDDING_LEGEND_FONTSIZE,
                "legend_fontweight": legend_fontweight,
            }
            figure_size = self.EMBEDDING_FIGSIZE
        else:
            legend_options = self._embedding_legend_options(adata, resolved_color)
            legend_options["legend_fontweight"] = legend_fontweight
            if legend_fontsize is not None:
                legend_options["legend_fontsize"] = legend_fontsize
            figure_size = self._compute_embedding_figsize(adata, resolved_color, legend_options)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "sort_order": sort_order,
            "na_color": na_color,
            "na_in_legend": na_in_legend,
            "add_outline": add_outline,
            "edges": edges,
            "edges_width": edges_width,
            "edges_color": edges_color,
            "ncols": ncols,
            "colorbar_loc": colorbar_loc,
            **legend_options,
        }
        if resolved_color:
            sc_kwargs["color"] = resolved_color
        if gene_symbols_column:
            sc_kwargs["gene_symbols"] = gene_symbols_column
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if size is not None:
            sc_kwargs["size"] = float(size)
        if frameon is not None:
            sc_kwargs["frameon"] = bool(frameon)
        if color_map:
            sc_kwargs["color_map"] = color_map
        if palette:
            sc_kwargs["palette"] = palette
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if groups:
            sc_kwargs["groups"] = groups
        if title:
            sc_kwargs["title"] = title

        # Use sc.pl.umap/tsne only for canonical keys; fall back to sc.pl.embedding for
        # non-standard names (e.g. X_umap_harmony, X_tsne_perplexity30)
        bare_basis = basis[2:] if basis.startswith("X_") else basis
        is_canonical_umap = bare_basis == "umap"
        is_canonical_tsne = bare_basis == "tsne"
        with rc_context({"figure.figsize": figure_size}):
            if is_canonical_umap:
                sc.pl.umap(adata, **sc_kwargs)
            elif is_canonical_tsne:
                sc.pl.tsne(adata, **sc_kwargs)
            else:
                sc.pl.embedding(adata, basis=bare_basis, **sc_kwargs)

        fig = plt.gcf()
        fig.tight_layout()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"umap_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        actual_plot_type = "umap" if is_canonical_umap else "tsne"
        coloring_label = self._build_coloring_label(
            adata, resolved_color, self._resolve_cell_type_column(adata)
        )
        return PlotResult(
            plot_type=actual_plot_type,
            output_file=out_file,
            embedding_basis=basis,
            color_columns=resolved_color or None,
            resolved_coloring_label=coloring_label,
            display_plot_type="UMAP embedding",
        )

    def run_diffmap(
        self,
        adata: ad.AnnData,
        *,
        color: list[str] | None = None,
        components: str | None = None,
        use_raw: bool | None = None,
        layer: str | None = None,
        gene_symbols: str | None = None,
        legend_loc: str = "right margin",
        legend_fontsize: float | None = None,
        legend_fontweight: str = "bold",
        colorbar_loc: str = "right",
        color_map: str | None = None,
        palette: str | None = None,
        na_color: str = "lightgray",
        na_in_legend: bool = True,
        size: float | None = None,
        frameon: bool | None = None,
        vmin: str | float | None = None,
        vmax: str | float | None = None,
        vcenter: float | None = None,
        add_outline: bool = False,
        sort_order: bool = True,
        edges: bool = False,
        edges_width: float = 0.1,
        edges_color: str = "grey",
        groups: list[str] | None = None,
        ncols: int = 4,
        title: str | None = None,
    ) -> PlotResult:
        """Fully-parameterized diffusion map plot — requires sc.tl.diffmap to have been run."""
        plt.close("all")

        diffmap_key = self._find_embedding_basis(adata, "diffmap")
        if diffmap_key is None:
            available = self._available_embeddings(adata)
            raise ValueError(
                f"Diffusion map embedding not found in adata.obsm. "
                f"Please run sc.tl.diffmap(adata) first — this requires sc.pp.neighbors(adata) to have been computed. "
                f"Available embeddings: {available or ['none']}."
            )

        resolved_color, gene_symbols_column = self._resolve_embedding_color(
            adata, color or [], preferred_gene_symbols_column=gene_symbols
        )
        if not resolved_color:
            inferred = self._resolve_cell_type_column(adata)
            if inferred:
                resolved_color = [inferred]

        override_legend = legend_loc != "right margin"
        if override_legend:
            legend_options: dict[str, object] = {
                "legend_loc": legend_loc,
                "legend_fontsize": legend_fontsize or self.EMBEDDING_LEGEND_FONTSIZE,
                "legend_fontweight": legend_fontweight,
            }
            figure_size = self.EMBEDDING_FIGSIZE
        else:
            legend_options = self._embedding_legend_options(adata, resolved_color)
            legend_options["legend_fontweight"] = legend_fontweight
            if legend_fontsize is not None:
                legend_options["legend_fontsize"] = legend_fontsize
            figure_size = self._compute_embedding_figsize(adata, resolved_color, legend_options)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "sort_order": sort_order,
            "na_color": na_color,
            "na_in_legend": na_in_legend,
            "add_outline": add_outline,
            "edges": edges,
            "edges_width": edges_width,
            "edges_color": edges_color,
            "ncols": ncols,
            "colorbar_loc": colorbar_loc,
            **legend_options,
        }
        if resolved_color:
            sc_kwargs["color"] = resolved_color
        if components:
            sc_kwargs["components"] = components
        if gene_symbols_column:
            sc_kwargs["gene_symbols"] = gene_symbols_column
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if size is not None:
            sc_kwargs["size"] = float(size)
        if frameon is not None:
            sc_kwargs["frameon"] = bool(frameon)
        if color_map:
            sc_kwargs["color_map"] = color_map
        if palette:
            sc_kwargs["palette"] = palette
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if groups:
            sc_kwargs["groups"] = groups
        if title:
            sc_kwargs["title"] = title

        bare_diffmap = diffmap_key[2:] if diffmap_key.startswith("X_") else diffmap_key
        is_canonical_diffmap = bare_diffmap == "diffmap"
        with rc_context({"figure.figsize": figure_size}):
            if is_canonical_diffmap:
                sc.pl.diffmap(adata, **sc_kwargs)
            else:
                sc.pl.embedding(adata, basis=bare_diffmap, **sc_kwargs)

        fig = plt.gcf()
        fig.tight_layout()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"diffmap_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        coloring_label = self._build_coloring_label(
            adata, resolved_color, self._resolve_cell_type_column(adata)
        )
        return PlotResult(
            plot_type="diffmap",
            output_file=out_file,
            embedding_basis=diffmap_key,
            color_columns=resolved_color or None,
            resolved_coloring_label=coloring_label,
            display_plot_type="Diffusion map",
        )

    def run_embedding(
        self,
        adata: ad.AnnData,
        basis: str,
        *,
        color: list[str] | None = None,
        components: str | None = None,
        use_raw: bool | None = None,
        layer: str | None = None,
        gene_symbols: str | None = None,
        legend_loc: str = "right margin",
        legend_fontsize: float | None = None,
        legend_fontweight: str = "bold",
        colorbar_loc: str = "right",
        color_map: str | None = None,
        palette: str | None = None,
        na_color: str = "lightgray",
        na_in_legend: bool = True,
        size: float | None = None,
        frameon: bool | None = None,
        vmin: str | float | None = None,
        vmax: str | float | None = None,
        vcenter: float | None = None,
        add_outline: bool = False,
        sort_order: bool = True,
        edges: bool = False,
        edges_width: float = 0.1,
        edges_color: str = "grey",
        groups: list[str] | None = None,
        projection: str = "2d",
        ncols: int = 4,
        title: str | None = None,
    ) -> PlotResult:
        """Generic embedding plot — user supplies the basis name (e.g. 'umap', 'pca', 'tsne')."""
        plt.close("all")

        # Try partial matching first (handles long names like X_umap_harmony, X_tsne_perplexity30)
        matched_key = self._find_embedding_basis(adata, basis)
        if matched_key is None:
            available = self._available_embeddings(adata)
            raise ValueError(
                f"Embedding '{basis}' not found in adata.obsm. "
                f"Available embeddings: {available or ['none']}."
            )
        # sc.pl.embedding expects the bare name (without X_ prefix)
        normalised_basis = matched_key[2:] if matched_key.startswith("X_") else matched_key
        obsm_key = matched_key if matched_key.startswith("X_") else f"X_{matched_key}"

        resolved_color, gene_symbols_column = self._resolve_embedding_color(
            adata, color or [], preferred_gene_symbols_column=gene_symbols
        )
        if not resolved_color:
            inferred = self._resolve_cell_type_column(adata)
            if inferred:
                resolved_color = [inferred]

        override_legend = legend_loc != "right margin"
        if override_legend:
            legend_options: dict[str, object] = {
                "legend_loc": legend_loc,
                "legend_fontsize": legend_fontsize or self.EMBEDDING_LEGEND_FONTSIZE,
                "legend_fontweight": legend_fontweight,
            }
            figure_size = self.EMBEDDING_FIGSIZE
        else:
            legend_options = self._embedding_legend_options(adata, resolved_color)
            legend_options["legend_fontweight"] = legend_fontweight
            if legend_fontsize is not None:
                legend_options["legend_fontsize"] = legend_fontsize
            figure_size = self._compute_embedding_figsize(adata, resolved_color, legend_options)

        sc_kwargs: dict[str, object] = {
            "basis": normalised_basis,
            "show": False,
            "sort_order": sort_order,
            "na_color": na_color,
            "na_in_legend": na_in_legend,
            "add_outline": add_outline,
            "edges": edges,
            "edges_width": edges_width,
            "edges_color": edges_color,
            "projection": projection,
            "ncols": ncols,
            "colorbar_loc": colorbar_loc,
            **legend_options,
        }
        if resolved_color:
            sc_kwargs["color"] = resolved_color
        if components:
            sc_kwargs["components"] = components
        if gene_symbols_column:
            sc_kwargs["gene_symbols"] = gene_symbols_column
            sc_kwargs["use_raw"] = False
        else:
            if gene_symbols:
                sc_kwargs["gene_symbols"] = gene_symbols
            sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if size is not None:
            sc_kwargs["size"] = float(size)
        if frameon is not None:
            sc_kwargs["frameon"] = bool(frameon)
        if color_map:
            sc_kwargs["color_map"] = color_map
        if palette:
            sc_kwargs["palette"] = palette
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if groups:
            sc_kwargs["groups"] = groups
        if title:
            sc_kwargs["title"] = title

        with rc_context({"figure.figsize": figure_size}):
            sc.pl.embedding(adata, **sc_kwargs)

        fig = plt.gcf()
        fig.tight_layout()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"embedding_{normalised_basis}_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        coloring_label = self._build_coloring_label(
            adata, resolved_color, self._resolve_cell_type_column(adata)
        )
        display_name = f"{normalised_basis.upper()} embedding"
        return PlotResult(
            plot_type="embedding",
            output_file=out_file,
            embedding_basis=obsm_key,
            color_columns=resolved_color or None,
            resolved_coloring_label=coloring_label,
            display_plot_type=display_name,
        )

    def run_rank_genes_groups(
        self,
        adata: ad.AnnData,
        *,
        groups: list[str] | None = None,
        n_genes: int = 5,
        gene_symbols: str | None = None,
        key: str = "rank_genes_groups",
        fontsize: int = 13,
        ncols: int = 3,
        sharey: bool = True,
        title: str | None = None,
    ) -> PlotResult:
        """Plot sc.pl.rank_genes_groups — raises if results are not yet computed."""
        plt.close("all")

        # Check that rank_genes_groups results exist — never auto-compute here.
        rgg = adata.uns.get(key)
        if not isinstance(rgg, dict) or "names" not in rgg:
            computed_keys = [k for k, v in adata.uns.items() if isinstance(v, dict) and "names" in v]
            hint = (
                f"Found related keys: {computed_keys}. Use the 'key' parameter to target one of them."
                if computed_keys
                else "No rank_genes_groups results found in adata.uns at all."
            )
            raise ValueError(
                f"rank_genes_groups results not found under adata.uns['{key}']. "
                f"{hint} "
                f"Please run sc.tl.rank_genes_groups(adata, groupby=...) first."
            )

        # Determine number of groups to size the figure
        try:
            group_names = list(rgg["names"].dtype.names or [])
        except Exception:
            group_names = []
        n_groups = len(groups) if groups else (len(group_names) if group_names else ncols)
        n_rows = max(1, math.ceil(n_groups / ncols))
        fig_width = min(ncols * 4.5, 28.0)
        fig_height = max(2.5, min(n_rows * (2.0 + n_genes * 0.12), 24.0))

        sc_kwargs: dict[str, object] = {
            "n_genes": n_genes,
            "key": key,
            "fontsize": fontsize,
            "ncols": ncols,
            "sharey": sharey,
            "show": False,
        }
        if groups:
            sc_kwargs["groups"] = groups
        sym_col = self._find_gene_symbol_column(adata)
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols
        elif sym_col:
            sc_kwargs["gene_symbols"] = sym_col

        with rc_context({
            "figure.figsize": (fig_width, fig_height),
            "axes.labelsize": 13,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
        }):
            sc.pl.rank_genes_groups(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))
        fig.tight_layout()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"rank_genes_groups_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        # Recover groupby used for this result
        params = rgg.get("params") if isinstance(rgg.get("params"), dict) else {}
        groupby_used = str(params.get("groupby", "")).strip() or None

        return PlotResult(
            plot_type="rank_genes_groups",
            output_file=out_file,
            resolved_groupby=groupby_used,
            display_plot_type="Rank genes groups",
        )

    def run_rank_genes_groups_violin(
        self,
        adata: ad.AnnData,
        *,
        groups: list[str] | None = None,
        n_genes: int = 5,
        gene_names: list[str] | None = None,
        gene_symbols: str | None = None,
        use_raw: bool | None = None,
        key: str = "rank_genes_groups",
        split: bool = True,
        density_norm: str = "width",
        strip: bool = True,
        jitter: bool | float = True,
        size: int = 1,
        title: str | None = None,
    ) -> PlotResult:
        """Plot rank_genes_groups_violin — raises if results are not yet computed."""
        plt.close("all")

        rgg = adata.uns.get(key)
        if not isinstance(rgg, dict) or "names" not in rgg:
            computed_keys = [k for k, v in adata.uns.items() if isinstance(v, dict) and "names" in v]
            hint = (
                f"Found related keys: {computed_keys}. Use the 'key' parameter to target one of them."
                if computed_keys
                else "No rank_genes_groups results found in adata.uns at all."
            )
            raise ValueError(
                f"rank_genes_groups results not found under adata.uns['{key}']. "
                f"{hint} "
                f"Please run sc.tl.rank_genes_groups(adata, groupby=...) first."
            )

        # Auto-size: one panel per group requested (or all groups)
        try:
            all_group_names = list(rgg["names"].dtype.names or [])
        except Exception:
            all_group_names = []
        display_groups = groups if groups else all_group_names
        n_panels = max(1, len(display_groups))
        fig_width = min(max(8.0, n_panels * 4.0), 32.0)
        fig_height = max(4.0, min(2.5 + n_genes * 0.18, 24.0))

        sc_kwargs: dict[str, object] = {
            "n_genes": n_genes,
            "key": key,
            "split": split,
            "density_norm": density_norm,
            "strip": strip,
            "jitter": jitter,
            "size": size,
            "show": False,
        }
        if groups:
            sc_kwargs["groups"] = groups
        if gene_names:
            sc_kwargs["gene_names"] = gene_names
        sym_col = self._find_gene_symbol_column(adata)
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols
        elif sym_col:
            sc_kwargs["gene_symbols"] = sym_col
        sc_kwargs["use_raw"] = use_raw if use_raw is not None else False

        with rc_context({"figure.figsize": (fig_width, fig_height)}):
            sc.pl.rank_genes_groups_violin(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()
        for ax in fig.axes:
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))
        fig.tight_layout()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"rank_genes_groups_violin_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        params = rgg.get("params") if isinstance(rgg.get("params"), dict) else {}
        groupby_used = str(params.get("groupby", "")).strip() or None

        return PlotResult(
            plot_type="rank_genes_groups_violin",
            output_file=out_file,
            resolved_groupby=groupby_used,
            display_plot_type="Rank genes groups violin",
        )

    def run_rank_genes_groups_stacked_violin(
        self,
        adata: ad.AnnData,
        *,
        groups: list[str] | None = None,
        n_genes: int | None = None,
        groupby: str | None = None,
        gene_symbols: str | None = None,
        var_names: list[str] | None = None,
        min_logfoldchange: float | None = None,
        key: str = "rank_genes_groups",
        swap_axes: bool = False,
        cmap: str = "Blues",
        stripplot: bool = False,
        jitter: bool | float = False,
        size: int = 1,
        row_palette: str | None = None,
        yticklabels: bool = False,
        standard_scale: str | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        colorbar_title: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Plot rank_genes_groups_stacked_violin — raises if results are not yet computed."""
        plt.close("all")

        rgg = adata.uns.get(key)
        if not isinstance(rgg, dict) or "names" not in rgg:
            computed_keys = [k for k, v in adata.uns.items() if isinstance(v, dict) and "names" in v]
            hint = (
                f"Found related keys: {computed_keys}. Use the 'key' parameter to target one of them."
                if computed_keys
                else "No rank_genes_groups results found in adata.uns at all."
            )
            raise ValueError(
                f"rank_genes_groups results not found under adata.uns['{key}']. "
                f"{hint} "
                f"Please run sc.tl.rank_genes_groups(adata, groupby=...) first."
            )

        params = rgg.get("params") if isinstance(rgg.get("params"), dict) else {}
        groupby_used = str(params.get("groupby", "")).strip() or None
        resolved_groupby = groupby or groupby_used

        # Auto-size: x-axis has n_genes * n_groups columns (one gene block per group)
        if figsize is None:
            try:
                all_group_names = list(rgg["names"].dtype.names or [])
            except Exception:
                all_group_names = []
            display_groups = groups if groups else all_group_names
            n_g = len(display_groups) if display_groups else 8
            n_genes_display = abs(n_genes) if n_genes else 5
            if swap_axes:
                width = max(6.0, n_g * 0.5 + 2)
                height = max(4.0, n_genes_display * n_g * 0.3 + 2)
            else:
                width = max(8.0, n_genes_display * n_g * 0.25 + 2)
                height = max(4.0, n_g * 0.45 + 2)
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "swap_axes": swap_axes,
            "cmap": cmap,
            "stripplot": stripplot,
            "jitter": jitter,
            "size": size,
            "yticklabels": yticklabels,
        }
        if key:
            sc_kwargs["key"] = key
        if groups:
            sc_kwargs["groups"] = groups
        if n_genes is not None:
            sc_kwargs["n_genes"] = n_genes
        if resolved_groupby:
            sc_kwargs["groupby"] = resolved_groupby
        if min_logfoldchange is not None:
            sc_kwargs["min_logfoldchange"] = min_logfoldchange
        if var_names:
            sc_kwargs["var_names"] = var_names
        use_gene_symbols = self._has_feature_name_column(adata)
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if row_palette:
            sc_kwargs["row_palette"] = row_palette
        if colorbar_title:
            sc_kwargs["colorbar_title"] = colorbar_title
        if figsize:
            sc_kwargs["figsize"] = figsize

        sc.pl.rank_genes_groups_stacked_violin(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()
        if swap_axes:
            for ax in fig.axes:
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"rank_genes_groups_stacked_violin_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="rank_genes_groups_stacked_violin",
            output_file=out_file,
            resolved_groupby=groupby_used,
            display_plot_type="Rank-genes stacked violin",
        )

    def run_rank_genes_groups_heatmap(
        self,
        adata: ad.AnnData,
        *,
        groups: list[str] | None = None,
        n_genes: int | None = None,
        groupby: str | None = None,
        gene_symbols: str | None = None,
        var_names: list[str] | None = None,
        min_logfoldchange: float | None = None,
        key: str = "rank_genes_groups",
        standard_scale: str | None = None,
        swap_axes: bool = False,
        show_gene_labels: bool | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        cmap: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Plot rank_genes_groups_heatmap — raises if results are not yet computed."""
        plt.close("all")

        rgg = adata.uns.get(key)
        if not isinstance(rgg, dict) or "names" not in rgg:
            computed_keys = [k for k, v in adata.uns.items() if isinstance(v, dict) and "names" in v]
            hint = (
                f"Found related keys: {computed_keys}. Use the 'key' parameter to target one of them."
                if computed_keys
                else "No rank_genes_groups results found in adata.uns at all."
            )
            raise ValueError(
                f"rank_genes_groups results not found under adata.uns['{key}']. "
                f"{hint} "
                f"Please run sc.tl.rank_genes_groups(adata, groupby=...) first."
            )

        params = rgg.get("params") if isinstance(rgg.get("params"), dict) else {}
        groupby_used = str(params.get("groupby", "")).strip() or None
        resolved_groupby = groupby or groupby_used

        # Auto-size: heatmap — genes × groups rows (or cols when swapped)
        if figsize is None:
            try:
                all_group_names = list(rgg["names"].dtype.names or [])
            except Exception:
                all_group_names = []
            display_groups = groups if groups else all_group_names
            n_g = len(display_groups) if display_groups else 8
            n_genes_display = abs(n_genes) if n_genes else 5
            if swap_axes:
                # genes on x-axis: total cols = n_genes * n_groups
                width = max(8.0, n_genes_display * n_g * 0.25 + 2)
                height = max(4.0, n_g * 0.45 + 2)
            else:
                # genes on y-axis, groups on x-axis
                width = max(6.0, n_g * 0.5 + 2)
                height = max(4.0, n_genes_display * n_g * 0.2 + 2)
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "swap_axes": swap_axes,
        }
        if key:
            sc_kwargs["key"] = key
        if groups:
            sc_kwargs["groups"] = groups
        if n_genes is not None:
            sc_kwargs["n_genes"] = n_genes
        if resolved_groupby:
            sc_kwargs["groupby"] = resolved_groupby
        if min_logfoldchange is not None:
            sc_kwargs["min_logfoldchange"] = min_logfoldchange
        if var_names:
            sc_kwargs["var_names"] = var_names
        use_gene_symbols = self._has_feature_name_column(adata)
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if show_gene_labels is not None:
            sc_kwargs["show_gene_labels"] = show_gene_labels
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if cmap:
            sc_kwargs["cmap"] = cmap
        if figsize:
            sc_kwargs["figsize"] = figsize

        with rc_context({"figure.figsize": figsize}):
            sc.pl.rank_genes_groups_heatmap(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()
        if not swap_axes:
            for ax in fig.axes:
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"rank_genes_groups_heatmap_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="rank_genes_groups_heatmap",
            output_file=out_file,
            resolved_groupby=groupby_used,
            display_plot_type="Rank-genes heatmap",
        )

    def run_rank_genes_groups_dotplot(
        self,
        adata: ad.AnnData,
        *,
        groups: list[str] | None = None,
        n_genes: int | None = None,
        groupby: str | None = None,
        gene_symbols: str | None = None,
        var_names: list[str] | None = None,
        min_logfoldchange: float | None = None,
        key: str = "rank_genes_groups",
        values_to_plot: str | None = None,
        standard_scale: str | None = None,
        dendrogram: bool = False,
        swap_axes: bool = False,
        cmap: str | None = None,
        dot_max: float | None = None,
        dot_min: float | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Plot rank_genes_groups_dotplot — raises if results are not yet computed."""
        plt.close("all")

        rgg = adata.uns.get(key)
        if not isinstance(rgg, dict) or "names" not in rgg:
            computed_keys = [k for k, v in adata.uns.items() if isinstance(v, dict) and "names" in v]
            hint = (
                f"Found related keys: {computed_keys}. Use the 'key' parameter to target one of them."
                if computed_keys
                else "No rank_genes_groups results found in adata.uns at all."
            )
            raise ValueError(
                f"rank_genes_groups results not found under adata.uns['{key}']. "
                f"{hint} "
                f"Please run sc.tl.rank_genes_groups(adata, groupby=...) first."
            )

        params = rgg.get("params") if isinstance(rgg.get("params"), dict) else {}
        groupby_used = str(params.get("groupby", "")).strip() or None
        resolved_groupby = groupby or groupby_used

        # Auto-size: x-axis has n_genes * n_groups columns (one gene block per group)
        if figsize is None:
            try:
                all_group_names = list(rgg["names"].dtype.names or [])
            except Exception:
                all_group_names = []
            display_groups = groups if groups else all_group_names
            n_g = len(display_groups) if display_groups else 8
            n_genes_display = abs(n_genes) if n_genes else 5
            if swap_axes:
                width = max(6.0, n_g * 0.5 + 2)
                height = max(4.0, n_genes_display * n_g * 0.3 + 2)
            else:
                width = max(8.0, n_genes_display * n_g * 0.25 + 2)
                height = max(4.0, n_g * 0.5 + 2)
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "swap_axes": swap_axes,
            "dendrogram": dendrogram,
        }
        if key:
            sc_kwargs["key"] = key
        if groups:
            sc_kwargs["groups"] = groups
        if n_genes is not None:
            sc_kwargs["n_genes"] = n_genes
        if resolved_groupby:
            sc_kwargs["groupby"] = resolved_groupby
        if min_logfoldchange is not None:
            sc_kwargs["min_logfoldchange"] = min_logfoldchange
        if var_names:
            sc_kwargs["var_names"] = var_names
        use_gene_symbols = self._has_feature_name_column(adata)
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols
        if values_to_plot:
            sc_kwargs["values_to_plot"] = values_to_plot
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if cmap:
            sc_kwargs["cmap"] = cmap
        if dot_max is not None:
            sc_kwargs["dot_max"] = dot_max
        if dot_min is not None:
            sc_kwargs["dot_min"] = dot_min
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if figsize:
            sc_kwargs["figsize"] = figsize

        with rc_context({"figure.figsize": figsize}):
            sc.pl.rank_genes_groups_dotplot(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"rank_genes_groups_dotplot_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="rank_genes_groups_dotplot",
            output_file=out_file,
            resolved_groupby=groupby_used,
            display_plot_type="Rank-genes dot plot",
        )

    def run_rank_genes_groups_matrixplot(
        self,
        adata: ad.AnnData,
        *,
        groups: list[str] | None = None,
        n_genes: int | None = None,
        groupby: str | None = None,
        gene_symbols: str | None = None,
        var_names: list[str] | None = None,
        min_logfoldchange: float | None = None,
        key: str = "rank_genes_groups",
        values_to_plot: str | None = None,
        standard_scale: str | None = None,
        dendrogram: bool = False,
        swap_axes: bool = False,
        cmap: str | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        colorbar_title: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Plot rank_genes_groups_matrixplot — raises if results are not yet computed."""
        plt.close("all")

        rgg = adata.uns.get(key)
        if not isinstance(rgg, dict) or "names" not in rgg:
            computed_keys = [k for k, v in adata.uns.items() if isinstance(v, dict) and "names" in v]
            hint = (
                f"Found related keys: {computed_keys}. Use the 'key' parameter to target one of them."
                if computed_keys
                else "No rank_genes_groups results found in adata.uns at all."
            )
            raise ValueError(
                f"rank_genes_groups results not found under adata.uns['{key}']. "
                f"{hint} "
                f"Please run sc.tl.rank_genes_groups(adata, groupby=...) first."
            )

        params = rgg.get("params") if isinstance(rgg.get("params"), dict) else {}
        groupby_used = str(params.get("groupby", "")).strip() or None
        resolved_groupby = groupby or groupby_used

        # Auto-size: x-axis has n_genes * n_groups columns (one gene block per group)
        if figsize is None:
            try:
                all_group_names = list(rgg["names"].dtype.names or [])
            except Exception:
                all_group_names = []
            display_groups = groups if groups else all_group_names
            n_g = len(display_groups) if display_groups else 8
            n_genes_display = abs(n_genes) if n_genes else 5
            if swap_axes:
                width = max(6.0, n_g * 0.5 + 2)
                height = max(4.0, n_genes_display * n_g * 0.3 + 2)
            else:
                width = max(8.0, n_genes_display * n_g * 0.25 + 2)
                height = max(4.0, n_g * 0.5 + 2)
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "swap_axes": swap_axes,
            "dendrogram": dendrogram,
        }
        if key:
            sc_kwargs["key"] = key
        if groups:
            sc_kwargs["groups"] = groups
        if n_genes is not None:
            sc_kwargs["n_genes"] = n_genes
        if resolved_groupby:
            sc_kwargs["groupby"] = resolved_groupby
        if min_logfoldchange is not None:
            sc_kwargs["min_logfoldchange"] = min_logfoldchange
        if var_names:
            sc_kwargs["var_names"] = var_names
        use_gene_symbols = self._has_feature_name_column(adata)
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols
        if values_to_plot:
            sc_kwargs["values_to_plot"] = values_to_plot
        if standard_scale:
            sc_kwargs["standard_scale"] = standard_scale
        if cmap:
            sc_kwargs["cmap"] = cmap
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if colorbar_title:
            sc_kwargs["colorbar_title"] = colorbar_title
        if figsize:
            sc_kwargs["figsize"] = figsize

        with rc_context({"figure.figsize": figsize}):
            sc.pl.rank_genes_groups_matrixplot(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"rank_genes_groups_matrixplot_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="rank_genes_groups_matrixplot",
            output_file=out_file,
            resolved_groupby=groupby_used,
            display_plot_type="Rank-genes matrix plot",
        )

    def run_rank_genes_groups_tracksplot(
        self,
        adata: ad.AnnData,
        *,
        groups: list[str] | None = None,
        n_genes: int | None = None,
        groupby: str | None = None,
        gene_symbols: str | None = None,
        min_logfoldchange: float | None = None,
        key: str = "rank_genes_groups",
        dendrogram: bool = False,
        use_raw: bool | None = None,
        log: bool = False,
        layer: str | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Plot rank_genes_groups_tracksplot — raises if results are not yet computed."""
        plt.close("all")

        rgg = adata.uns.get(key)
        if not isinstance(rgg, dict) or "names" not in rgg:
            computed_keys = [k for k, v in adata.uns.items() if isinstance(v, dict) and "names" in v]
            hint = (
                f"Found related keys: {computed_keys}. Use the 'key' parameter to target one of them."
                if computed_keys
                else "No rank_genes_groups results found in adata.uns at all."
            )
            raise ValueError(
                f"rank_genes_groups results not found under adata.uns['{key}']. "
                f"{hint} "
                f"Please run sc.tl.rank_genes_groups(adata, groupby=...) first."
            )

        params = rgg.get("params") if isinstance(rgg.get("params"), dict) else {}
        groupby_used = str(params.get("groupby", "")).strip() or None
        resolved_groupby = groupby or groupby_used

        # Auto-size: tracksplot — x-axis has n_genes * n_groups columns, height scales with groups
        if figsize is None:
            try:
                all_group_names = list(rgg["names"].dtype.names or [])
            except Exception:
                all_group_names = []
            display_groups = groups if groups else all_group_names
            n_g = len(display_groups) if display_groups else 8
            n_genes_display = abs(n_genes) if n_genes else 5
            width = max(10.0, n_genes_display * n_g * 0.25 + 2)
            height = max(4.0, n_g * 0.8 + 2)
            figsize = (width, height)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "dendrogram": dendrogram,
            "log": log,
        }
        if key:
            sc_kwargs["key"] = key
        if groups:
            sc_kwargs["groups"] = groups
        if n_genes is not None:
            sc_kwargs["n_genes"] = n_genes
        if resolved_groupby:
            sc_kwargs["groupby"] = resolved_groupby
        if min_logfoldchange is not None:
            sc_kwargs["min_logfoldchange"] = min_logfoldchange
        sym_col = self._find_gene_symbol_column(adata)
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols
        elif sym_col:
            sc_kwargs["gene_symbols"] = sym_col
        sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if figsize:
            sc_kwargs["figsize"] = figsize

        with rc_context({"figure.figsize": figsize}):
            sc.pl.rank_genes_groups_tracksplot(adata, **sc_kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

        fig = plt.gcf()
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"rank_genes_groups_tracksplot_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="rank_genes_groups_tracksplot",
            output_file=out_file,
            resolved_groupby=groupby_used,
            display_plot_type="Rank-genes tracks plot",
        )

    def run_correlation_matrix(
        self,
        adata: ad.AnnData,
        *,
        groupby: str | None = None,
        show_correlation_numbers: bool = False,
        dendrogram: bool | None = None,
        cmap: str | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        vcenter: float | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Plot correlation_matrix — requires groupby."""
        plt.close("all")

        resolved_groupby = self._resolve_cell_type_column(adata, groupby)
        if not resolved_groupby:
            raise ValueError("correlation_matrix requires a groupby column and no inferred group column was found.")

        # Auto-size: square matrix — size scales with number of categories
        if figsize is None:
            try:
                n_cats = adata.obs[resolved_groupby].nunique()
            except Exception:
                n_cats = 8
            side = max(6.0, min(n_cats * 0.6 + 2, 20.0))
            figsize = (side, side)

        sc_kwargs: dict[str, object] = {
            "show": False,
            "show_correlation_numbers": show_correlation_numbers,
        }
        if dendrogram is not None:
            sc_kwargs["dendrogram"] = dendrogram
        if cmap:
            sc_kwargs["cmap"] = cmap
        if vmin is not None:
            sc_kwargs["vmin"] = vmin
        if vmax is not None:
            sc_kwargs["vmax"] = vmax
        if vcenter is not None:
            sc_kwargs["vcenter"] = vcenter
        if figsize:
            sc_kwargs["figsize"] = figsize

        with rc_context({"figure.figsize": figsize}):
            sc.pl.correlation_matrix(adata, resolved_groupby, **sc_kwargs)

        fig = plt.gcf()
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"correlation_matrix_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="correlation_matrix",
            output_file=out_file,
            resolved_groupby=resolved_groupby,
            display_plot_type="Correlation matrix",
        )

    def run_highest_expr_genes(
        self,
        adata: ad.AnnData,
        *,
        n_top: int = 30,
        layer: str | None = None,
        gene_symbols: str | None = None,
        log: bool = False,
        title: str | None = None,
    ) -> PlotResult:
        """Plot sc.pl.highest_expr_genes — top genes by mean fraction of counts per cell."""
        plt.close("all")

        # Auto-detect gene_symbols column if not supplied
        if not gene_symbols:
            gene_symbols = self._find_gene_symbol_column(adata)

        # Height scales with n_top so all gene labels are readable
        fig_height = max(5.0, n_top * 0.28 + 1.5)
        fig_width = 8.0

        sc_kwargs: dict[str, object] = {
            "n_top": n_top,
            "log": log,
            "show": False,
        }
        if layer:
            sc_kwargs["layer"] = layer
        if gene_symbols:
            sc_kwargs["gene_symbols"] = gene_symbols

        with rc_context({"figure.figsize": (fig_width, fig_height)}):
            sc.pl.highest_expr_genes(adata, **sc_kwargs)

        fig = plt.gcf()
        if title:
            if fig.axes:
                fig.axes[0].set_title(str(title))
        fig.tight_layout()

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"highest_expr_genes_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="highest_expr_genes",
            output_file=out_file,
            display_plot_type="Highest expressed genes",
        )

    def run_scatter(
        self,
        adata: ad.AnnData,
        *,
        x: str | None = None,
        y: str | None = None,
        color: list[str] | None = None,
        basis: str | None = None,
        use_raw: bool | None = None,
        layer: str | None = None,
        groups: list[str] | None = None,
        components: str | None = None,
        sort_order: bool = True,
        legend_loc: str = "right margin",
        size: float | None = None,
        color_map: str | None = None,
        frameon: bool | None = None,
        title: str | None = None,
        figsize_width: float | None = None,
        figsize_height: float | None = None,
    ) -> PlotResult:
        """Plot sc.pl.scatter — flexible scatter plot over obs columns or embeddings."""
        plt.close("all")

        n_panels = max(1, len(color)) if color else 1
        fw = figsize_width if figsize_width else max(5.0, n_panels * 5.0)
        fh = figsize_height if figsize_height else 4.0

        sc_kwargs: dict[str, object] = {
            "sort_order": sort_order,
            "legend_loc": legend_loc,
            "show": False,
        }
        if x:
            sc_kwargs["x"] = x
        if y:
            sc_kwargs["y"] = y
        if color:
            sc_kwargs["color"] = color
        if basis:
            sc_kwargs["basis"] = basis
        sc_kwargs["use_raw"] = use_raw if use_raw is not None else False
        if layer:
            sc_kwargs["layer"] = layer
        if groups:
            sc_kwargs["groups"] = groups
        if components:
            sc_kwargs["components"] = components
        if size is not None:
            sc_kwargs["size"] = size
        if color_map:
            sc_kwargs["color_map"] = color_map
        if frameon is not None:
            sc_kwargs["frameon"] = frameon
        if title:
            sc_kwargs["title"] = title

        with rc_context({"figure.figsize": (fw, fh)}):
            sc.pl.scatter(adata, **sc_kwargs)

        fig = plt.gcf()
        fig.tight_layout()

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"scatter_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="scatter",
            output_file=out_file,
            display_plot_type="Scatter plot",
        )

    def run_spatial_scatter(
        self,
        adata: ad.AnnData,
        *,
        color: list[str] | None = None,
        figsize_width: float | None = None,
        figsize_height: float | None = None,
        title: str | None = None,
        # Layout
        wspace: float | None = None,
        hspace: float | None = None,
        ncols: int | None = None,
        # Shape / size
        shape: str | None = None,
        size: float | None = None,
        alpha: float | None = None,
        # Color / style
        cmap: str | None = None,
        palette: str | None = None,
        na_color: str | None = None,
        # Image
        img: bool | None = None,
        img_alpha: float | None = None,
        # Crop
        crop_coord: tuple[int, int, int, int] | None = None,
        # Legend / colorbar
        legend_loc: str | None = None,
        legend_fontsize: float | None = None,
        colorbar: bool | None = None,
        frameon: bool | None = None,
        # Outline
        outline: bool | None = None,
        # Groups
        groups: list[str] | None = None,
        # Layer
        layer: str | None = None,
        use_raw: bool | None = None,
    ) -> PlotResult:
        """Spatial scatter plot using squidpy sq.pl.spatial_scatter."""
        plt.close("all")

        if "spatial" not in adata.uns or "spatial" not in adata.obsm:
            raise ValueError(
                "Dataset is not a spatial transcriptomics dataset. "
                "Requires adata.uns['spatial'] and adata.obsm['spatial']."
            )

        # Remove is_single key if present (squidpy doesn't expect it)
        if "is_single" in adata.uns["spatial"]:
            adata.uns["spatial"].pop("is_single")

        # Resolve color columns — default to cell type column
        if not color:
            resolved_ct = self._resolve_cell_type_column(adata)
            color = [resolved_ct] if resolved_ct else []
        else:
            obs_lower = {str(c).lower(): str(c) for c in adata.obs.columns}
            var_names_set = set(adata.var_names)
            resolved_color = []
            # Separate tokens into obs columns vs potential gene names
            obs_tokens: list[str] = []
            gene_tokens: list[str] = []
            for c in color:
                c_str = str(c).strip()
                if c_str in adata.obs.columns:
                    obs_tokens.append(c_str)
                elif c_str.lower() in obs_lower:
                    obs_tokens.append(obs_lower[c_str.lower()])
                else:
                    gene_tokens.append(c_str)
            # Resolve gene tokens via the same ENS ID / gene symbol logic used by dotplot etc.
            resolved_genes = self._resolve_gene_names(adata, gene_tokens) if gene_tokens else []
            # Filter to genes actually present in var_names (resolved_genes may contain ENS IDs)
            valid_genes = [g for g in resolved_genes if g in var_names_set]
            # Any original token that produced no entry in valid_genes is truly unresolvable;
            # pass it through so squidpy raises a meaningful error (don't re-add symbols that
            # were already resolved to ENS IDs — that's what caused the original failure).
            n_resolved = len(valid_genes)
            n_expected = len(gene_tokens)
            truly_unresolved: list[str] = []
            if n_resolved < n_expected:
                # Re-resolve one-by-one to find which original tokens had no match
                for tok in gene_tokens:
                    r = self._resolve_gene_names(adata, [tok])
                    if not r or r[0] not in var_names_set:
                        truly_unresolved.append(tok)
            resolved_color = obs_tokens + valid_genes + truly_unresolved
            color = resolved_color

        # Figure size: scale width by number of panels
        n_panels = max(1, len(color))
        fw = figsize_width if figsize_width else max(10.0, n_panels * 8.0)
        fh = figsize_height if figsize_height else 8.0
        figsize = (fw, fh)

        # Build per-panel titles: always one title per color panel (gene symbol or obs column name).
        # A combined string title from the agent (e.g. "Spatial scatter by EPCAM and MUC1") is
        # only used when there is exactly one panel — otherwise always use individual panel titles.
        if color:
            label_map = self._resolve_feature_labels_for_var_keys(adata, color)
            panel_titles: list[str] = [label_map.get(c, c) for c in color]
            if title and len(panel_titles) == 1:
                panel_titles = [title]
        else:
            panel_titles = [title] if title else []

        sq_kwargs: dict[str, object] = {
            "color": color if color else None,
            "figsize": figsize,
        }
        if panel_titles:
            sq_kwargs["title"] = panel_titles
        if wspace is not None:
            sq_kwargs["wspace"] = wspace
        if hspace is not None:
            sq_kwargs["hspace"] = hspace
        if ncols is not None:
            sq_kwargs["ncols"] = ncols
        if alpha is not None:
            sq_kwargs["alpha"] = alpha
        if cmap is not None:
            sq_kwargs["cmap"] = cmap
        if palette is not None:
            sq_kwargs["palette"] = palette
        if na_color is not None:
            sq_kwargs["na_color"] = na_color
        if img_alpha is not None:
            sq_kwargs["img_alpha"] = img_alpha
        if crop_coord is not None:
            sq_kwargs["crop_coord"] = crop_coord
        if legend_loc is not None:
            sq_kwargs["legend_loc"] = legend_loc
        if legend_fontsize is not None:
            sq_kwargs["legend_fontsize"] = legend_fontsize
        if colorbar is not None:
            sq_kwargs["colorbar"] = colorbar
        if frameon is not None:
            sq_kwargs["frameon"] = frameon
        if outline is not None:
            sq_kwargs["outline"] = outline
        if groups is not None:
            sq_kwargs["groups"] = groups
        if layer is not None:
            sq_kwargs["layer"] = layer
        if use_raw is not None:
            sq_kwargs["use_raw"] = use_raw

        # Determine if spatial uns has library entries (image + coordinates info)
        has_spatial_info = len(adata.uns["spatial"]) > 0

        if has_spatial_info:
            if shape is not None:
                sq_kwargs["shape"] = shape
            if size is not None:
                sq_kwargs["size"] = size
            if img is not None:
                sq_kwargs["img"] = img
            sq.pl.spatial_scatter(adata, **sq_kwargs)
        else:
            sq_kwargs["shape"] = shape if shape is not None else None
            sq_kwargs["size"] = size if size is not None else 3
            sq_kwargs["img"] = img if img is not None else False
            sq.pl.spatial_scatter(adata, **sq_kwargs)

        fig = plt.gcf()

        # Resize each colorbar axes to match the height of its paired image axes.
        # Squidpy places colorbars as narrow axes that are shorter than the plot panel;
        # detect them by their narrow width relative to height (width < 10% of height).
        def _is_cbar(ax: plt.Axes) -> bool:
            pos = ax.get_position()
            return pos.width > 0 and (pos.width / max(pos.height, 1e-6)) < 0.15

        all_axes = fig.get_axes()
        image_axes = [ax for ax in all_axes if not _is_cbar(ax)]
        colorbar_axes = [ax for ax in all_axes if _is_cbar(ax)]
        for cax in colorbar_axes:
            # Find the image axes closest horizontally to this colorbar
            if image_axes:
                paired = min(image_axes, key=lambda ax: abs(ax.get_position().x1 - cax.get_position().x0))
                ipos = paired.get_position()
                cpos = cax.get_position()
                cax.set_position([cpos.x0, ipos.y0, cpos.width, ipos.height])

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"spatial_scatter_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="spatial_scatter",
            output_file=out_file,
            color_columns=color or None,
            display_plot_type="Spatial scatter plot",
        )

    def run_nhood_enrichment(
        self,
        adata: ad.AnnData,
        *,
        cluster_key: str | None = None,
        mode: str = "zscore",
        figsize_width: float | None = None,
        figsize_height: float | None = None,
        title: str | None = None,
    ) -> PlotResult:
        """Neighborhood enrichment analysis and heatmap using squidpy."""
        plt.close("all")

        if "spatial" not in adata.uns or "spatial" not in adata.obsm:
            raise ValueError(
                "Dataset is not a spatial transcriptomics dataset. "
                "Requires adata.uns['spatial'] and adata.obsm['spatial']."
            )

        # Remove is_single key if present
        if "is_single" in adata.uns["spatial"]:
            adata.uns["spatial"].pop("is_single")

        # Resolve cluster_key — default to cell type column
        resolved_cluster_key = self._resolve_cell_type_column(adata, cluster_key)
        if not resolved_cluster_key:
            raise ValueError(
                "Could not resolve a cluster key. Pass an explicit adata.obs column via cluster_key."
            )

        # Figure size: nhood_enrichment is a square matrix; scale with category count
        n_cats = int(adata.obs[resolved_cluster_key].nunique(dropna=True)) if resolved_cluster_key in adata.obs.columns else 1
        fw = figsize_width if figsize_width else max(10.0, n_cats * 0.6 + 4.0)
        fh = figsize_height if figsize_height else max(8.0, n_cats * 0.6 + 2.0)
        figsize = (fw, fh)

        # Fill NA values in the cluster column before spatial graph computation —
        # squidpy does not handle NaN categories and will error or produce wrong results.
        col = adata.obs[resolved_cluster_key]
        if hasattr(col, "cat"):
            if "Unknown" not in col.cat.categories:
                adata.obs[resolved_cluster_key] = col.cat.add_categories("Unknown")
                # Drop the pre-existing color palette — it has one fewer entry than the new
                # category list, which causes matplotlib BoundaryNorm to raise
                # "ncolors must equal or exceed the number of bins".
                colors_key = f"{resolved_cluster_key}_colors"
                adata.uns.pop(colors_key, None)
            adata.obs[resolved_cluster_key] = adata.obs[resolved_cluster_key].fillna("Unknown")
        else:
            adata.obs[resolved_cluster_key] = col.fillna("Unknown")

        # Compute spatial graph and neighborhood enrichment
        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, cluster_key=resolved_cluster_key)

        # Use a diverging colormap with enough discrete levels to cover all category pairs.
        # The heatmap has n_cats rows/cols; matplotlib BoundaryNorm needs ncolors >= n_bins.
        import matplotlib.colors as mcolors
        n_levels = max(n_cats + 2, 256)
        enrichment_cmap = plt.get_cmap("RdBu_r", n_levels)

        sq_kwargs: dict[str, object] = {
            "cluster_key": resolved_cluster_key,
            "mode": mode,
            "figsize": figsize,
            "show": False,
            "cmap": enrichment_cmap,
        }
        if title:
            sq_kwargs["title"] = title

        sq.pl.nhood_enrichment(adata, **sq_kwargs)

        fig = plt.gcf()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = self.output_dir / f"nhood_enrichment_{timestamp}.png"
        fig.savefig(out_file, dpi=150, bbox_inches="tight")
        plt.close("all")

        return PlotResult(
            plot_type="nhood_enrichment",
            output_file=out_file,
            resolved_groupby=resolved_cluster_key,
            display_plot_type="Neighborhood enrichment",
        )

    def _derive_plot_metadata(
        self,
        *,
        adata: ad.AnnData,
        request: PlotRequest,
        requested_plot_type: str,
        actual_plot_type: str,
        embedding_basis: str | None,
        color_columns: list[str] | None,
    ) -> dict[str, object]:
        resolved_genes: list[str] | None = None
        resolved_groupby: str | None = None
        resolved_coloring_label: str | None = None
        resolved_color_columns = color_columns or None

        if actual_plot_type in {"umap", "tsne", "gene_cell_embedding"} or embedding_basis:
            requested_colors = list(color_columns or [])
            if not requested_colors:
                requested_colors = list(request.color or request.genes or [])
            resolved_color_columns, _ = self._resolve_embedding_color(
                adata,
                requested_colors,
                preferred_gene_symbols_column=request.gene_symbols_column,
            )
            if not resolved_color_columns:
                inferred_color = self._resolve_cell_type_column(adata, request.groupby)
                if inferred_color:
                    resolved_color_columns = [inferred_color]
            resolved_genes = self._resolve_embedding_genes(adata, request, resolved_color_columns or [])
            resolved_groupby = self._resolve_cell_type_column(adata, request.groupby)
            resolved_coloring_label = self._build_coloring_label(adata, resolved_color_columns or [], resolved_groupby)
        elif actual_plot_type in {"violin", "dotplot", "matrixplot", "heatmap"}:
            resolved_genes = self._resolve_distribution_genes(adata, request)
            resolved_groupby = self._resolve_distribution_groupby(adata, request)
        elif actual_plot_type.startswith("rank_genes_groups_"):
            resolved_groupby = self._rank_genes_groups_used_groupby or self._resolve_cell_type_column(adata, request.groupby)
        elif actual_plot_type == "correlation_matrix":
            resolved_groupby = self._resolve_cell_type_column(adata, request.groupby)
        elif actual_plot_type == "cell_counts_barplot":
            resolved_groupby = self._resolve_cell_type_column(adata, request.groupby)
        elif actual_plot_type == "cell_type_proportion_barplot":
            resolved_groupby = self._resolve_cell_type_column(adata, request.groupby)

        return {
            "color_columns": resolved_color_columns,
            "resolved_genes": resolved_genes,
            "resolved_groupby": resolved_groupby,
            "resolved_coloring_label": resolved_coloring_label,
            "display_plot_type": self._display_plot_type(actual_plot_type),
        }

    @staticmethod
    def _display_plot_type(plot_type: str) -> str:
        display_map = {
            "umap": "UMAP embedding",
            "tsne": "tSNE embedding",
            "gene_cell_embedding": "Gene/cell embedding",
            "violin": "Violin plot",
            "dotplot": "Dot plot",
            "matrixplot": "Matrix plot",
            "heatmap": "Heatmap",
            "rank_genes_groups_dotplot": "Rank-genes dot plot",
            "rank_genes_groups_matrixplot": "Rank-genes matrix plot",
            "rank_genes_groups_stacked_violin": "Rank-genes stacked violin",
            "rank_genes_groups_heatmap": "Rank-genes heatmap",
            "correlation_matrix": "Correlation matrix",
            "cell_counts_barplot": "Cell counts bar plot",
            "cell_type_proportion_barplot": "Cell type proportion bar plot",
            "tracksplot": "Tracks plot",
            "stacked_violin": "Stacked violin plot",
            "clustermap": "Cluster map",
            "dendrogram": "Dendrogram",
            "rank_genes_groups": "Rank genes groups",
            "rank_genes_groups_violin": "Rank genes groups violin",
        }
        return display_map.get(plot_type, plot_type.replace("_", " ").title())

    def _resolve_distribution_genes(self, adata: ad.AnnData, request: PlotRequest) -> list[str] | None:
        genes = request.genes or request.color
        genes = self._filter_placeholder_markers([str(gene) for gene in genes])
        if not genes:
            genes = list(self._session_markers)
        if not genes:
            return None
        return self._resolve_gene_names(adata, genes)

    def _resolve_distribution_groupby(self, adata: ad.AnnData, request: PlotRequest) -> str | None:
        return self._resolve_cell_type_column(adata, request.groupby)

    def _resolve_embedding_genes(
        self,
        adata: ad.AnnData,
        request: PlotRequest,
        resolved_color_columns: list[str],
    ) -> list[str] | None:
        requested_genes = [str(gene).strip() for gene in request.genes if str(gene).strip()]
        if requested_genes:
            return self._resolve_gene_names(adata, requested_genes)

        gene_like = [
            token
            for token in resolved_color_columns
            if token not in adata.obs.columns
        ]
        return gene_like or None

    def _build_coloring_label(
        self,
        adata: ad.AnnData,
        color_columns: list[str],
        resolved_groupby: str | None,
    ) -> str | None:
        if not color_columns:
            return "Uncolored embedding"

        if len(color_columns) == 1 and color_columns[0] in adata.obs.columns:
            column = color_columns[0]
            if resolved_groupby and column == resolved_groupby:
                return f"Cell types ({column})"
            return f"Colored by {column}"

        return f"Colored by {', '.join(color_columns)}"

    def _ensure_embedding(self, adata: ad.AnnData, embedding: str) -> None:
        target_key = f"X_{embedding}"
        if target_key in adata.obsm.keys():
            return

        if adata.isbacked:
            raise ValueError(
                f"AnnData is backed and missing '{target_key}'. Re-run without backed mode to compute embeddings."
            )

        if "X_pca" not in adata.obsm.keys():
            sc.pp.pca(adata)
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata)

        if embedding == "umap":
            sc.tl.umap(adata)
        elif embedding == "tsne":
            sc.tl.tsne(adata)

    @staticmethod
    def _normalize_embedding_key(key: str) -> str:
        return "".join(ch for ch in key.lower() if ch.isalnum())

    def _available_embeddings(self, adata: ad.AnnData) -> list[str]:
        """Return bare embedding names (without X_ prefix) sorted alphabetically."""
        return sorted(
            k[2:] if k.startswith("X_") else k
            for k in adata.obsm.keys()
            if k.startswith("X_") or not k.startswith("_")
        )

    def _find_embedding_basis(self, adata: ad.AnnData, embedding: str) -> str | None:
        target = self._normalize_embedding_key(embedding)
        obsm_keys = list(adata.obsm.keys())
        if not obsm_keys:
            return None

        priority_keys = [
            f"X_{embedding}",
            embedding,
        ]
        for key in priority_keys:
            if key in adata.obsm.keys():
                return key

        candidates: list[str] = []
        for key in obsm_keys:
            normalized = self._normalize_embedding_key(key)
            if target and target in normalized:
                candidates.append(key)

        if candidates:
            candidates.sort(key=lambda k: (0 if k.startswith("X_") else 1, len(k)))
            return candidates[0]
        return None

    def _find_first_matching_embedding_basis(self, adata: ad.AnnData, embeddings: list[str]) -> str | None:
        for embedding in embeddings:
            basis = self._find_embedding_basis(adata, embedding)
            if basis is not None:
                return basis
        return None

    @staticmethod
    def _basis_aliases(basis: str) -> list[str]:
        aliases = [basis]
        if basis.startswith("X_"):
            aliases.append(basis[2:])
        else:
            aliases.append(f"X_{basis}")
        deduped: list[str] = []
        for item in aliases:
            if item not in deduped:
                deduped.append(item)
        return deduped

    def _embedding_legend_options(self, adata: ad.AnnData, color: list[str]) -> dict[str, object]:
        options: dict[str, object] = {"legend_fontsize": self.EMBEDDING_LEGEND_FONTSIZE}
        if not color:
            return options

        primary_color = color[0]
        if primary_color not in adata.obs.columns:
            return options

        series = adata.obs[primary_color]
        if series is None:
            return options

        category_count = int(series.nunique(dropna=True))
        non_null = series.dropna().astype(str)
        max_label_len = int(non_null.str.len().max()) if not non_null.empty else 0

        normalized_primary_color = self._normalize_column_name(primary_color)
        is_cell_type_color = "celltype" in normalized_primary_color

        options["legend_loc"] = "right margin"

        if options.get("legend_loc") == "right margin":
            font_size = float(self.EMBEDDING_LEGEND_FONTSIZE)
            # Reduce for long labels
            if max_label_len > 50:
                font_size -= 1.0
            elif max_label_len > 35:
                font_size -= 0.5
            # Reduce further for many categories
            if category_count > 80:
                font_size -= 1.0
            elif category_count > 50:
                font_size -= 0.5
            options["legend_fontsize"] = max(self.EMBEDDING_LEGEND_MIN_FONTSIZE, font_size)

        return options

    # Scanpy places ~19 legend entries per column in right-margin mode.
    _LEGEND_ENTRIES_PER_COL = 19

    def _compute_embedding_figsize(
        self,
        adata: ad.AnnData,
        color: list[str],
        legend_options: dict[str, object],
    ) -> tuple[float, float]:
        """Adapt figure width and height to the legend content.

        Goals:
        - Scatter axes fills its area as a square.
        - Legend panel width matches the actual label lengths x number of columns.
        - Figure height is snug around the legend rows so there is no blank space.
        """
        if legend_options.get("legend_loc") != "right margin":
            return self.EMBEDDING_FIGSIZE

        primary_color = color[0] if color else None
        if not primary_color or primary_color not in adata.obs.columns:
            return self.EMBEDDING_FIGSIZE_WITH_RIGHT_LEGEND

        non_null = adata.obs[primary_color].dropna().astype(str)
        category_count = int(non_null.nunique())
        max_label_len = int(non_null.str.len().max()) if not non_null.empty else 10

        font_size = float(legend_options.get("legend_fontsize", self.EMBEDDING_LEGEND_FONTSIZE))
        scale = font_size / self.EMBEDDING_LEGEND_FONTSIZE

        # Estimate the number of columns scanpy will use (~19 entries/col).
        ncols = max(1, math.ceil(category_count / self._LEGEND_ENTRIES_PER_COL))
        n_rows = math.ceil(category_count / ncols)

        # Figure height: snug around legend rows + top/bottom margins.
        # Each row is ~0.22" at 7pt; scale with actual font size.
        row_height_in = 0.22 * scale
        legend_content_h = n_rows * row_height_in
        fig_height = max(6.0, min(legend_content_h + 1.4, 22.0))

        # Legend panel width: per-column label width x ncols + outer margin.
        # At 7pt, each character is ~0.068" wide (proportional font average).
        col_width_in = max_label_len * 0.068 * scale + 0.55
        legend_panel_w = ncols * col_width_in + 0.8

        # Figure width: square scatter area (== fig_height) + legend panel.
        fig_width = max(10.0, min(fig_height + legend_panel_w, 32.0))

        return (fig_width, fig_height)

    def _plot_embedding(
        self,
        adata: ad.AnnData,
        basis: str,
        color: list[str],
        preferred_gene_symbols_column: str | None = None,
        frameon: bool | None = None,
        ncols: int | None = None,
        point_size: float | None = None,
        palette: str | None = None,
    ) -> None:
        last_error: Exception | None = None
        plot_color, gene_symbols_column = self._resolve_embedding_color(
            adata,
            color,
            preferred_gene_symbols_column=preferred_gene_symbols_column,
        )
        legend_options = self._embedding_legend_options(adata, color)
        figure_size = self._compute_embedding_figsize(adata, color, legend_options)
        for alias in self._basis_aliases(basis):
            try:
                with rc_context({"figure.figsize": figure_size}):
                    if plot_color:
                        kwargs: dict[str, object] = {
                            "basis": alias,
                            "color": plot_color,
                            "show": False,
                            **legend_options,
                        }
                        if point_size is not None:
                            kwargs["size"] = float(point_size)
                        if palette:
                            kwargs["palette"] = str(palette)
                        if frameon is not None:
                            kwargs["frameon"] = bool(frameon)
                        if ncols is not None and ncols > 0:
                            kwargs["ncols"] = int(ncols)
                        if gene_symbols_column:
                            kwargs["gene_symbols"] = gene_symbols_column
                            kwargs["use_raw"] = False
                        sc.pl.embedding(
                            adata,
                            **kwargs,
                        )
                    else:
                        kwargs_no_color: dict[str, object] = {
                            "basis": alias,
                            **legend_options,
                            "show": False,
                        }
                        if point_size is not None:
                            kwargs_no_color["size"] = float(point_size)
                        if palette:
                            kwargs_no_color["palette"] = str(palette)
                        if frameon is not None:
                            kwargs_no_color["frameon"] = bool(frameon)
                        if ncols is not None and ncols > 0:
                            kwargs_no_color["ncols"] = int(ncols)
                        sc.pl.embedding(adata, **kwargs_no_color)
                return
            except Exception as exc:
                last_error = exc
        if last_error is not None:
            raise last_error

    @staticmethod
    def _normalize_column_name(name: str) -> str:
        return "".join(ch for ch in str(name).lower() if ch.isalnum())

    def _filter_placeholder_markers(self, markers: list[str]) -> list[str]:
        filtered: list[str] = []
        for marker in markers:
            normalized = self._normalize_column_name(marker)
            if normalized in self._PLACEHOLDER_MARKER_TOKENS:
                continue
            filtered.append(marker)
        return filtered

    @staticmethod
    def _is_usable_categorical_series(series: pd.Series) -> bool:
        if series is None:
            return False
        if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_categorical_dtype(series):
            return False
        unique_count = int(series.nunique(dropna=True))
        return unique_count > 1

    @staticmethod
    def _series_looks_like_ontology_ids(series: pd.Series) -> bool:
        non_null = series.dropna()
        if non_null.empty:
            return False

        sample = non_null.astype(str).head(300)
        if sample.empty:
            return False

        id_pattern = re.compile(r"^(?:[A-Za-z][A-Za-z0-9_]+:\d+|CL:\d+)$")
        matched = sum(1 for value in sample if id_pattern.match(value.strip()))
        return (matched / len(sample)) >= 0.6

    @staticmethod
    def _series_looks_like_numeric_labels(series: pd.Series) -> bool:
        non_null = series.dropna()
        if non_null.empty:
            return False

        sample = non_null.astype(str).head(300)
        if sample.empty:
            return False

        numeric_pattern = re.compile(r"^[+-]?(?:\d+\.\d+|\d+)$")
        matched = sum(1 for value in sample if numeric_pattern.match(value.strip()))
        return (matched / len(sample)) >= 0.8

    def _infer_cell_type_color_column(self, adata: ad.AnnData) -> str | None:
        columns = [str(column) for column in adata.obs.columns]
        lower_to_original = {column.lower(): column for column in columns}

        human_readable_candidates: dict[str, int] = {}
        fallback_candidates: dict[str, int] = {}

        for column in columns:
            lowered = column.lower()
            if "_ontology_term_id" not in lowered:
                continue

            base_name = lowered.replace("_ontology_term_id", "")
            base_normalized = self._normalize_column_name(base_name)
            if "celltype" not in base_normalized:
                continue

            readable_column = lower_to_original.get(base_name)
            if not readable_column:
                continue

            readable_series = adata.obs[readable_column]
            if not self._is_usable_categorical_series(readable_series):
                continue

            unique_count = int(readable_series.nunique(dropna=True))
            is_ontology_like = self._series_looks_like_ontology_ids(readable_series)
            is_numeric_like = self._series_looks_like_numeric_labels(readable_series)
            if not is_ontology_like and not is_numeric_like:
                human_readable_candidates[readable_column] = unique_count
            fallback_candidates[readable_column] = unique_count

        for column in columns:
            normalized = self._normalize_column_name(column)
            looks_like_cell_annotation = (
                "cell" in normalized
                or "cluster" in normalized
                or "annotation" in normalized
                or "label" in normalized
                or "subtype" in normalized
                or "lineage" in normalized
                or "identity" in normalized
                or "class" in normalized
            )
            if not looks_like_cell_annotation:
                continue

            series = adata.obs[column]
            if not self._is_usable_categorical_series(series):
                continue

            unique_count = int(series.nunique(dropna=True))
            is_ontology_like = self._series_looks_like_ontology_ids(series)
            is_numeric_like = self._series_looks_like_numeric_labels(series)
            if not is_ontology_like and not is_numeric_like:
                human_readable_candidates[str(column)] = unique_count
            fallback_candidates[str(column)] = unique_count

        if human_readable_candidates:
            return max(human_readable_candidates.items(), key=lambda item: (item[1], item[0].lower()))[0]
        if fallback_candidates:
            return max(fallback_candidates.items(), key=lambda item: (item[1], item[0].lower()))[0]
        return None

    def _remember_cell_type_column_if_valid(self, adata: ad.AnnData, column: str | None) -> None:
        if not column:
            return
        if column not in adata.obs.columns:
            return
        series = adata.obs[column]
        if not self._is_usable_categorical_series(series):
            return
        if self._series_looks_like_ontology_ids(series):
            return
        if self._series_looks_like_numeric_labels(series):
            return
        self._session_cell_type_column = str(column)

    def _resolve_cell_type_column(self, adata: ad.AnnData, explicit_column: str | None = None) -> str | None:
        explicit = (explicit_column or "").strip()
        if explicit:
            normalized_explicit = self._normalize_column_name(explicit)
            placeholder_tokens = {"obscolumn", "groupby", "column", "none", "null", "na"}
            if normalized_explicit not in placeholder_tokens:
                # Exact match first.
                if explicit in adata.obs.columns:
                    series = adata.obs[explicit]
                    if self._is_usable_categorical_series(series):
                        self._session_cell_type_column = str(explicit)
                        return str(explicit)
                # Then case-insensitive exact match.
                lower_to_original = {str(col).lower(): str(col) for col in adata.obs.columns}
                matched = lower_to_original.get(explicit.lower())
                if matched:
                    series = adata.obs[matched]
                    if self._is_usable_categorical_series(series):
                        self._session_cell_type_column = matched
                        return matched

        remembered = self._session_cell_type_column
        if remembered and remembered in adata.obs.columns:
            series = adata.obs[remembered]
            if (
                self._is_usable_categorical_series(series)
                and not self._series_looks_like_ontology_ids(series)
                and not self._series_looks_like_numeric_labels(series)
            ):
                return remembered

        inferred = self._infer_cell_type_color_column(adata)
        if inferred:
            self._session_cell_type_column = inferred
        return inferred

    _GENE_SYMBOL_COLUMNS = (
        "feature_name",
        "gene_name",
        "gene_names",
        "gene_symbols",
        "symbol",
        "hgnc_symbol",
    )

    @staticmethod
    def _find_gene_symbol_column(adata: ad.AnnData) -> str | None:
        """Return the first var column that stores gene symbols, or None."""
        for col in ScanpyPlotExecutor._GENE_SYMBOL_COLUMNS:
            if col in adata.var.columns:
                return col
        return None

    @staticmethod
    def _has_feature_name_column(adata: ad.AnnData) -> bool:
        return ScanpyPlotExecutor._find_gene_symbol_column(adata) is not None

    @staticmethod
    def _gene_to_ens(df: pd.DataFrame, genes: list[str], column: str = "feature_name") -> list[str]:
        if column not in df.columns:
            raise KeyError(f"Column '{column}' not found in DataFrame")

        normalized_genes = [str(gene).strip().lower() for gene in genes if str(gene).strip()]
        if not normalized_genes:
            return []

        normalized_gene_set = set(normalized_genes)
        column_values = df[column].astype(str).str.strip().str.lower()
        return [str(item) for item in df.index[column_values.isin(normalized_gene_set)].tolist()]

    @staticmethod
    def _canonicalize_gene_tokens(genes: list[str]) -> list[str]:
        canonical_genes: list[str] = []
        seen_canonical: set[str] = set()
        for gene in genes:
            token = str(gene).strip()
            if not token:
                continue
            mapped = GENE_SYNONYM_MAP.get(token.upper(), token)
            lowered = mapped.lower()
            if lowered not in seen_canonical:
                canonical_genes.append(mapped)
                seen_canonical.add(lowered)
        return canonical_genes

    def _resolve_feature_name_genes(self, adata: ad.AnnData, genes: list[str]) -> list[str]:
        canonical_genes = self._canonicalize_gene_tokens(genes)
        if not canonical_genes:
            return []
        gene_sym_col = self._find_gene_symbol_column(adata)
        if not gene_sym_col:
            return canonical_genes

        feature_series = adata.var[gene_sym_col].astype(str)
        var_names = [str(name) for name in adata.var_names]

        feature_lookup: dict[str, str] = {}
        var_to_feature: dict[str, str] = {}
        for var_name, feature_name in zip(var_names, feature_series, strict=False):
            canonical_feature = str(feature_name).strip()
            if not canonical_feature or canonical_feature.lower() == "nan":
                continue
            var_to_feature[str(var_name).lower()] = canonical_feature
            lowered = canonical_feature.lower()
            if lowered not in feature_lookup:
                feature_lookup[lowered] = canonical_feature

        resolved: list[str] = []
        seen: set[str] = set()
        for gene in canonical_genes:
            gene_lower = gene.lower()
            resolved_name = feature_lookup.get(gene_lower) or var_to_feature.get(gene_lower) or gene
            resolved_key = resolved_name.lower()
            if resolved_key not in seen:
                resolved.append(resolved_name)
                seen.add(resolved_key)
        return resolved

    def _resolve_gene_names(self, adata: ad.AnnData, genes: list[str]) -> list[str]:
        if not genes:
            return []

        canonical_genes = self._canonicalize_gene_tokens(genes)

        if not self._has_feature_name_column(adata):
            resolved_without_feature_name: list[str] = []
            seen_without_feature_name: set[str] = set()
            for mapped in canonical_genes:
                if mapped not in seen_without_feature_name:
                    resolved_without_feature_name.append(mapped)
                    seen_without_feature_name.add(mapped)
            return resolved_without_feature_name

        resolved: list[str] = []
        seen: set[str] = set()
        resolved_gene_ids = self._gene_to_ens(adata.var, canonical_genes, column=self._find_gene_symbol_column(adata) or "feature_name")
        for gene_id in resolved_gene_ids:
            if gene_id not in seen:
                resolved.append(gene_id)
                seen.add(gene_id)

        # Keep compatibility when input already uses var_names or no feature_name match exists.
        available_var_names = {str(name).lower() for name in adata.var_names}
        for gene in canonical_genes:
            if gene.lower() in available_var_names and gene not in seen:
                resolved.append(gene)
                seen.add(gene)
            elif gene.lower() not in available_var_names and gene not in seen and not resolved_gene_ids:
                resolved.append(gene)
                seen.add(gene)

        return resolved

    def _resolve_feature_labels_for_var_keys(self, adata: ad.AnnData, var_keys: list[str]) -> dict[str, str]:
        gene_sym_col = self._find_gene_symbol_column(adata)
        if not var_keys or not gene_sym_col:
            return {}

        label_map: dict[str, str] = {}
        var_df = adata.var
        if var_df is None or var_df.empty:
            return label_map

        for key in var_keys:
            token = str(key).strip()
            if not token or token in label_map:
                continue
            if token not in var_df.index:
                continue
            feature_name = str(var_df.at[token, gene_sym_col]).strip()
            if feature_name and feature_name.lower() != "nan":
                label_map[token] = feature_name
        return label_map

    @staticmethod
    @staticmethod
    def _build_full_label_map(adata: ad.AnnData) -> dict[str, str]:
        """Return a var_name → symbol map for all genes (used for post-hoc axis relabeling)."""
        gene_sym_col = ScanpyPlotExecutor._find_gene_symbol_column(adata)
        if not gene_sym_col:
            return {}
        result: dict[str, str] = {}
        for var_name, symbol in zip(adata.var_names, adata.var[gene_sym_col].astype(str)):
            sym = symbol.strip()
            if sym and sym.lower() != "nan":
                result[str(var_name)] = sym
        return result

    @staticmethod
    def _relabel_var_axes(label_map: dict[str, str]) -> None:
        """Rename Ensembl-ID tick labels to gene symbols on all axes of the current figure."""
        if not label_map:
            return
        fig = plt.gcf()
        for ax in fig.axes:
            x_labels = [t.get_text() for t in ax.get_xticklabels()]
            new_x = [label_map.get(t, t) for t in x_labels]
            if x_labels and x_labels != new_x:
                ax.set_xticklabels(new_x)
            y_labels = [t.get_text() for t in ax.get_yticklabels()]
            new_y = [label_map.get(t, t) for t in y_labels]
            if y_labels and y_labels != new_y:
                ax.set_yticklabels(new_y)

    @staticmethod
    def _relabel_violin_axes(label_map: dict[str, str]) -> None:
        if not label_map:
            return

        fig = plt.gcf()
        for ax in fig.axes:
            ylabel = str(ax.get_ylabel()).strip()
            if ylabel in label_map:
                ax.set_ylabel(label_map[ylabel])
            title = str(ax.get_title()).strip()
            if title in label_map:
                ax.set_title(label_map[title])

    def _resolve_embedding_color(
        self,
        adata: ad.AnnData,
        color: list[str],
        preferred_gene_symbols_column: str | None = None,
    ) -> tuple[list[str], str | None]:
        if not color:
            return color, None

        obs_columns = {str(column) for column in adata.obs.columns}
        requested_symbol_column = (preferred_gene_symbols_column or "").strip()
        detected_sym_col = self._find_gene_symbol_column(adata)
        symbol_columns: list[str] = []
        if requested_symbol_column and requested_symbol_column in adata.var.columns:
            symbol_columns.append(requested_symbol_column)
        if detected_sym_col and detected_sym_col not in symbol_columns:
            symbol_columns.append(detected_sym_col)

        available_var_names = {str(name).lower() for name in adata.var_names}
        available_symbol_names_by_column: dict[str, set[str]] = {}
        for column in symbol_columns:
            symbol_series = adata.var[column].astype(str)
            available_symbol_names_by_column[column] = {
                str(name).strip().lower()
                for name in symbol_series
                if str(name).strip() and str(name).strip().lower() != "nan"
            }

        gene_like_tokens: list[str] = []
        for token in color:
            token_str = str(token)
            if token_str in obs_columns:
                continue
            gene_like_tokens.append(token_str)

        replacement_map: dict[str, str] = {}
        if gene_like_tokens:
            if requested_symbol_column and requested_symbol_column in adata.var.columns and requested_symbol_column != detected_sym_col:
                resolved_gene_like_tokens = self._canonicalize_gene_tokens(gene_like_tokens)
            else:
                resolved_gene_like_tokens = self._resolve_gene_names(adata, gene_like_tokens)
            replacement_map = dict(zip(gene_like_tokens, resolved_gene_like_tokens, strict=False))

        resolved_color: list[str] = []
        resolved_gene_symbols_column: str | None = None
        for token in color:
            token_str = str(token)
            if token_str in obs_columns:
                resolved_color.append(token_str)
                continue

            candidate = replacement_map.get(token_str, token_str)
            candidate_lower = candidate.lower()
            if candidate_lower in available_var_names:
                resolved_color.append(candidate)
                continue
            for symbol_column in symbol_columns:
                available_symbol_names = available_symbol_names_by_column.get(symbol_column, set())
                if candidate_lower in available_symbol_names:
                    resolved_color.append(candidate)
                    if resolved_gene_symbols_column is None:
                        resolved_gene_symbols_column = symbol_column
                    break

        return resolved_color, resolved_gene_symbols_column

    def _compute_dotplot_figsize(self, adata: ad.AnnData, groupby: str, gene_count: int) -> tuple[float, float]:
        category_count = 0
        if groupby in adata.obs.columns:
            category_count = int(adata.obs[groupby].nunique(dropna=True))

        raw_width = (
            self.DOTPLOT_BASE_WIDTH
            + self.DOTPLOT_WIDTH_PER_GENE * max(1, gene_count)
            + self.DOTPLOT_WIDTH_PER_GROUP * max(1, category_count)
        )
        bounded_width = min(self.DOTPLOT_MAX_WIDTH, max(self.DOTPLOT_MIN_WIDTH, raw_width))
        return float(bounded_width), float(self.DOTPLOT_HEIGHT)

    def _plot_umap(self, adata: ad.AnnData, request: PlotRequest) -> tuple[str, str, list[str] | None]:
        basis = self._find_first_matching_embedding_basis(adata, ["umap", "tsne"])

        if basis is None:
            try:
                self._ensure_embedding(adata, "umap")
                basis = self._find_first_matching_embedding_basis(adata, ["umap", "tsne"])
            except ValueError:
                tsne_basis = self._find_embedding_basis(adata, "tsne")
                if tsne_basis is not None:
                    basis = tsne_basis
                else:
                    self._ensure_embedding(adata, "tsne")
                    basis = self._find_embedding_basis(adata, "tsne")

        if basis is None:
            raise ValueError("No UMAP/tSNE-like embedding basis found or computed.")

        color = request.color or request.genes or []
        if request.color:
            self._remember_cell_type_column_if_valid(adata, request.color[0])
        resolved_color_preview, _ = self._resolve_embedding_color(adata, color) if color else ([], None)
        if not resolved_color_preview:
            inferred_color = self._resolve_cell_type_column(adata)
            if inferred_color:
                color = [inferred_color]
        self._plot_embedding(
            adata,
            basis=basis,
            color=color,
            point_size=self.UMAP_TSNE_POINT_SIZE,
        )

        if "tsne" in self._normalize_embedding_key(basis):
            return "tsne", basis, color or None
        return basis, basis, color or None

    def _plot_tsne(self, adata: ad.AnnData, request: PlotRequest) -> tuple[str, str, list[str] | None]:
        basis = self._find_embedding_basis(adata, "tsne")
        if basis is None:
            self._ensure_embedding(adata, "tsne")
            basis = self._find_embedding_basis(adata, "tsne")
        if basis is None:
            raise ValueError("No tSNE-like embedding basis found or computed.")

        color = request.color or request.genes or []
        if request.color:
            self._remember_cell_type_column_if_valid(adata, request.color[0])
        resolved_color_preview, _ = self._resolve_embedding_color(adata, color) if color else ([], None)
        if not resolved_color_preview:
            inferred_color = self._resolve_cell_type_column(adata)
            if inferred_color:
                color = [inferred_color]
        self._plot_embedding(
            adata,
            basis=basis,
            color=color,
            point_size=self.UMAP_TSNE_POINT_SIZE,
        )
        return "tsne", basis, color or None

    def _plot_gene_cell_embedding(self, adata: ad.AnnData, request: PlotRequest) -> tuple[str, str, list[str] | None]:
        basis = self._find_first_matching_embedding_basis(adata, ["umap", "tsne"])

        if basis is None:
            try:
                self._ensure_embedding(adata, "umap")
                basis = self._find_first_matching_embedding_basis(adata, ["umap", "tsne"])
            except ValueError:
                tsne_basis = self._find_embedding_basis(adata, "tsne")
                if tsne_basis is not None:
                    basis = tsne_basis
                else:
                    self._ensure_embedding(adata, "tsne")
                    basis = self._find_embedding_basis(adata, "tsne")

        if basis is None:
            raise ValueError("No UMAP/tSNE-like embedding basis found or computed.")

        color: list[str] = []
        seen_color_tokens: set[str] = set()
        for token in request.color:
            token_str = str(token).strip()
            if not token_str:
                continue
            token_lower = token_str.lower()
            if token_lower not in seen_color_tokens:
                color.append(token_str)
                seen_color_tokens.add(token_lower)

        for gene in request.genes:
            gene_str = str(gene).strip()
            if not gene_str:
                continue
            gene_lower = gene_str.lower()
            if gene_lower not in seen_color_tokens:
                color.append(gene_str)
                seen_color_tokens.add(gene_lower)

        requested_groupby = (request.groupby or "").strip()
        if requested_groupby:
            self._remember_cell_type_column_if_valid(adata, requested_groupby)
            groupby_lower = requested_groupby.lower()
            if groupby_lower not in seen_color_tokens:
                color.append(requested_groupby)
                seen_color_tokens.add(groupby_lower)

        resolved_color_preview, _ = self._resolve_embedding_color(
            adata,
            color,
            preferred_gene_symbols_column=request.gene_symbols_column,
        )
        has_obs_color = any(token in adata.obs.columns for token in resolved_color_preview)
        if not has_obs_color:
            inferred_groupby = self._resolve_cell_type_column(adata, request.groupby)
            if inferred_groupby:
                inferred_lower = inferred_groupby.lower()
                if inferred_lower not in seen_color_tokens:
                    color.append(inferred_groupby)

        self._plot_embedding(
            adata,
            basis=basis,
            color=color,
            preferred_gene_symbols_column=request.gene_symbols_column,
            frameon=False,
            ncols=3,
        )

        if "tsne" in self._normalize_embedding_key(basis):
            return "gene_cell_embedding", basis, color or None
        return "gene_cell_embedding", basis, color or None

    def _plot_cell_counts_barplot(self, adata: ad.AnnData, request: PlotRequest) -> None:
        groupby = self._resolve_cell_type_column(adata, request.groupby)
        if not groupby or groupby not in adata.obs.columns:
            raise ValueError(
                f"Column '{groupby}' not found in adata.obs. "
                f"Available columns: {', '.join(str(c) for c in adata.obs.columns[:20])}"
            )

        cell_counts = adata.obs[groupby].value_counts().reset_index()
        cell_counts.columns = [groupby, "count"]
        cell_counts = cell_counts.sort_values("count", ascending=False).reset_index(drop=True)

        n_categories = len(cell_counts)
        max_label_len = int(cell_counts[groupby].astype(str).str.len().max()) if n_categories else 0

        # Width: 0.55 in per bar, min 8, max 36
        fig_width = max(8.0, min(n_categories * 0.55, 36.0))
        fig_height = 6.0

        # X-label font size scales down for many/long labels
        if n_categories > 50 or max_label_len > 25:
            xlabel_fontsize = 6
        elif n_categories > 25 or max_label_len > 15:
            xlabel_fontsize = 8
        else:
            xlabel_fontsize = 10

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.barplot(
            data=cell_counts,
            x=groupby,
            y="count",
            order=cell_counts[groupby].tolist(),
            ax=ax,
        )

        title = str(request.title or "").strip() or f"Cell Counts by {groupby.replace('_', ' ').title()}"
        ax.set_title(title)
        ax.set_xlabel(groupby.replace("_", " ").title())
        ax.set_ylabel("Cell Count")
        # ha='right' anchors each label's right edge at its tick mark, preventing overlap
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=xlabel_fontsize)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

        # Count labels on top of each bar
        count_fontsize = max(5, min(xlabel_fontsize, 9))
        for patch, count in zip(ax.patches, cell_counts["count"]):
            ax.text(
                patch.get_x() + patch.get_width() / 2,
                patch.get_height(),
                f"{int(count):,}",
                ha="center",
                va="bottom",
                fontsize=count_fontsize,
                clip_on=False,
            )

        # Bottom margin proportional to label length so rotated labels aren't clipped
        bottom_margin = min(0.12 + max_label_len * 0.009, 0.45)
        fig.subplots_adjust(bottom=bottom_margin)

    def _resolve_sample_column(self, adata: ad.AnnData, explicit_column: str | None = None) -> str | None:
        """Find a sample/donor ID column in adata.obs.

        Priority:
        1. Explicit column name (exact then case-insensitive match).
        2. Columns whose lowercased name matches a known sample-ID keyword pattern.
        3. Any column whose name contains 'sample'.
        Returns None if nothing suitable is found.
        """
        obs_cols = list(adata.obs.columns)
        lower_to_col = {str(c).lower(): str(c) for c in obs_cols}

        # 1. Explicit column provided by user
        if explicit_column and explicit_column.strip():
            explicit = explicit_column.strip()
            if explicit in obs_cols:
                return explicit
            matched = lower_to_col.get(explicit.lower())
            if matched:
                return matched

        # 2. Known sample-ID keyword patterns (ordered by priority)
        _SAMPLE_KEYWORDS = [
            "sample_id",
            "sampleid",
            "sample",
            "sample_name",
            "donor_id",
            "donor",
            "patient_id",
            "patient",
            "subject_id",
            "subject",
            "batch",
        ]
        for keyword in _SAMPLE_KEYWORDS:
            matched = lower_to_col.get(keyword)
            if matched:
                return matched

        # 3. Fallback: any column whose name contains 'sample'
        for lower_name, original in lower_to_col.items():
            if "sample" in lower_name:
                return original

        return None

    def _plot_cell_type_proportion_barplot(self, adata: ad.AnnData, request: PlotRequest) -> None:
        """Stacked bar plot: percentage of each cell type per sample."""
        # Resolve cell type column (default: author_cell_type, user-customisable via groupby)
        cell_type_col = self._resolve_cell_type_column(adata, request.groupby)
        if not cell_type_col or cell_type_col not in adata.obs.columns:
            raise ValueError(
                f"Cell type column '{cell_type_col}' not found in adata.obs. "
                f"Available columns: {', '.join(str(c) for c in adata.obs.columns[:20])}"
            )

        # Resolve sample column – try color list first, then auto-detect
        sample_col_hint = request.color[0] if request.color else None
        sample_col = self._resolve_sample_column(adata, sample_col_hint)
        if not sample_col or sample_col not in adata.obs.columns:
            raise ValueError(
                "Could not find a sample/donor ID column in adata.obs. "
                "Please specify the column name explicitly. "
                f"Available columns: {', '.join(str(c) for c in adata.obs.columns[:20])}"
            )

        # Build percentage cross-tab: rows = cell types, columns = samples
        cross_tab = pd.crosstab(
            adata.obs[cell_type_col],
            adata.obs[sample_col],
            normalize="columns",
        ) * 100.0

        n_samples = cross_tab.shape[1]
        n_cell_types = cross_tab.shape[0]

        fig_width = max(8.0, min(n_samples * 0.8 + 3.0, 36.0))
        fig_height = 6.0
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        cross_tab.T.plot(kind="bar", stacked=True, ax=ax, width=0.8)

        # Legend outside the plot area; adjust right margin to make room
        legend_cols = max(1, n_cell_types // 25 + 1)
        ax.legend(
            title=cell_type_col.replace("_", " ").title(),
            bbox_to_anchor=(1.02, 1.0),
            loc="upper left",
            ncol=legend_cols,
            fontsize=7,
            title_fontsize=8,
        )

        title = str(request.title or "").strip() or (
            f"Cell Type Proportions per {sample_col.replace('_', ' ').title()}"
        )
        ax.set_title(title)
        ax.set_xlabel(sample_col.replace("_", " ").title())
        ax.set_ylabel("Percentage (%)")

        max_label_len = int(
            pd.Series(cross_tab.columns.astype(str)).str.len().max()
        ) if n_samples else 0
        if n_samples > 40 or max_label_len > 20:
            xlabel_fontsize = 6
        elif n_samples > 20 or max_label_len > 12:
            xlabel_fontsize = 8
        else:
            xlabel_fontsize = 10
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=xlabel_fontsize)

        bottom_margin = min(0.12 + max_label_len * 0.009, 0.45)
        fig.subplots_adjust(bottom=bottom_margin)

    def _plot_violin(self, adata: ad.AnnData, request: PlotRequest) -> None:
        keys = request.genes or request.color
        keys = self._filter_placeholder_markers([str(key) for key in keys])
        if not keys:
            keys = list(self._session_markers)
        if not keys:
            self._plot_rank_genes_groups_stacked_violin(adata, request)
            return
        if request.color:
            self._remember_cell_type_column_if_valid(adata, request.color[0])

        plot_keys = self._resolve_gene_names(adata, keys)
        label_map = self._resolve_feature_labels_for_var_keys(adata, plot_keys)

        # Always default to author_cell_type groupby unless user explicitly specified something else
        groupby = self._resolve_cell_type_column(adata, request.groupby or self.DEFAULT_ACTIVE_CELL_TYPE_COLUMN)
        kwargs: dict[str, object] = {"keys": plot_keys, "show": False}

        n_categories = 0
        if groupby and groupby in adata.obs.columns:
            n_categories = int(adata.obs[groupby].nunique(dropna=True))
            max_label_len = int(adata.obs[groupby].dropna().astype(str).str.len().max()) if n_categories else 0
        else:
            max_label_len = 0

        # Dynamic figure width: at least 0.6 in per category, min 10, max 36
        fig_width = max(10.0, min(n_categories * 0.6, 36.0)) if n_categories else 12.0
        # Dynamic x-label font size: shrink for many/long labels
        if n_categories > 40 or max_label_len > 30:
            xlabel_fontsize = 6
        elif n_categories > 20 or max_label_len > 20:
            xlabel_fontsize = 8
        else:
            xlabel_fontsize = 10

        if groupby:
            kwargs["groupby"] = groupby
            kwargs["use_raw"] = False
        elif len(plot_keys) > 1:
            kwargs["multi_panel"] = True
            kwargs["jitter"] = 0.4

        with rc_context({"figure.figsize": (fig_width, 6.0)}):
            sc.pl.violin(adata, **kwargs)

        # Apply ha='right' so rotated labels fan out without overlapping
        fig = plt.gcf()
        for ax in fig.axes:
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=xlabel_fontsize)
        # Reserve bottom space proportional to label length and rotation
        bottom_margin = min(0.08 + max_label_len * 0.008, 0.45)
        fig.subplots_adjust(bottom=bottom_margin)

        self._relabel_violin_axes(label_map)

    def _ensure_rank_genes_groups(self, adata: ad.AnnData, request: PlotRequest) -> str:
        requested_groupby = self._resolve_cell_type_column(adata, request.groupby)
        if not requested_groupby:
            raise ValueError("rank_genes_groups plots require groupby, and no inferred group column was found.")

        cache_state = adata.uns.get(self._RANK_GENES_CACHE_KEY)
        if isinstance(cache_state, dict):
            cached_groupby = str(cache_state.get("groupby", "")).strip()
            if cached_groupby and isinstance(adata.uns.get("rank_genes_groups"), dict):
                self._rank_genes_groups_used_groupby = cached_groupby
                self._rank_genes_groups_notice = (
                    f"Reusing cached rank_genes_groups computed with groupby '{cached_groupby}'."
                )
                if cached_groupby != requested_groupby:
                    self._rank_genes_groups_notice += (
                        f" Requested groupby '{requested_groupby}' was ignored to avoid recomputation."
                    )
                return cached_groupby

        existing = adata.uns.get("rank_genes_groups")
        if isinstance(existing, dict):
            params = existing.get("params") if isinstance(existing.get("params"), dict) else {}
            existing_groupby = str(params.get("groupby", "")).strip()
            if existing_groupby:
                adata.uns[self._RANK_GENES_CACHE_KEY] = {
                    "groupby": existing_groupby,
                    "computed_by": "existing",
                }
                self._rank_genes_groups_used_groupby = existing_groupby
                self._rank_genes_groups_notice = (
                    f"Reusing existing rank_genes_groups computed with groupby '{existing_groupby}'."
                )
                if existing_groupby != requested_groupby:
                    self._rank_genes_groups_notice += (
                        f" Requested groupby '{requested_groupby}' was ignored to avoid recomputation."
                    )
                return existing_groupby

        if adata.isbacked:
            raise ValueError("rank_genes_groups computation requires non-backed AnnData. Re-run with backed=false.")

        print(
            f"[genopixel] {self._RANK_GENES_PREP_WARNING} groupby='{requested_groupby}'",
            flush=True,
        )
        self._rank_genes_groups_notice = (
            f"rank_genes_groups was computed for this active dataset using groupby '{requested_groupby}'. "
            "Subsequent rank-gene plots will reuse this cached result."
        )
        self._rank_genes_groups_computed_this_run = True

        sc.tl.rank_genes_groups(
            adata,
            groupby=requested_groupby,
            use_raw=False,
            method="wilcoxon",
            log=True,
        )

        # Prepare group ordering used by rank_genes_groups_* plotting functions.
        sc.tl.dendrogram(
            adata,
            groupby=requested_groupby,
            use_raw=False,
        )
        adata.uns[self._RANK_GENES_CACHE_KEY] = {
            "groupby": requested_groupby,
            "computed_by": "genopixel",
        }
        self._rank_genes_groups_used_groupby = requested_groupby
        return requested_groupby

    @staticmethod
    def _rgg_figsize(rgg: dict, n_genes: int, swap_axes: bool = False) -> tuple[float, float]:
        """Compute figure size for rank_genes_groups plots where total x-axis cols = n_genes × n_groups."""
        try:
            all_group_names = list(rgg["names"].dtype.names or [])
        except Exception:
            all_group_names = []
        n_g = len(all_group_names) if all_group_names else 8
        if swap_axes:
            return (max(6.0, n_g * 0.5 + 2), max(4.0, n_genes * n_g * 0.3 + 2))
        return (max(8.0, n_genes * n_g * 0.25 + 2), max(4.0, n_g * 0.5 + 2))

    def _plot_rank_genes_groups_dotplot(self, adata: ad.AnnData, request: PlotRequest) -> None:
        groupby = self._ensure_rank_genes_groups(adata, request)
        n_genes = 5
        figsize = self._rgg_figsize(adata.uns["rank_genes_groups"], n_genes)
        sym_col = self._find_gene_symbol_column(adata)
        kwargs: dict[str, object] = {
            "groupby": groupby,
            "standard_scale": "var",
            "n_genes": n_genes,
            "dendrogram": True,
            "figsize": figsize,
            "show": False,
        }
        if sym_col:
            kwargs["gene_symbols"] = sym_col
        sc.pl.rank_genes_groups_dotplot(adata, **kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

    def _plot_rank_genes_groups_matrixplot(self, adata: ad.AnnData, request: PlotRequest) -> None:
        groupby = self._ensure_rank_genes_groups(adata, request)
        n_genes = 5
        figsize = self._rgg_figsize(adata.uns["rank_genes_groups"], n_genes)
        sym_col = self._find_gene_symbol_column(adata)
        kwargs: dict[str, object] = {
            "groupby": groupby,
            "n_genes": n_genes,
            "use_raw": False,
            "vmin": -3,
            "vmax": 3,
            "cmap": "bwr",
            "figsize": figsize,
            "show": False,
        }
        if sym_col:
            kwargs["gene_symbols"] = sym_col
        sc.pl.rank_genes_groups_matrixplot(adata, **kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

    def _plot_rank_genes_groups_stacked_violin(self, adata: ad.AnnData, request: PlotRequest) -> None:
        groupby = self._ensure_rank_genes_groups(adata, request)
        n_genes = 5
        figsize = self._rgg_figsize(adata.uns["rank_genes_groups"], n_genes)
        sym_col = self._find_gene_symbol_column(adata)
        kwargs: dict[str, object] = {
            "groupby": groupby,
            "n_genes": n_genes,
            "cmap": "viridis_r",
            "figsize": figsize,
            "show": False,
        }
        if sym_col:
            kwargs["gene_symbols"] = sym_col
        sc.pl.rank_genes_groups_stacked_violin(adata, **kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

    def _plot_rank_genes_groups_heatmap(self, adata: ad.AnnData, request: PlotRequest) -> None:
        groupby = self._ensure_rank_genes_groups(adata, request)
        n_genes = 5
        # swap_axes=True: genes on x-axis, groups on y-axis
        figsize = self._rgg_figsize(adata.uns["rank_genes_groups"], n_genes, swap_axes=True)
        kwargs: dict[str, object] = {
            "groupby": groupby,
            "n_genes": n_genes,
            "use_raw": False,
            "swap_axes": True,
            "show_gene_labels": False,
            "vmin": -3,
            "vmax": 3,
            "cmap": "bwr",
            "figsize": figsize,
            "show": False,
        }
        sc.pl.rank_genes_groups_heatmap(adata, **kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

    def _plot_correlation_matrix(self, adata: ad.AnnData, request: PlotRequest) -> None:
        groupby = self._resolve_cell_type_column(adata, request.groupby)
        if not groupby:
            raise ValueError("correlation_matrix requires groupby, and no inferred group column was found.")
        sc.pl.correlation_matrix(adata, groupby, show=False)

    def _plot_dotplot(self, adata: ad.AnnData, request: PlotRequest) -> None:
        genes = request.genes or request.color
        genes = self._filter_placeholder_markers([str(gene) for gene in genes])
        if not genes:
            # Notebook-aligned marker workflow:
            # 1) compute rank_genes_groups by inferred/provided cell-type grouping
            # 2) render grouped marker dotplot
            self._plot_rank_genes_groups_dotplot(adata, request)
            return
        if request.color:
            self._remember_cell_type_column_if_valid(adata, request.color[0])

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = self._resolve_gene_names(adata, genes)

        groupby = self._resolve_cell_type_column(adata, request.groupby)
        if not groupby:
            raise ValueError("dotplot requires groupby column, and no inferred group column was found.")
        figure_size = self._compute_dotplot_figsize(adata, groupby=groupby, gene_count=len(plot_genes))
        kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": groupby,
            "dendrogram": False,
            "figsize": figure_size,
            "use_raw": False,
            "show": False,
        }
        sc.pl.dotplot(adata, **kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

    def _plot_matrixplot(self, adata: ad.AnnData, request: PlotRequest) -> None:
        genes = request.genes or request.color
        genes = self._filter_placeholder_markers([str(gene) for gene in genes])
        if not genes:
            self._plot_rank_genes_groups_matrixplot(adata, request)
            return
        if request.color:
            self._remember_cell_type_column_if_valid(adata, request.color[0])

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = self._resolve_gene_names(adata, genes)

        groupby = self._resolve_cell_type_column(adata, request.groupby)
        if not groupby:
            raise ValueError("matrixplot requires groupby column, and no inferred group column was found.")
        kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": groupby,
            "dendrogram": False,
            "use_raw": False,
            "show": False,
        }
        sc.pl.matrixplot(adata, **kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))

    def _plot_heatmap(self, adata: ad.AnnData, request: PlotRequest) -> None:
        genes = request.genes or request.color
        genes = self._filter_placeholder_markers([str(gene) for gene in genes])
        if not genes:
            self._plot_rank_genes_groups_heatmap(adata, request)
            return
        if request.color:
            self._remember_cell_type_column_if_valid(adata, request.color[0])

        use_gene_symbols = self._has_feature_name_column(adata)
        plot_genes = self._resolve_gene_names(adata, genes)

        groupby = self._resolve_cell_type_column(adata, request.groupby)
        if not groupby:
            raise ValueError("heatmap requires groupby column, and no inferred group column was found.")
        kwargs: dict[str, object] = {
            "var_names": plot_genes,
            "groupby": groupby,
            "use_raw": False,
            "show": False,
        }
        sc.pl.heatmap(adata, **kwargs)
        self._relabel_var_axes(self._build_full_label_map(adata))
