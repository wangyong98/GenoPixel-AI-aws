from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


ALL_SHEET_SEARCH_COLUMNS = [
    "doi",
    "cellxgene_doi",
    "title",
    "author",
    "year",
    "journal",
    "tissue",
    "tissue_type",
    "disease",
    "organism",
]


@dataclass
class SearchFilters:
    keywords: list[str] = field(default_factory=list)
    column_filters: dict[str, str] = field(default_factory=dict)


@dataclass
class PlotRequest:
    plot_type: str = "umap"
    color: list[str] = field(default_factory=list)
    genes: list[str] = field(default_factory=list)
    groupby: str | None = None
    gene_symbols_column: str | None = None
    title: str | None = None


@dataclass
class UserIntent:
    raw_query: str
    search: SearchFilters
    plot: PlotRequest


@dataclass
class MatchCandidate:
    sheet: str
    row_number: int
    row_data: dict
    exact_hits: int
    fuzzy_score: float


@dataclass
class ResolvedDataset:
    selected_all_row: MatchCandidate
    selected_multiple_row: MatchCandidate | None
    h5ad_file_value: str
    h5ad_path: Path


@dataclass
class PlotResult:
    plot_type: str
    output_file: Path
    embedding_basis: str | None = None
    color_columns: list[str] | None = None
    resolved_genes: list[str] | None = None
    resolved_groupby: str | None = None
    resolved_coloring_label: str | None = None
    display_plot_type: str | None = None
    rank_genes_groups_computed: bool | None = None
    rank_genes_groups_notice: str | None = None
