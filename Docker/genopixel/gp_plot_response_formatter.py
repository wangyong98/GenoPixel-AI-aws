from __future__ import annotations

from typing import Any


EMBEDDING_PLOT_TYPES = {"umap", "tsne", "gene_cell_embedding"}
DISTRIBUTION_PLOT_TYPES = {"violin", "dotplot", "matrixplot", "heatmap", "tracksplot", "stacked_violin"}
RANK_GENES_PLOT_TYPES = {
    "rank_genes_groups_dotplot",
    "rank_genes_groups_matrixplot",
    "rank_genes_groups_stacked_violin",
    "rank_genes_groups_heatmap",
}

DISPLAY_NAMES = {
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
    "cell_type_proportion_barplot": "Cell type proportions",
    "tracksplot": "Tracks plot",
    "stacked_violin": "Stacked violin plot",
    "clustermap": "Cluster map",
    "dendrogram": "Dendrogram",
    "tsne": "tSNE embedding",
    "diffmap": "Diffusion map",
    "embedding": "Embedding",
    "rank_genes_groups": "Rank genes groups",
    "rank_genes_groups_violin": "Rank genes groups violin",
}


def _format_int(value: Any) -> str:
    try:
        return f"{int(value):,}"
    except Exception:
        return str(value)


def _join_tokens(tokens: list[str] | None) -> str:
    return ", ".join(str(token).strip() for token in (tokens or []) if str(token).strip())


def _display_name(plot_type: str, explicit: str | None) -> str:
    if explicit and explicit.strip():
        return explicit.strip()
    return DISPLAY_NAMES.get(plot_type, plot_type.replace("_", " ").title())


def build_canonical_response_markdown(
    active_dataset: dict[str, Any],
    plot_payload: dict[str, Any],
    output_markdown: str,
) -> str | None:
    plot_type = str(plot_payload.get("plot_type") or "").strip().lower()
    if not plot_type or not output_markdown:
        return None

    display_plot_type = _display_name(plot_type, plot_payload.get("display_plot_type"))
    resolved_genes = [str(item).strip() for item in (plot_payload.get("resolved_genes") or []) if str(item).strip()]
    resolved_groupby = str(plot_payload.get("resolved_groupby") or "").strip()
    resolved_coloring_label = str(plot_payload.get("resolved_coloring_label") or "").strip()
    rank_genes_groups_notice = str(plot_payload.get("rank_genes_groups_notice") or "").strip()

    parts: list[str] = [display_plot_type]
    if resolved_coloring_label:
        parts.append(resolved_coloring_label)
    elif resolved_groupby:
        parts.append(resolved_groupby)
    if resolved_genes:
        parts.append(_join_tokens(resolved_genes[:5]))
    if rank_genes_groups_notice:
        parts.append(rank_genes_groups_notice)

    summary = " • ".join(p for p in parts if p)
    return f"**{summary}**\n\n{output_markdown}"
