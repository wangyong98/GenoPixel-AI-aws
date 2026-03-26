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

    dataset_title = str(active_dataset.get("title") or "Loaded dataset").strip()
    total_cells = active_dataset.get("total_cells")
    display_plot_type = _display_name(plot_type, plot_payload.get("display_plot_type"))
    resolved_genes = [str(item).strip() for item in (plot_payload.get("resolved_genes") or []) if str(item).strip()]
    resolved_groupby = str(plot_payload.get("resolved_groupby") or "").strip()
    resolved_coloring_label = str(plot_payload.get("resolved_coloring_label") or "").strip()
    rank_genes_groups_notice = str(plot_payload.get("rank_genes_groups_notice") or "").strip()

    lines = [f"I've successfully generated a {display_plot_type.lower()} for your dataset!", "", f"**Dataset:** {dataset_title}"]
    if total_cells is not None:
        lines.append(f"- **Total cells:** {_format_int(total_cells)} cells")

    if plot_type in EMBEDDING_PLOT_TYPES:
        lines.append(f"- **Visualization type:** {display_plot_type}")
        if resolved_coloring_label:
            lines.append(f"- **Coloring:** {resolved_coloring_label}")
        if resolved_genes:
            lines.append(f"- **Genes:** {_join_tokens(resolved_genes)}")
    elif plot_type in DISTRIBUTION_PLOT_TYPES:
        lines.append(f"- **Plot type:** {display_plot_type}")
        if resolved_genes:
            lines.append(f"- **Genes:** {_join_tokens(resolved_genes)}")
        if resolved_groupby:
            lines.append(f"- **Grouping column:** {resolved_groupby}")
    elif plot_type in RANK_GENES_PLOT_TYPES:
        lines.append(f"- **Plot type:** {display_plot_type}")
        if resolved_groupby:
            lines.append(f"- **Ranking basis:** {resolved_groupby}")
        if rank_genes_groups_notice:
            lines.append(f"- **Rank-genes status:** {rank_genes_groups_notice}")
    elif plot_type == "correlation_matrix":
        lines.append(f"- **Plot type:** {display_plot_type}")
        if resolved_groupby:
            lines.append(f"- **Grouping column:** {resolved_groupby}")
        if resolved_genes:
            lines.append(f"- **Features:** {_join_tokens(resolved_genes)}")
    elif plot_type == "cell_counts_barplot":
        lines.append(f"- **Plot type:** {display_plot_type}")
        if resolved_groupby:
            lines.append(f"- **Category column:** {resolved_groupby}")
    else:
        lines.append(f"- **Plot type:** {display_plot_type}")
        if resolved_genes:
            lines.append(f"- **Features:** {_join_tokens(resolved_genes)}")
        if resolved_groupby:
            lines.append(f"- **Grouping column:** {resolved_groupby}")
        if resolved_coloring_label:
            lines.append(f"- **Coloring:** {resolved_coloring_label}")

    lines.extend(["", "**View the plot here:**", "", output_markdown])
    return "\n".join(lines)
