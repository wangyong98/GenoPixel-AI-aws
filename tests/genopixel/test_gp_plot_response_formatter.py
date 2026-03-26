from __future__ import annotations

from gp_plot_response_formatter import build_canonical_response_markdown


def test_build_canonical_response_markdown_for_umap() -> None:
    markdown = build_canonical_response_markdown(
        {
            "title": "Dataset A",
            "total_cells": 32458,
        },
        {
            "plot_type": "umap",
            "display_plot_type": "UMAP embedding",
            "resolved_coloring_label": "Cell types (author_cell_type1)",
            "resolved_genes": [],
        },
        "![Plot](http://localhost/assets/genopixel/plot.png)",
    )

    assert markdown is not None
    assert "I've successfully generated a umap embedding" in markdown.lower()
    assert "- **Total cells:** 32,458 cells" in markdown
    assert "- **Visualization type:** UMAP embedding" in markdown
    assert "- **Coloring:** Cell types (author_cell_type1)" in markdown


def test_build_canonical_response_markdown_for_violin() -> None:
    markdown = build_canonical_response_markdown(
        {
            "title": "Dataset B",
            "total_cells": 6789,
        },
        {
            "plot_type": "violin",
            "display_plot_type": "Violin plot",
            "resolved_genes": ["OSMR"],
            "resolved_groupby": "author_cell_type1",
        },
        "![Plot](http://localhost/assets/genopixel/violin.png)",
    )

    assert markdown is not None
    assert "- **Plot type:** Violin plot" in markdown
    assert "- **Genes:** OSMR" in markdown
    assert "- **Grouping column:** author_cell_type1" in markdown


def test_build_canonical_response_markdown_for_rank_genes_includes_cache_status() -> None:
    markdown = build_canonical_response_markdown(
        {
            "title": "Dataset C",
            "total_cells": 3210,
        },
        {
            "plot_type": "rank_genes_groups_dotplot",
            "display_plot_type": "Rank-genes dot plot",
            "resolved_groupby": "author_cell_type1",
            "rank_genes_groups_notice": "Reusing cached rank_genes_groups computed with groupby 'author_cell_type1'.",
        },
        "![Plot](http://localhost/assets/genopixel/rank.png)",
    )

    assert markdown is not None
    assert "- **Plot type:** Rank-genes dot plot" in markdown
    assert "- **Ranking basis:** author_cell_type1" in markdown
    assert "- **Rank-genes status:** Reusing cached rank_genes_groups computed with groupby 'author_cell_type1'." in markdown
