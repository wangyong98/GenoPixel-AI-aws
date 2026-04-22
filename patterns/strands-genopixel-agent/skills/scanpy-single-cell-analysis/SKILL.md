---
name: scanpy-single-cell-analysis
description: Optional secondary skill for GenoPixels that adds brief Scanpy biological interpretation hints. Must not override GenoPixel tool availability rules or plot formatting rules.
---

# Scanpy Interpretation Hints

Use this only as a secondary skill after the GenoPixel core skills.

- Add brief biological context when it helps the user understand a generated plot (e.g. what a cluster separation in UMAP means, or why a gene shows high expression in a cell type).
- Do not suggest analysis steps (QC, normalization, clustering, trajectory) that are not exposed by the current GenoPixel tool manifest.
- Do not override tool-selection rules from genopixel-tool-usage.
- Do not imply access to raw Python execution, notebook environments, or file system operations.
