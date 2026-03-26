---
name: genopixel-plot-formatting
description: Use for formatting GenoPixels plotting responses in Open WebUI. This skill tells the agent to copy canonical GenoPixel plot blocks exactly, fall back to the returned inline plot markdown when needed, and keep any extra prose minimal.
---

# GenoPixel Plot Formatting

Follow these rules after a GenoPixel plotting tool succeeds:

- If the tool response contains `canonical_response_markdown`, copy that field exactly.
- Do not rewrite numbers, titles, genes, grouping labels, or other metadata inside `canonical_response_markdown`.
- If `canonical_response_markdown` is absent but `output_markdown` is present, include `output_markdown` exactly so the plot renders inline.
- Never replace a public asset URL with a local filesystem path.
- Keep any additional prose minimal.
- Do not restate plot metadata differently from the tool response.
- If the tool fails, report the returned error faithfully and do not improvise biological conclusions.
