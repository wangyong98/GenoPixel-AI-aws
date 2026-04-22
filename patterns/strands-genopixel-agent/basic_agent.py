"""
GenoPixel AI — single-cell genomics analysis agent.

On every invocation the agent checks DynamoDB for the user's active dataset
selection (set via POST /api/catalog/{row}/analyze) and pre-loads the h5ad file
if it differs from what is currently in RUNTIME_STATE.  This gives users a
seamless experience: select a dataset in the browser → open chat → ask questions
immediately without a "load dataset" step.
"""

from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path

import boto3
from bedrock_agentcore.memory.integrations.strands.config import AgentCoreMemoryConfig
from bedrock_agentcore.memory.integrations.strands.session_manager import (
    AgentCoreMemorySessionManager,
)
from bedrock_agentcore.runtime import BedrockAgentCoreApp, RequestContext
from strands import Agent
from strands.models import BedrockModel
from tools.gateway import create_gateway_mcp_client
from utils.auth import extract_user_id_from_context

from gp_runtime_state import RUNTIME_STATE  # type: ignore[import-untyped]
from tools.gp_tools import ALL_GENOPIXEL_TOOLS  # noqa: F401 — imported for side-effect (tool registration)

logger = logging.getLogger(__name__)

_SKILLS_DIR = Path(__file__).parent / "skills"
_FRONTMATTER_RE = re.compile(r"^---\s*\n.*?\n---\s*\n", re.DOTALL)


def _load_skills() -> str:
    """Load all SKILL.md files in priority order and return concatenated content."""
    order = ["genopixel-tool-usage", "genopixel-plot-formatting", "scanpy-single-cell-analysis"]
    chunks: list[str] = []
    for name in order:
        path = _SKILLS_DIR / name / "SKILL.md"
        if path.exists():
            raw = path.read_text(encoding="utf-8")
            content = _FRONTMATTER_RE.sub("", raw).strip()
            if content:
                chunks.append(content)
    return "\n\n".join(chunks)


app = BedrockAgentCoreApp()
DEFAULT_BEDROCK_MODEL_ID = "us.anthropic.claude-haiku-4-5-20251001-v1:0"

SYSTEM_PROMPT = """You are GenoPixel AI, an expert assistant for single-cell RNA-seq (scRNA-seq) analysis.
You help researchers explore and visualize h5ad datasets (AnnData format) using tools backed by Scanpy.

## Your capabilities
- **Dataset status**: call get_active_dataset_info() at the start of every conversation to know what is loaded.
- **Visualizations**: use GenoPixel-compatible tools first (generate_umap_plot, generate_tsne_plot, generate_violin_plot, generate_dotplot_plot, generate_heatmap_plot, cell_counts_barplot, cell_type_proportion_barplot).
- **Exploration**: list obs columns, get unique values per column (cell types, diseases, tissues, etc.).
- **Marker genes**: set default session markers with set_markers(), then omit genes in subsequent plot calls.

## Workflow rules
1. **Always call get_active_dataset_info() first** so you know whether a dataset is already loaded.
2. If get_active_dataset_info returns `loaded=false` AND includes a `pending_selection` field, **immediately call load_dataset** using the `primary_file`, `all_excel_row`, and `title` from that field — do NOT ask the user to reselect.
3. If no dataset is loaded and there is no pending selection, tell the user to select one in the **Datasets** tab of the browser, then come back to chat.
3. When the user asks for a plot, call the appropriate tool and include inline_markdown (or markdown image) from the tool response verbatim.
4. After generating a plot, briefly interpret it: notable clusters, top cell types, expression patterns.
5. If a tool returns an error (e.g. gene not found), explain what happened and suggest alternatives.

## Gene name tips
- Use official HGNC symbols (e.g. CD3E, not "CD3 epsilon").
- The dataset may use Ensembl IDs in adata.var — check get_obs_columns() if genes are not found.

## Plot parameter guidance
- `color_by` for UMAP/tSNE: obs column such as "author_cell_type", "disease", "tissue_type"
- `groupby` for violin/dotplot/heatmap: same obs columns
- `genes`: comma-separated, e.g. "CD3E,CD4,CD8A"
"""


def _create_session_manager(user_id: str, session_id: str) -> AgentCoreMemorySessionManager:
    memory_id = os.environ.get("MEMORY_ID")
    if not memory_id:
        raise ValueError("MEMORY_ID environment variable is required")
    config = AgentCoreMemoryConfig(
        memory_id=memory_id, session_id=session_id, actor_id=user_id
    )
    return AgentCoreMemorySessionManager(
        agentcore_memory_config=config,
        region_name=os.environ.get("AWS_DEFAULT_REGION", "us-east-1"),
    )


def _is_truthy(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes", "on"}


def _is_under_base(path: Path, base: Path) -> bool:
    try:
        path.resolve().relative_to(base.resolve())
        return True
    except Exception:
        return False


def _parse_int_env(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError:
        logger.warning("Invalid %s=%r; using default=%d", name, raw, default)
        return default


def _try_preload_active_dataset(user_id: str) -> str:
    """
    Look up the user's active-dataset selection in DynamoDB and pre-load the h5ad
    if it isn't already in RUNTIME_STATE.

    Returns a short status string appended to the system prompt so the agent
    knows the current dataset context without calling get_active_dataset_info().
    """
    table_name = os.environ.get("ACTIVE_DATASET_TABLE", "")
    if not table_name:
        return ""

    try:
        ddb = boto3.resource("dynamodb", region_name=os.environ.get("AWS_DEFAULT_REGION", "us-east-1"))
        table = ddb.Table(table_name)
        resp = table.get_item(Key={"sessionId": user_id})
        item = resp.get("Item")
    except Exception as exc:
        logger.warning("Could not read active-dataset from DynamoDB: %s", exc)
        return ""

    if not item:
        return ""

    all_excel_row = item.get("all_excel_row")
    primary_file = item.get("primary_file", "")
    title = item.get("title", f"Dataset {all_excel_row}")
    multiple_excel_row = item.get("multiple_excel_row")

    # Always record the selection so get_active_dataset_info() can return it
    # even when the h5ad hasn't been pre-loaded into memory yet.
    RUNTIME_STATE.set_pending_selection(
        all_excel_row=all_excel_row,
        multiple_excel_row=multiple_excel_row,
        title=title,
        primary_file=primary_file,
    )

    # Check if RUNTIME_STATE already has this dataset loaded
    payload = RUNTIME_STATE.get_active_dataset_payload()
    already_loaded = (
        payload.get("loaded")
        and payload.get("all_excel_row") == all_excel_row
        and payload.get("multiple_excel_row") == (multiple_excel_row or None)
    )
    if already_loaded:
        return f"\n\n[Context: Dataset already loaded — {title}, {payload.get('total_cells'):,} cells]"

    if not primary_file:
        return (
            f"\n\n[Context: User has selected '{title}' (row {all_excel_row}) "
            "but the h5ad filename is not recorded. Call load_dataset with the correct filename.]"
        )

    # Attempt pre-load
    logger.info("Pre-loading active dataset: row=%s file=%s", all_excel_row, primary_file)
    try:
        from tools.gp_tools import _resolve_h5ad_path, H5AD_BASE_DIR
        path = _resolve_h5ad_path(primary_file)
        resolved_path = Path(path).resolve()
        is_efs_backed = _is_under_base(resolved_path, H5AD_BASE_DIR)

        # Keep chat responsive: avoid auto-preload from S3 fallback unless explicitly enabled.
        if not is_efs_backed and not _is_truthy(os.environ.get("ALLOW_S3_PRELOAD")):
            logger.info("Skipping auto-preload for non-EFS path: %s", resolved_path)
            return (
                f"\n\n[Context: User selected '{title}' (row {all_excel_row}, file '{primary_file}') "
                "but auto-preload is skipped on S3 fallback storage to keep chat responsive. "
                "Call load_dataset if analysis is needed now.]"
            )

        # Optional safeguard for very large files.
        max_preload_bytes = _parse_int_env("PRELOAD_MAX_BYTES", 1_000_000_000)
        try:
            file_size = resolved_path.stat().st_size
        except OSError:
            file_size = None
        if file_size is not None and file_size > max_preload_bytes:
            logger.info(
                "Skipping auto-preload for large file (%s bytes > %s): %s",
                file_size,
                max_preload_bytes,
                resolved_path,
            )
            return (
                f"\n\n[Context: User selected '{title}' (row {all_excel_row}, file '{primary_file}') "
                f"but auto-preload is skipped because file size is {file_size:,} bytes. "
                "Call load_dataset to load it on demand.]"
            )

        backed = is_efs_backed or _is_truthy(os.environ.get("S3_PRELOAD_BACKED", "true"))
        RUNTIME_STATE.load_active_dataset(
            h5ad_path=str(resolved_path),
            all_excel_row=all_excel_row,
            multiple_excel_row=multiple_excel_row,
            title=title,
            backed=backed,
            force_reload=False,
        )
        payload = RUNTIME_STATE.get_active_dataset_payload()
        return (
            f"\n\n[Context: Dataset pre-loaded — '{title}', "
            f"{payload.get('total_cells'):,} cells, backed={backed}]"
        )
    except FileNotFoundError:
        return (
            f"\n\n[Context: User selected '{title}' (row {all_excel_row}, file '{primary_file}') "
            "but the file is not yet available on EFS or S3. "
            "Tell the user the data files need to be uploaded first.]"
        )
    except Exception as exc:
        logger.warning("Pre-load failed for %s: %s", primary_file, exc)
        return (
            f"\n\n[Context: Pre-load of '{title}' failed ({exc}). "
            "Call load_dataset to try again manually.]"
        )


_SKILL_CONTENT = _load_skills()


def _create_agent(user_id: str, session_id: str, dataset_context: str) -> Agent:
    model_id = os.environ.get("BEDROCK_MODEL_ID", DEFAULT_BEDROCK_MODEL_ID)
    bedrock_model = BedrockModel(
        model_id=model_id, temperature=0.1
    )

    session_manager = _create_session_manager(user_id, session_id)
    gateway_client = create_gateway_mcp_client()

    full_prompt = SYSTEM_PROMPT
    if _SKILL_CONTENT:
        full_prompt += "\n\n" + _SKILL_CONTENT
    full_prompt += dataset_context

    return Agent(
        name="genopixel_agent",
        system_prompt=full_prompt,
        tools=[gateway_client, *ALL_GENOPIXEL_TOOLS],
        model=bedrock_model,
        session_manager=session_manager,
        trace_attributes={"user.id": user_id, "session.id": session_id},
    )


@app.entrypoint
async def invocations(payload: dict, context: RequestContext):
    user_query = payload.get("prompt")
    session_id = payload.get("runtimeSessionId")

    if not all([user_query, session_id]):
        yield {"status": "error", "error": "Missing required fields: prompt or runtimeSessionId"}
        return

    try:
        user_id = extract_user_id_from_context(context)
    except Exception as exc:
        logger.warning("Could not extract user_id from JWT: %s", exc)
        user_id = session_id  # fallback

    # Pre-load the user's active dataset (non-blocking if unavailable)
    dataset_context = _try_preload_active_dataset(user_id)

    try:
        agent = _create_agent(user_id, session_id, dataset_context)
        async for event in agent.stream_async(user_query):
            yield json.loads(json.dumps(dict(event), default=str))
    except Exception as exc:
        logger.exception("Agent run failed")
        yield {"status": "error", "error": str(exc)}


if __name__ == "__main__":
    app.run()
