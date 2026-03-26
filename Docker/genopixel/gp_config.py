from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from dotenv import load_dotenv


@dataclass
class Settings:
    ollama_base_url: str
    ollama_model: str
    ollama_timeout_seconds: int
    metadata_xlsx: Path
    h5ad_base_dir: Path
    output_dir: Path
    default_backed: bool


def _as_bool(value: str, default: bool = False) -> bool:
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "y"}


def _as_int(value: str | None, default: int) -> int:
    if value is None:
        return default
    try:
        parsed = int(value)
        return parsed if parsed > 0 else default
    except ValueError:
        return default


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_output_dir(raw_output_dir: str) -> Path:
    output_dir = Path(raw_output_dir).expanduser()
    if output_dir.is_absolute():
        return output_dir
    return (_project_root() / output_dir).resolve()


def load_settings() -> Settings:
    load_dotenv()
    metadata_xlsx = Path(
        os.getenv(
            "METADATA_XLSX",
            "/Volumes/cx10/Single_cell_data_0117_2026/Final_metadata/cellxgene_HCA_final_02182026.xlsx",
        )
    )
    h5ad_base_dir = Path(
        os.getenv(
            "H5AD_BASE_DIR",
            "/Volumes/cx10/Single_cell_data_0117_2026/cellxgene_final",
        )
    )
    output_dir = _resolve_output_dir(os.getenv("OUTPUT_DIR", "outputs"))

    return Settings(
        ollama_base_url=os.getenv("OLLAMA_BASE_URL", "http://localhost:11434"),
        ollama_model=os.getenv("OLLAMA_MODEL", "qwen3-coder:30b"),
        ollama_timeout_seconds=_as_int(os.getenv("OLLAMA_TIMEOUT_SECONDS"), default=180),
        metadata_xlsx=metadata_xlsx,
        h5ad_base_dir=h5ad_base_dir,
        output_dir=output_dir,
        default_backed=_as_bool(os.getenv("DEFAULT_BACKED", "false"), default=False),
    )
