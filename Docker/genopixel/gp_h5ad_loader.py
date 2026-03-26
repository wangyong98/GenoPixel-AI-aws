from __future__ import annotations

from pathlib import Path

import anndata as ad


def load_h5ad(path: Path, backed: bool = False) -> ad.AnnData:
    if not path.exists():
        raise FileNotFoundError(f"h5ad file does not exist: {path}")

    if backed:
        return ad.read_h5ad(path, backed="r")
    return ad.read_h5ad(path)
