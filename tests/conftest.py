from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
GENOPIXEL_DIR = ROOT / "Docker" / "genopixel"

if str(GENOPIXEL_DIR) not in sys.path:
    sys.path.insert(0, str(GENOPIXEL_DIR))

from gp_runtime_state import RUNTIME_STATE


@pytest.fixture(autouse=True)
def clear_runtime_state():
    RUNTIME_STATE.clear_active_dataset()
    yield
    RUNTIME_STATE.clear_active_dataset()
