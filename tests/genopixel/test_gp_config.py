from __future__ import annotations

import gp_config


def test_load_settings_defaults_to_non_backed(monkeypatch) -> None:
    monkeypatch.delenv('DEFAULT_BACKED', raising=False)

    settings = gp_config.load_settings()

    assert settings.default_backed is False
