from __future__ import annotations

from typing import Any
from urllib.parse import quote_plus

try:
    from starlette.responses import HTMLResponse
except Exception:  # pragma: no cover - Open WebUI runtime should provide Starlette.
    HTMLResponse = None  # type: ignore[assignment]


class Action:
    name = "Browse Datasets"
    description = "Open the GenoPixel dataset browser inside chat."
    icon = "table"
    priority = 20

    async def action(
        self,
        body: dict[str, Any],
        __event_call__=None,
        __event_emitter__=None,
        __user__=None,
    ) -> Any:
        query = self._extract_latest_user_text(body)
        iframe_url = "/apps/genopixel-datasets/"
        if query:
            iframe_url = f"{iframe_url}?q={quote_plus(query)}"

        html = self._build_embed(iframe_url)

        if __event_call__ is not None:
            try:
                return await __event_call__(
                    {
                        "type": "html",
                        "data": {
                            "content": html,
                        },
                    }
                )
            except Exception:
                pass

        if __event_emitter__ is not None:
            try:
                await __event_emitter__(
                    {
                        "type": "message",
                        "data": {
                            "content": html,
                        },
                    }
                )
                return None
            except Exception:
                pass

        if HTMLResponse is not None:
            return HTMLResponse(content=html)
        return {"type": "html", "data": {"content": html}}

    @staticmethod
    def _extract_latest_user_text(body: dict[str, Any]) -> str:
        messages = body.get("messages") or []
        for message in reversed(messages):
            if message.get("role") != "user":
                continue
            content = message.get("content")
            if isinstance(content, str) and content.strip():
                return content.strip()
            if isinstance(content, list):
                text_parts = []
                for entry in content:
                    if isinstance(entry, dict) and entry.get("type") == "text":
                        text = str(entry.get("text") or "").strip()
                        if text:
                            text_parts.append(text)
                if text_parts:
                    return " ".join(text_parts)
        return ""

    @staticmethod
    def _build_embed(iframe_url: str) -> str:
        return f"""
<div style=\"display:grid;gap:12px;margin:8px 0 4px;\">
  <div style=\"display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap;\">
    <div>
      <div style=\"font-size:12px;letter-spacing:0.18em;text-transform:uppercase;opacity:0.65;\">GenoPixel Catalog</div>
      <div style=\"font-size:16px;font-weight:600;\">Interactive dataset browser</div>
    </div>
    <a href=\"{iframe_url}\" target=\"_blank\" rel=\"noreferrer\" style=\"padding:8px 12px;border:1px solid rgba(0,0,0,0.14);border-radius:999px;text-decoration:none;color:inherit;\">Open full page</a>
  </div>
  <iframe
    src=\"{iframe_url}\"
    title=\"GenoPixel dataset browser\"
    style=\"width:100%;min-height:780px;border:1px solid rgba(0,0,0,0.08);border-radius:18px;background:#fff;\"
  ></iframe>
</div>
        """.strip()
