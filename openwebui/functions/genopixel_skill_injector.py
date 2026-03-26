"""
title: GenoPixel Skill Injector
author: OpenAI
description: Inject mirrored GenoPixel skill artifacts and the live GenoPixel tool manifest into GenoPixels chats.
required_open_webui_version: 0.6.0
version: 1.0.0
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from open_webui.models.skills import Skills
from open_webui.utils.misc import is_string_allowed
from open_webui.utils.tools import get_tool_servers


class Valves(BaseModel):
    priority: int = Field(default=100, description="Run after generic filters and before the model call.")
    target_model_names: str = Field(default="GenoPixels", description="Comma-separated model display names.")
    target_model_ids: str = Field(default="", description="Comma-separated model ids.")
    required_skill_ids: str = Field(
        default="genopixel-tool-usage,genopixel-plot-formatting",
        description="Comma-separated skill ids injected for every GenoPixels chat.",
    )
    optional_skill_ids: str = Field(
        default="",
        description="Comma-separated secondary skill ids injected after the required skills.",
    )
    tool_server_name_match: str = Field(
        default="genopixel",
        description="Only include attached tool servers whose id, title, or description matches this token.",
    )


valves = Valves()


class Filter:
    name = "GenoPixel Skill Injector"

    async def inlet(
        self,
        body: dict[str, Any],
        __model__=None,
        __request__=None,
        __user__=None,
        __metadata__=None,
    ) -> dict[str, Any]:
        if not isinstance(body, dict):
            return body
        if not self._is_target_model(body, __model__):
            return body
        if __request__ is None:
            return body

        ordered_skill_ids = self._ordered_skill_ids()
        skill_entries = self._load_skill_entries(ordered_skill_ids)
        tool_manifest = await self._build_tool_manifest(__model__, __request__)
        runtime_context = self._build_runtime_context(skill_entries, tool_manifest)
        if not runtime_context:
            return body

        body["messages"] = self._upsert_runtime_context_message(body.get("messages"), runtime_context)
        return body

    @staticmethod
    def _csv_values(value: str) -> list[str]:
        return [item.strip() for item in str(value or "").split(",") if item.strip()]

    def _ordered_skill_ids(self) -> list[str]:
        ordered: list[str] = []
        for source in (self._csv_values(valves.required_skill_ids), self._csv_values(valves.optional_skill_ids)):
            for skill_id in source:
                if skill_id not in ordered:
                    ordered.append(skill_id)
        return ordered

    def _is_target_model(self, body: dict[str, Any], model: Any) -> bool:
        model_id = str(body.get("model") or self._model_value(model, "id") or "").strip()
        model_name = str(self._model_value(model, "name") or "").strip()

        target_ids = set(self._csv_values(valves.target_model_ids))
        target_names = {value.lower() for value in self._csv_values(valves.target_model_names)}

        if model_id and model_id in target_ids:
            return True
        if model_name and model_name.lower() in target_names:
            return True
        return False

    @staticmethod
    def _model_meta(model: Any) -> dict[str, Any]:
        if not isinstance(model, dict):
            return {}
        if isinstance(model.get("info"), dict) and isinstance(model["info"].get("meta"), dict):
            return model["info"]["meta"]
        if isinstance(model.get("meta"), dict):
            return model["meta"]
        return {}

    def _model_value(self, model: Any, key: str) -> Any:
        if not isinstance(model, dict):
            return None
        if key in model:
            return model.get(key)
        if isinstance(model.get("info"), dict) and key in model["info"]:
            return model["info"].get(key)
        meta = self._model_meta(model)
        return meta.get(key)

    def _load_skill_entries(self, ordered_skill_ids: list[str]) -> list[dict[str, str]]:
        entries: list[dict[str, str]] = []
        for skill_id in ordered_skill_ids:
            skill = Skills.get_skill_by_id(skill_id)
            if skill is None or not skill.is_active:
                continue
            content = str(skill.content or "").strip()
            if not content:
                continue
            entries.append(
                {
                    "id": skill.id,
                    "name": skill.name,
                    "description": str(skill.description or "").strip(),
                    "content": content,
                }
            )
        return entries

    async def _build_tool_manifest(self, model: Any, request: Any) -> str:
        meta = self._model_meta(model)
        tool_ids = meta.get("toolIds") or []
        if not isinstance(tool_ids, list):
            return ""

        attached_server_ids: list[tuple[str, set[str] | None]] = []
        for tool_id in tool_ids:
            token = str(tool_id or "").strip()
            if not token.startswith("server:"):
                continue
            _, _, remainder = token.partition(":")
            server_id, _, function_blob = remainder.partition("|")
            function_names = {name.strip() for name in function_blob.split(",") if name.strip()} if function_blob else None
            attached_server_ids.append((server_id.strip(), function_names))

        if not attached_server_ids:
            return ""

        servers = await get_tool_servers(request)
        servers_by_id = {str(server.get("id")): server for server in servers}
        manifest_lines = ["<genopixel_tool_manifest>"]
        matched_any_server = False

        for server_id, explicit_function_names in attached_server_ids:
            server = servers_by_id.get(server_id)
            if not isinstance(server, dict):
                continue

            title = str(server.get("openapi", {}).get("info", {}).get("title") or server.get("info", {}).get("title") or server_id)
            description = str(server.get("openapi", {}).get("info", {}).get("description") or "").strip()
            match_tokens = [title.lower(), description.lower(), str(server_id).lower()]
            matcher = str(valves.tool_server_name_match or "").strip().lower()
            if matcher and not any(matcher in token for token in match_tokens if token):
                continue

            matched_any_server = True
            manifest_lines.append(f"server: {title}")

            connection = request.app.state.config.TOOL_SERVER_CONNECTIONS[server.get("idx", 0)]
            connection_filters = connection.get("config", {}).get("function_name_filter_list", "") or ""
            if isinstance(connection_filters, str):
                connection_filters = [item.strip() for item in connection_filters.split(",") if item.strip()]

            visible_specs: list[dict[str, Any]] = []
            for spec in server.get("specs", []):
                function_name = str(spec.get("name") or "").strip()
                if not function_name:
                    continue
                if explicit_function_names and function_name not in explicit_function_names:
                    continue
                if connection_filters and not is_string_allowed(function_name, connection_filters):
                    continue
                visible_specs.append(spec)

            if not visible_specs:
                manifest_lines.append("- (no available tools)")
                continue

            for spec in visible_specs:
                function_name = str(spec.get("name") or "").strip()
                function_description = " ".join(str(spec.get("description") or "").split())
                if not function_description:
                    function_description = "No description provided."
                manifest_lines.append(f"- {function_name}: {function_description}")

        if not matched_any_server:
            manifest_lines.append("server: genopixel-tools")
            manifest_lines.append("- (no available tools)")

        manifest_lines.append("</genopixel_tool_manifest>")
        return "\n".join(manifest_lines)

    @staticmethod
    def _build_runtime_context(skill_entries: list[dict[str, str]], tool_manifest: str) -> str:
        lines = ["<genopixel_runtime_context>"]
        if tool_manifest:
            lines.append(tool_manifest)
        if skill_entries:
            lines.append("<genopixel_runtime_skills>")
            for skill in skill_entries:
                lines.append(f'<skill id="{skill["id"]}" name="{skill["name"]}">')
                lines.append(skill["content"])
                lines.append("</skill>")
            lines.append("</genopixel_runtime_skills>")
        lines.append("</genopixel_runtime_context>")
        content = "\n".join(lines).strip()
        return content if skill_entries or tool_manifest else ""

    @staticmethod
    def _upsert_runtime_context_message(messages: Any, runtime_context: str) -> list[dict[str, Any]]:
        marker = "<genopixel_runtime_context>"
        normalized: list[dict[str, Any]] = []
        if isinstance(messages, list):
            for message in messages:
                if isinstance(message, dict):
                    normalized.append(dict(message))

        for message in normalized:
            if message.get("role") != "system":
                continue
            content = message.get("content")
            if isinstance(content, str) and marker in content:
                message["content"] = runtime_context
                return normalized

        normalized.insert(0, {"role": "system", "content": runtime_context})
        return normalized
