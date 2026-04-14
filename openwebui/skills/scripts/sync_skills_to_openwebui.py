#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SKILLS_ROOT = REPO_ROOT / "openwebui" / "skills"
FILTER_FILE = REPO_ROOT / "openwebui" / "functions" / "genopixel_skill_injector.py"
CORE_SKILL_IDS = ["genopixel-tool-usage", "genopixel-plot-formatting"]


def _parse_frontmatter(text: str) -> tuple[dict[str, str], str]:
    if not text.startswith("---\n"):
        raise ValueError("SKILL.md is missing frontmatter")
    _, _, remainder = text.partition("---\n")
    frontmatter_text, sep, body = remainder.partition("\n---\n")
    if not sep:
        raise ValueError("SKILL.md frontmatter terminator not found")

    data: dict[str, str] = {}
    for raw_line in frontmatter_text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, sep, value = line.partition(":")
        if not sep:
            continue
        data[key.strip()] = value.strip()
    return data, body.lstrip()


def _skill_tags(skill_dir: Path) -> list[str]:
    relative = skill_dir.relative_to(SKILLS_ROOT)
    tags = ["genopixel", "skill"]
    if relative.parts and relative.parts[0] == "public":
        tags.append("public")
    else:
        tags.append("core")
    return tags


def _load_skill_payload(skill_dir: Path) -> dict[str, Any]:
    skill_md = skill_dir / "SKILL.md"
    parsed, body = _parse_frontmatter(skill_md.read_text(encoding="utf-8"))
    skill_id = skill_dir.name
    name = parsed.get("name", skill_id)
    description = parsed.get("description", "")
    if not description:
        raise ValueError(f"Skill '{skill_id}' is missing a description")
    return {
        "id": skill_id,
        "name": name,
        "description": description,
        "content": body.strip(),
        "meta": {"tags": _skill_tags(skill_dir)},
        "is_active": True,
        "access_grants": [
            {
                "principal_type": "user",
                "principal_id": "*",
                "permission": "read",
            }
        ],
    }


def _discover_skills() -> list[dict[str, Any]]:
    skill_dirs = [path.parent for path in SKILLS_ROOT.rglob("SKILL.md") if path.parent != SKILLS_ROOT]
    seen: set[Path] = set()
    ordered_dirs: list[Path] = []
    for skill_dir in sorted(skill_dirs, key=lambda path: ("/public/" in str(path), str(path))):
        if skill_dir in seen:
            continue
        seen.add(skill_dir)
        ordered_dirs.append(skill_dir)

    skills = [_load_skill_payload(skill_dir) for skill_dir in ordered_dirs]
    missing_core = [skill_id for skill_id in CORE_SKILL_IDS if skill_id not in {skill["id"] for skill in skills}]
    if missing_core:
        raise ValueError(f"Missing required core skills: {', '.join(missing_core)}")
    return skills


def _build_payload(model_name: str, filter_id: str) -> dict[str, Any]:
    filter_content = FILTER_FILE.read_text(encoding="utf-8")
    all_skills = _discover_skills()
    optional_ids = [s["id"] for s in all_skills if s["id"] not in CORE_SKILL_IDS]
    return {
        "skills": all_skills,
        "function": {
            "id": filter_id,
            "name": "GenoPixel Skill Injector",
            "type": "filter",
            "content": filter_content,
            "meta": {
                "description": "Inject mirrored GenoPixel skill artifacts and the live GenoPixel tool manifest into GenoPixels chats.",
                "manifest": {
                    "title": "GenoPixel Skill Injector",
                    "description": "Inject mirrored GenoPixel skill artifacts and the live GenoPixel tool manifest into GenoPixels chats.",
                },
            },
            "valves": {
                "priority": 100,
                "target_model_names": model_name,
                "target_model_ids": "",
                "required_skill_ids": ",".join(CORE_SKILL_IDS),
                "optional_skill_ids": ",".join(optional_ids),
                "tool_server_name_match": "genopixel",
            },
            "is_active": True,
            "is_global": False,
        },
        "model": {
            "name": model_name,
            "filter_id": filter_id,
        },
    }


REMOTE_SYNC_SCRIPT = r'''
import json
import sqlite3
import sys
import time
import uuid

payload = json.load(sys.stdin)
conn = sqlite3.connect('/app/backend/data/webui.db')
conn.row_factory = sqlite3.Row
cur = conn.cursor()
now = int(time.time())

admin_row = cur.execute("select id from user where role = 'admin' order by created_at asc limit 1").fetchone()
if admin_row is None:
    raise SystemExit('No admin user found in Open WebUI database')
admin_user_id = admin_row['id']

for skill in payload['skills']:
    existing = cur.execute('select id, created_at from skill where id = ?', (skill['id'],)).fetchone()
    meta_json = json.dumps(skill['meta'])
    if existing is None:
        cur.execute(
            'insert into skill (id, user_id, name, description, content, meta, is_active, updated_at, created_at) values (?, ?, ?, ?, ?, ?, ?, ?, ?)',
            (
                skill['id'],
                admin_user_id,
                skill['name'],
                skill['description'],
                skill['content'],
                meta_json,
                1 if skill.get('is_active', True) else 0,
                now,
                now,
            ),
        )
    else:
        cur.execute(
            'update skill set user_id = ?, name = ?, description = ?, content = ?, meta = ?, is_active = ?, updated_at = ? where id = ?',
            (
                admin_user_id,
                skill['name'],
                skill['description'],
                skill['content'],
                meta_json,
                1 if skill.get('is_active', True) else 0,
                now,
                skill['id'],
            ),
        )

    cur.execute('delete from access_grant where resource_type = ? and resource_id = ?', ('skill', skill['id']))
    for grant in skill.get('access_grants', []):
        cur.execute(
            'insert into access_grant (id, resource_type, resource_id, principal_type, principal_id, permission, created_at) values (?, ?, ?, ?, ?, ?, ?)',
            (
                str(uuid.uuid4()),
                'skill',
                skill['id'],
                grant['principal_type'],
                grant['principal_id'],
                grant['permission'],
                now,
            ),
        )

function = payload['function']
function_existing = cur.execute('select id from function where id = ?', (function['id'],)).fetchone()
function_meta_json = json.dumps(function['meta'])
function_valves_json = json.dumps(function['valves'])
if function_existing is None:
    cur.execute(
        'insert into function (id, user_id, name, type, content, meta, created_at, updated_at, valves, is_active, is_global) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
        (
            function['id'],
            admin_user_id,
            function['name'],
            function['type'],
            function['content'],
            function_meta_json,
            now,
            now,
            function_valves_json,
            1 if function.get('is_active', True) else 0,
            1 if function.get('is_global', False) else 0,
        ),
    )
else:
    cur.execute(
        'update function set user_id = ?, name = ?, type = ?, content = ?, meta = ?, updated_at = ?, valves = ?, is_active = ?, is_global = ? where id = ?',
        (
            admin_user_id,
            function['name'],
            function['type'],
            function['content'],
            function_meta_json,
            now,
            function_valves_json,
            1 if function.get('is_active', True) else 0,
            1 if function.get('is_global', False) else 0,
            function['id'],
        ),
    )

model = payload['model']
model_row = cur.execute('select id, meta, params from model where name = ? or id = ? order by updated_at desc limit 1', (model['name'], model['name'])).fetchone()
if model_row is not None:
    meta = json.loads(model_row['meta'] or '{}')
    params = json.loads(model_row['params'] or '{}')
    filter_ids = meta.get('filterIds') or []
    if not isinstance(filter_ids, list):
        filter_ids = []
    if model['filter_id'] not in filter_ids:
        filter_ids.append(model['filter_id'])
    meta['filterIds'] = filter_ids
    meta['skillIds'] = [s['id'] for s in payload['skills']]
    params['system'] = ''
    cur.execute(
        'update model set meta = ?, params = ?, updated_at = ? where id = ?',
        (
            json.dumps(meta),
            json.dumps(params),
            now,
            model_row['id'],
        ),
    )

conn.commit()
print(json.dumps({
    'skills': [skill['id'] for skill in payload['skills']],
    'function_id': function['id'],
    'model_name': model['name'],
}))
'''

def _sync_via_docker(payload: dict[str, Any], service: str) -> None:
    cmd = ["docker", "compose", "exec", "-T", service, "python", "-c", REMOTE_SYNC_SCRIPT]
    completed = subprocess.run(
        cmd,
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=False,
        cwd=REPO_ROOT,
    )
    if completed.returncode != 0:
        if completed.stderr:
            sys.stderr.write(completed.stderr)
        raise SystemExit(completed.returncode)
    if completed.stdout.strip():
        print(completed.stdout.strip())


def main() -> None:
    parser = argparse.ArgumentParser(description="Sync GenoPixel skills and runtime filter into Open WebUI.")
    parser.add_argument("--service", default="openwebui", help="Docker Compose service name for Open WebUI.")
    parser.add_argument("--model-name", default="GenoPixel", help="Open WebUI model name to attach the filter to.")
    parser.add_argument("--filter-id", default="genopixel_skill_injector", help="Function id for the runtime injector filter.")
    args = parser.parse_args()

    payload = _build_payload(model_name=args.model_name, filter_id=args.filter_id)
    _sync_via_docker(payload, service=args.service)


if __name__ == "__main__":
    main()
