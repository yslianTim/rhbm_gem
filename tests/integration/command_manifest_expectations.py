from __future__ import annotations

import re
from pathlib import Path


_COMMAND_PATTERN = re.compile(
    r"RHBM_GEM_COMMAND\(\s*"
    r"(?P<id>[A-Za-z_][A-Za-z0-9_]*)\s*,\s*"
    r"(?P<command>[A-Za-z_][A-Za-z0-9_]*)\s*,"
)


def _manifest_entries() -> list[tuple[str, str]]:
    manifest_path = Path(__file__).resolve().parents[2] / "src" / "core" / "internal" / "CommandList.def"
    return [
        (match.group("id"), match.group("command"))
        for match in _COMMAND_PATTERN.finditer(manifest_path.read_text(encoding="utf-8"))
    ]


def expected_request_names() -> list[str]:
    return ["CommonCommandRequest"] + [
        f"{command.removesuffix('Command')}Request"
        for _, command in _manifest_entries()
    ]


def expected_run_functions() -> list[str]:
    return [
        f"Run{command.removesuffix('Command')}"
        for _, command in _manifest_entries()
    ]
