#!/usr/bin/env python3
"""Shared command manifest helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class CommandEntry:
    command_id: str
    command_type: str
    cli_name: str
    description: str
    profile: str


COMMAND_PATTERN = re.compile(
    r"RHBM_GEM_COMMAND\(\s*"
    r"(?P<id>[A-Za-z_][A-Za-z0-9_]*)\s*,\s*"
    r"(?P<command>[A-Za-z_][A-Za-z0-9_]*)\s*,\s*"
    r"\"(?P<cli>[^\"]+)\"\s*,\s*"
    r"\"(?P<desc>[^\"]*)\"\s*,\s*"
    r"(?P<profile>[A-Za-z_][A-Za-z0-9_]*)\s*"
    r"\)",
    re.MULTILINE,
)


def parse_commands(path: Path) -> list[CommandEntry]:
    text = path.read_text(encoding="utf-8")
    return [
        CommandEntry(
            command_id=match.group("id"),
            command_type=match.group("command"),
            cli_name=match.group("cli"),
            description=match.group("desc"),
            profile=match.group("profile"),
        )
        for match in COMMAND_PATTERN.finditer(text)
    ]


def command_surface_paths(root: Path, command_type: str) -> list[Path]:
    return [
        root / "src" / "core" / "command" / f"{command_type}.hpp",
        root / "src" / "core" / "command" / f"{command_type}.cpp",
        root / "tests" / "core" / "command" / f"{command_type}_test.cpp",
    ]


def collect_drift_paths(outputs: list[tuple[Path, str]]) -> list[Path]:
    drift_paths: list[Path] = []
    for path, expected_content in outputs:
        current_content = path.read_text(encoding="utf-8") if path.exists() else ""
        if current_content != expected_content:
            drift_paths.append(path)
    return drift_paths
