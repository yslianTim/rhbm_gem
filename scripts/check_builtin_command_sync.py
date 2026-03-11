#!/usr/bin/env python3
"""Validate built-in command manifest synchronization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys


@dataclass(frozen=True)
class CommandEntry:
    command_type: str
    cli_name: str
    description: str
    python_binding_name: str


def parse_builtins(path: Path) -> list[CommandEntry]:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r"RHBM_GEM_BUILTIN_COMMAND\(\s*"
        r"(?P<command>[A-Za-z_][A-Za-z0-9_]*)\s*,\s*"
        r"\"(?P<cli>[^\"]+)\"\s*,\s*"
        r"\"(?P<desc>[^\"]*)\"\s*,\s*"
        r"\"(?P<py>[^\"]+)\"\s*"
        r"\)",
        re.MULTILINE,
    )
    return [
        CommandEntry(
            command_type=m.group("command"),
            cli_name=m.group("cli"),
            description=m.group("desc"),
            python_binding_name=m.group("py"),
        )
        for m in pattern.finditer(text)
    ]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    builtins = parse_builtins(root / "src" / "core" / "internal" / "BuiltInCommandList.def")
    if not builtins:
        print("No built-in commands parsed from BuiltInCommandList.def")
        return 1

    catalog_text = (root / "src" / "core" / "command" / "BuiltInCommandCatalog.cpp").read_text(encoding="utf-8")
    sources_text = (root / "src" / "rhbm_gem_sources.cmake").read_text(encoding="utf-8")
    core_tests_text = (root / "tests" / "cmake" / "core_tests.cmake").read_text(encoding="utf-8")

    errors: list[str] = []
    for entry in builtins:
        command = entry.command_type
        stem = command.removesuffix("Command")
        expected_paths = [
            root / "include" / "rhbm_gem" / "core" / "command" / f"{command}.hpp",
            root / "src" / "core" / "command" / f"{command}.cpp",
            root / "bindings" / f"{stem}Bindings.cpp",
            root / "tests" / "core" / "command" / f"{command}_test.cpp",
        ]
        for p in expected_paths:
            if not p.exists():
                errors.append(f"missing file for {command}: {p.relative_to(root)}")

        header_include = f"#include <rhbm_gem/core/command/{command}.hpp>"
        if header_include not in catalog_text:
            errors.append(f"BuiltInCommandCatalog.cpp missing include for {command}")

        source_ref = f"core/command/{command}.cpp"
        if source_ref not in sources_text:
            errors.append(f"rhbm_gem_sources.cmake missing source entry: {source_ref}")

        test_ref = f"core/command/{command}_test.cpp"
        if test_ref not in core_tests_text:
            errors.append(f"core_tests.cmake missing test entry: {test_ref}")

    if errors:
        print("Built-in command sync check failed:")
        for err in errors:
            print(f" - {err}")
        return 1

    print(f"Built-in command sync check passed ({len(builtins)} commands).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
