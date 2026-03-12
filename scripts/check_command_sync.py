#!/usr/bin/env python3
"""Validate command manifest synchronization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
import sys


@dataclass(frozen=True)
class CommandEntry:
    command_type: str
    cli_name: str
    description: str


def parse_commands(path: Path) -> list[CommandEntry]:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r"RHBM_GEM_COMMAND\(\s*"
        r"(?P<command>[A-Za-z_][A-Za-z0-9_]*)\s*,\s*"
        r"\"(?P<cli>[^\"]+)\"\s*,\s*"
        r"\"(?P<desc>[^\"]*)\"\s*"
        r"\)",
        re.MULTILINE,
    )
    return [
        CommandEntry(
            command_type=m.group("command"),
            cli_name=m.group("cli"),
            description=m.group("desc"),
        )
        for m in pattern.finditer(text)
    ]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    commands = parse_commands(root / "src" / "core" / "internal" / "CommandList.def")
    if not commands:
        print("No commands parsed from CommandList.def")
        return 1

    core_bindings_text = (root / "bindings" / "CoreBindings.cpp").read_text(encoding="utf-8")

    errors: list[str] = []
    generator = root / "scripts" / "generate_command_artifacts.py"
    generated = subprocess.run(
        [sys.executable, str(generator), "--check"],
        check=False,
        capture_output=True,
        text=True,
    )
    if generated.returncode != 0:
        message = generated.stdout.strip()
        if generated.stderr.strip():
            message += ("\n" if message else "") + generated.stderr.strip()
        errors.append(f"generated artifacts drift detected:\n{message}")

    for entry in commands:
        command = entry.command_type
        stem = command.removesuffix("Command")
        binding_unit = f"{stem}Bindings.cpp"
        expected_paths = [
            root / "include" / "rhbm_gem" / "core" / "command" / f"{command}.hpp",
            root / "src" / "core" / "command" / f"{command}.cpp",
            root / "bindings" / binding_unit,
            root / "tests" / "core" / "command" / f"{command}_test.cpp",
        ]
        for p in expected_paths:
            if not p.exists():
                errors.append(f"missing file for {command}: {p.relative_to(root)}")

        binding_text = (root / "bindings" / binding_unit).read_text(encoding="utf-8")
        bind_template_ref = f"BindCommandClass<{command}>"
        if bind_template_ref not in binding_text:
            errors.append(f"{binding_unit} missing command bind template reference: {bind_template_ref}")
        bind_registration_ref = f"BindCommand<{command}>"
        if bind_registration_ref not in binding_text:
            errors.append(
                f"{binding_unit} missing command registration specialization: {bind_registration_ref}"
            )

    if core_bindings_text.count('#include "internal/CommandList.def"') < 1:
        errors.append("CoreBindings.cpp must include internal/CommandList.def")
    if "BindCommand<rhbm_gem::COMMAND_TYPE>" not in core_bindings_text:
        errors.append("CoreBindings.cpp must register commands via manifest macro expansion")

    if errors:
        print("Command sync check failed:")
        for err in errors:
            print(f" - {err}")
        return 1

    print(f"Command sync check passed ({len(commands)} commands).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
