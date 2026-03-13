#!/usr/bin/env python3
"""Validate command manifest synchronization."""

from __future__ import annotations

from pathlib import Path

import generate_command_artifacts
from command_manifest import command_surface_paths, collect_drift_paths, parse_commands


ALLOWED_PROFILES = {"FileWorkflow", "DatabaseWorkflow"}


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    commands = parse_commands(root / "src" / "core" / "internal" / "CommandList.def")
    if not commands:
        print("No commands parsed from CommandList.def")
        return 1

    errors: list[str] = []
    generated_outputs = generate_command_artifacts.compute_expected_outputs(root, commands)
    generated_drift = collect_drift_paths(generated_outputs)
    if generated_drift:
        drift_lines = "\n".join(
            f" - {path.relative_to(root)}" for path in generated_drift
        )
        errors.append(
            "generated artifacts drift detected:\n"
            f"{drift_lines}\n"
            "Run: python3 scripts/developer/generate_command_artifacts.py"
        )

    for entry in commands:
        command = entry.command_type
        expected_paths = command_surface_paths(root, command)
        for p in expected_paths:
            if not p.exists():
                errors.append(f"missing file for {command}: {p.relative_to(root)}")

    command_ids = [entry.command_id for entry in commands]
    if len(command_ids) != len(set(command_ids)):
        errors.append("CommandList.def has duplicated command ids.")
    cli_names = [entry.cli_name for entry in commands]
    if len(cli_names) != len(set(cli_names)):
        errors.append("CommandList.def has duplicated CLI names.")
    for entry in commands:
        if entry.profile not in ALLOWED_PROFILES:
            errors.append(
                f"{entry.command_type} has unsupported profile '{entry.profile}'. "
                "Expected FileWorkflow or DatabaseWorkflow."
            )

    if errors:
        print("Command sync check failed:")
        for err in errors:
            print(f" - {err}")
        return 1

    print(f"Command sync check passed ({len(commands)} commands).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
