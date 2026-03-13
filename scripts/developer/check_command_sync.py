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

    core_bindings_text = (root / "bindings" / "CoreBindings.cpp").read_text(encoding="utf-8")

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
        binding_unit = f"{command.removesuffix('Command')}Bindings.cpp"
        expected_paths = command_surface_paths(root, command)
        for p in expected_paths:
            if not p.exists():
                errors.append(f"missing file for {command}: {p.relative_to(root)}")

        binding_path = root / "bindings" / binding_unit
        if not binding_path.exists():
            continue
        binding_text = binding_path.read_text(encoding="utf-8")
        bind_template_ref = f"BindCommandClass<{command}>"
        if bind_template_ref not in binding_text:
            errors.append(f"{binding_unit} missing command bind template reference: {bind_template_ref}")
        bind_registration_ref = f"BindCommand<{command}>"
        if bind_registration_ref not in binding_text:
            errors.append(
                f"{binding_unit} missing command registration specialization: {bind_registration_ref}"
            )

    command_ids = [entry.command_id for entry in commands]
    if len(command_ids) != len(set(command_ids)):
        errors.append("CommandList.def has duplicated command ids.")
    cli_names = [entry.cli_name for entry in commands]
    if len(cli_names) != len(set(cli_names)):
        errors.append("CommandList.def has duplicated CLI names.")
    python_names = [entry.python_binding_name for entry in commands]
    if len(python_names) != len(set(python_names)):
        errors.append("CommandList.def has duplicated Python binding names.")
    for entry in commands:
        if entry.profile not in ALLOWED_PROFILES:
            errors.append(
                f"{entry.command_type} has unsupported profile '{entry.profile}'. "
                "Expected FileWorkflow or DatabaseWorkflow."
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
    raise SystemExit(main())
