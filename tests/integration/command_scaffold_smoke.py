#!/usr/bin/env python3

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    project_root = Path(sys.argv[1]).resolve()
    scaffold_script = project_root / "resources" / "tools" / "developer" / "command_scaffold.py"

    result = subprocess.run(
        [
            sys.executable,
            str(scaffold_script),
            "--name",
            "Example",
            "--dry-run",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    stdout = result.stdout
    expected_paths = [
        project_root / "src" / "core" / "command" / "ExampleCommand.cpp",
        project_root / "docs" / "developer" / "commands" / "example.md",
    ]

    for expected_path in expected_paths:
        if str(expected_path) not in stdout:
            raise AssertionError(f"Expected dry-run output to mention {expected_path}.\n{stdout}")

    unexpected_fragments = [
        str(project_root / "resources" / "src"),
        str(project_root / "resources" / "tests"),
        str(project_root / "resources" / "docs"),
        "src/core/" + "internal" + "/command",
        "internal" + "/command/ExampleCommand.hpp",
        str(project_root / "src" / "core" / "command" / "ExampleCommand.hpp"),
    ]
    for unexpected_fragment in unexpected_fragments:
        if unexpected_fragment in stdout:
            raise AssertionError(f"Dry-run output used resources/ as repo root.\n{stdout}")

    guidance_fragments = [
        str(project_root / "tests" / "core" / "command"),
        str(project_root / "tests" / "core" / "command" / "CommandValidationScenarios_test.cpp"),
        str(project_root / "tests" / "core" / "command" / "CommandWorkflowScenarios_test.cpp"),
        "extend the grouped command tests",
    ]
    for fragment in guidance_fragments:
        if fragment not in stdout:
            raise AssertionError(f"Expected dry-run output to mention {fragment}.\n{stdout}")

    scaffold_source = scaffold_script.read_text(encoding="utf-8")
    if "NormalizeAndValidateRequest" not in scaffold_source:
        raise AssertionError("Scaffold template should use NormalizeAndValidateRequest().")
    if "ValidatePreparedRequest" not in scaffold_source:
        raise AssertionError("Scaffold template should use ValidatePreparedRequest().")
    if "CommandExecutor.hpp" not in scaffold_source:
        raise AssertionError("Scaffold template should use CommandExecutor.hpp.")
    removed_runtime_hook_name = "Reset" + "Runtime" + "State"
    if removed_runtime_hook_name in scaffold_source:
        raise AssertionError(
            f"Scaffold template still references {removed_runtime_hook_name}().")
    old_hook_name = "Normalize" + "Request"
    if old_hook_name in scaffold_source:
        raise AssertionError(f"Scaffold template still references {old_hook_name}().")
    old_validation_hook_name = "Validate" + "Options"
    if old_validation_hook_name in scaffold_source:
        raise AssertionError(
            f"Scaffold template still references {old_validation_hook_name}().")
    removed_header_suffix = "Command" + ".hpp"
    if removed_header_suffix in scaffold_source:
        raise AssertionError(f"Scaffold template still references per-command {removed_header_suffix}.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
