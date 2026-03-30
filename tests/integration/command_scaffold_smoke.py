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
        project_root / "src" / "core" / "internal" / "command" / "ExampleCommand.hpp",
        project_root / "src" / "core" / "command" / "ExampleCommand.cpp",
        project_root / "tests" / "core" / "command" / "ExampleCommand_test.cpp",
        project_root / "docs" / "developer" / "commands" / "example.md",
    ]

    for expected_path in expected_paths:
        if str(expected_path) not in stdout:
            raise AssertionError(f"Expected dry-run output to mention {expected_path}.\n{stdout}")

    unexpected_fragments = [
        str(project_root / "resources" / "src"),
        str(project_root / "resources" / "tests"),
        str(project_root / "resources" / "docs"),
    ]
    for unexpected_fragment in unexpected_fragments:
        if unexpected_fragment in stdout:
            raise AssertionError(f"Dry-run output used resources/ as repo root.\n{stdout}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
