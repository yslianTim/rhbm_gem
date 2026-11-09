"""Check the effective CMake compile commands for test instrumentation leakage."""

import json
from pathlib import Path
import re
import sys


def main():
    instrumented = {
        "rhbm_gem": False,
        "rhbm_gem_cli": False,
        "rhbm_gem_test": True,
        "rhbm_tests": True,
    }
    target_path = re.compile(r"CMakeFiles/([^/]+)\.dir/")
    instrumentation = re.compile(
        r"""(?:^|\s)(?:-D|/D)\s*["']?RHBM_GEM_TEST_INSTRUMENTATION(?=[\s="']|$)"""
    )
    entries = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    sources = {target: set() for target in instrumented}
    errors = []
    for entry in entries:
        command = entry["command"]
        # Older CMake versions only provide the object path in the command.
        match = target_path.search((entry.get("output") or command).replace("\\", "/"))
        if not match or match[1] not in instrumented:
            continue
        target = match[1]
        source = entry["file"].replace("\\", "/")
        sources[target].add(source)
        if bool(instrumentation.search(command)) != instrumented[target]:
            errors.append(f"{target}: incorrect instrumentation definition for {source}")
        if re.search(r"(?:^|/)tests/support/", source) and target != "rhbm_tests":
            errors.append(f"{target}: test support source compiled outside rhbm_tests: {source}")

    for target, files in sources.items():
        if not files:
            errors.append(f"{target}: no compile commands found")
        print(f"{target}: {len(files)} sources")
    if sources["rhbm_gem"] != sources["rhbm_gem_test"]:
        errors.append("rhbm_gem and rhbm_gem_test compile different source sets")
    if errors:
        sys.exit("\n".join(errors))


if __name__ == "__main__":
    main()
