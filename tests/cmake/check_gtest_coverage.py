#!/usr/bin/env python3

"""Check actual grouped CTest filters against the compiled GoogleTest registry."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def run(command: list[str], **kwargs) -> str:
    try:
        return subprocess.run(
            command, check=True, capture_output=True, text=True, timeout=30, **kwargs
        ).stdout
    except (OSError, subprocess.SubprocessError) as error:
        detail = getattr(error, "stderr", "")
        raise ValueError(f"Cannot list registration with {command!r}: {error}\n{detail}") from error


def parse_tests(output: str) -> set[str]:
    tests = set()
    suite = None
    suite_has_tests = False
    for raw in output.splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line:
            continue
        if suite is None and line.startswith("Running main() from "):
            continue
        if re.fullmatch(r"[A-Za-z0-9_/]+\.", line):
            if suite is not None and not suite_has_tests:
                raise ValueError(f"Incomplete GoogleTest suite: {suite}")
            suite = line
            suite_has_tests = False
        elif suite is not None and re.fullmatch(r"  [A-Za-z0-9_/]+", line):
            name = suite + line.strip()
            if name in tests:
                raise ValueError(f"Duplicate name in GoogleTest listing: {name}")
            tests.add(name)
            suite_has_tests = True
        else:
            raise ValueError(f"Cannot parse GoogleTest listing line: {raw!r}")
    if not tests or not suite_has_tests:
        raise ValueError("Empty or incomplete GoogleTest listing")
    return tests


def is_disabled(name: str) -> bool:
    return any(part.startswith("DISABLED_") for part in re.split(r"[./]", name))


def audit(binary: Path, registrations: dict) -> tuple[set[str], set[str], dict]:
    if not binary.is_file():
        raise ValueError(f"GoogleTest executable does not exist: {binary}")
    compiled = parse_tests(run([
        str(binary), "--gtest_list_tests", "--gtest_color=no", "--gtest_filter=*"
    ]))
    disabled = {name for name in compiled if is_disabled(name)}
    required = compiled - disabled
    if not required:
        raise ValueError("No enabled GoogleTests registered")

    entries = registrations.get("tests")
    if not isinstance(entries, list) or not entries:
        raise ValueError("Missing CTest registration information")
    coverage = defaultdict(list)
    groups = set()
    for entry in entries:
        command = entry.get("command")
        if not isinstance(command, list) or not command or not all(isinstance(x, str) for x in command):
            raise ValueError(f"Missing CTest command: {entry!r}")
        if Path(command[0]).resolve() != binary.resolve():
            continue
        name = entry["name"]
        if name in groups:
            raise ValueError(f"Duplicate CTest group: {name}")
        groups.add(name)
        filters = [arg for arg in command if arg.startswith("--gtest_filter=")]
        if len(filters) != 1 or not filters[0].removeprefix("--gtest_filter="):
            raise ValueError(f"Expected one explicit nonempty GoogleTest filter: {name}")
        properties = {prop["name"]: prop["value"] for prop in entry.get("properties", [])}
        if properties.get("DISABLED", False):
            print(f"{name}: excluded disabled CTest group")
            continue
        environment = dict(os.environ)
        for assignment in properties.get("ENVIRONMENT", []):
            key, value = assignment.split("=", 1)
            environment[key] = value
        selected = parse_tests(run(
            [*command, "--gtest_list_tests", "--gtest_color=no"],
            cwd=properties.get("WORKING_DIRECTORY"), env=environment,
        ))
        if not selected <= compiled:
            raise ValueError(f"Group {name} lists tests outside the full registry: {sorted(selected - compiled)}")
        for test in sorted(selected - disabled):
            coverage[test].append(name)
        print(f"{name}: {len(selected - disabled)} enabled GoogleTests")
    if not groups:
        raise ValueError("No CTest groups found for the GoogleTest executable")
    return compiled, disabled, dict(coverage)


def report(compiled: set[str], disabled: set[str], coverage: dict) -> int:
    missing = compiled - disabled - coverage.keys()
    duplicates = {name: owners for name, owners in coverage.items() if len(owners) > 1}
    print(f"Registered: {len(compiled)}; enabled: {len(compiled - disabled)}; "
          f"covered unique: {len(coverage)}; missing: {len(missing)}; "
          f"duplicate coverage: {len(duplicates)}; disabled excluded: {len(disabled)}")
    for name in sorted(disabled):
        print(f"DISABLED (excluded; not run by default): {name}")
    for name in sorted(missing):
        print(f"MISSING: {name}")
    for name, owners in sorted(duplicates.items()):
        intentional = name.startswith("DataObjectPersistenceTest.") and set(owners) == {
            "joint_persistence_contract_test", "rhbm_tests_data_schema"
        }
        print(f"DUPLICATE ({'intentional' if intentional else 'unexpected'}): {name}: {', '.join(owners)}")
    return int(bool(missing))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--ctest", required=True)
    parser.add_argument("--config", default="")
    args = parser.parse_args()
    command = [args.ctest, "--test-dir", str(args.build_dir), "--show-only=json-v1"]
    if args.config:
        command.extend(["-C", args.config])
    try:
        registrations = json.loads(run(command))
        return report(*audit(args.binary, registrations))
    except (ValueError, KeyError, TypeError, AttributeError) as error:
        print(f"GoogleTest registration coverage FAILED: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
