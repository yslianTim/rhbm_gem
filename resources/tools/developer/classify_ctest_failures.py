#!/usr/bin/env python3
"""Classify failed CTest cases into contract/runtime buckets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


def parse_failed_tests(ctest_log: str) -> list[str]:
    failed: list[str] = []

    block_pattern = re.compile(r"The following tests FAILED:\n(?P<body>(?:\s+.+\n)+)", re.MULTILINE)
    block_match = block_pattern.search(ctest_log)
    if block_match:
        for line in block_match.group("body").splitlines():
            m = re.search(r"\s+\d+\s*-\s*(.+?)\s+\((?:Failed|Timeout|SEGFAULT)\)", line)
            if m:
                failed.append(m.group(1).strip())

    if failed:
        return sorted(set(failed))

    inline_pattern = re.compile(r"Test #\d+:\s+(.+?)\s+\.*\*{3}(?:Failed|Timeout|Exception)")
    return sorted({m.group(1).strip() for m in inline_pattern.finditer(ctest_log)})


def load_test_labels(build_dir: Path) -> dict[str, set[str]]:
    proc = subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "--show-only=json-v1"],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return {}
    payload = json.loads(proc.stdout)
    label_map: dict[str, set[str]] = {}
    for test_entry in payload.get("tests", []):
        name = test_entry.get("name")
        labels: set[str] = set()
        for prop in test_entry.get("properties", []):
            if prop.get("name") == "LABELS":
                labels.update(prop.get("value", []))
        if name:
            label_map[name] = labels
    return label_map


def main() -> int:
    parser = argparse.ArgumentParser(description="Classify ctest failures by intent label.")
    parser.add_argument("--build-dir", required=True, help="CTest build directory.")
    parser.add_argument("--ctest-log", required=True, help="Path to captured ctest stdout/stderr.")
    args = parser.parse_args()

    ctest_log = Path(args.ctest_log).read_text(encoding="utf-8", errors="ignore")
    failed_tests = parse_failed_tests(ctest_log)
    if not failed_tests:
        print("Could not parse failed test names from ctest output.")
        return 0

    label_map = load_test_labels(Path(args.build_dir))
    categories = {"contract": [], "runtime": []}
    for name in failed_tests:
        labels = label_map.get(name, set())
        if "intent:contract" in labels:
            categories["contract"].append(name)
        else:
            categories["runtime"].append(name)

    print("CTest failure classification:")
    print(f" - contract: {len(categories['contract'])}")
    for name in categories["contract"]:
        print(f"   * {name}")
    print(f" - runtime: {len(categories['runtime'])}")
    for name in categories["runtime"]:
        print(f"   * {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
