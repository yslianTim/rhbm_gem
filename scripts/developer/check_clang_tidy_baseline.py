#!/usr/bin/env python3
"""Check clang-tidy warning counts against a stored baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


WARNING_PATTERN = re.compile(r"^(?P<path>[^:\n]+):\d+:\d+:\s+warning:", re.MULTILINE)


def parse_report(report_text: str) -> tuple[int, dict[str, int]]:
    per_file: dict[str, int] = {}
    for match in WARNING_PATTERN.finditer(report_text):
        path = match.group("path")
        per_file[path] = per_file.get(path, 0) + 1
    total = sum(per_file.values())
    return total, per_file


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare clang-tidy warnings to a baseline.")
    parser.add_argument("--report", required=True, help="Path to captured clang-tidy output.")
    parser.add_argument("--baseline", required=True, help="Path to baseline json file.")
    parser.add_argument("--update-baseline", action="store_true", help="Rewrite baseline from report.")
    args = parser.parse_args()

    report_path = Path(args.report)
    baseline_path = Path(args.baseline)
    report_text = report_path.read_text(encoding="utf-8", errors="ignore")
    total, per_file = parse_report(report_text)

    if args.update_baseline:
        payload = {
            "total_warnings_max": total,
            "per_file_max": dict(sorted(per_file.items())),
        }
        baseline_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Updated clang-tidy baseline: total={total}, files={len(per_file)}")
        return 0

    if not baseline_path.exists():
        print(f"Missing baseline file: {baseline_path}")
        return 1

    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    total_max = int(baseline.get("total_warnings_max", 0))
    per_file_max = baseline.get("per_file_max", {})

    errors: list[str] = []
    if total > total_max:
        errors.append(f"total warnings exceeded baseline: current={total}, baseline={total_max}")

    for file_path, count in per_file.items():
        allowed = int(per_file_max.get(file_path, total_max))
        if count > allowed:
            errors.append(
                f"{file_path} warnings exceeded baseline: current={count}, baseline={allowed}"
            )

    if errors:
        print("clang-tidy baseline check failed:")
        for err in errors:
            print(f" - {err}")
        return 1

    print(
        "clang-tidy baseline check passed: "
        f"current_total={total}, baseline_total={total_max}, files={len(per_file)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
