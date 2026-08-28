#!/usr/bin/env python3
"""Analyze developer-only frozen-IRLS trust-model experiment logs."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
from typing import Any, Sequence


FUNNEL_MARKER = "Trust-model funnel:"
TRIAL_MARKER = "Trust-model shadow:"
FIELD_PATTERN = re.compile(
    r"(?:^|, )(?P<name>[a-z][a-z0-9-]*)=(?P<value>[^,]+)")


def _fields(line: str, marker: str) -> dict[str, str] | None:
    position = line.find(marker)
    if position < 0:
        return None
    payload = line[position + len(marker):].strip()
    return {
        match.group("name"): match.group("value").strip().rstrip(".")
        for match in FIELD_PATTERN.finditer(payload)
    }


def parse_log(text: str) -> dict[str, list[dict[str, str]]]:
    funnels: list[dict[str, str]] = []
    trials: list[dict[str, str]] = []
    for line in text.splitlines():
        if fields := _fields(line, FUNNEL_MARKER):
            if fields.get("schema") == "1":
                funnels.append(fields)
        if fields := _fields(line, TRIAL_MARKER):
            if fields.get("schema") == "2":
                trials.append(fields)
    return {"funnels": funnels, "trials": trials}


def analyze(parsed: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    trials = parsed["trials"]
    status_counts = Counter(record["status"] for record in trials)
    rho_values = [
        float(record["rho"]) for record in trials if record.get("rho") != "-"
    ]
    rho_bins = Counter(
        "low" if rho < 0.25 else "mid" if rho <= 0.75 else "high"
        for rho in rho_values)
    action_confusion = Counter(
        f'{record["current-action"]}->{record["shadow-action"]}'
        for record in trials
        if record.get("readiness-eligible") == "1" and
        record.get("shadow-action") != "suppressed")
    return {
        "schema_version": 1,
        "diagnostic_only": True,
        "funnel_record_count": len(parsed["funnels"]),
        "trial_record_count": len(trials),
        "status_counts": dict(sorted(status_counts.items())),
        "rho_bins": dict(sorted(rho_bins.items())),
        "action_confusion": dict(sorted(action_confusion.items())),
        "elapsed_milliseconds": sum(
            float(record.get("elapsed-ms", "0")) for record in trials),
    }


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path)
    parser.add_argument("--json", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    parsed = {"funnels": [], "trials": []}
    for path in args.logs:
        current = parse_log(path.read_text(encoding="utf-8"))
        parsed["funnels"].extend(current["funnels"])
        parsed["trials"].extend(current["trials"])
    payload = json.dumps(analyze(parsed), indent=2, sort_keys=True) + "\n"
    if args.json is None:
        print(payload, end="")
    else:
        args.json.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
