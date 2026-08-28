#!/usr/bin/env python3
"""Parse production-only second-stage convergence audit records."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
from typing import Any, Iterable, Sequence


TRAJECTORY_MARKER = "Convergence safeguard audit:"
TERMINAL_MARKER = "Second-stage audit terminal:"
TERMINAL_ATOM_MARKER = "Second-stage audit terminal atom:"
FIELD_PATTERN = re.compile(
    r"(?:^|, )(?P<name>[a-z][a-z0-9-]*)(?:\[[^]]+\])?=(?P<value>[^,]+)"
)
CURRENT_SCHEMA = "9"
COMPATIBILITY_SCHEMA = "8"
CERTIFICATE_DEFINITION = "1"
PRODUCTION_FIELDS = (
    "try",
    "acc",
    "accepted-active-population",
    "operator-nominal-population",
    "certificate",
    "accepted-active-p99",
    "accepted-active-max",
    "operator-nominal-residual-p99",
    "operator-nominal-residual-max",
    "operator-nominal-unavailable",
    "operator-nominal-unavailable-reasons",
    "operator-nominal-tail",
    "residual-state",
    "qualification",
    "unified-search",
    "path",
    "limiters",
    "fixed",
    "blockers",
    "accepted-equals-operator",
)
CURRENT_REQUIRED_FIELDS = {
    "schema",
    "certificate-definition",
    "atoms",
    "quarantine",
    "operator-shadow-refit",
    "joint-status",
    "local-status",
    "offset-groups",
    "ratios",
    "certificate-blockers",
    *PRODUCTION_FIELDS,
}
FORBIDDEN_CURRENT_FIELDS = {
    "comparator-set",
    "stop-candidates",
    "exposures",
    "production-maximum-predicates",
}


def _fields(line: str, marker: str) -> dict[str, str] | None:
    position = line.find(marker)
    if position < 0:
        return None
    payload = line[position + len(marker):].strip()
    return {
        match.group("name"): match.group("value").strip().rstrip(".")
        for match in FIELD_PATTERN.finditer(payload)
    }


def _integers(value: str) -> list[int]:
    return [int(float(item)) for item in value.split("/")]


def _validate_certificate(record: dict[str, str]) -> None:
    certificate = _integers(record["certificate"])
    if len(certificate) != 7:
        raise ValueError("Certificate vector must contain seven values")
    if any(value not in (0, 1) for value in certificate):
        raise ValueError("Certificate vector must contain boolean values")
    if certificate[6] != int(all(certificate[:6])):
        raise ValueError("Serialized production decision is not the conjunction")


def parse_record(line: str) -> dict[str, str] | None:
    fields = _fields(line, TRAJECTORY_MARKER)
    if fields is None:
        return None
    schema = fields.get("schema")
    if schema not in (CURRENT_SCHEMA, COMPATIBILITY_SCHEMA):
        return None
    missing = CURRENT_REQUIRED_FIELDS - fields.keys()
    if missing:
        raise ValueError(
            "Convergence audit record is missing: " + ", ".join(sorted(missing)))
    if fields["certificate-definition"] != CERTIFICATE_DEFINITION:
        raise ValueError("Unsupported convergence certificate definition")
    if schema == CURRENT_SCHEMA:
        forbidden = FORBIDDEN_CURRENT_FIELDS & fields.keys()
        historical = {name for name in fields if name.startswith("historical-")}
        if forbidden or historical:
            names = ", ".join(sorted(forbidden | historical))
            raise ValueError("Current convergence record contains retired fields: " + names)
    elif fields.get("comparator-set") != "1":
        raise ValueError("Schema-8 compatibility record has an unknown comparator set")
    _validate_certificate(fields)
    return fields


def normalize_record(record: dict[str, str]) -> dict[str, str]:
    if record.get("schema") not in (CURRENT_SCHEMA, COMPATIBILITY_SCHEMA):
        raise ValueError("Unsupported convergence audit schema")
    return {
        "source-schema": record["schema"],
        **{name: record.get(name, "-") for name in PRODUCTION_FIELDS},
    }


def parse_log(text: str) -> dict[str, Any]:
    trajectory_records: list[dict[str, str]] = []
    terminal: dict[str, str] | None = None
    terminal_atoms: list[dict[str, str]] = []
    for line in text.splitlines():
        if record := parse_record(line):
            trajectory_records.append(record)
        if fields := _fields(line, TERMINAL_MARKER):
            if fields.get("schema") not in ("1", "2"):
                continue
            required = {"reason", "try", "acc", "fixed-domain", "objective"}
            if not required.issubset(fields):
                raise ValueError("Terminal schema-2 record is incomplete")
            terminal = fields
        if fields := _fields(line, TERMINAL_ATOM_MARKER):
            if fields.get("schema") != "1":
                continue
            required = {"serial", "group", "amplitude", "width", "offset"}
            if not required.issubset(fields):
                raise ValueError("Terminal atom schema-1 record is incomplete")
            terminal_atoms.append(fields)
    return {
        "trajectory_records": trajectory_records,
        "terminal": terminal,
        "terminal_atoms": terminal_atoms,
    }


def analyze_records(records: Iterable[dict[str, str]]) -> dict[str, Any]:
    rows = list(records)
    normalized = [normalize_record(record) for record in rows]
    return {
        "schema_version": 1,
        "record_count": len(rows),
        "production_convergence_count": sum(
            _integers(record["certificate"])[6] for record in rows),
        "residual_state_counts": dict(sorted(Counter(
            record["residual-state"] for record in rows).items())),
        "records": normalized,
    }


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path)
    parser.add_argument("--json", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    records = [
        record
        for path in args.logs
        for line in path.read_text(encoding="utf-8").splitlines()
        if (record := parse_record(line)) is not None
    ]
    report = analyze_records(records)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json is None:
        print(payload, end="")
    else:
        args.json.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
