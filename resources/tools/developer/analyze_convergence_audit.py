#!/usr/bin/env python3
"""Aggregate second-stage convergence safeguard Debug records."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import re
from typing import Iterable


MARKER = "Convergence safeguard audit:"
FIELD_PATTERN = re.compile(
    r"(?:^|, )(?P<name>[a-z][a-z0-9-]*)(?:\[[^]]+\])?=(?P<value>[^,]+)"
)
PREDICATE_NAMES = (
    "qualification",
    "accepted_p99",
    "residual_p99",
)
LEGACY_MAXIMUM_PREDICATE_NAMES = (
    "qualification",
    "accepted_p99",
    "accepted_max",
    "guarded_p99",
    "guarded_max",
)
TRACK_FIELDS = {
    "production": "certificate",
    "historical_all_selected": "historical-all-selected-predicates",
    "historical_active_proposal": "historical-active-proposal-predicates",
}
POPULATION_FIELDS = {
    "production": "accepted-active-population",
    "historical_all_selected": "historical-all-selected-population",
    "historical_active_proposal": "accepted-active-population",
}
SUMMARY_FIELDS = {
    "production": (
        "accepted-active-p99",
        "accepted-active-max",
        "operator-nominal-residual-p99",
        "operator-nominal-residual-max",
    ),
    "historical_all_selected": (
        "historical-all-selected-accepted-p99",
        "historical-all-selected-accepted-max",
        "historical-active-proposal-p99",
        "historical-active-proposal-max",
    ),
    "historical_active_proposal": (
        "accepted-active-p99",
        "accepted-active-max",
        "historical-active-proposal-p99",
        "historical-active-proposal-max",
    ),
}
EXPOSURE_NAMES = (
    "historical-all-selected",
    "historical-cluster-active-proposal-maximum",
    "historical-active-proposal",
    "production-maximum",
)
SCHEMA_8_CERTIFICATE_LABEL = (
    "certificate[solver/accepted-p99/operator-complete/operator-p99/"
    "invariants/orthogonal/production]="
)
SCHEMA_8_STOP_LABEL = (
    "stop-candidates[production/historical-all-selected/"
    "historical-cluster-active-proposal-maximum/historical-active-proposal/"
    "production-maximum]="
)
SCHEMA_8_EXPOSURE_LABEL = (
    "exposures[historical-all-selected/"
    "historical-cluster-active-proposal-maximum/historical-active-proposal/"
    "production-maximum]="
)


def _numbers(value: str) -> list[float]:
    return [float(item) for item in value.rstrip(".").split("/")]


def _integers(value: str) -> list[int]:
    return [int(float(item)) for item in value.rstrip(".").split("/")]


def summary_field_names(record: dict[str, str], track: str) -> tuple[str, ...]:
    return SUMMARY_FIELDS[track]


def _join(values: list[int] | list[str]) -> str:
    return "/".join(str(value) for value in values)


def _normalize_legacy(fields: dict[str, str]) -> dict[str, str]:
    schema = fields["schema"]
    operator = _integers(fields["operator-predicates"])
    stops = _integers(fields["stop-candidates"])
    qualification = _integers(fields["qualification"])
    legacy_maximum = _integers(fields["legacy-maximum-predicates"])
    exposures = _integers(fields["exposures"])
    if len(operator) != 5 or len(stops) != 7 or len(exposures) != 5:
        raise ValueError("Legacy convergence audit vectors have invalid sizes")
    guarded_p99 = fields.get(
        "guarded-proposal-p99", fields["fixed-point-residual-p99"])
    guarded_max = fields.get(
        "guarded-proposal-max", fields["fixed-point-residual-max"])
    blockers = _integers(fields["blockers"])
    return {
        **fields,
        "source-schema": schema,
        "accepted-active-population": fields["production-population"],
        "historical-all-selected-population": fields["legacy-population"],
        "operator-nominal-population": fields["operator-population"],
        "certificate": _join([
            qualification[1], operator[1], operator[3], operator[2],
            int(len(qualification) < 17 or qualification[16] == 0),
            stops[0], stops[1],
        ]),
        "historical-all-selected-predicates": fields[
            "legacy-population-predicates"],
        "historical-cluster-active-proposal-maximum-predicates": fields[
            "legacy-maximum-predicates"],
        "historical-active-proposal-predicates": fields[
            "solver-qualified-predicates"],
        "production-maximum-predicates": _join([
            stops[1], legacy_maximum[2], operator[4]]),
        "accepted-active-p99": fields["production-accepted-p99"],
        "accepted-active-max": fields["production-accepted-max"],
        "historical-all-selected-accepted-p99": fields["legacy-accepted-p99"],
        "historical-all-selected-accepted-max": fields["legacy-accepted-max"],
        "historical-active-proposal-p99": guarded_p99,
        "historical-active-proposal-max": guarded_max,
        "operator-nominal-residual-p99": fields["fixed-point-residual-p99"],
        "operator-nominal-residual-max": fields["fixed-point-residual-max"],
        "operator-nominal-unavailable": fields["operator-unavailable"],
        "operator-nominal-unavailable-reasons": fields[
            "operator-unavailable-reasons"],
        "operator-nominal-tail": fields["operator-tail"],
        "certificate-blockers": _join([
            blockers[3], blockers[2], blockers[0], blockers[1]]),
        "stop-candidates": _join([
            stops[1], stops[2], stops[3], stops[4], stops[6]]),
        "exposures": _join([
            exposures[0], exposures[1], exposures[2], exposures[4]]),
    }


def parse_record(line: str) -> dict[str, str] | None:
    marker_position = line.find(MARKER)
    if marker_position < 0:
        return None
    payload = line[marker_position + len(MARKER):].strip()
    fields = {
        match.group("name"): match.group("value").strip()
        for match in FIELD_PATTERN.finditer(payload)
    }
    schema = fields.get("schema")
    if schema not in {"6", "7", "8"}:
        return None
    if schema in {"6", "7"}:
        legacy_required = {
            "production-population", "legacy-population", "operator-population",
            "production-predicates", "legacy-population-predicates",
            "legacy-maximum-predicates", "solver-qualified-predicates",
            "operator-predicates", "production-accepted-p99",
            "production-accepted-max", "legacy-accepted-p99",
            "legacy-accepted-max", "fixed-point-residual-p99",
            "fixed-point-residual-max", "operator-unavailable",
            "operator-unavailable-reasons", "operator-tail", "blockers",
            "qualification", "stop-candidates", "exposures", "path", "ratios",
            "offset-groups", "joint-status", "residual-state", "limiters",
        }
        if schema == "7":
            legacy_required.add("unified-search")
        if not legacy_required.issubset(fields):
            missing = ", ".join(sorted(legacy_required - fields.keys()))
            raise ValueError(f"Convergence audit record is missing: {missing}")
        if {"certificate-definition", "comparator-set", "certificate"} & fields.keys():
            raise ValueError("Convergence audit record mixes schema definitions")
        return _normalize_legacy(fields)

    required = set(TRACK_FIELDS.values()) | set(POPULATION_FIELDS.values()) | {
        "certificate-definition", "comparator-set",
        "qualification",
        "stop-candidates",
        "exposures",
        "historical-cluster-active-proposal-maximum-predicates",
        "production-maximum-predicates",
        "path",
        "ratios",
        "offset-groups",
        "joint-status",
        "operator-nominal-unavailable",
        "operator-nominal-unavailable-reasons",
        "operator-nominal-tail",
        "residual-state",
        "limiters",
        "certificate-blockers",
        "unified-search",
    }
    required.update(field for names in SUMMARY_FIELDS.values() for field in names)
    if not required.issubset(fields):
        missing = ", ".join(sorted(required - fields.keys()))
        raise ValueError(f"Convergence audit record is missing: {missing}")
    legacy_only = {
        "production-population", "legacy-population", "operator-population",
        "production-predicates", "legacy-population-predicates",
        "legacy-maximum-predicates", "solver-qualified-predicates",
        "operator-predicates", "production-accepted-p99",
        "fixed-point-residual-p99",
    }
    if legacy_only & fields.keys():
        raise ValueError("Convergence audit record mixes schema definitions")
    if fields["certificate-definition"] != "1" or fields["comparator-set"] != "1":
        raise ValueError("Unsupported convergence certificate definition")
    if not all(label in payload for label in (
            SCHEMA_8_CERTIFICATE_LABEL,
            SCHEMA_8_STOP_LABEL,
            SCHEMA_8_EXPOSURE_LABEL)):
        raise ValueError("Schema-8 comparator names or order do not match")
    return fields


def normalize_record(record: dict[str, str]) -> dict[str, str]:
    """Normalize a previously parsed schema-6/7/8 record."""
    if record.get("schema") in {"6", "7"}:
        return _normalize_legacy(dict(record))
    if record.get("schema") == "8":
        return dict(record)
    raise ValueError("Unsupported convergence audit schema")


def _truth_tables(vectors: list[list[int]]) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    for left in range(len(PREDICATE_NAMES)):
        for right in range(left + 1, len(PREDICATE_NAMES)):
            counts = Counter(
                f"{vector[left]}{vector[right]}" for vector in vectors
            )
            result[f"{PREDICATE_NAMES[left]}|{PREDICATE_NAMES[right]}"] = {
                key: counts.get(key, 0) for key in ("00", "01", "10", "11")
            }
    return result


def _bin_ratio(value: float) -> str:
    if value == 0.0:
        return "0"
    if value <= 0.25:
        return "(0,0.25]"
    if value <= 0.5:
        return "(0.25,0.5]"
    if value <= 0.9:
        return "(0.5,0.9]"
    return "(0.9,1]"


def _bin_group_size(value: float) -> str:
    if value <= 1:
        return "<=1"
    if value <= 10:
        return "2-10"
    if value <= 91:
        return "11-91"
    return ">91"


def _record_stratum(
    table: dict[str, Counter[str]],
    name: str,
    value: str,
    actual_exposure: bool,
) -> None:
    table[name][f"{value}:records"] += 1
    if actual_exposure:
        table[name][f"{value}:actual_exposures"] += 1


def analyze_records(records: Iterable[dict[str, str]]) -> dict[str, object]:
    record_list = list(records)
    vectors_by_track: dict[str, list[list[int]]] = defaultdict(list)
    blocker_count: dict[str, Counter[str]] = defaultdict(Counter)
    unique_blocker_count: dict[str, Counter[str]] = defaultdict(Counter)
    implication_counterexamples: Counter[str] = Counter()
    implication_examples: dict[str, list[str]] = defaultdict(list)
    p99_implies_max: dict[str, Counter[str]] = defaultdict(Counter)
    p99_implies_max_examples: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: defaultdict(list)
    )
    strata: dict[str, Counter[str]] = defaultdict(Counter)
    exposure_count: Counter[str] = Counter()
    actual_exposure_count = 0
    actual_exposure_examples: list[str] = []
    diagnostic_mismatch_count = 0

    for record in record_list:
        record_ref = f"try={record.get('try', '?')},acc={record.get('acc', '?')}"
        vector_by_track: dict[str, list[int]] = {}
        for track, field in TRACK_FIELDS.items():
            vector = _integers(record[field])
            if track == "production":
                if len(vector) != 7:
                    raise ValueError(
                        "certificate must contain seven predicates")
                expected_production = int(all(vector[:6]))
                if (record.get("schema") == "8" and
                        vector[6] != expected_production):
                    raise ValueError(
                        "Serialized production decision contradicts certificate")
                vector = [vector[0], vector[1], vector[3]]
            elif len(vector) != len(PREDICATE_NAMES):
                raise ValueError(f"{field} must contain three predicates")
            vector_by_track[track] = vector
            vectors_by_track[track].append(vector)
            failed = [index for index, passed in enumerate(vector) if not passed]
            for index in failed:
                blocker_count[track][PREDICATE_NAMES[index]] += 1
            if len(failed) == 1:
                unique_blocker_count[track][PREDICATE_NAMES[failed[0]]] += 1

        legacy_maximum_vector = _integers(
            record["historical-cluster-active-proposal-maximum-predicates"])
        if len(legacy_maximum_vector) != len(LEGACY_MAXIMUM_PREDICATE_NAMES):
            raise ValueError(
                "historical maximum predicates must contain five predicates")
        legacy_maximum_failed = [
            index for index, passed in enumerate(legacy_maximum_vector)
            if not passed
        ]
        for index in legacy_maximum_failed:
            blocker_count["historical_cluster_active_proposal_maximum"][
                LEGACY_MAXIMUM_PREDICATE_NAMES[index]
            ] += 1
        if len(legacy_maximum_failed) == 1:
            unique_blocker_count["historical_cluster_active_proposal_maximum"][
                LEGACY_MAXIMUM_PREDICATE_NAMES[legacy_maximum_failed[0]]
            ] += 1

        qualification = _integers(record["qualification"])
        stops = _integers(record["stop-candidates"])
        exposures = _integers(record["exposures"])
        if len(stops) != 5:
            raise ValueError("stop-candidates must contain five values")
        if len(exposures) != len(EXPOSURE_NAMES):
            raise ValueError("exposures must contain four values")
        for name, exposed in zip(EXPOSURE_NAMES, exposures):
            exposure_count[name] += int(bool(exposed))

        if vector_by_track["production"][0] and not qualification[1]:
            implication_name = "production_qualification=>solver_qualified"
            implication_counterexamples[implication_name] += 1
            implication_examples[implication_name].append(record_ref)
        for comparison_track in (
                "historical_all_selected", "historical_active_proposal"):
            for index, name in enumerate(PREDICATE_NAMES[1:], start=1):
                if (vector_by_track["production"][index] and
                        not vector_by_track[comparison_track][index]):
                    implication_name = (
                        f"production_{name}=>{comparison_track}_{name}")
                    implication_counterexamples[implication_name] += 1
                    implication_examples[implication_name].append(record_ref)

        actual_exposure = any(bool(value) for value in exposures)
        actual_exposure_count += int(actual_exposure)
        if actual_exposure:
            actual_exposure_examples.append(record_ref)
        production_predicate_pass = all(vector_by_track["production"])
        comparison_mismatch = (
            not all(vector_by_track["historical_all_selected"]) or
            not all(legacy_maximum_vector) or
            not all(vector_by_track["historical_active_proposal"])
        )
        diagnostic_mismatch_count += int(
            production_predicate_pass and comparison_mismatch and not stops[0])

        for track, population_field in POPULATION_FIELDS.items():
            populations = _integers(record[population_field])
            accepted_p99, accepted_max, residual_p99, residual_max = (
                _numbers(record[field]) for field in summary_field_names(record, track)
            )
            for coordinate, population_size in enumerate(populations):
                size_bin = "N<=91" if population_size <= 91 else "N>91"
                for state, p99_values, max_values in (
                    ("accepted", accepted_p99, accepted_max),
                    ("residual", residual_p99, residual_max),
                ):
                    p99_pass = math.isfinite(p99_values[coordinate]) and p99_values[coordinate] < 1.0e-4
                    max_pass = math.isfinite(max_values[coordinate]) and max_values[coordinate] < 1.0e-3
                    key = f"{state}:{size_bin}"
                    p99_implies_max[track][f"{key}:samples"] += 1
                    if p99_pass and not max_pass:
                        p99_implies_max[track][f"{key}:counterexamples"] += 1
                        p99_implies_max_examples[track][key].append(
                            f"{record_ref},coordinate={coordinate},N={population_size}"
                        )

        ratios = _numbers(record["ratios"])
        group_stats = _numbers(record["offset-groups"])
        joint_status = _integers(record["joint-status"])
        path = _integers(record["path"])
        _record_stratum(strata, "shape_active_ratio", _bin_ratio(ratios[0]), actual_exposure)
        _record_stratum(strata, "offset_active_ratio", _bin_ratio(ratios[1]), actual_exposure)
        _record_stratum(strata, "quarantine_ratio", _bin_ratio(ratios[2]), actual_exposure)
        _record_stratum(strata, "maximum_group_size", _bin_group_size(group_stats[-1]), actual_exposure)
        joint_status_names = (
            "converged",
            "system_build_failed",
            "empty_system",
            "initial_solve_failed",
            "irls_solve_failed",
            "objective_deteriorated",
            "maximum_iterations",
        )
        solver_label = "+".join(
            name for name, count in zip(joint_status_names, joint_status) if count
        ) or "none"
        _record_stratum(strata, "solver_status", solver_label, actual_exposure)
        path_names = (
            ("limited", "polish", "boundary", "rescue")
            if record.get("schema") in {"7", "8"}
            else ("trust", "backtrack", "polish", "boundary", "rescue")
        )
        active_paths = [name for name, count in zip(path_names, path) if count]
        _record_stratum(
            strata,
            "proposal_path",
            "+".join(active_paths) if active_paths else "normal",
            actual_exposure,
        )
        _record_stratum(
            strata,
            "fixed_point_interpretation",
            record["residual-state"],
            actual_exposure,
        )
        unavailable = _integers(record["operator-nominal-unavailable"])
        _record_stratum(
            strata,
            "operator_availability",
            "complete" if not any(unavailable) else "restricted",
            actual_exposure,
        )

    normalized_blockers = {
        track: {
            predicate: blocker_count[track].get(predicate, 0)
            for predicate in PREDICATE_NAMES
        }
        for track in TRACK_FIELDS
    }
    normalized_blockers["historical_cluster_active_proposal_maximum"] = {
        predicate: blocker_count[
            "historical_cluster_active_proposal_maximum"].get(predicate, 0)
        for predicate in LEGACY_MAXIMUM_PREDICATE_NAMES
    }
    normalized_unique = {
        track: {
            predicate: unique_blocker_count[track].get(predicate, 0)
            for predicate in PREDICATE_NAMES
        }
        for track in TRACK_FIELDS
    }
    normalized_unique["historical_cluster_active_proposal_maximum"] = {
        predicate: unique_blocker_count[
            "historical_cluster_active_proposal_maximum"].get(predicate, 0)
        for predicate in LEGACY_MAXIMUM_PREDICATE_NAMES
    }
    return {
        "record_count": len(record_list),
        "actual_stop_exposure_count": actual_exposure_count,
        "actual_stop_exposure_examples": actual_exposure_examples,
        "exposure_count": {
            name: exposure_count.get(name, 0) for name in EXPOSURE_NAMES
        },
        "diagnostic_mismatch_count": diagnostic_mismatch_count,
        "blocker_count": normalized_blockers,
        "unique_blocker_count": normalized_unique,
        "truth_tables": {
            track: _truth_tables(vectors)
            for track, vectors in vectors_by_track.items()
        },
        "implication_counterexamples": dict(sorted(implication_counterexamples.items())),
        "implication_counterexample_examples": {
            name: examples for name, examples in sorted(implication_examples.items())
        },
        "p99_implies_max": {
            track: dict(sorted(counts.items()))
            for track, counts in p99_implies_max.items()
        },
        "p99_implies_max_examples": {
            track: dict(sorted(examples.items()))
            for track, examples in p99_implies_max_examples.items()
        },
        "strata": {
            name: dict(sorted(counts.items()))
            for name, counts in strata.items()
        },
    }


def format_markdown(report: dict[str, object]) -> str:
    lines = [
        "# Convergence safeguard audit aggregation",
        "",
        f"- Records: {report['record_count']}",
        f"- Actual stop exposures: {report['actual_stop_exposure_count']}",
        f"- Exposure counts: {json.dumps(report['exposure_count'], sort_keys=True)}",
        f"- Diagnostic-only mismatches: {report['diagnostic_mismatch_count']}",
        "",
        "## Blockers",
        "",
        "| Track | Predicate | Blockers | Unique blockers |",
        "| --- | --- | ---: | ---: |",
    ]
    blockers = report["blocker_count"]
    unique = report["unique_blocker_count"]
    for track in TRACK_FIELDS:
        for predicate in PREDICATE_NAMES:
            lines.append(
                f"| {track} | {predicate} | "
                f"{blockers.get(track, {}).get(predicate, 0)} | "
                f"{unique.get(track, {}).get(predicate, 0)} |"
            )
    for predicate in LEGACY_MAXIMUM_PREDICATE_NAMES:
        lines.append(
            f"| legacy_maximum | {predicate} | "
            f"{blockers.get('legacy_maximum', {}).get(predicate, 0)} | "
            f"{unique.get('legacy_maximum', {}).get(predicate, 0)} |"
        )
    lines.extend(("", "## Implication counterexamples", ""))
    counterexamples = report["implication_counterexamples"]
    if counterexamples:
        for name, count in counterexamples.items():
            lines.append(f"- `{name}`: {count}")
    else:
        lines.append("- None observed.")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "log", type=Path,
        help="Debug log containing schema=6, schema=7, or schema=8 audit records")
    parser.add_argument("--format", choices=("json", "markdown"), default="markdown")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    records = []
    for line in args.log.read_text(encoding="utf-8", errors="replace").splitlines():
        record = parse_record(line)
        if record is not None:
            records.append(record)
    report = analyze_records(records)
    output = (
        json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else format_markdown(report)
    )
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
