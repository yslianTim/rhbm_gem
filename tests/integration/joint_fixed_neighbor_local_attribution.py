"""Summarize aggregate local FixedNeighbor SearchProfile attribution."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


PHASES = (
    ("profile_basis_seconds", "Basis construction"),
    ("linear_matrix_preparation_seconds", "Linear solve setup"),
    ("linear_symbolic_seconds", "Symbolic QR"),
    ("linear_numeric_seconds", "Numeric QR"),
    ("linear_rhs_solve_seconds", "RHS / solve"),
    ("linear_certificate_seconds", "Certificate"),
    ("derivative_prepare_seconds", "Derivative preparation"),
    ("derivative_reduce_seconds", "Derivative reduction"),
    ("replay_trust_seconds", "Replay / trust"),
    ("lm_overhead_seconds", "LM / control overhead"),
)


def _finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def _fixed(case):
    return case.get("fixed_neighbor", case)


def _profile_work(fixed):
    return fixed.get("local_profile_work") or (fixed.get("fixed_neighbor_work") or {}).get("local_profile_work", {})


def analyze_case(case):
    fixed = _fixed(case)
    work = _profile_work(fixed)
    total = work.get("total") or {}
    local_seconds = work.get("total_seconds")
    search_seconds = fixed.get("search_seconds")
    rows = []
    for field, label in PHASES:
        seconds = work.get("lm_overhead_seconds") if field == "lm_overhead_seconds" else total.get(field)
        seconds = seconds if _finite(seconds) else 0.0
        rows.append({
            "phase": label,
            "seconds": seconds,
            "percent_of_local_search": seconds / local_seconds * 100.0
            if _finite(local_seconds) and local_seconds > 0 else None,
            "percent_of_total_search": seconds / search_seconds * 100.0
            if _finite(search_seconds) and search_seconds > 0 else None,
        })
    return {
        "topology": case.get("topology"),
        "atoms": case.get("atoms"),
        "profile_evaluations": total.get("evaluations", 0),
        "initial_profile_evaluations": (work.get("initial_profile") or {}).get("evaluations", 0),
        "trial_profile_evaluations": (work.get("trial_profile") or {}).get("evaluations", 0),
        "accepted_endpoint_evaluations": (work.get("accepted_endpoint") or {}).get("evaluations", 0),
        "reference_evaluations": (work.get("reference_evaluation") or {}).get("evaluations", 0),
        "local_search_seconds": local_seconds,
        "search_seconds": search_seconds,
        "phases": rows,
    }


def analyze(cases):
    return {"phase": "FixedNeighbor local SearchProfile attribution",
            "cases": sorted((analyze_case(case) for case in cases),
                            key=lambda row: (row.get("topology") or "", row.get("atoms") or 0))}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args(argv)
    cases = []
    for path in args.cases:
        value = json.loads(Path(path).read_text())
        cases.append(value.get("result", value))
    args.output.write_text(json.dumps(analyze(cases), indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
