"""Summarize exact numeric factor reuse opportunities for FixedNeighbor searches."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


COUNTERS = (
    "numeric_factor_requests",
    "numeric_factor_exact_reuse_opportunities",
    "numeric_factor_pattern_only_matches",
    "numeric_factor_value_mismatches",
    "numeric_factor_column_mismatches",
    "numeric_factor_policy_mismatches",
    "numeric_factor_pattern_mismatches",
    "initial_profile_exact_reuse_opportunities",
    "trial_profile_exact_reuse_opportunities",
    "reference_exact_reuse_opportunities",
    "accepted_endpoint_exact_reuse_opportunities",
)


def _fixed(case):
    return case.get("fixed_neighbor", case)


def analyze_case(case):
    fixed = _fixed(case)
    values = {name: fixed.get(name, 0) for name in COUNTERS}
    requests = values["numeric_factor_requests"]
    exact = values["numeric_factor_exact_reuse_opportunities"]
    values.update({
        "topology": case.get("topology"),
        "atoms": case.get("atoms"),
        "exact_reuse_rate": exact / requests if requests else None,
    })
    return values


def analyze(cases):
    return {
        "phase": "FixedNeighbor exact numeric factor reuse opportunity census",
        "cases": sorted(
            (analyze_case(case) for case in cases),
            key=lambda row: (row.get("topology") or "", row.get("atoms") or 0),
        ),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args(argv)
    cases = [json.loads(Path(path).read_text()) for path in args.cases]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(analyze(cases), indent=2) + "\n")


if __name__ == "__main__":
    main()
