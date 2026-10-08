"""Qualify the selected fixed-neighbor local route on 256/512 cases."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from joint_fixed_neighbor_hybrid import _correct, summarize


EXPECTED = {(topology, atoms) for topology in ("chain", "cube") for atoms in (256, 512)}
SELECTED_POLICY = "+2 Legacy polish"


def analyze(rows):
    groups = {}
    for row in rows:
        groups.setdefault((row["topology"], row["atoms"]), []).append(row)
    cases = []
    for key, entries in sorted(groups.items()):
        cases.append({"topology": key[0], "atoms": key[1], "policies": entries,
                      "correct_policies": [row["policy"] for row in entries if _correct(row)]})
    seen = {(case["topology"], case["atoms"]) for case in cases}
    missing = sorted(EXPECTED - seen)
    selected_rows = [row for row in rows if row["policy"] == SELECTED_POLICY]
    route_passed = not missing and len(selected_rows) == len(EXPECTED) and all(
        _correct(row) for row in selected_rows)
    return {
        "phase": "P6/7 fixed-neighbor selected local route qualification",
        "selected_route": {"local_search": "OperatorPcg", "local_schwarz_core_atoms": 32,
                           "local_schwarz_overlap_hops": 0, "local_schwarz_max_block_atoms": 64,
                           "local_work": SELECTED_POLICY, "polish_solver": "LegacyCompact",
                           "maximum_polish_sweeps": 2},
        "expected_cases": [{"topology": topology, "atoms": atoms}
                           for topology, atoms in sorted(EXPECTED)],
        "missing_cases": [{"topology": topology, "atoms": atoms} for topology, atoms in missing],
        "cases": cases,
        "qualification": "passed" if route_passed else "failed",
        "promotion_blocked": not route_passed,
        "selection_basis": "all 256/512 endpoint, confirmation and RuntimeConvergence gates are mandatory",
        "thresholds": {"assessment_local": 1e-10, "global_ac_kkt": 1e-10,
                       "width_gradient_inf_norm": 1e-12, "eta_confirmation": 1e-10},
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case", action="append", required=True,
                        help="topology:atoms:raw-json; raw JSON must be selected route")
    args = parser.parse_args(argv)
    rows = []
    for specification in args.case:
        topology, atoms, path = specification.split(":", 2)
        case = json.loads(Path(path).read_text())
        if case.get("topology") != topology or case.get("atoms") != int(atoms):
            parser.error(f"invalid case specification: {specification}")
        rows.append(summarize(case, SELECTED_POLICY))
    report = analyze(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "qualification": report["qualification"],
                      "promotion_blocked": report["promotion_blocked"]}))
    return 0 if report["qualification"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
