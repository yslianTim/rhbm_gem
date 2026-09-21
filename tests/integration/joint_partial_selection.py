"""Measure the complete selected-domain estimator in independent serial processes."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import statistics
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("repetitions must be positive")
    root = Path(__file__).resolve().parents[2]
    fixture = root / "tests/support/JointPartialSelection.hpp"
    cases = {}
    for case in ("all", "partial", "bridge", "weak"):
        runs = []
        for _ in range(args.repetitions):
            result = subprocess.run([str(args.executable.resolve()), case], capture_output=True, text=True, check=True)
            run = json.loads(result.stdout.splitlines()[-1])
            outcome = run.pop("outcome")
            run["costs"] = outcome["costs"]
            run["initialization"] = outcome["initialization"]
            run["software"] = outcome["metadata"]["software"]
            run["runtime_convergence"] = outcome["runtime_convergence"]
            run["search_completed"] = outcome["search_completed"]
            run["objective"] = outcome["objective"]
            run["components"] = [{key: c[key] for key in ("id", "stop_reason", "runtime_convergence", "state", "ranks", "evidence")} for c in outcome["components"]]
            runs.append(run)
        cases[case] = {"runs": runs, "median_total_seconds": statistics.median(r["total_seconds"] for r in runs),
                       "median_peak_rss_bytes": statistics.median(r["process_peak_rss_bytes"] for r in runs)}
    report = {"schema_version": 1, "measurement": "serial-independent-processes",
              "platform": platform.platform(), "fixture_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
              "scope": "Map/Model builder, model copy, contributor initialization, fitting and assembly; fixture generation precedes timing but contributes to process peak RSS; matched/omitted controls run after timing and RSS capture",
              "limits": "Small synthetic cases only; no maximum supported component size or total memory/time bound established.",
              "cases": cases}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
