"""Qualify fixed-neighbor coordinate-confirmation candidates against endpoint checks."""

import argparse
import csv
import json
from pathlib import Path


EXPECTED_CASES = {(topology, atoms) for topology in ("chain", "cube") for atoms in (32, 128, 256)}
ETA_THRESHOLD = 1e-10
BETA_THRESHOLD = 1e-10
KKT_THRESHOLD = 1e-10
WIDTH_THRESHOLD = 1e-12
RULES = ("eta_only", "beta_only", "eta_beta")


def _first_sweep(sweeps, predicate):
    return next((sweep["sweep"] for sweep in sweeps if predicate(sweep)), None)


def analyze(cases, expected=EXPECTED_CASES):
    rows = []
    seen = set()
    for case in cases:
        topology, atoms = case["topology"], case["atoms"]
        seen.add((topology, atoms))
        fixed = case["fixed_neighbor"]
        sweeps = [dict(item, sweep=item.get("sweep", index))
            for index, item in enumerate(fixed["sweep_telemetry"], start=1)]
        assessments = {item["sweep"]: item for item in case["endpoint_assessment_by_sweep"]}
        cheap = _first_sweep(sweeps, lambda item:
            item["global_ac_kkt"] <= KKT_THRESHOLD and
            item["global_width_gradient_inf_norm"] <= WIDTH_THRESHOLD)
        certified = _first_sweep(sweeps, lambda item:
            item["sweep"] in assessments and
            all(assessments[item["sweep"]][name] is True
                for name in ("inner", "gradient", "local", "identified", "endpoint_trust")) and
            assessments[item["sweep"]]["runtime_convergence"] == "Passed")

        def confirmed(item, require_eta, require_beta):
            if item["global_ac_kkt"] > KKT_THRESHOLD or item["global_width_gradient_inf_norm"] > WIDTH_THRESHOLD:
                return False
            if not item["coordinate_confirmation_available"]:
                return False
            return ((not require_eta or item["eta_change_inf"] <= ETA_THRESHOLD) and
                    (not require_beta or item["beta_scaled_change"] <= BETA_THRESHOLD))

        eta_only = _first_sweep(sweeps, lambda item: confirmed(item, True, False))
        beta_only = _first_sweep(sweeps, lambda item: confirmed(item, False, True))
        eta_beta = _first_sweep(sweeps, lambda item: confirmed(item, True, True))
        candidate_sweeps = {"eta_only": eta_only, "beta_only": beta_only, "eta_beta": eta_beta}
        candidate_results = {}
        for rule, stop_sweep in candidate_sweeps.items():
            precedes = stop_sweep < certified if stop_sweep is not None and certified is not None else None
            delay = stop_sweep - certified if stop_sweep is not None and certified is not None else None
            candidate_results[rule] = {
                "confirmation_sweep": stop_sweep,
                "confirmation_precedes_certification": precedes,
                "confirmation_delay": delay,
            }

        endpoint_by_sweep = []
        for sweep in sweeps:
            endpoint = assessments.get(sweep["sweep"])
            if endpoint is None:
                endpoint_by_sweep.append({"sweep": sweep["sweep"], "endpoint_certified": False,
                    "runtime_convergence": None})
            else:
                endpoint_by_sweep.append({"sweep": sweep["sweep"],
                    "inner": endpoint["inner"], "gradient": endpoint["gradient"],
                    "local": endpoint["local"], "identified": endpoint["identified"],
                    "endpoint_trust": endpoint["endpoint_trust"],
                    "endpoint_certified": all(endpoint[name] is True for name in
                        ("inner", "gradient", "local", "identified", "endpoint_trust")),
                    "fully_certified": all(endpoint[name] is True for name in
                        ("inner", "gradient", "local", "identified", "endpoint_trust")) and
                        endpoint["runtime_convergence"] == "Passed",
                    "runtime_convergence": endpoint["runtime_convergence"]})

        rows.append({"topology": topology, "atoms": atoms,
            "cheap_stationarity_sweep": cheap,
            "first_certified_sweep": certified,
            "eta_only_confirmation_sweep": eta_only,
            "beta_only_confirmation_sweep": beta_only,
            "eta_beta_confirmation_sweep": eta_beta,
            "candidate_results": candidate_results,
            "search_sweeps_executed": fixed["sweeps"],
            "endpoint_certification_by_sweep": endpoint_by_sweep,
            "runtime_convergence_by_sweep": [
                {"sweep": item["sweep"], "runtime_convergence": item["runtime_convergence"]}
                for item in endpoint_by_sweep],
            "final_runtime_convergence": fixed["runtime_convergence"]})

    missing = sorted(expected - seen)
    rule_safety = {}
    for rule in RULES:
        checks = [row["candidate_results"][rule] for row in rows]
        rule_safety[rule] = {
            "safe_across_cases": bool(rows) and not missing and all(
                check["confirmation_sweep"] is not None and
                check["confirmation_precedes_certification"] is False for check in checks),
            "unsafe_cases": [
                {"topology": row["topology"], "atoms": row["atoms"]}
                for row in rows
                if row["candidate_results"][rule]["confirmation_precedes_certification"] is True],
            "missing_confirmation_cases": [
                {"topology": row["topology"], "atoms": row["atoms"]}
                for row in rows if row["candidate_results"][rule]["confirmation_sweep"] is None],
        }
    return {
        "qualification": "passed" if any(item["safe_across_cases"] for item in rule_safety.values()) else "failed",
        "expected_cases": [{"topology": topology, "atoms": atoms} for topology, atoms in sorted(expected)],
        "missing_cases": [{"topology": topology, "atoms": atoms} for topology, atoms in missing],
        "thresholds": {"global_ac_kkt": KKT_THRESHOLD, "width_gradient_inf_norm": WIDTH_THRESHOLD,
            "eta_change_inf": ETA_THRESHOLD, "beta_scaled_change": BETA_THRESHOLD},
        "candidate_rule_safety": rule_safety,
        "cases": sorted(rows, key=lambda row: (row["topology"], row["atoms"])),
    }


def write_report(report, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    columns = ["topology", "atoms", "cheap_stationarity_sweep", "first_certified_sweep",
        "eta_only_confirmation_sweep", "beta_only_confirmation_sweep", "eta_beta_confirmation_sweep",
        "search_sweeps_executed", "final_runtime_convergence"]
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in report["cases"]:
            writer.writerow({column: row[column] for column in columns})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="qualification artifact directory")
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    cases = [json.loads(Path(path).read_text()) for path in args.cases]
    report = analyze(cases)
    write_report(report, Path(args.output))
    print(report["qualification"])
    return 0 if report["qualification"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
