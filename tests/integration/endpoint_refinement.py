#!/usr/bin/env python3
"""Gated testing-only endpoint refinement experiments; no production gate changes."""
from __future__ import annotations

import argparse
import collections
import json
import os
from pathlib import Path
import subprocess
import time
from types import SimpleNamespace

import mdpde_experiment as experiment
import fold_168_regression as fold

ROOT = Path(__file__).resolve().parents[2]
POLICIES = ("failed-only", "fresh-residual")


def read(path):
    return json.loads(Path(path).read_text())


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fresh(path):
    require(not path.exists(), f"Use a fresh output directory: {path}")
    path.mkdir(parents=True)


def distribution(values):
    values = sorted(values)
    return {"count": len(values), "min": values[0], "median": values[len(values)//2],
            "p99": values[min(len(values)-1, int(.99*len(values)))], "max": values[-1], "total": sum(values)}


def choose_budget(maximum):
    return next((n for n in (32, 64, 128, 256) if n >= 2*maximum), None)


def refine_run(args, directory, budget):
    fresh(directory)
    before = experiment.provenance(args.executable.resolve(), args.source.resolve())
    hashes = {str(p.resolve()): fold.sha256_file(p) for p in sorted(args.captures.glob("*.txt"))}
    command = [str(args.executable.resolve()), "refine", str(args.captures.resolve()), str(directory.resolve()), str(budget)]
    start = time.perf_counter()
    with (directory/"run.log").open("w") as log:
        code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    experiment.write(directory/"execution.json", {"command": command, "returncode": code, "seconds": time.perf_counter()-start})
    experiment.write(directory/"provenance.json", before)
    experiment.write(directory/"input-hashes.json", hashes)
    require(code == 0 and before == experiment.provenance(args.executable.resolve(),args.source.resolve()),
            "Refinement failed or source/executable changed.")
    require(all(fold.sha256_file(Path(p)) == h for p,h in hashes.items()), "Capture inputs changed.")
    index = read(directory/"index.json")
    require(index["shape_count"] == 337 and index["independent_count"] == 4 and len(index["files"]) == 341,
            "Stage A requires all 337 captures and four independent fixtures.")
    results = [read(directory/name) for name in index["files"]]
    require(all(r["exact_replay"] for r in results if r["fixture"].startswith("shape-")), "Inexact replay.")
    require(all(r["candidate_equation_evaluations"] <= budget for r in results), "Equation budget exceeded.")
    return results


def calibrate(args):
    fresh(args.output)
    results = refine_run(args, args.output/"profile", 2000)
    work = distribution([r["candidate_equation_evaluations"] for r in results])
    budget = choose_budget(work["max"])
    failures = [{"fixture": r["fixture"], "reason": r["reason"]} for r in results if not r["accepted"]]
    passed = not failures and budget is not None
    frozen = refine_run(args,args.output/"frozen",budget) if passed else []
    passed = passed and all(r["accepted"] and r["branch"]["pass"] for r in frozen)
    report = {"schema_version":1, "stage":"A", "passed":passed, "budget":budget,
              "budget_rule":"smallest of 32,64,128,256 >= twice maximum observed equation work",
              "equation_tolerance":1e-8,"reference_tolerance":1e-10,"branch_tolerance":1e-6,
              "captures":337,"independent_fixtures":4,"profile_work":work,"profile_failures":failures,
              "frozen_accepted":sum(r["accepted"] for r in frozen),
              "reference_updates":distribution([r["reference_updates"] for r in results]),
              "profile_provenance_sha256":fold.sha256_file(args.output/"profile/provenance.json")}
    experiment.write(args.output/"stage-a.json",report)
    print(json.dumps(report))
    require(passed,"Stage A did not pass; operator and closed-loop experiments are not authorized by this gate.")


def validate_refinements(records, budget):
    failures = collections.Counter()
    accepted = 0
    for record in records:
        if not record.get("triggered"):
            continue
        r = record["refinement"]
        require(r["candidate_equation_evaluations"] <= budget, "Live refinement exceeded the frozen budget.")
        if r["accepted"]:
            require(r["branch"]["pass"] and r["root"]["equation_pass"] and r["reference"]["reference_pass"],
                    "Accepted refinement has incomplete equation/branch evidence.")
            require(record["effective_status"] == 0,"Accepted endpoint is not qualified.")
            accepted += 1
        else:
            require(record["effective_status"] != 0,"Rejected refinement inherited SUCCESS.")
            failures[r["reason"]] += 1
    return {"accepted":accepted,"rejections":dict(failures)}


def compare_report(directory, budget):
    checkpoints = {}
    passed = {p:True for p in POLICIES}
    for phase in ("final", "recovery-current"):
        path = directory/"endpoint"/phase
        data = read(path/"comparison.json")
        before, after = fold.sha256_file(path/"input-before.json"), fold.sha256_file(path/"input-after.json")
        require(data["input_unchanged"] and before == after,"Operator comparison changed its input.")
        variants = {v["policy"]:v for v in data["variants"]}
        require(set(variants) == {"legacy",*POLICIES},"Missing operator policies.")
        require(variants["legacy"]["legacy_exact"],"Legacy operator replay differs from the original evaluation.")
        summary = {"input_sha256":before,"legacy_exact":True,"variants":{}}
        for name, v in variants.items():
            validation = validate_refinements(v["solves"],budget)
            require(v["offset_solves_equal"] and all(row[2] == 0 for row in v["delta_abc_from_legacy"]),
                    "Same-state joint-offset solves changed.")
            if name in passed:
                passed[name] &= v["qualified"] and v["complete"] and not validation["rejections"]
            summary["variants"][name] = {k:v[k] for k in ("qualified","complete","original_solver_qualified",
                "nominal_p99","nominal_max","recovery_residual_mean_square","candidate_and_diagnostic_equations","reference_updates")}
            summary["variants"][name].update(validation)
            summary["variants"][name]["delta_abc_max"] = [max(abs(row[k]) for row in v["delta_abc_from_legacy"]) for k in range(3)]
        checkpoints[phase] = summary
    return {"schema_version":1,"stage":"B","passed_policies":[p for p in POLICIES if passed[p]],"checkpoints":checkpoints}


def run(args):
    stage_a = read(args.stage_a) if args.stage_a else None
    if args.mode != "legacy":
        require(stage_a and stage_a["passed"],"Stage A must pass first.")
    if args.mode in POLICIES:
        require(args.stage_b and args.mode in read(args.stage_b)["passed_policies"],"This policy has not passed stage B.")
    budget = stage_a["budget"] if stage_a else 0
    fresh(args.output)
    variables = {"RHBM_TEST_ENDPOINT_DIR":str((args.output/"endpoint").resolve()),
                 "RHBM_TEST_ENDPOINT_POLICY":"legacy" if args.mode == "compare" else args.mode,
                 "RHBM_TEST_ENDPOINT_BUDGET":str(budget),"RHBM_TEST_ENDPOINT_COMPARE":"1" if args.mode == "compare" else "0"}
    saved = {k:os.environ.get(k) for k in variables}
    try:
        os.environ.update(variables)
        experiment.capture(SimpleNamespace(executable=args.executable,source=args.source,model=args.model,map=args.map,
            output=args.output,verbosity=args.verbosity,jobs=args.jobs,no_capture=False,score_existing=False))
    finally:
        for k,v in saved.items():
            if v is None: os.environ.pop(k,None)
            else: os.environ[k]=v
    experiment.write(args.output/"policy.json",{"mode":args.mode,"budget":budget,"environment":variables,
        "stage_a_sha256":fold.sha256_file(args.stage_a) if args.stage_a else None,
        "stage_b_sha256":fold.sha256_file(args.stage_b) if args.stage_b else None})
    session = read(args.output/"endpoint/session.json")
    require(session["complete"],"Endpoint session did not complete.")
    records = [json.loads(line) for line in (args.output/"endpoint/solves.jsonl").read_text().splitlines()]
    validation = validate_refinements(records,budget)
    if args.reference:
        if args.mode in ("legacy","compare"):
            experiment.verify(SimpleNamespace(reference=args.reference,runs=[args.output],output=args.output/"neutrality.json"))
        if (args.reference/"endpoint/initial.json").exists():
            initial, expected = read(args.output/"endpoint/initial.json"), read(args.reference/"endpoint/initial.json")
            initial.pop("options"); expected.pop("options")
            require(initial == expected,"Second-stage inputs differ before applying the policy.")
    if args.mode == "compare":
        require(args.reference,"Operator comparison requires a baseline run.")
        report = compare_report(args.output,budget)
        experiment.write(args.output/"stage-b.json",report)
        print(json.dumps(report))
    else:
        actual = read(args.output/"actual.json")
        experiment.write(args.output/"endpoint-summary.json",{"session":session,"validation":validation,
            "summary":actual["second_stage_summary"],"certificate":actual["production_fitting"],"quality":actual["quality_metrics"]})


def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest="operation",required=True)
    p=sub.add_parser("calibrate")
    for name in ("executable","captures","output"): p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--source",type=Path,default=ROOT)
    p=sub.add_parser("run")
    for name in ("executable","model","map","output"): p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--source",type=Path,default=ROOT)
    p.add_argument("--mode",choices=("legacy","compare",*POLICIES),required=True)
    p.add_argument("--stage-a",type=Path); p.add_argument("--stage-b",type=Path); p.add_argument("--reference",type=Path)
    p.add_argument("--jobs",type=int,default=4); p.add_argument("--verbosity",type=int,default=3)
    args=parser.parse_args(); globals()[args.operation](args)


if __name__ == "__main__":
    main()
