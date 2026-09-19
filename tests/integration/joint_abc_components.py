#!/usr/bin/env python3
"""Frozen-snapshot exact-component experiment and explicit scientific gates."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

NUMERIC_ENV = {name: "1" for name in
               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
os.environ.update(NUMERIC_ENV)

import numpy as np

import joint_abc_certification as certification
import joint_abc_component_records as records
from joint_abc_profile import scientific
from observation_matching import read, require

ROOT = Path(__file__).resolve().parents[2]
VARIANTS = ("legacy", "guarded", "guarded-log")
COMPOSITES = {
    "regular-two": ("baseline", "weak-1e-2"),
    "regular-three": ("baseline", "baseline", "near-0.10"),
    "regular-near": ("baseline", "near-0.02"),
    "regular-weak": ("baseline", "weak-1e-4"),
    "regular-active": ("baseline", "active-a"),
    "regular-zero": ("baseline", "zero-signal"),
    "regular-duplicate": ("baseline", "duplicate"),
}


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance(executable):
    files = [*sorted((ROOT/"tests/support").glob("*.cpp")),
             *sorted((ROOT/"tests/support").glob("*.hpp")),
             *sorted((ROOT/"tests/utils/hrl").glob("JointABC*.cpp")),
             ROOT/"tests/experiments/mdpde_experiment.cpp", *sorted((ROOT/"tests/integration").glob("*.py")),
             ROOT/"tests/CMakeLists.txt", ROOT/"docs/developer/joint_abc_components_contract.md"]
    return {"executable_sha256": sha(executable),
            "shared_libraries": {p.name: sha(p) for pattern in ("*.dylib", "*.so") for p in sorted((executable.parent.parent/"src").glob(pattern)) if p.is_file()},
            "source_hashes": {str(p.relative_to(ROOT)): sha(p) for p in files}}


def execute(command, log):
    env = dict(os.environ, **NUMERIC_ENV)
    start = time.monotonic()
    with log.open("w") as stream:
        subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT, env=env, check=True)
    return time.monotonic()-start


def regression(args):
    baseline, output, executable = args.baseline.resolve(), args.output.resolve(), args.executable.resolve()
    require(not output.exists(), "Use a fresh regression directory.")
    datasets = sorted(p for p in (baseline/"datasets").iterdir() if p.is_dir())
    require(len(datasets) == 9, "Regression requires all nine frozen datasets.")
    output.mkdir(parents=True); before = provenance(executable); write(output/"provenance.json", before)

    def run_one(source):
        seconds = execute([executable, "joint-abc-component-regression", source, output/"datasets"/source.name], output/(source.name+".log"))
        print(f"Completed regression: {source.name}", flush=True)
        return {"dataset": source.name, "seconds": seconds}

    with ThreadPoolExecutor(max_workers=3) as pool:
        costs = list(pool.map(run_one, datasets))
    require(before == provenance(executable), "Regression sources or executable changed.")
    differences, files, counts, census = [], 0, {v: 0 for v in VARIANTS}, []
    for source in datasets:
        target = output/"datasets"/source.name
        census.append(read(target/"census.json"))
        for variant in VARIANTS:
            for fit_path in sorted((target/variant/"fits").glob("*.json")):
                relative = fit_path.relative_to(target)
                fit, frozen = read(fit_path), read(source/relative)
                audit_path = target/variant/"audits"/fit_path.name
                audit = read(audit_path)
                for p, a, b in ((relative, fit, frozen), (audit_path.relative_to(target), audit, read(source/audit_path.relative_to(target)))):
                    files += 1
                    if scientific(a) != scientific(b): differences.append(f"{source.name}/{p}")
                old_certificate = read(source/variant/"audits"/(fit_path.stem+"-certificate.json"))
                certificate = certification.certificate(fit, audit, old_certificate["checks"]["endpoint_replay"])
                write(audit_path.with_name(fit_path.stem+"-certificate.json"), certificate)
                if scientific(certificate) != scientific(old_certificate): differences.append(f"{source.name}/{variant}/{fit_path.stem}-certificate")
                counts[variant] += certificate["regular_qualified"]
    result = {"schema_version": 1, "passed": not differences and counts == {"legacy": 38, "guarded": 40, "guarded-log": 40},
              "fit_count": files//2, "compared_fit_and_audit_files": files, "regular_certificates": counts,
              "differences": differences, "costs": costs,
              "replay_evidence": "Frozen independent raw-row replay; exact trial and endpoint identity required before reuse."}
    write(output/"regression.json", result); write(output/"census.json", census)
    print(json.dumps({k: v for k, v in result.items() if k != "costs"}), flush=True)
    require(result["passed"], "Monolithic regression differs; inspect regression.json.")


def inputs(baseline, output):
    """Copy immutable observations and recorded starts; compose only identities."""
    root = output/"inputs"
    require({p.name for p in (baseline/"datasets").iterdir()} == set(certification.coverage.DATASETS), "Expected the nine original frozen datasets.")
    for source in sorted((baseline/"datasets").iterdir()):
        target = root/source.name
        target.mkdir(parents=True)
        for name in ("snapshot.json", "dataset.json", "voxels.csv", "contributors.csv"):
            shutil.copyfile(source/name, target/name)
        shutil.copytree(source/"guarded/fits", target/"guarded/fits")
        for certificate in sorted((source/"guarded/audits").glob("*-certificate.json")):
            destination = target/"guarded/certificates"/(certificate.stem.removesuffix("-certificate")+".json")
            destination.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(certificate, destination)
    for name, sources in COMPOSITES.items():
        write(root/name/"composition.json", {
            "schema_version": 1, "name": name, "kind": "disjoint-frozen-snapshot-composition",
            "sources": [{"path": "../"+source, "identity_prefix": f"{i:02d}:{source}:",
                         "snapshot_sha256": sha(root/source/"snapshot.json")}
                        for i, source in enumerate(sources)]})
    return sorted(root.iterdir())


def same_state(args):
    output, executable = args.output.resolve(), args.executable.resolve()
    require(not output.exists(), "Use a fresh same-state directory.")
    output.mkdir(parents=True)
    datasets = inputs(args.baseline.resolve(), output)
    before = provenance(executable); write(output/"provenance.json", before)

    def run_one(source):
        seconds = execute([executable, "joint-abc-component-same-state", source, output/"datasets"/source.name], output/(source.name+".log"))
        print(f"Completed same-state: {source.name}", flush=True)
        return {"dataset": source.name, "seconds": seconds}

    with ThreadPoolExecutor(max_workers=3) as pool:
        costs = list(pool.map(run_one, datasets))
    require(before == provenance(executable), "Same-state sources or executable changed.")
    results = {p.name: read(output/"datasets"/p.name/"completion.json") for p in datasets}
    for p in datasets:
        records.verify_census(records.load(p), read(output/"datasets"/p.name/"census.json"))
    failures = [f"{name}/{row['id']}" for name, result in results.items() for row in result["states"] if not row["passed"]]
    result = {"schema_version": 1, "passed": not failures, "failures": failures,
              "states": sum(len(r["states"]) for r in results.values()),
              "full_equivalence": sum(row["full_equivalence"] for r in results.values() for row in r["states"]),
              "costs": costs}
    write(output/"same-state.json", result); print(json.dumps(result), flush=True)
    require(result["passed"], "Same-state gate differs; inspect scientific records.")


def inventory(root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


def run(args):
    output, executable = args.output.resolve(), args.executable.resolve()
    require(not output.exists(), "Use a fresh experiment directory.")
    gate = read(args.regression_gate)
    require(gate["passed"] and gate["fit_count"] == 216 and gate["regular_certificates"]["guarded"] == 40,
            "The 216-case monolithic regression gate must pass first.")
    output.mkdir(parents=True); datasets = inputs(args.baseline.resolve(), output)
    census = read(args.regression_gate.parent/"census.json")
    require({row["snapshot_sha256"] for row in census} ==
            {sha(p/"snapshot.json") for p in datasets if p.name not in COMPOSITES}, "Regression gate belongs to different observations.")
    before = provenance(executable)
    write(output/"provenance.json", dict(before, regression_gate_sha256=sha(args.regression_gate),
          numeric_threads=1, independent_processes=3, audit_cache="exact-input precision references within one case; no cross-run reuse"))
    write(output/"input-hashes.json", inventory(output/"inputs"))
    for relative in before["source_hashes"]:
        target = output/"source"/relative; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/relative, target)

    def run_one(source):
        seconds = execute([executable, "joint-abc-components-run", source, output/"datasets"/source.name], output/(source.name+"-run.log"))
        print(f"Completed searches: {source.name}", flush=True)
        audit_seconds = execute([executable, "joint-abc-components-audit", source, output/"datasets"/source.name,
                                 output/"audits"/source.name], output/(source.name+"-audit.log"))
        print(f"Completed fresh audits: {source.name}", flush=True)
        return {"dataset": source.name, "run_seconds": seconds, "audit_seconds": audit_seconds}

    with ThreadPoolExecutor(max_workers=3) as pool:
        costs = list(pool.map(run_one, datasets))
    require(before == provenance(executable), "Experiment sources or executable changed during execution.")
    write(output/"completion.json", {"complete": True, "datasets": [p.name for p in datasets], "costs": costs})
    summarize(output)


def compatible_search_kernels(run_root, current):
    previous = read(run_root/"provenance.json")["source_hashes"]
    keys = [key for key in previous if key.startswith("tests/support/") and
            key not in ("tests/support/JointABCComponentExperiment.cpp", "tests/support/JointABCComponentExperiment.hpp")]
    require(all(previous[key] == current["source_hashes"][key] for key in keys),
            "Search numerical kernels changed; rerun the searches before comparing components.")
    return keys


def audit(args):
    run_root, output, executable = args.run.resolve(), args.output.resolve(), args.executable.resolve()
    require(not output.exists(), "Use a fresh audit directory.")
    require(inventory(run_root/"inputs") == read(run_root/"input-hashes.json"), "Changed frozen inputs.")
    start = time.monotonic(); before = provenance(executable); kernels = compatible_search_kernels(run_root, before)
    output.mkdir(parents=True)
    write(output/"provenance.json", dict(before, search_provenance_sha256=sha(run_root/"provenance.json"),
          unchanged_search_kernels=kernels, numeric_threads=1, independent_processes=3,
          audit_cache="fresh case-local references; no reuse from previous audits or another run"))
    for relative in before["source_hashes"]:
        target = output/"source"/relative; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/relative, target)

    datasets = sorted((run_root/"inputs").iterdir())
    cases = [(source, case) for source in datasets for case in read(run_root/"datasets"/source.name/"completion.json")["cases"]]
    # Every process owns a case and starts with empty reference caches. This
    # also prevents one expensive boundary case from serializing an entire run.
    def run_one(item):
        source, case = item
        seconds = execute([executable, "joint-abc-components-audit", source, run_root/"datasets"/source.name, output/source.name, case],
                          output/(source.name+"-"+case+".log"))
        print(f"Completed endpoint audit: {source.name}/{case}", flush=True)
        return {"dataset": source.name, "case": case, "audit_seconds": seconds}

    with ThreadPoolExecutor(max_workers=3) as pool:
        costs = list(pool.map(run_one, cases))
    for source in datasets:
        write(output/source.name/"completion.json", {"complete": True,
              "cases": [case for path, case in cases if path == source],
              "cache_policy": "fresh process per case; no cross-case or cross-run cache"})
    require(before == provenance(executable), "Audit sources or executable changed during execution.")
    write(output/"completion.json", {"complete": True, "costs": costs, "seconds": time.monotonic()-start})
    summarize(run_root, output)


def subset(data, component):
    atoms, rows = component["atoms"], component["rows"]
    table = data["table"][np.isin(data["table"]["atom"], atoms)].copy()
    atom_map = np.full(len(data["ids"]), -1); atom_map[atoms] = np.arange(len(atoms))
    row_map = np.full(len(data["y64"]), -1); row_map[rows] = np.arange(len(rows))
    table["atom"] = atom_map[table["atom"]]; table["row"] = row_map[table["row"]]
    return dict(data, table=table, ids=[data["ids"][a] for a in atoms], y64=data["y64"][rows], y32=data["y32"][rows])


def replay(data, y, endpoint, scale):
    result = certification.snapshot_replay(data["table"], len(data["ids"]), y, endpoint)
    ratio = max(1., np.linalg.norm(y))/scale
    result["projected_kkt"] *= ratio; result["b_gradient"] *= ratio*ratio
    return result


def replay_passed(data, y, endpoint, scale):
    if not endpoint or not endpoint.get("valid"): return None
    raw = replay(data, y, endpoint, scale)
    return bool(abs(raw["projected_kkt"]-endpoint["projected_kkt"]) <= 1e-13 and
                np.allclose(raw["b_gradient"], endpoint["b_gradient"], rtol=2e-9, atol=1e-13) and
                abs(np.linalg.norm(raw["residual"])-np.sqrt(endpoint["rss"])) <= 512*np.finfo(float).eps*scale)


def certify(data, y, fit, audit_record, scope):
    certification.validate_audit(audit_record, fit)
    scale = fit["residual_scale"]
    for trial in fit.get("trials", []):
        if trial["accepted"]:
            require(trial.get("trust", {}).get("passed") is True, "Accepted an untrusted component trial.")
            require(replay_passed(data, y, trial, scale), "Accepted trial failed independent parent-scale replay.")
    result = certification.certificate(fit, audit_record, replay_passed(data, y, fit.get("primary"), scale))
    if scope == "assembled-global":
        result["checks"]["assembled_profile"] = fit.get("assembled_profile_agrees")
        result["checks"]["usable_global_state"] = fit["prediction_available"]
        result["failures"] = [k for k, v in result["checks"].items() if v is False]
        result["unavailable"] = [k for k, v in result["checks"].items() if v is None]
        result["regular_qualified"] = all(v is True for v in result["checks"].values())
    result["scope"] = scope
    return result


def endpoint_parity(data, y, monolithic, assembled, certificates, required):
    both = all(c["regular_qualified"] for c in certificates)
    row = {"required_regular": required, "both_regular": both, "available": False,
           "passed": not (required or both), "reason": "nonregular-or-unavailable; see certificates"}
    a, b = monolithic.get("primary"), assembled.get("assembled_state")
    if a and b and a["valid"] and b["valid"]:
        def difference(x, z):
            x, z = np.asarray(x), np.asarray(z)
            return float(np.max(np.abs(x-z)/(1+np.maximum(np.abs(x), np.abs(z)))))
        scale = monolithic["residual_scale"]
        pa, pb = replay(data, y, a, scale), replay(data, y, b, scale)
        values = {"scaled_ac_difference": difference(a["beta"], b["beta"]),
                  "scaled_b_difference": difference(a["b"], b["b"]),
                  "maximum_log_b_difference": float(np.max(np.abs(np.array(a["eta"])-b["eta"]))),
                  "normalized_prediction_inf_difference": float(np.max(np.abs(pa["prediction"]-pb["prediction"]))/scale),
                  "normalized_objective_difference": float(abs(pa["residual"]@pa["residual"]-pb["residual"]@pb["residual"])/2/scale**2)}
        numeric = all(values[k] <= (1e-12 if k == "normalized_objective_difference" else 1e-8) for k in values)
        row.update(values, available=True, numeric_passed=numeric,
                   passed=(both and numeric) if required or both else True,
                   reason="regular-parity" if both and numeric else "missing-regular-certificate" if required and not both else
                          "search-endpoint-difference" if both and not numeric else "nonregular-control")
    return row


def csv_write(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys); writer.writeheader(); writer.writerows(rows)


def validate_endpoint_scope(fit, context, global_context=None, atoms=None):
    directions = np.asarray(context["audit"]["directions"])
    if global_context is not None:
        expected = np.asarray(global_context["audit"]["directions"])[:, atoms]
        require(np.array_equal(directions, expected), "Component audit directions were changed or renormalized.")
    elif "width_spectrum" in fit:
        require(np.array_equal(directions[2], fit["width_spectrum"]["weak_directions"][0]),
                "Final audit does not test the endpoint weakest direction.")
    return True


def setup_cost(dataset, case, target, initial_state):
    setup = read(target/"setup.json"); primary = references = 0
    if setup["direction_reference"].startswith("frozen-initial-global-profile"):
        context = read(target/"context.json"); context["audit"]["directions"] = []
        require(context == initial_state["context"], "Setup count evidence uses a different initial context.")
        primary = references = 1
        if all(row["monolithic"]["valid"] for row in initial_state["profile_solves"]):
            derivative = initial_state.get("derivative", {})
            require("monolithic_reason" in derivative, "Initial differential evidence is unavailable for setup counts.")
            if derivative["monolithic_reason"] == "full-profile-derivative": references += 12
    return {"dataset": dataset, "case": case, "primary_profile_calls": primary,
            "reference_profile_calls": references, "total_profile_calls": primary+references,
            "seconds": setup["seconds"],
            "count_evidence": "Assess control flow replayed from the identical saved initial-state validity and differential evidence; separate from search budgets"}


def summarize(run_root, audit_root=None):
    run_root = Path(run_root).resolve(); audit_root = Path(audit_root or run_root/"audits")
    require(read(run_root/"completion.json")["complete"], "Incomplete component run.")
    require(inventory(run_root/"inputs") == read(run_root/"input-hashes.json"), "Changed frozen inputs.")
    output = run_root/"summary" if audit_root.resolve() == (run_root/"audits").resolve() else audit_root/"summary"
    output.mkdir(exist_ok=True)
    pairs, failures, costs, states, regressions, census_rows, endpoint_states, setup_costs = [], [], [], [], [], [], [], []
    regular = {"original_monolithic": 0, "original_assembled": 0, "composite_monolithic": 0, "composite_assembled": 0}
    for source in sorted((run_root/"inputs").iterdir()):
        data = records.load(source); root = run_root/"datasets"/source.name
        census = read(root/"census.json"); records.verify_census(data, census)
        census_rows.append({k: v for k, v in census.items() if k in ("atoms", "rows", "memberships", "component_count", "largest_atom_count", "largest_row_count")}|{"dataset": source.name, "constant_row_count": len(census["constant_rows"]),
                       "unobserved_atom_count": len(census["unobserved_atoms"]),
                       "component_sizes": ";".join(str(c["atom_count"])+" atoms/"+str(c["row_count"])+" rows" for c in census["components"])})
        require(read(root/"same-state/completion.json")["passed"], "Same-state prerequisite failed.")
        for path in sorted((root/"same-state").glob("*-initial.json"))+sorted((root/"same-state").glob("*-endpoint.json")):
            row = read(path)
            correction = row.get("derivative", {}).get("local_correction_scaled_difference")
            require(correction is None or correction <= 1e-10, "Same-state local correction disagreement.")
            states.append({"dataset": source.name, "state": path.stem, "passed": row["passed"], "full_equivalence": row["full_equivalence"]})
        for case in read(root/"completion.json")["cases"]:
            target = root/"cases"/case; y = data["y64" if case.endswith("double") else "y32"]
            setup_costs.append(setup_cost(source.name, case, target, read(root/"same-state"/(case+"-initial.json"))))
            fits = [read(target/(v+"-fit.json")) for v in ("monolithic", "assembled")]
            required = all(read(p/"guarded/certificates"/(case+".json"))["regular_qualified"] for p in data["sources"])
            certs = []
            for variant, fit in zip(("monolithic", "assembled"), fits):
                audit_dir = audit_root/source.name/case/variant; audit_record = read(audit_dir/"audit.json")
                scoped_fit = read(audit_dir/"fit.json")
                # The direction was selected once from this frozen reference;
                # a repeated singular value need not retain a second SVD basis.
                validate_endpoint_scope(fit, read(audit_dir/"context.json"))
                certificate = certify(data, y, scoped_fit, audit_record, "assembled-global" if variant == "assembled" else "monolithic-global")
                write(audit_dir/"certificate.json", certificate); certs.append(certificate)
                kind = "composite" if source.name in COMPOSITES else "original"
                regular[kind+"_"+variant] += certificate["regular_qualified"]
                failures.append({"dataset": source.name, "case": case, "scope": variant,
                                 "regular": certificate["regular_qualified"], "search_completed": not fit.get("search_stopped_without_convergence", False),
                                 "usable_state": fit.get("usable_state", fit.get("prediction_available", fit.get("primary", {}).get("valid", False))),
                                 "failures": ";".join(certificate["failures"]),
                                 "unavailable": ";".join(certificate["unavailable"]), "derivative_status": certificate["derivative_status"]})
                parity = read(target/(variant+"-same-state.json"))
                if parity.get("available", True): require(parity["passed"], f"Endpoint same-state disagreement: {source.name}/{case}/{variant}")
                endpoint_states.append({"dataset": source.name, "case": case, "variant": variant,
                                        "available": parity.get("available", True), "passed": parity.get("passed"),
                                        "full_equivalence": parity.get("full_equivalence", False),
                                        "limitation": parity.get("reason", parity.get("derivative", {}).get("limitation", ""))})
                costs.append(cost_row(source.name, case, variant, fit, audit_record, len(data["ids"]), len(y),
                                      scoped_fit, read(audit_dir/"boundary.json") if (audit_dir/"boundary.json").exists() else {}))
            if source.name not in COMPOSITES:
                frozen = read(source/"guarded/fits"/(case+".json"))
                differences = [key for key in frozen.keys() & fits[0].keys() if scientific(frozen[key]) != scientific(fits[0][key]) and
                               not key.endswith("seconds") and key != "process_peak_rss_bytes"]
                if differences: regressions.append({"dataset": source.name, "case": case, "fields": differences})
            pairs.append({"dataset": source.name, "case": case, **endpoint_parity(data, y, *fits, certs, required)})
            for k, component in enumerate(census["components"]):
                local = subset(data, component); local_y = local["y64" if case.endswith("double") else "y32"]
                fit = read(target/"components"/(str(k)+"-fit.json")); audit_dir = audit_root/source.name/case/"components"/str(k)
                audit_record = read(audit_dir/"audit.json"); scoped_fit = read(audit_dir/"fit.json")
                validate_endpoint_scope(scoped_fit, read(audit_dir/"context.json"),
                                        read(audit_root/source.name/case/"assembled/context.json"), component["atoms"])
                certificate = certify(local, local_y, scoped_fit, audit_record, "component-search-policy-supplement")
                write(audit_dir/"certificate.json", certificate)
                failures.append({"dataset": source.name, "case": case, "scope": "component:"+component["id"],
                                 "regular": certificate["regular_qualified"], "search_completed": not fit.get("search_stopped_without_convergence", False),
                                 "usable_state": fit.get("usable_state", fit.get("prediction_available", fit.get("primary", {}).get("valid", False))),
                                 "failures": ";".join(certificate["failures"]),
                                 "unavailable": ";".join(certificate["unavailable"]), "derivative_status": certificate["derivative_status"]})
                costs.append(cost_row(source.name, case, component["id"], fit, audit_record, len(local["ids"]), len(local_y),
                                      scoped_fit, read(audit_dir/"boundary.json") if (audit_dir/"boundary.json").exists() else {}))
    summary = {"schema_version": 1, "passed": not regressions and all(p["passed"] for p in pairs) and
               regular["original_monolithic"] == 40 and regular["original_assembled"] == 40,
               "regular_certificates": regular, "same_state_records": len(states),
               "full_same_state_equivalence": sum(s["full_equivalence"] for s in states),
               "endpoint_same_state_records": len(endpoint_states),
               "endpoint_full_same_state_equivalence": sum(s["full_equivalence"] for s in endpoint_states),
               "paired_branches": len(pairs), "required_regular_pairs": sum(p["required_regular"] for p in pairs),
               "failed_pairs": [p for p in pairs if not p["passed"]], "historical_differences": regressions}
    write(output/"summary.json", summary)
    for name, rows in (("census", census_rows), ("same-state-parity", states), ("endpoint-parity", pairs), ("same-state-at-fit-endpoints", endpoint_states), ("failure-matrix", failures), ("costs", costs), ("context-setup-costs", setup_costs)):
        csv_write(output/(name+".csv"), rows)
    print(json.dumps(summary), flush=True)
    require(summary["passed"], "Component endpoint experiment has unresolved gates; inspect summary.")
    return summary


def cost_row(dataset, case, scope, fit, audit_record, atoms, rows, scoped_fit, boundary):
    resources = fit.get("resources", {})
    return {"dataset": dataset, "case": case, "scope": scope, "atoms": atoms, "rows": rows,
            "search_counter_scope": "sum-of-components" if scope == "assembled" else "this-search",
            "profile_evaluations": fit.get("profile_evaluations", 0), "accepted_updates": fit.get("accepted_updates", 0),
            "search_reference_solves": fit.get("search_reference_evaluations", 0),
            "endpoint_reference_solves": (1+fit.get("directional_evaluations", 0)) if "primary" in fit else 0,
            "audit_reference_solves": (1+2*sum(len(l["samples"]) for l in audit_record.get("ladders", []))) if "endpoint_trust" in audit_record else 0,
            "post_search_trust_reference_solves": int("endpoint_trust" in fit),
            "refreshed_assessment_reference_solves": (1+scoped_fit.get("directional_evaluations", 0)) if scoped_fit.get("audit_assessment_refreshed") and fit.get("primary", {}).get("valid") else 0,
            "precision_reference_pairs_requested": int("precision" in audit_record),
            "boundary_reference_pairs_requested": int(bool(boundary)),
            "boundary_seconds": boundary.get("seconds", 0),
            "refreshed_assessment_seconds": scoped_fit.get("audit_assessment_seconds", 0),
            "search_seconds": resources.get("search", {}).get("seconds", 0),
            "search_reference_seconds_nested": fit.get("search_reference_seconds", 0),
            "endpoint_audit_seconds": resources.get("endpoint_audit", {}).get("seconds", fit.get("assembly_seconds", 0)),
            "audit_seconds": audit_record.get("seconds", 0),
            "precision_seconds_nested": audit_record.get("precision", {}).get("seconds", 0),
            "process_peak_rss_bytes": max(resources.get("search", {}).get("process_peak_rss_bytes", 0), audit_record.get("process_peak_rss_bytes", 0))}


def compare(args):
    left, right = args.left.resolve(), args.right.resolve()
    left_audits = (args.left_audits or left/"audits").resolve()
    right_audits = (args.right_audits or right/"audits").resolve()
    summarize(left, left_audits); summarize(right, right_audits)
    differences, count = [], 0
    for label, aroot, broot in (("search", left/"datasets", right/"datasets"), ("audit", left_audits, right_audits)):
        # Audit roots additionally contain source/provenance/summary manifests;
        # the scientific dataset records are the objects being paired here.
        paths = lambda root: {str(p.relative_to(root)) for p in root.rglob("*.json")
                              if label == "search" or p.relative_to(root).parts[0] in certification.coverage.DATASETS or
                              p.relative_to(root).parts[0] in COMPOSITES}
        a, b = paths(aroot), paths(broot); count += len(a)
        differences.extend(label+"/"+name for name in sorted(a ^ b))
        differences.extend(label+"/"+name for name in sorted(a & b) if scientific(read(aroot/name)) != scientific(read(broot/name)))
    if read(left/"input-hashes.json") != read(right/"input-hashes.json"): differences.append("input-hashes.json")
    result = {"passed": not differences, "scientific_records": count, "differences": differences,
              "excluded": ["elapsed seconds", "process peak RSS", "execution paths and logs"]}
    write(args.output, result); print(json.dumps(result), flush=True)
    require(result["passed"], "Independent scientific records differ.")


def rerun_component(args):
    root, output = args.run.resolve(), args.output.resolve()
    data = root/"inputs"/args.dataset; target = root/"datasets"/args.dataset/"cases"/args.case
    require(inventory(root/"inputs") == read(root/"input-hashes.json"), "Changed frozen inputs.")
    before = provenance(args.executable.resolve()); compatible_search_kernels(root, before)
    census = read(root/"datasets"/args.dataset/"census.json")
    found = [(k, c) for k, c in enumerate(census["components"]) if c["id"] == args.component]
    require(len(found) == 1, "Use the exact stable component ID from census.json.")
    k, component = found[0]
    require(not output.exists(), "Use a fresh rerun directory."); output.parent.mkdir(parents=True, exist_ok=True)
    audit_root = (args.audits or root/"audits").resolve()
    original_audit = audit_root/args.dataset/args.case/"components"/str(k)
    bundle = output.with_suffix(".context.json")
    write(bundle, {"search_context": read(target/"context.json"),
                   "audit_context": read(audit_root/args.dataset/args.case/"assembled/context.json")})
    execute([args.executable.resolve(), "joint-abc-rerun-component", data, args.case, component["id"], bundle, output], output.with_suffix(".log"))
    require(before == provenance(args.executable.resolve()), "Rerun sources or executable changed.")
    write(output/"provenance.json", before); shutil.copyfile(bundle, output/"context-bundle.json")
    differences = []
    for a, b in ((output/"fit.json", target/"components"/(str(k)+"-fit.json")),
                 (output/"audit/audit.json", original_audit/"audit.json"),
                 (output/"audit/fit.json", original_audit/"fit.json"),
                 (output/"audit/context.json", original_audit/"context.json")):
        if scientific(read(a)) != scientific(read(b)): differences.append(a.name)
    if (original_audit/"boundary.json").exists() and scientific(read(original_audit/"boundary.json")) != scientific(read(output/"audit/boundary.json")):
        differences.append("boundary.json")
    result = {"passed": not differences, "dataset": args.dataset, "case": args.case,
              "component_id": component["id"], "differences": differences}
    write(output/"comparison.json", result); print(json.dumps(result), flush=True)
    require(result["passed"], "Isolated component rerun differs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ("regression", "same-state"):
        sub = commands.add_parser(command)
        for name in ("baseline", "executable", "output"): sub.add_argument("--"+name, type=Path, required=True)
    sub = commands.add_parser("run")
    for name in ("baseline", "executable", "output", "regression-gate"): sub.add_argument("--"+name, type=Path, required=True)
    sub = commands.add_parser("audit")
    for name in ("run", "executable", "output"): sub.add_argument("--"+name, type=Path, required=True)
    sub = commands.add_parser("summarize"); sub.add_argument("--run", type=Path, required=True)
    sub.add_argument("--audits", type=Path)
    sub = commands.add_parser("compare")
    for name in ("left", "right", "output"): sub.add_argument("--"+name, type=Path, required=True)
    for name in ("left-audits", "right-audits"): sub.add_argument("--"+name, type=Path)
    sub = commands.add_parser("rerun-component")
    for name in ("run", "executable", "output"): sub.add_argument("--"+name, type=Path, required=True)
    for name in ("dataset", "case", "component"): sub.add_argument("--"+name, required=True)
    sub.add_argument("--audits", type=Path)
    args = parser.parse_args()
    {"regression": regression, "same-state": same_state, "run": run, "audit": audit,
     "summarize": lambda a: summarize(a.run, a.audits), "compare": compare, "rerun-component": rerun_component}[args.command](args)


if __name__ == "__main__":
    main()
