#!/usr/bin/env python3
"""Paired legacy/guarded certification, with frozen observations and honest failures."""
from __future__ import annotations

import argparse
from pathlib import Path
import os
import math
import resource
import subprocess
import sys
import time

import joint_abc_coverage as coverage
import simulation_contract
import numpy as np
from observation_matching import require, read

fixed = coverage.fixed
write = coverage.write
ROOT = coverage.ROOT
VARIANTS = ("legacy", "guarded", "guarded-log")
BASELINE = ROOT/"docs/developer/figures/joint-abc-coverage"


def paths_and_manifest(directory):
    paths = {k: Path(v) for k, v in read(directory/"inputs.json").items()}
    manifest, contract = simulation_contract.resolve(paths["manifest"], paths)
    fixture = read(coverage.FIXTURE)
    require(all(fixed.fold.sha256_file(paths[k]) == fixture["input_hashes"][k] for k in ("model", "map")), "Changed certification main inputs.")
    require(contract["widths"] == fixture["width_contract"]["b"], "Wrong main width vector.")
    require(manifest["settings"] == fixture["generation_contract"]["settings"], "Changed generation settings.")
    return paths, manifest, contract, fixture


def snapshot(directory, points, atoms):
    value = read(directory/"snapshot.json")
    require(value["support_policy"] == "sphere-fma-v1" and value["rows"] == len(points) and value["atoms"] == len(atoms), "Changed snapshot identity.")
    for key, name in (("voxels_sha256", "voxels.csv"), ("contributors_sha256", "contributors.csv")):
        require(fixed.fold.sha256_file(directory/name) == value[key], "Snapshot hash mismatch.")
    table = fixed.load_table(directory/"contributors.csv")
    row, atom = table["row"].astype(int), table["atom"].astype(int)
    require(np.array_equal(row, table["row"]) and np.array_equal(atom, table["atom"]), "Noninteger support indices.")
    require(len(table) == value["memberships"] and np.all((row >= 0) & (row < len(points))) and
            np.all((atom >= 0) & (atom < len(atoms))), "Invalid support population.")
    pairs = np.column_stack((row, atom))
    require(np.array_equal(np.lexsort((atom, row)), np.arange(len(table))) and len(np.unique(pairs, axis=0)) == len(table), "Duplicate/unordered support.")
    offsets = np.r_[0, np.cumsum(np.bincount(row, minlength=len(points)))]
    require(offsets.tolist() == value["row_offsets"], "Changed contributor CSR offsets.")
    for k, a in enumerate(atoms):
        square = coverage.squared_distances(points-a["position"]); expected = np.flatnonzero(square <= 6.25)
        selected = atom == k
        require(np.array_equal(row[selected], expected) and np.allclose(table["square"][selected], square[expected], rtol=0, atol=2e-15), "Support membership/distance differs.")
    return value


def legacy_science(value):
    ignored = {"trust", "lm", "variant", "initial_accepted", "search_reference_evaluations", "search_reference_seconds", "search_stopped_without_convergence", "observation_snapshot_sha256"}
    if isinstance(value, dict):
        return {k: legacy_science(v) for k, v in value.items() if k not in ignored and not k.endswith("seconds") and k != "process_peak_rss_bytes"}
    if isinstance(value, list): return [legacy_science(v) for v in value]
    return value


def certificate(fit, audit, replay):
    raw = fit.get("qualification_checks", {})
    checks = {k: raw.get(k) for k in ("inner", "b_gradient", "local_correction", "identified")}
    local_source = "legacy-double-QR"
    precision = audit.get("precision")
    if precision and precision.get("agreement_passed") and precision.get("fixed_face_feasible"):
        require(precision["local_correction_passed"] == (precision["local_correction_inf"] <= 1e-10), "Incorrect precision local correction flag.")
        checks["local_correction"] = precision["local_correction_passed"]; local_source = "independent-100-digit-QR"
    checks["derivative"] = None if audit["derivative_status"] == "unavailable" else audit["derivative_verified"]
    checks["endpoint_replay"] = replay
    checks["endpoint_trust"] = audit.get("endpoint_trust", {}).get("passed")
    checks["search_completed"] = fit.get("execution_complete", False) and not fit.get("search_stopped_without_convergence", False)
    if precision:
        from decimal import Decimal
        checks["precision_agreement"] = precision["agreement_passed"]
        checks["precision_feasibility"] = precision.get("fixed_face_feasible")
        checks["precision_kkt"] = Decimal(precision["projected_kkt"]) <= Decimal("1e-10") if "projected_kkt" in precision else None
        checks["precision_b_gradient"] = max(map(abs, precision["b_gradient"])) <= 1e-12 if "b_gradient" in precision else None
    # Strict complementarity is required for a regular boundary certificate.
    constraints = audit.get("constraint_evidence", [])
    weak_boundary = any(c["exactly_active"] and (c["dual"] is None or c["dual"] <= 1e-10) for c in constraints)
    checks["regular_active_face"] = not weak_boundary if constraints else None
    failures = [k for k, v in checks.items() if v is False]
    unavailable = [k for k, v in checks.items() if v is None]
    regular = all(v is True for v in checks.values())
    return {"schema_version": 2, "checks": checks, "failures": failures, "unavailable": unavailable,
            "regular_qualified": regular, "derivative_status": audit["derivative_status"],
            "local_correction_source": local_source,
            "identifiability_scope": "endpoint local numerical rank; generating-model evidence is separate",
            "truth_used": False, "boundary_status": "weakly-active-nonregular" if weak_boundary else "interior-or-strict-face", "evidence": evidence(fit, audit),
            "unavailable_reasons": {k: audit.get("reason", "not executed") for k in unavailable}}


def evidence(fit, audit):
    entries = []
    def add(name, value, threshold, source, comparison="<="):
        entries.append({"check": name, "value": value, "threshold": threshold,
                        "comparison": comparison, "source": source,
                        "status": "unavailable" if value is None else "pass" if
                        (value == threshold if comparison == "==" else value <= threshold) else "fail"})
    for label in ("primary", "reference"):
        endpoint = fit.get(label, {})
        for name, threshold, comparison in (("feasible", True, "=="), ("projected_kkt", 1e-10, "<="), ("b_gradient_inf", 1e-12, "<=")):
            add(label+"_"+name, endpoint.get(name), threshold, "fit."+label+"."+name, comparison)
    add("coefficient_agreement", fit.get("scaled_reference_difference"), 1e-10, "fit.scaled_reference_difference")
    for name, expected in (("design_spectrum", len(fit.get("primary", {}).get("beta", []))),
                           ("width_spectrum", len(fit.get("primary", {}).get("b", []))),
                           ("profile_jacobian_spectrum", len(fit.get("primary", {}).get("b", [])))):
        add(name+"_rank", fit.get(name, {}).get("rank"), expected, "fit."+name+".rank", "==")
    add("legacy_local_correction", fit.get("local_correction_inf"), 1e-10, "fit.local_correction_inf")
    if "precision" in audit:
        add("precision_local_correction", audit["precision"].get("local_correction_inf"), 1e-10, "audit.precision.local_correction_inf")
        for name in ("valid50", "valid100", "agreement_passed", "fixed_face_feasible", "derivative_passed"):
            add("fixed_face_reference_"+name, audit["precision"].get(name), True, "audit.precision."+name, "==")
        p = audit["precision"]
        add("precision_projected_kkt", float(p["projected_kkt"]) if "projected_kkt" in p else None, 1e-10, "audit.precision.projected_kkt")
        add("precision_b_gradient_inf", max(map(abs, p["b_gradient"])) if "b_gradient" in p else None, 1e-12, "audit.precision.b_gradient")
    for ladder in audit.get("ladders", []):
        for key, threshold in (("estimated_relative_error", 1e-7), ("first_relative_error", 1e-6), ("second_relative_error", 1e-6)):
            add("direction_"+str(ladder["direction"])+"_"+key, ladder[key], threshold, "audit.ladders."+str(ladder["direction"])+"."+key)
    return entries


def snapshot_replay(table, atom_count, y, endpoint):
    """Independent scalar kernels on immutable CSR memberships and distances."""
    beta, widths = np.asarray(endpoint["beta"]), np.asarray(endpoint["b"])
    basis, derivative = [], []
    for a, b in enumerate(widths):
        cells = table[table["atom"] == a]; ids = cells["row"].astype(int); square = cells["square"]
        radius = np.sqrt(square); exponent = np.exp(-square/(2*b*b))
        gaussian = (2*math.pi*b*b)**(-1.5)*exponent; center = math.sqrt(2/math.pi)/b
        charge = np.fromiter((center if r < 1e-5 else math.erf(r/b/math.sqrt(2))/r for r in radius), float)
        basis.extend(((ids, gaussian), (ids, charge)))
        derivative.append((ids, beta[2*a]*gaussian*(square/(b*b)-3)-beta[2*a+1]*center*np.where(square < 1e-10, 1., exponent)))
    require(len(widths) == atom_count, "Wrong endpoint atom count.")
    prediction = fixed.predict(basis, beta, len(y)); residual = prediction-y
    scale = max(1., np.linalg.norm(y)); norms = np.array([np.linalg.norm(v) for _, v in basis])
    require(np.all(norms > 0), "Zero design column.")
    u = norms*beta/scale; gradient = np.array([v@residual[ids] for ids, v in basis])/norms/scale
    projected = u-gradient; projected[::2] = np.maximum(0, projected[::2])
    return {"prediction": prediction, "residual": residual, "projected_kkt": float(np.max(np.abs(u-projected))),
            "b_gradient": np.array([v@residual[ids]/scale/scale for ids, v in derivative])}


def validate_audit(audit, fit):
    require(audit["schema_version"] == 2 and audit["case"] == fit["case"] and audit["dataset"] == fit["dataset"] and
            audit["legacy_joint_qualified"] == fit["joint_qualified"], "Wrong audit identity.")
    ladders = audit.get("ladders", [])
    require(not audit.get("precision") or len(ladders) == 3, "Precision target is missing its full derivative ladder.")
    if ladders:
        require(len(ladders) == 3, "Missing derivative direction.")
        for ladder in ladders:
            require(len(ladder["samples"]) == 17 and [s["h"] for s in ladder["samples"]] == [.01*2**-k for k in range(17)], "Changed step ladder.")
            candidates = [c for c in ladder["candidates"] if c["estimated_relative_error"] is not None]
            selected = min(candidates, key=lambda c: (c["estimated_relative_error"], c["index"])) if candidates else None
            require(ladder["selected"] == (selected["index"] if selected else None), "Step selection is not independent minimum estimated error.")
            passed = selected is not None and selected["estimated_relative_error"] <= 1e-7 and ladder["first_relative_error"] <= 1e-6 and ladder["second_relative_error"] <= 1e-6
            require(ladder["passed"] == passed, "Incorrect Richardson gate.")
        precision = audit.get("precision", {})
        passed = all(d["passed"] for d in ladders) and precision.get("agreement_passed", False) and precision.get("derivative_passed", False)
        require(audit["derivative_verified"] == passed, "False derivative certificate.")
    elif audit["derivative_verified"]:
        require(fit.get("derivative_verified") is True, "Missing derivative evidence.")
    if audit.get("precision", {}).get("valid50") and audit["precision"]["valid100"]:
        from decimal import Decimal
        p = audit["precision"]
        require(p["agreement_passed"] == (Decimal(p["maximum_scaled_precision_difference"]) <= Decimal("1e-20")), "False precision agreement.")
        require(p["derivative_passed"] == (p["fixed_face_feasible"] and max(p["derivative_relative_errors"]) <= 1e-6), "False high precision derivative flag.")


def validate_accepted(fit, table, atom_count, y):
    failures = []
    if fit.get("variant") == "legacy": return failures
    for trial in fit.get("trials", []):
        if not trial["accepted"]: continue
        require(trial.get("trust", {}).get("passed") is True, "Guard accepted an untrusted trial.")
        raw = snapshot_replay(table, atom_count, y, trial)
        if abs(raw["projected_kkt"]-trial["projected_kkt"]) > 1e-13: failures.append([trial["evaluation"], "KKT"])
        if not np.allclose(raw["b_gradient"], trial["b_gradient"], rtol=2e-9, atol=1e-13): failures.append([trial["evaluation"], "B-gradient"])
        if abs(np.linalg.norm(raw["residual"])-math.sqrt(trial["rss"])) > 512*np.finfo(float).eps*max(1., np.linalg.norm(y)):
            failures.append([trial["evaluation"], "residual"])
    return failures


def summarize(directory):
    directory = directory.resolve(); paths, manifest, contract, fixture = paths_and_manifest(directory)
    require(read(directory/"completion.json") == {"complete": True, "datasets": coverage.DATASETS}, "Incomplete certification datasets.")
    results = []; flat = []; estimates = []; parity = []; regressions = []; accepted_failures = []
    for name in coverage.DATASETS:
        root = directory/"datasets"/name
        points, atoms, truth, y64, y32 = coverage.validate_dataset(root, name, manifest, paths, fixture)
        domain = snapshot(root, points, atoms)
        contributor_table = fixed.load_table(root/"contributors.csv")
        for variant in VARIANTS:
            path = root/variant
            require(read(path/"snapshot.json") == domain and fixed.fold.sha256_file(path/"voxels.csv") == domain["voxels_sha256"], "Variant changed observations/domain.")
            result, rows = coverage.summarize_dataset(path, name, manifest, paths, fixture)
            result["variant"] = variant
            for row in rows: row["variant"] = variant
            estimates.extend(rows)
            for case in result["cases"]:
                label = case["case"]; fit = read(path/"fits"/f"{label}.json"); audit = read(path/"audits"/f"{label}.json")
                require(fit["observation_snapshot_sha256"] == fixed.fold.sha256_file(root/"snapshot.json"), "Modified observation snapshot.")
                validate_audit(audit, fit)
                y = y64 if label.endswith("double") else y32
                replay = {}
                for endpoint_name in ("primary", "reference"):
                    endpoint = fit.get(endpoint_name, {})
                    if not endpoint.get("valid"): replay[endpoint_name] = None; continue
                    raw = snapshot_replay(contributor_table, len(atoms), y, endpoint)
                    checks = {"kkt": abs(raw["projected_kkt"]-endpoint["projected_kkt"]) <= 1e-13,
                        "gradient": bool(np.allclose(raw["b_gradient"], endpoint["b_gradient"], rtol=2e-9, atol=1e-13)),
                        "residual": abs(np.linalg.norm(raw["residual"])-math.sqrt(endpoint["rss"])) <= 512*np.finfo(float).eps*max(1., np.linalg.norm(y))}
                    if endpoint_name == "primary":
                        residuals = fixed.load_table(path/"residuals"/f"{label}.csv")
                        checks["prediction"] = bool(np.allclose(raw["prediction"], residuals["prediction"], rtol=2e-13, atol=2e-12))
                    replay[endpoint_name] = {"checks": {k: bool(v) for k, v in checks.items()}, "passed": all(checks.values()),
                        "projected_kkt": raw["projected_kkt"], "b_gradient": raw["b_gradient"].tolist()}
                raw_pass = all(v is not None and v["passed"] for v in replay.values())
                write(path/"audits"/f"{label}-raw-replay.json", replay)
                cert = certificate(fit, audit, case["independent_endpoint_replay_passed"] and raw_pass)
                case["certificate"] = cert
                write(path/"audits"/f"{label}-certificate.json", cert)
                errors = validate_accepted(fit, contributor_table, len(atoms), y)
                if errors: accepted_failures.append({"dataset": name, "variant": variant, "case": label, "failures": errors})
                old = read(BASELINE/"datasets"/name/"fits"/f"{label}.json")
                if variant == "legacy":
                    same = legacy_science(old) == legacy_science(fit)
                    if not same: parity.append({"dataset": name, "case": label})
                elif old["joint_qualified"] and not (case["joint_qualified"] and
                        (case["oracle_recovered"] if label.endswith("double") else case["float32_accuracy_passed"])):
                    regressions.append({"dataset": name, "variant": variant, "case": label})
                flat.append({"dataset": name, "variant": variant, "case": label, "legacy_joint_qualified": case["joint_qualified"],
                    "regular_qualified": cert["regular_qualified"], "failures": ";".join(cert["failures"]), "unavailable": ";".join(cert["unavailable"]),
                    "derivative_status": cert["derivative_status"], "local_correction_source": cert["local_correction_source"],
                    "replay_passed": case["independent_endpoint_replay_passed"], "oracle_recovered": case["oracle_recovered"],
                    "float32_accuracy_passed": case["float32_accuracy_passed"], "rss": case.get("rss"), "stop_reason": case["stop_reason"]})
            result["regular_representatives"] = {}
            for precision in ("double", "float32"):
                eligible = [c for c in result["cases"] if c["case"].endswith(precision) and c["certificate"]["regular_qualified"]]
                result["regular_representatives"][precision] = min(eligible, key=lambda c: c["rss"])["case"] if eligible else None
            results.append(result)
    summary = {"experiment": "joint-abc-certification", "case_count": len(flat), "matrix_complete": len(flat) == 216,
        "variants": VARIANTS, "results": results, "manifest_contract": contract, "legacy_parity_failures": parity,
        "qualified_baseline_regressions": regressions, "untrusted_accepted_states": accepted_failures,
        "counts": {v: {"legacy_qualified": sum(c["legacy_joint_qualified"] for c in flat if c["variant"] == v),
                       "regular_qualified": sum(c["regular_qualified"] for c in flat if c["variant"] == v)} for v in VARIANTS}}
    write(directory/"results.json", summary); fixed.write_csv(directory/"failure-matrix.csv", flat)
    fixed.write_csv(directory/"estimates.csv", estimates)
    write(directory/"artifact-index.json", {str(p.relative_to(directory)): fixed.fold.sha256_file(p) for p in sorted(directory.rglob('*'))
         if p.is_file() and p.name not in ("artifact-index.json", "reproducibility.json")})
    print({k: summary[k] for k in ("case_count", "counts", "legacy_parity_failures", "qualified_baseline_regressions", "untrusted_accepted_states")}, flush=True)
    return summary


def execute(command, directory, phase):
    start = time.perf_counter()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    with (directory/f"{phase}.log").open('w') as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    write(directory/f"{phase}-execution.json", {"returncode": result.returncode, "seconds": time.perf_counter()-start,
        "command": command, "process_peak_rss_bytes": peak if sys.platform == 'darwin' else peak*1024})
    require(result.returncode == 0, f"{phase} failed; retained artifacts and log.")


def run(args):
    output = args.output.resolve(); require(not output.exists(), "Use a fresh output directory.")
    paths = {k: getattr(args, k).resolve() for k in ("model", "map", "manifest")}
    simulation_contract.resolve(paths["manifest"], paths)
    executable = args.executable.resolve(); before = fixed.experiment.provenance(executable, ROOT)
    output.mkdir(parents=True)
    write(output/"inputs.json", {k: str(v) for k, v in paths.items()}); paths_and_manifest(output)
    write(output/"input-hashes.json", {k: fixed.fold.sha256_file(v) for k, v in paths.items()})
    write(output/"provenance.json", before); write(output/"fixture.json", read(coverage.FIXTURE))
    write(output/"environment.json", {"python": sys.version, "jobs": 1})
    execute([str(executable), "joint-abc-certification", *(str(paths[k]) for k in ("model", "map", "manifest")), str(output)], output, "search")
    require(before == fixed.experiment.provenance(executable, ROOT), "Sources changed during search.")
    audit(output, executable)
    require(before == fixed.experiment.provenance(executable, ROOT), "Sources changed during audit.")
    write(output/"execution-status.json", {"complete": True, "inputs_and_sources_stable": True})
    require(read(output/"input-hashes.json") == {k: fixed.fold.sha256_file(v) for k, v in paths.items()}, "Inputs changed during run.")
    summarize(output)


def audit(directory, executable):
    paths_and_manifest(directory)
    before = fixed.experiment.provenance(executable.resolve(), ROOT)
    execute([str(executable.resolve()), "joint-abc-certification-audit", str(directory.resolve())], directory, "audit")
    require(before == fixed.experiment.provenance(executable.resolve(), ROOT), "Sources changed during audit.")
    paths_and_manifest(directory)
    write(directory/"audit-provenance.json", before)


def compare(left, right, output):
    for p in (left, right): summarize(p)
    excluded = {"inputs.json", "environment.json", "artifact-index.json", "reproducibility.json", "search-execution.json", "audit-execution.json"}
    def inventory(p): return {str(f.relative_to(p)) for f in p.rglob('*') if f.suffix in ('.json', '.csv') and f.name not in excluded}
    a, b = inventory(left), inventory(right); differences = sorted(a ^ b)
    for name in sorted(a & b):
        same = coverage.joint.scientific(read(left/name)) == coverage.joint.scientific(read(right/name)) if name.endswith('.json') else fixed.fold.sha256_file(left/name) == fixed.fold.sha256_file(right/name)
        if not same: differences.append(name)
    write(output, {"passed": not differences, "case_count": 216, "differences": differences,
                   "comparison": "exact scientific JSON and CSV bytes; timing, RSS, execution paths excluded"})
    require(not differences, "Scientific repeat differs.")


def plots(directory, output):
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    rows = read(directory/"results.json")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    for i, key in enumerate(("legacy_qualified", "regular_qualified")):
        axes[i].bar(VARIANTS, [rows["counts"][v][key] for v in VARIANTS]); axes[i].set(ylim=(0,72), title=key, ylabel="Branches / 72")
    for suffix in ('png', 'pdf'): fig.savefig(output/f"qualification.{suffix}", dpi=160)
    plt.close(fig)
    data = []
    for name in ("weak-1e-4", "near-0.02"):
        p = directory/'datasets'/name/'legacy/audits/first-stage-double.json'
        if p.exists(): data.append((name, read(p)))
    fig, axes = plt.subplots(1, len(data), figsize=(12, 4), squeeze=False, layout="constrained")
    for axis, (name, value) in zip(axes[0], data):
        for d in value.get("ladders", []):
            axis.loglog([s['h'] for s in d['samples']], [s['relative_error'] for s in d['samples']], marker='.', label=f"direction {d['direction']}")
        axis.axhline(1e-6, color='black', linestyle='--'); axis.set(title=name, xlabel='h', ylabel='Relative derivative discrepancy'); axis.legend()
    for suffix in ('png', 'pdf'): fig.savefig(output/f"derivative-ladders.{suffix}", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
    for axis, precision in zip(axes, ("double", "float32")):
        for variant in VARIANTS:
            fit = read(directory/"datasets/near-0.02"/variant/"fits"/f"narrower-{precision}.json")
            trials = fit["trials"]
            axis.semilogy([t["evaluation"] for t in trials], [t["b"][0] for t in trials], '.-', label=variant)
            accepted = [t for t in trials if t["accepted"]]
            axis.scatter([t["evaluation"] for t in accepted], [t["b"][0] for t in accepted], marker='o', facecolors='none', edgecolors='black')
        axis.set(title="near-0.02/narrower/"+precision, xlabel="Profile trial", ylabel="First atom B (angstrom)"); axis.legend()
    for suffix in ('png', 'pdf'): fig.savefig(output/f"pathological-trials.{suffix}", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
    for axis, location in zip(axes, ("truth-boundary-double.json", "legacy/audits/first-stage-double-boundary.json")):
        value = read(directory/"datasets/active-a"/location)
        for atom in (1,5,9):
            rows = [r for r in value["precision100"]["rows"] if r["atom"] == atom]
            axis.loglog([float(r['step']) for r in rows], [float(r['structure_loss']) for r in rows], label=f"atom {atom+1}")
        axis.set(title="Generating model" if location.startswith("truth") else "Actual endpoint", xlabel="Positive log-B step", ylabel="Prediction change squared / 2"); axis.legend()
    for suffix in ('png', 'pdf'): fig.savefig(output/f"boundary-compensation.{suffix}", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__); commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('run')
    for name in ('executable', 'model', 'map', 'manifest', 'output'): p.add_argument('--'+name, required=True, type=Path)
    p = commands.add_parser('audit'); p.add_argument('directory', type=Path); p.add_argument('--executable', type=Path, required=True)
    p = commands.add_parser('summarize'); p.add_argument('directory', type=Path)
    p = commands.add_parser('compare')
    for name in ('left', 'right', 'output'): p.add_argument('--'+name, type=Path, required=True)
    p = commands.add_parser('plots'); p.add_argument('directory', type=Path); p.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'run': run(args)
    elif args.command == 'audit': audit(args.directory, args.executable)
    elif args.command == 'summarize': summarize(args.directory)
    elif args.command == 'compare': compare(args.left, args.right, args.output)
    else: plots(args.directory, args.output)


if __name__ == '__main__': main()
