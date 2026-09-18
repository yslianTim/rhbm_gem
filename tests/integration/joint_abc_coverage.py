#!/usr/bin/env python3
"""Frozen joint-LS coverage matrix; failures are evidence, never relaxed gates."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import math
import os
from pathlib import Path
import platform
import resource
import signal
import struct
import subprocess
import sys
import time

import joint_abc_profile as joint  # Sets single-thread environment before NumPy.
import numpy as np
from observation_matching import read, require

fixed = joint.fixed
write = fixed.experiment.write
ROOT = fixed.ROOT
FIXTURE = ROOT / "tests/benchmarks/joint_abc_coverage.json"
DATASETS = ["heterogeneous-168", "baseline", "weak-1e-2", "weak-1e-4", "active-a",
            "zero-signal", "near-0.10", "near-0.02", "duplicate"]
STARTS = ["first-stage", "narrower", "wider", "mixed"]
CASES = [f"{s}-{p}" for p in ("double", "float32") for s in STARTS]


def element_width(element):
    return .4 if element == 8 else .45 if element == 7 else .5


def grid_axes(dataset):
    # Match the recorded Release generator's fused coordinate arithmetic exactly.
    return [np.array([math.fma(float(i), spacing, origin) for i in range(n)]) for n, spacing, origin in
            zip(dataset["grid_size"], dataset["generation_spacing"], dataset["generation_origin"])]


def squared_distances(delta):
    square = delta[:, 0]**2+delta[:, 1]**2+delta[:, 2]**2
    # Only the hard support decision needs bit-exact arithmetic. No radius slack.
    boundary = np.flatnonzero(np.abs(square-6.25) <= 64*np.finfo(float).eps*6.25)
    for i in boundary:
        x, y, z = delta[i]; square[i] = math.fma(z, z, math.fma(y, y, x*x))
    return square


def sphere_membership(dataset):
    size = np.array(dataset["grid_size"]); axes = grid_axes(dataset)
    origin = np.array(dataset["generation_origin"]); spacing = np.array(dataset["generation_spacing"])
    coverage = np.zeros(int(np.prod(size)), dtype=np.uint16); nearest = np.full(len(coverage), np.inf)
    for atom in dataset["atoms"]:
        center = np.asarray(atom["position"])
        lo = np.maximum(np.floor((center-2.5-origin)/spacing).astype(int), 0)
        hi = np.minimum(np.floor((center+2.5-origin)/spacing).astype(int), size-1)
        z, y, x = np.meshgrid(*(np.arange(lo[k], hi[k]+1) for k in (2, 1, 0)), indexing="ij")
        ids = (x+size[0]*(y+size[1]*z)).ravel()
        points = np.column_stack((axes[0][x.ravel()], axes[1][y.ravel()], axes[2][z.ravel()]))
        square = squared_distances(points-center); keep = square <= 6.25; selected = ids[keep]
        coverage[selected] += 1; nearest[selected] = np.minimum(nearest[selected], square[keep])
    ids = np.flatnonzero(coverage)
    return ids, coverage[ids], np.sqrt(nearest[ids])


def columns(points, atoms, widths):
    out = []
    for atom, b in zip(atoms, widths):
        square = squared_distances(points-atom["position"]); ids = np.flatnonzero(square <= 6.25)
        r = np.sqrt(square[ids]); g = (2*math.pi*b*b)**(-1.5)*np.exp(-square[ids]/(2*b*b))
        k = np.fromiter((math.sqrt(2/math.pi)/b if v < 1e-5 else math.erf(v/b/math.sqrt(2))/v for v in r), float, count=len(r))
        out.extend(((ids, g), (ids, k)))
    return out


def raw_certificate(points, atoms, y, endpoint):
    return joint.raw_certificate(points, atoms, y, endpoint, column_builder=columns)


def numerical_check(qualified, failures, passed, message):
    if not passed:
        # Keep failed-endpoint evidence without weakening any tolerance or gate.
        require(not qualified, message)
        failures.append(message)


def validate_inputs(paths, fixture):
    require(fixture == read(FIXTURE), "Changed frozen coverage fixture.")
    hashes = fixed.fold.validate_input_hashes(paths, fixture["input_hashes"])
    manifest = fixed.fold.load_simulation_manifest(paths["manifest"], hashes)
    require(all(manifest[k] == v for k, v in fixture["generation_contract"].items()), "Changed generation contract.")
    require([a["preparation_index"] for a in manifest["atoms"]] == list(range(168)), "Changed preparation order.")
    require(len({a["serial_id"] for a in manifest["atoms"]}) == 168, "Duplicate identity.")
    require(fixture["width_contract"]["b"] == [element_width(a["element"]) for a in manifest["atoms"]], "Wrong width rule.")
    return manifest


def synthetic_truth(name):
    require(name in DATASETS[1:], "Unknown synthetic dataset.")
    atoms = []; truth = []
    for z in range(2):
        for y in range(2):
            for x in range(3):
                serial = 1+x+3*y+6*z; element = 6+(serial-1) % 3
                atoms.append({"serial_id": serial, "chain_id": "T", "sequence_id": serial, "component_id": "ALA",
                              "atom_id": {6: "C", 7: "N", 8: "O"}[element], "alternate_indicator": ".",
                              "position": [1.2*x, 1.2*y, 1.2*z], "element": element})
                truth.append([float(element), element_width(element), .2 if serial % 2 else -.2])
    truth = np.asarray(truth)
    for k in (1, 5, 9):
        if name in ("weak-1e-2", "weak-1e-4"):
            truth[k, [0, 2]] *= .01 if name == "weak-1e-2" else .0001
        if name == "active-a": truth[k, [0, 2]] = [0, .2]
    if name == "zero-signal": truth[5, [0, 2]] = 0
    if name in ("near-0.10", "near-0.02", "duplicate"):
        atoms[3]["position"] = [.1 if name == "near-0.10" else .02 if name == "near-0.02" else 0., 0., 0.]
    return atoms, truth


def initial_widths(b0, atoms, start):
    b0 = np.asarray(b0, dtype=float)
    require(b0.shape == (len(atoms),) and np.isfinite(b0).all() and np.all(b0 > 0), "Invalid first-stage B.")
    require(start in STARTS, "Unknown coverage start.")
    return joint.initial_widths(b0, atoms, "checkpoint" if start == "first-stage" else start)


def validate_initialization(value, atoms):
    require(value["source"] == "first-stage-float32-map" and value["sampling_method"] == "FibonacciDeterministic" and
            value["jobs"] == 1 and value["seed_abc"] == [0, 1, 0] and
            value["uses_peeling"] is False and value["uses_truth"] is False, "Changed initialization contract.")
    b0 = np.asarray(value["b0"], dtype=float)
    valid = b0.shape == (len(atoms),) and np.isfinite(b0).all() and bool(np.all(b0 > 0))
    require(value["valid"] == valid, "Incorrect initialization validity.")
    if not value["atoms"]:
        require(not valid and value["reason"] not in ("valid-widths", "invalid-widths"), "Missing initialization records.")
        return None
    require([a["identity"] for a in value["atoms"]] == atoms, "Initialization identity mismatch.")
    widths = np.asarray([a["mdpde"][1] for a in value["atoms"]], dtype=float)
    require(np.array_equal(widths, b0, equal_nan=True), "Initialization B does not come from first-stage MDPDE.")
    warnings = 0
    for atom in value["atoms"]:
        require(atom["ols"][2] == atom["mdpde"][2] == 0 and atom["raw_sample_count"] >= atom["sample_count"],
                "Changed fixed-offset first-stage contract.")
        require(atom["alpha"] is not None and 0 <= atom["alpha"] <= 1, "Invalid trained alpha.")
        diagnostics = atom["diagnostics"]
        warnings += diagnostics is None or diagnostics["qualification"] == 0
    require(value["unqualified_local_count"] == warnings, "Incorrect initialization warning count.")
    return b0 if valid else None


def validate_dataset(directory, name, manifest, paths, fixture):
    dataset = read(directory/"dataset.json")
    if name == DATASETS[0]:
        keys = fixed.IDENTITY + ["element"]
        atoms = [{k: a[k] for k in keys} for a in manifest["atoms"]]
        truth = np.asarray([[a["element"], element_width(a["element"]), a["charge_used"]] for a in manifest["atoms"]])
        geometry = {"grid_size": manifest["settings"]["grid_size"], "generation_origin": manifest["settings"]["origin"],
                    "generation_spacing": manifest["settings"]["grid_spacing"]}
    else:
        atoms, truth = synthetic_truth(name)
        geometry = {"grid_size": [33, 29, 29], "generation_origin": [-3.6]*3, "generation_spacing": [.3]*3}
    require(dataset["experiment"] == "joint-abc-coverage" and dataset["name"] == name and
            dataset["atoms"] == atoms and dataset["radius"] == 2.5 and dataset["membership_geometry"] == "generation" and
            all(dataset[k] == v for k, v in geometry.items()), "Changed dataset geometry/identity.")
    require(read(directory/"scoring-truth.json") == {"beta": truth[:, [0, 2]].ravel().tolist(), "b": truth[:, 1].tolist(),
            "known_unidentified": name in ("zero-signal", "duplicate")}, "Wrong width rule or scoring truth.")
    ids, multiplicity, nearest = sphere_membership(dataset)
    require(len(ids) == dataset["row_count"] and int(multiplicity.sum()) == dataset["memberships"], "Wrong ROI population.")
    if name == DATASETS[0]:
        require(len(ids) == fixture["expected_rows"] and int(multiplicity.sum()) == fixture["expected_memberships"] and
                hashlib.sha256(ids.astype("<u8").tobytes()).hexdigest() == fixture["voxel_ids_sha256_le_u64"], "Wrong main ROI.")
    require(fixed.fold.sha256_file(directory/"voxels.csv") == dataset["voxel_table_sha256"], "Voxel hash mismatch.")
    table = fixed.load_table(directory/"voxels.csv")
    require(len(table) == len(ids) and np.array_equal(table["row"], np.arange(len(ids))) and
            np.array_equal(table["index"], ids) and np.array_equal(table["multiplicity"], multiplicity) and
            np.allclose(table["nearest_distance"], nearest, rtol=0, atol=2e-13), "Wrong voxel rows.")
    nx, ny, _ = dataset["grid_size"]
    axes = grid_axes(dataset)
    points = np.column_stack((axes[0][ids % nx], axes[1][ids//nx % ny], axes[2][ids//(nx*ny)]))
    require(np.array_equal(np.column_stack([table[k] for k in ("x", "y", "z")]), points), "Wrong voxel coordinates/arithmetic contract.")
    y64, y32, q = (table[k] for k in ("reference_double", "observed", "quantization_delta"))
    fixed.validate_quantization(y64, y32, q)
    basis = columns(points, atoms, truth[:, 1]); beta = truth[:, [0, 2]].ravel()
    absolute = fixed.predict([(ids, np.abs(v)) for ids, v in basis], np.abs(beta), len(y64))
    require(np.all(np.abs(fixed.predict(basis, beta, len(y64))-y64) <= 512*np.finfo(float).eps*np.maximum(1, absolute)),
            "Independent analytic forward mismatch.")
    if name == DATASETS[0]:
        with paths["map"].open("rb") as stream:
            header = stream.read(1024)
            require(struct.unpack_from("<3i", header) == tuple(dataset["grid_size"]) and
                    struct.unpack_from("<i", header, 12)[0] == 2 and struct.unpack_from("<3i", header, 64) == (1, 2, 3), "Invalid map header.")
            stream.seek(1024+struct.unpack_from("<i", header, 92)[0]); stored = np.fromfile(stream, dtype="<f4")
        require(len(stored) == math.prod(dataset["grid_size"]) and np.array_equal(stored[ids].astype(float), y32), "Map bytes differ.")
        outside = np.ones(len(stored), dtype=bool); outside[ids] = False
        require(np.all(stored[outside] == 0), "Nonzero contributions outside fixed support.")
    forward = read(directory/"forward-status.json")
    require(forward["passed"] and forward["generator_jobs"] == 1 and forward["maximum_quantized_difference"] == 0 and
            forward["checked_map_voxels"] == math.prod(dataset["grid_size"]), "Missing full-map forward evidence.")
    return points, atoms, truth, y64, y32


def truth_identifiability(points, atoms, truth, y):
    """Independent dense SVD for the small controls, including exact zero columns."""
    basis = columns(points, atoms, truth[:, 1]); n = len(points); m = len(atoms)
    x = np.zeros((n, 2*m)); d = np.zeros((n, m))
    for k, (ids, values) in enumerate(basis): x[ids, k] = values
    for k, (atom, b) in enumerate(zip(atoms, truth[:, 1])):
        ids, gaussian = basis[2*k]; r2 = np.sum((points[ids]-atom["position"])**2, axis=1)
        dk = -math.sqrt(2/math.pi)/b*np.where(r2 < 1e-10, 1, np.exp(-r2/(2*b*b)))
        d[ids, k] = truth[k, 0]*gaussian*(r2/b**2-3)+truth[k, 2]*dk
    def decompose(a):
        u, s, v = np.linalg.svd(a, full_matrices=False)
        threshold = np.finfo(float).eps*max(a.shape)*s[0]
        return u, s, v, int(np.count_nonzero(s > threshold)), float(threshold)
    _, sx, _, rx, tx = decompose(x/np.linalg.norm(x, axis=0))
    free = np.array([k for k in range(2*m) if k % 2 or truth[k//2, 0] > 0])
    z = x[:, free]/np.linalg.norm(x[:, free], axis=0)
    u, _, _, rank, _ = decompose(z)
    projected = (d-u[:, :rank]@(u[:, :rank].T@d))/max(1., np.linalg.norm(y))
    _, sw, vw, rw, tw = decompose(projected)
    norms = np.linalg.norm(projected, axis=0)
    normalized = projected/np.where(norms > 0, norms, 1)
    _, sn, _, rn, tn = decompose(normalized)
    return {"design_rank": rx, "design_singular_values": sx.tolist(), "design_rank_threshold": tx,
            "free_design_rank": rank, "width_rank": rw, "width_singular_values": sw.tolist(), "width_rank_threshold": tw,
            "width_column_norms": norms.tolist(), "normalized_width_singular_values": sn.tolist(),
            "normalized_width_rank": rn, "normalized_width_rank_threshold": tn,
            "weak_directions": vw[-3:].tolist(), "active_face_only": bool(np.any(truth[:, 0] == 0)),
            "location": "generating truth; diagnostic only, never optimizer initialization"}


def compare_starts(fits, estimates):
    pairs = []; representatives = {}; consistent = True
    for precision in ("double", "float32"):
        eligible = [c for c in CASES if c.endswith(precision) and fits[c]["joint_qualified"]]
        representatives[precision] = min(eligible, key=lambda c: fits[c]["primary"]["rss"]) if eligible else None
        consistent &= len(eligible) == 4
        # Retain endpoint disagreements even when a branch does not qualify.
        for left, right in itertools.combinations([c for c in CASES if c.endswith(precision) and c in estimates], 2):
            a, b = estimates[left], estimates[right]
            diff = fixed.scaled_difference(a, b); log_diff = float(np.max(np.abs(np.log(a[:, 1])-np.log(b[:, 1]))))
            passed = diff <= 1e-8 and log_diff <= 1e-8
            if left in eligible and right in eligible: consistent &= passed
            pairs.append({"left": left, "right": right, "scaled_abc_difference": diff, "maximum_log_b_difference": log_diff,
                          "both_qualified": left in eligible and right in eligible, "passed": passed})
    return bool(consistent), pairs, representatives


def validate_controls(directory, points, atoms, truth, y64, y32, b0, fits, estimates):
    rows = []
    for precision, y in (("double", y64), ("float32", y32)):
        for kind in ("true-b", "first-stage-b"):
            if kind == "first-stage-b" and b0 is None: continue
            control = read(directory/"controls"/f"{kind}-{precision}.json")
            expected_b = truth[:, 1] if kind == "true-b" else b0
            require(control["b"] == expected_b.tolist() and control["alpha"] == 0, "Wrong fixed-B control widths.")
            for label in ("primary", "reference"):
                endpoint = control[label]
                if not endpoint["valid"]: continue
                raw = raw_certificate(points, atoms, y, dict(endpoint, b=control["b"]))
                require(abs(raw["projected_kkt"]-endpoint["projected_kkt"]) <= 1e-13 and
                        raw["feasible"] == endpoint["feasible"] and
                        endpoint["kkt_passed"] == (endpoint["feasible"] and endpoint["projected_kkt"] <= 1e-10),
                        "Independent control KKT differs.")
                require(abs(np.linalg.norm(raw["residual"])-math.sqrt(endpoint["rss"])) <=
                        512*np.finfo(float).eps*max(1., np.linalg.norm(y)), "Independent control residual differs.")
            p, r = (control[k] for k in ("primary", "reference"))
            difference = fixed.scaled_difference(np.asarray(p["beta"]), np.asarray(r["beta"]))
            qualified = all(e["valid"] and e["kkt_passed"] for e in (p, r)) and difference <= 1e-10 and control["design_spectrum"]["rank"] == 2*len(atoms)
            require(control["qualified"] == qualified, "False fixed-B qualification.")
            control_abc = np.column_stack((np.asarray(p["beta"])[::2], control["b"], np.asarray(p["beta"])[1::2]))
            if kind == "first-stage-b":
                initial = fits[f"first-stage-{precision}"]["initial"]
                if initial["valid"] and p["valid"]:
                    require(fixed.scaled_difference(np.asarray(initial["beta"]), np.asarray(p["beta"])) <= 1e-10,
                            "Initial profile does not reproduce fixed-B0 control.")
            for case, abc in estimates.items():
                if not case.endswith(precision): continue
                rows.append({"case": case, "control": f"{kind}-{precision}", "control_qualified": qualified,
                             "joint_rss": fits[case]["primary"]["rss"], "control_rss": p["rss"],
                             **{f"control_{key}_rmse": fixed.statistics((control_abc-truth)[:, k])["rmse"] for k, key in enumerate(("A", "B", "C"))}})
    return rows


def summarize_dataset(directory, name, manifest, paths, fixture):
    completion = read(directory/"completion.json")
    require(completion == {"complete": True, "cases": CASES}, "Missing completed coverage case.")
    points, atoms, truth, y64, y32 = validate_dataset(directory, name, manifest, paths, fixture)
    initialization = read(directory/"initialization.json"); b0 = validate_initialization(initialization, atoms)
    known_unidentified = name in ("zero-signal", "duplicate")
    evidence = None
    if name != DATASETS[0]:
        evidence = truth_identifiability(points, atoms, truth, y64)
        if name == "zero-signal": require(evidence["width_rank"] < 12 and evidence["width_column_norms"][5] == 0, "Missing zero-width evidence.")
        if name == "duplicate": require(evidence["design_rank"] < 24, "Missing duplicate-rank evidence.")
        write(directory/"truth-identifiability.json", evidence)
    fits = {}; estimates = {}; summaries = []; atom_rows = []
    for case in CASES:
        fit = read(directory/"fits"/f"{case}.json"); fits[case] = fit
        require(fit["case"] == case and fit["dataset"] == name, "Wrong case identity.")
        start, precision = case.rsplit("-", 1); y = y64 if precision == "double" else y32
        values = None; numeric_failures = []
        if b0 is not None:
            require(fit["initialization_status"] == "valid", "Skipped valid initialization.")
            values = joint.validate_fit(directory, case, fit, points, atoms, y, initial_widths(b0, atoms, start),
                certificate_fn=raw_certificate,
                numeric_check=lambda passed, message: numerical_check(fit["joint_qualified"], numeric_failures, passed, message))
        else:
            require(fit["initialization_status"] == "invalid" and fit["execution_complete"] is False and
                    fit["joint_qualified"] is False and fit["qualification_failure"] == "initialization-failed",
                    "Invalid initialization was reported successful.")
        checks = fit.get("qualification_checks", {})
        row = {"dataset": name, "case": case, "initialization_valid": b0 is not None,
               "execution_complete": fit["execution_complete"], "joint_qualified": fit["joint_qualified"],
               "qualification_failure": fit["qualification_failure"], "stop_reason": fit.get("stop_reason"),
               "numerical_checks": {k: checks.get(k) for k in ("inner", "b_gradient", "local_correction", "derivative")},
               "endpoint_identified": checks.get("identified"), "known_unidentified": known_unidentified,
               "independent_endpoint_replay_passed": (not numeric_failures) if values is not None else None,
               "independent_replay_failures": numeric_failures,
               "oracle_recovered": None if known_unidentified else False,
               "float32_accuracy_passed": None if known_unidentified else False}
        if values is not None:
            estimates[case] = values; error = values-truth
            row.update({key: fixed.statistics(error[:, k]) for k, key in enumerate(("A", "B", "C"))})
            row.update({k: fit["primary"][k] for k in ("rss", "relative_residual", "residual_rmse", "projected_kkt", "b_gradient_inf")})
            row.update({k: fit.get(k) for k in ("local_correction_inf", "profile_evaluations", "accepted_updates")})
            if not known_unidentified:
                row["oracle_recovered"] = precision == "double" and fixed.scaled_difference(values, truth) <= 1e-10 and row["relative_residual"] <= 1e-12
                row["float32_accuracy_passed"] = precision == "float32" and bool(np.all(np.abs(error) < .01))
            for k, atom in enumerate(atoms):
                atom_rows.append({"dataset": name, "case": case, "serial_id": atom["serial_id"],
                                  **{key: values[k, q] for q, key in enumerate(("A", "B", "C"))},
                                  **{key+"_error": error[k, q] for q, key in enumerate(("A", "B", "C"))},
                                  "truth_comparison_identifiable": not known_unidentified})
        summaries.append(row)
    consistent, pairs, representatives = compare_starts(fits, estimates)
    controls = validate_controls(directory, points, atoms, truth, y64, y32, b0, fits, estimates) if name == DATASETS[0] else []
    passed = (all(s["joint_qualified"] for s in summaries) and consistent and
              all(s["oracle_recovered"] if s["case"].endswith("double") else s["float32_accuracy_passed"] for s in summaries)) if not known_unidentified else None
    results = {"dataset": name, "matrix_complete": True, "initialization_valid": b0 is not None,
               "known_unidentified": known_unidentified, "recovery_passed": passed,
               "nonidentifiability_evidence_verified": bool(known_unidentified and evidence is not None),
               "multistart_consistent": consistent, "representatives": representatives, "cases": summaries,
               "pairwise_comparisons": pairs, "fixed_b_comparisons": controls}
    if b0 is not None:
        initial_abc = np.asarray([a["mdpde"] for a in initialization["atoms"]], dtype=float)
        results["first_stage_errors"] = {key: fixed.statistics((initial_abc-truth)[:, k]) for k, key in enumerate(("A", "B", "C"))}
    write(directory/"results.json", results)
    if atom_rows: fixed.write_csv(directory/"estimates.csv", atom_rows)
    if pairs: fixed.write_csv(directory/"multistart.csv", pairs)
    if controls: fixed.write_csv(directory/"fixed-b-comparison.csv", controls)
    return results, atom_rows


def summarize(directory):
    directory = directory.resolve(); clock = time.perf_counter()
    fixture = read(directory/"fixture.json"); paths = {k: Path(v) for k, v in read(directory/"inputs.json").items()}
    manifest = validate_inputs(paths, fixture)
    require(all(fixed.fold.sha256_file(Path(p)) == h for p, h in read(directory/"input-hashes.json").items()), "Changed inputs.")
    require(read(directory/"completion.json") == {"complete": True, "datasets": DATASETS}, "Missing completed coverage dataset.")
    results = []; atoms = []
    for name in DATASETS:
        result, rows = summarize_dataset(directory/"datasets"/name, name, manifest, paths, fixture)
        results.append(result); atoms.extend(rows)
    summary = {"experiment": "joint-abc-coverage", "matrix_complete": True, "case_count": 72,
               "main_recovery_passed": results[0]["recovery_passed"], "datasets": results,
               "qualified_count": sum(c["joint_qualified"] for r in results for c in r["cases"]),
               "validation_passed": all(c["independent_endpoint_replay_passed"] is not False for r in results for c in r["cases"]),
               "qualified_endpoint_replay_passed": True, "data_contract_passed": True}
    write(directory/"results.json", summary)
    if atoms: fixed.write_csv(directory/"estimates.csv", atoms)
    flat = []
    for r in results:
        for row in r["cases"]:
            flat.append({**{k: row.get(k) for k in ("dataset", "case", "initialization_valid", "execution_complete", "joint_qualified",
                         "qualification_failure", "endpoint_identified", "known_unidentified", "oracle_recovered", "float32_accuracy_passed",
                         "independent_endpoint_replay_passed",
                         "rss", "residual_rmse", "projected_kkt", "b_gradient_inf", "local_correction_inf", "profile_evaluations", "accepted_updates")},
                         **{f"{p}_{k}": row.get(p, {}).get(k) for p in ("A", "B", "C") for k in ("bias", "rmse", "p99", "max", "below_0_01")}})
    fixed.write_csv(directory/"comparison.csv", flat)
    write(directory/"independent-validation.json", {"complete": True, "passed": summary["validation_passed"],
        "qualified_endpoints_passed": True, "data_contract_passed": True,
        "failed_endpoint_replays": [{"dataset": r["dataset"], "case": c["case"], "failures": c["independent_replay_failures"]}
                                    for r in results for c in r["cases"] if c["independent_endpoint_replay_passed"] is False],
        "seconds": time.perf_counter()-clock,
        "checks": ["frozen inputs, identities, heterogeneous widths", "complete 72-case matrix including failed branches",
                   "original-row forward, float32 map bytes, ROI and quantization", "first-stage provenance and legal widths",
                   "independent A/C KKT, B gradient, residual and qualification replay", "fixed-B0 initial profile parity",
                   "12-atom independent truth rank and width sensitivity", "truth-free RSS representative selection"]})
    write(directory/"raw-artifact-index.json", {str(p.relative_to(directory)): fixed.fold.sha256_file(p)
        for p in sorted(directory.rglob("*")) if p.is_file() and p.name not in ("raw-artifact-index.json", "reproducibility.json")})
    print({k: summary[k] for k in ("matrix_complete", "main_recovery_passed", "qualified_count", "validation_passed")})
    return summary


def run(args):
    output = args.output.resolve(); require(not output.exists(), "Use a fresh output directory.")
    paths = {k: getattr(args, k).resolve() for k in ("model", "map", "manifest")}; fixture = read(FIXTURE)
    validate_inputs(paths, fixture); executable = args.executable.resolve()
    before = fixed.experiment.provenance(executable, ROOT)
    hashes = {str(p): fixed.fold.sha256_file(p) for p in [*paths.values(), FIXTURE]}
    output.mkdir(parents=True)
    for name, value in (("inputs", {k: str(v) for k, v in paths.items()}), ("fixture", fixture), ("provenance", before), ("input-hashes", hashes),
                        ("environment", {"python": sys.version, "platform": platform.platform(), "jobs": 1,
                         "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()})):
        write(output/f"{name}.json", value)
    command = [str(executable), "joint-abc-coverage", *(str(paths[k]) for k in ("model", "map", "manifest")), str(output)]
    start = time.perf_counter(); interrupted = False
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with (output/"run.log").open("w") as log:
        child = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        try: status = child.wait()
        except KeyboardInterrupt:
            interrupted = True; child.send_signal(signal.SIGINT); status = child.wait()
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    write(output/"execution.json", {"command": command, "returncode": status, "seconds": time.perf_counter()-start,
          "jobs": 1, "user_interrupted": interrupted, "peak_child_rss_bytes": peak if sys.platform == "darwin" else peak*1024})
    stable = before == fixed.experiment.provenance(executable, ROOT) and all(fixed.fold.sha256_file(Path(p)) == h for p, h in hashes.items())
    write(output/"execution-status.json", {"inputs_and_sources_stable": stable, "complete": status == 0 and stable,
                                          "reason": "user-interrupted" if interrupted else "technical-failure" if status or not stable else "completed"})
    require(status == 0 and stable, "Incomplete execution or changed provenance; retained completed cases.")
    summarize(output)


def compare(left, right, output):
    for directory in (left, right): summarize(directory)
    excluded = {"inputs.json", "environment.json", "execution.json", "raw-artifact-index.json", "reproducibility.json"}
    def inventory(directory):
        return {str(p.relative_to(directory)) for p in directory.rglob("*") if p.suffix in (".json", ".csv") and p.name not in excluded}
    a, b = inventory(left), inventory(right); differences = sorted(a ^ b)
    for name in sorted(a & b):
        same = (joint.scientific(read(left/name)) == joint.scientific(read(right/name))) if name.endswith(".json") else fixed.fold.sha256_file(left/name) == fixed.fold.sha256_file(right/name)
        if not same: differences.append(name)
    write(output, {"passed": not differences, "case_count": 72, "differences": differences,
                   "comparison": "exact scientific JSON and CSV bytes", "excluded": ["timing", "process RSS high-water marks", "execution paths/environment"]})
    require(not differences, f"Independent repeat differs: {differences}")


def plots(directory, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    rows = fixed.load_table(directory/"estimates.csv")
    fig, axes = plt.subplots(3, 2, figsize=(12, 9), constrained_layout=True)
    for q, parameter in enumerate(("A", "B", "C")):
        for col, precision in enumerate(("double", "float32")):
            ax = axes[q, col]
            for start in STARTS:
                selected = rows[(rows["dataset"] == DATASETS[0]) & (rows["case"] == f"{start}-{precision}")]
                ax.plot(selected["serial_id"], selected[parameter+"_error"], label=start, linewidth=.8)
            ax.set(xlabel="Atom serial ID", ylabel=parameter+" error", title=precision); ax.grid(alpha=.2)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            if q == 0: ax.legend(fontsize=8)
    for suffix in ("png", "pdf"): fig.savefig(output/f"parameter-errors.{suffix}", dpi=160)
    plt.close(fig)
    results = read(directory/"results.json")["datasets"]
    matrix = np.array([[int(c["joint_qualified"]) for c in d["cases"]] for d in results])
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    axes[0].imshow(matrix, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
    axes[0].set(yticks=range(9), yticklabels=DATASETS, xticks=range(8), xticklabels=CASES, title="Joint qualification (green = passed)")
    axes[0].tick_params(axis="x", labelrotation=70)
    for name in DATASETS[1:]:
        evidence = read(directory/"datasets"/name/"truth-identifiability.json")
        axes[1].semilogy(range(1, 13), np.maximum(evidence["width_singular_values"], 1e-18), marker=".", label=name)
    axes[1].set(xlabel="Singular value index", ylabel="Projected width sensitivity at truth", title="Independent 12-atom diagnostics")
    axes[1].legend(fontsize=8); axes[1].grid(alpha=.2)
    for suffix in ("png", "pdf"): fig.savefig(output/f"coverage-spectrum.{suffix}", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__); commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("run")
    for name in ("executable", "model", "map", "manifest", "output"): p.add_argument("--"+name, type=Path, required=True)
    p = commands.add_parser("summarize"); p.add_argument("directory", type=Path)
    p = commands.add_parser("compare")
    for name in ("left", "right", "output"): p.add_argument("--"+name, type=Path, required=True)
    p = commands.add_parser("plots"); p.add_argument("directory", type=Path); p.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "run": run(args)
    elif args.command == "summarize": summarize(args.directory)
    elif args.command == "compare": compare(args.left, args.right, args.output)
    else: plots(args.directory, args.output)


if __name__ == "__main__": main()
