#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sqlite3
import subprocess
import tempfile
import time
from contextlib import closing
from pathlib import Path
from typing import Any, Sequence


SCHEMA_VERSION = 7
SAVED_KEY = "fold_gaus_charge1"
COMMAND_ARGUMENT_TEMPLATE = [
    "potential_analysis",
    "-d", "{database}",
    "-a", "{model}",
    "-m", "{map}",
    "-k", SAVED_KEY,
    "-j", "4",
    "-v", "3",
    "--simulation", "true",
    "-r", "0.50",
    "--exclude-hydrogen", "true",
]
EXPECTED_ATOM_COUNT = 168
MAXIMUM_ACCEPTED_ITERATIONS = 25
IDENTITY_FIELDS = (
    "serial_id", "chain_id", "sequence_id", "component_id", "atom_id",
    "alternate_indicator",
)
TRUTH_DEFINITION = {
    "amplitude": "element_atomic_number",
    "width": "settings.blurring_width",
    "offset": "atoms[].charge_used",
}
CHARGE_FAILURE_STATUSES = {
    "unsupported_residue", "unsupported_structure", "unsupported_spot",
    "table_data_mismatch",
}
SECOND_STAGE_SUMMARY_PATTERN = re.compile(
    r"Second-stage local fitting summary: "
    r"accepted_iterations=(?P<accepted_iterations>\d+), "
    r"best_iteration=(?P<best_iteration>initial|unavailable|\d+), "
    r"stop_reason=(?P<stop_reason>[a-z-]+), "
    r"best_audit_objective=(?P<best_audit_objective>unavailable|"
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?), "
    r"final_uses_polish=(?P<final_uses_polish>yes|no|unavailable), "
    r"final_state_source=(?P<final_state_source>"
    r"best-audit|latest-validated|unavailable)\."
)
SECOND_STAGE_MULTILINE_SUMMARY_PATTERN = re.compile(
    r"Second-Stage Local Fitting Summary\s*:\s*\n"
    r"\s*- accepted_iterations\s*=\s*(?P<accepted_iterations>\d+)\s*\n"
    r"\s*- best_iteration\s*=\s*(?P<best_iteration>initial|unavailable|\d+)\s*\n"
    r"\s*- stop_reason\s*=\s*(?P<stop_reason>[a-z-]+)\s*\n"
    r"\s*- best_audit_objective\s*=\s*(?P<best_audit_objective>unavailable|"
    r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)\s*\n"
    r"\s*- final_uses_polish\s*=\s*(?P<final_uses_polish>yes|no|unavailable)\s*\n"
    r"\s*- final_state_source\s*=\s*(?P<final_state_source>"
    r"best-audit|latest-validated|unavailable)",
    re.MULTILINE,
)
ATOM_CUTOFF_PATTERN = re.compile(
    r"Local-fitting atom cutoff: "
    r"atoms=(?P<atom_count>\d+), "
    r"limit=(?P<limit>\d+), "
    r"clusters=(?P<cluster_count>\d+), "
    r"max-atoms=(?P<maximum_atom_count>\d+), "
    r"cutoff-edges=(?P<cut_edge_count>\d+)\."
)


class RegressionError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RegressionError(message)


def require_fields(value: Any, fields: dict[str, type], context: str) -> None:
    require(isinstance(value, dict), f"{context}: expected an object.")
    for name, field_type in fields.items():
        require(type(value.get(name)) is field_type,
                f"{context}.{name}: missing or invalid {field_type.__name__}.")


def finite_number(value: Any, context: str) -> float:
    require(type(value) in (int, float), f"{context}: expected a finite number.")
    try:
        number = float(value)
    except OverflowError as error:
        raise RegressionError(f"{context}: number outside float64 range.") from error
    require(math.isfinite(number), f"{context}: non-finite number.")
    return number


def require_vector(value: Any, context: str, positive: bool = False) -> None:
    require(isinstance(value, list) and len(value) == 3,
            f"{context}: expected three coordinates.")
    for coordinate in value:
        number = finite_number(coordinate, context)
        require(not positive or number > 0.0, f"{context}: expected positive values.")


def read_json(path: Path) -> dict[str, Any]:
    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            require(key not in result, f"Duplicate JSON key: {key}.")
            result[key] = value
        return result

    def check_finite(value: Any) -> None:
        if isinstance(value, dict):
            for item in value.values():
                check_finite(item)
        elif isinstance(value, list):
            for item in value:
                check_finite(item)
        elif type(value) in (int, float):
            finite_number(value, str(path))

    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_object)
    except (OSError, ValueError) as error:
        raise RegressionError(f"Failed to read JSON {path}: {error}") from error
    require(isinstance(value, dict), f"{path}: expected a JSON object.")
    check_finite(value)
    return value


def validate_atom_identity(atom: dict[str, Any]) -> None:
    require_fields(atom, {
        "serial_id": int, "sequence_id": int, "chain_id": str,
        "component_id": str, "atom_id": str, "alternate_indicator": str,
        "element": int, "structure": int, "position": list,
    }, "atom")
    require(1 <= atom["element"] <= 118, "atom.element: invalid atomic number.")
    require_vector(atom["position"], "atom.position")


def load_simulation_manifest(path: Path, input_hashes: dict[str, str]) -> dict[str, Any]:
    manifest = read_json(path)
    require_fields(manifest, {
        "schema_version": int, "generator": dict, "source": dict, "output": dict,
        "settings": dict, "kernel": dict, "execution": dict,
        "charge_semantics": str, "atom_count": int, "fallback_charge_count": int,
        "atoms": list,
    }, "manifest")
    require(manifest["schema_version"] == 1, "Unsupported simulation manifest schema version.")
    require(manifest["charge_semantics"] == "argument_passed_to_electric_potential",
            "Unsupported charge semantics.")
    require_fields(manifest["generator"], {
        "version": str, "source_sha256": str, "configuration_sha256": str,
        "build_sha256": str,
    }, "generator")
    for name in ("source_sha256", "configuration_sha256", "build_sha256"):
        require(re.fullmatch(r"[0-9a-f]{64}", manifest["generator"][name]) is not None,
                f"generator.{name}: invalid SHA-256.")
    require_fields(manifest["source"], {
        "model_path": str, "pdb_id": str, "model_sha256": str,
    }, "source")
    require_fields(manifest["output"], {
        "map_file": str, "map_sha256": str, "format": str,
        "calculation_precision": str, "storage_precision": str,
    }, "output")
    # Original names are provenance; a renamed or relocated pair is identified by bytes.
    require(manifest["source"]["model_sha256"] == input_hashes["model"],
            "Manifest model SHA-256 mismatch.")
    require(manifest["output"]["map_sha256"] == input_hashes["map"],
            "Manifest map SHA-256 mismatch.")
    for name, expected in (("format", "ccp4"), ("calculation_precision", "float64"),
                           ("storage_precision", "float32")):
        require(manifest["output"][name] == expected, f"Unsupported output.{name}.")

    settings = manifest["settings"]
    require_fields(settings, {
        "potential_model": str, "potential_model_code": int,
        "charge_mode": str, "charge_mode_code": int, "blurring_width_list": list,
        "grid_spacing": list, "grid_size": list, "origin": list,
        "coordinate_unit": str, "missing_charge_policy": str,
        "exclude_hydrogen": bool, "only_backbone": bool,
        "occupancy_applied": bool, "temperature_factor_applied": bool,
        "normalization_applied": bool,
    }, "settings")
    require(settings["potential_model"] == "single_gaus" and settings["potential_model_code"] == 0,
            "Only single_gaus parameter truth is supported.")
    require((settings["charge_mode"], settings["charge_mode_code"]) in
            (("neutral", 0), ("partial", 1), ("amber", 2)), "Invalid charge mode.")
    require(settings["coordinate_unit"] == "angstrom" and
            settings["missing_charge_policy"] == "zero_with_status", "Unsupported generation policy.")
    require(not any(settings[name] for name in
                    ("occupancy_applied", "temperature_factor_applied", "normalization_applied")),
            "Transformed simulation values do not support direct parameter truth.")
    width = finite_number(settings.get("blurring_width"), "settings.blurring_width")
    require(width > 0.0, "Blurring width must be positive.")
    require(settings["blurring_width_list"] and all(
        finite_number(w, "settings.blurring_width_list") > 0.0
        for w in settings["blurring_width_list"]), "Invalid blurring width list.")
    require(width in settings["blurring_width_list"], "Blurring width absent from generation list.")
    require(finite_number(settings.get("cutoff_distance"), "settings.cutoff_distance") > 0.0,
            "Cutoff distance must be positive.")
    require_vector(settings["grid_spacing"], "settings.grid_spacing", positive=True)
    require_vector(settings["origin"], "settings.origin")
    require(len(settings["grid_size"]) == 3 and all(type(n) is int and n > 0
            for n in settings["grid_size"]), "Invalid grid size.")
    kernel = manifest["kernel"]
    for name in ("charge_term_cutoff", "near_zero_distance", "effective_charge_width"):
        require(finite_number(kernel.get(name), f"kernel.{name}") > 0.0,
                f"kernel.{name}: must be positive.")
    require("minimum_charge_width" in kernel and kernel["minimum_charge_width"] is None and
            kernel["effective_charge_width"] == width, "Inconsistent single_gaus kernel width.")
    require_fields(manifest["execution"], {
        "openmp_enabled": bool, "requested_job_count": int, "actual_job_count": int,
        "accumulation": str,
    }, "execution")
    require(manifest["execution"]["requested_job_count"] > 0 and
            manifest["execution"]["actual_job_count"] > 0 and
            manifest["execution"]["accumulation"] == "z_planes_in_preparation_order",
            "Unsupported execution settings.")

    atoms = manifest["atoms"]
    require(manifest["atom_count"] == len(atoms) and bool(atoms), "Invalid manifest atom count.")
    serial_ids: set[int] = set()
    indices: set[int] = set()
    fallback_count = 0
    for atom in atoms:
        validate_atom_identity(atom)
        require_fields(atom, {"preparation_index": int, "residue": int, "spot": int,
                              "lookup_status": str}, "atom")
        require(atom["serial_id"] not in serial_ids, "Duplicate atom serial ID in manifest.")
        require(atom["preparation_index"] not in indices, "Duplicate preparation index.")
        serial_ids.add(atom["serial_id"])
        indices.add(atom["preparation_index"])
        require("table" in atom and "lookup_charge" in atom, "Missing charge lookup evidence.")
        table, status = atom["table"], atom["lookup_status"]
        charge = finite_number(atom.get("charge_used"), "atom.charge_used")
        if settings["charge_mode"] == "neutral":
            require(status == "neutral_mode" and table is None and
                    atom["lookup_charge"] is None and charge == 0.0, "Invalid neutral charge evidence.")
        else:
            allowed_tables = ("amber95",) if settings["charge_mode"] == "amber" else (
                "buried", "helix", "sheet")
            require(table is None or table in allowed_tables, "Invalid charge table category.")
            if status == "found":
                require(table is not None and
                        finite_number(atom["lookup_charge"], "atom.lookup_charge") == charge,
                        "Successful lookup charge differs from charge_used.")
            else:
                require(status in CHARGE_FAILURE_STATUSES and atom["lookup_charge"] is None and
                        charge == 0.0, "Invalid zero fallback evidence.")
                fallback_count += 1
    require(indices == set(range(len(atoms))), "Preparation indices must cover 0..atom_count-1.")
    require(fallback_count == manifest["fallback_charge_count"], "Incorrect fallback charge count.")
    return manifest


def read_atom_results(database_path: Path) -> list[dict[str, Any]]:
    query = """
        SELECT local.serial_id, atom.chain_id, atom.sequence_id, atom.component_id,
               atom.atom_id, atom.indicator AS alternate_indicator,
               atom.element, atom.structure, atom.position_x, atom.position_y, atom.position_z,
               local.amplitude_estimate_mdpde_2nd AS amplitude_mdpde,
               local.width_estimate_mdpde_2nd AS width_mdpde,
               local.intercept_estimate_mdpde_2nd AS intercept_mdpde,
               local.alpha_r_2nd AS alpha_r
          FROM model_atom_local_potential AS local
          LEFT JOIN model_atom AS atom
            ON atom.key_tag = local.key_tag
           AND atom.serial_id = local.serial_id
         WHERE local.key_tag = ?
         ORDER BY local.serial_id
    """
    with closing(sqlite3.connect(database_path.resolve().as_uri() + "?mode=ro", uri=True)) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(query, (SAVED_KEY,)).fetchall()

    atoms: list[dict[str, Any]] = []
    seen_serial_ids: set[int] = set()
    for row in rows:
        atom = dict(row)
        atom["position"] = [atom.pop(f"position_{axis}") for axis in "xyz"]
        validate_atom_identity(atom)
        serial_id = atom["serial_id"]
        require(serial_id not in seen_serial_ids, f"Duplicate atom serial ID in output database: {serial_id}.")
        seen_serial_ids.add(serial_id)
        for name in ("amplitude_mdpde", "width_mdpde", "intercept_mdpde", "alpha_r"):
            atom[name] = finite_number(atom[name], f"atom {serial_id}.{name}")
        require(atom["amplitude_mdpde"] > 0.0 and atom["width_mdpde"] > 0.0,
                f"Atom {serial_id}: amplitude and width must be positive.")
        require(0.0 <= atom["alpha_r"] <= 1.0, f"Atom {serial_id}: alpha_r is outside [0, 1].")
        atoms.append(atom)
    return atoms


def pair_atom_truth(atoms: Sequence[dict[str, Any]], manifest: dict[str, Any]) -> list[dict[str, Any]]:
    def identity(atom: dict[str, Any]) -> tuple[Any, ...]:
        return tuple(atom[name] for name in IDENTITY_FIELDS)

    estimates = {identity(atom): atom for atom in atoms}
    truths = {identity(atom): atom for atom in manifest["atoms"]}
    require(len(estimates) == len(atoms) and len(truths) == len(manifest["atoms"]),
            "Duplicate atom identity in truth pairing.")
    require(estimates.keys() == truths.keys(),
            f"Atom identity mismatch: missing={list(truths.keys() - estimates.keys())}, "
            f"unexpected={list(estimates.keys() - truths.keys())}.")
    paired: list[dict[str, Any]] = []
    for source in sorted(manifest["atoms"], key=lambda atom: atom["preparation_index"]):
        estimate = estimates[identity(source)]
        for name in ("element", "structure", "position"):
            require(estimate[name] == source[name],
                    f"Atom {source['serial_id']}: {name} differs from simulation truth.")
        truth = {"amplitude": float(source["element"]),
                 "width": manifest["settings"]["blurring_width"], "offset": source["charge_used"]}
        errors = {name: finite_number(estimate[field] - truth[name], f"atom error.{name}")
                  for name, field in (("amplitude", "amplitude_mdpde"), ("width", "width_mdpde"),
                                      ("offset", "intercept_mdpde"))}
        paired.append({**source, **estimate, "truth": truth, "error": errors})
    return paired


def calculate_quality_metrics(atoms: Sequence[dict[str, Any]]) -> dict[str, float]:
    require(bool(atoms), "Cannot calculate fold-168 quality metrics without atoms.")
    metrics: dict[str, float] = {}
    for name in ("amplitude", "width", "offset"):
        errors = [finite_number(atom["error"][name], f"error.{name}") for atom in atoms]
        # Scaling avoids overflowing intermediate squares for finite errors.
        scale = max(abs(error) for error in errors)
        metrics[f"{name}_rmse"] = (scale * math.sqrt(
            math.fsum((error / scale) ** 2 for error in errors) / len(errors)) if scale else 0.0)
    metrics["offset_bias"] = math.fsum(atom["error"]["offset"] / len(atoms) for atom in atoms)
    metrics["offset_max_absolute_error"] = max(abs(atom["error"]["offset"]) for atom in atoms)
    for name, value in metrics.items():
        finite_number(value, name)
    return metrics


def parse_second_stage_summary(log_text: str) -> dict[str, Any]:
    matches = [
        *SECOND_STAGE_SUMMARY_PATTERN.finditer(log_text),
        *SECOND_STAGE_MULTILINE_SUMMARY_PATTERN.finditer(log_text),
    ]
    if len(matches) != 1:
        raise RegressionError(
            "Expected exactly one parseable second-stage final summary, "
            f"found {len(matches)}.")
    values = matches[0].groupdict()
    best_objective_text = values["best_audit_objective"]
    best_objective = (
        None if best_objective_text == "unavailable" else float(best_objective_text))
    if best_objective is not None and not math.isfinite(best_objective):
        raise RegressionError("Second-stage best audit objective is not finite.")
    final_uses_polish_text = values["final_uses_polish"]
    final_uses_polish = (
        None if final_uses_polish_text == "unavailable" else
        final_uses_polish_text == "yes")
    return {
        "accepted_iterations": int(values["accepted_iterations"]),
        "best_iteration": values["best_iteration"],
        "stop_reason": values["stop_reason"],
        "best_audit_objective": best_objective,
        "final_uses_polish": final_uses_polish,
        "final_state_source": values["final_state_source"],
    }


def parse_atom_cutoff_summary(log_text: str) -> dict[str, int]:
    matches = list(ATOM_CUTOFF_PATTERN.finditer(log_text))
    if len(matches) != 1:
        raise RegressionError(
            "Expected exactly one parseable atom cutoff summary, "
            f"found {len(matches)}.")
    return {
        name: int(value)
        for name, value in matches[0].groupdict().items()
    }


def evaluate_gates(baseline: dict[str, Any], actual: dict[str, Any]) -> dict[str, Any]:
    atom_differences: list[str] = []
    atoms = actual["atoms"]
    if len(atoms) != EXPECTED_ATOM_COUNT or sorted(atom["serial_id"] for atom in atoms) != baseline["serial_ids"]:
        atom_differences.append("atoms: expected the complete fixed 168-atom set")
    cutoff = actual["atom_cutoff_summary"]
    if cutoff is None:
        atom_differences.append("atom_cutoff_summary: missing")
    else:
        checks = {
            "atom_count": cutoff["atom_count"] == baseline["expected_atom_count"],
            "limit": cutoff["limit"] == baseline["maximum_atoms_per_cluster"],
            "cluster_count": cutoff["cluster_count"] >= baseline["minimum_cluster_count"],
            "maximum_atom_count": 0 < cutoff["maximum_atom_count"] <= baseline["maximum_atoms_per_cluster"],
        }
        atom_differences.extend(f"atom_cutoff_summary.{name}: violates baseline ({cutoff[name]})"
                                for name, passed in checks.items() if not passed)
    summary = actual["second_stage_summary"]
    iteration_differences: list[str] = []
    if summary is None:
        iteration_differences.append("second_stage_summary: missing")
    elif summary["accepted_iterations"] > MAXIMUM_ACCEPTED_ITERATIONS:
        iteration_differences.append(
            f"second_stage_summary.accepted_iterations: expected <= {MAXIMUM_ACCEPTED_ITERATIONS}, "
            f"got {summary['accepted_iterations']}")
    return {
        "quality_gate": {
            "status": "uncalibrated", "passed": False,
            "differences": ["quality_gate: uncalibrated; no parameter accuracy threshold has been established"],
        },
        "iteration_gate": {"passed": not iteration_differences, "differences": iteration_differences},
        "atom_cutoff_gate": {"passed": not atom_differences, "differences": atom_differences},
    }


def validate_input_hashes(paths: dict[str, Path], expected_hashes: dict[str, str]) -> dict[str, str]:
    actual_hashes: dict[str, str] = {}
    failures: list[str] = []
    for name, expected_hash in expected_hashes.items():
        path = paths[name]
        if not path.is_file():
            failures.append(f"{name} input does not exist: {path}")
            continue
        actual_hash = sha256_file(path)
        actual_hashes[name] = actual_hash
        if actual_hash != expected_hash:
            failures.append(f"{name} SHA-256 mismatch: expected {expected_hash}, got {actual_hash}")
    require(not failures, "\n".join(failures))
    return actual_hashes


def validate_fixture(manifest: dict[str, Any], baseline: dict[str, Any]) -> None:
    for name, expected in baseline["generation_contract"].items():
        require(manifest.get(name) == expected, f"Simulation {name} differs from the fixed benchmark contract.")
    require(manifest["atom_count"] == EXPECTED_ATOM_COUNT and
            sorted(atom["serial_id"] for atom in manifest["atoms"]) == baseline["serial_ids"],
            "Simulation atoms differ from the fixed 168-atom fixture.")


def build_command(
    executable: Path,
    database: Path,
    model: Path,
    map_path: Path,
) -> list[str]:
    substitutions = {
        "database": str(database),
        "model": str(model),
        "map": str(map_path),
    }
    return [
        str(executable),
        *(argument.format(**substitutions) for argument in COMMAND_ARGUMENT_TEMPLATE),
    ]


def make_empty_actual(input_hashes: dict[str, str]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "input_hashes": input_hashes,
        "command_arguments": COMMAND_ARGUMENT_TEMPLATE,
        "atoms": [],
        "quality_metrics": None,
        "diagnostics": None,
        "generation_record": None,
        "truth_definition": TRUTH_DEFINITION,
        "estimate_source": "sqlite_final_second_stage_mdpde",
        "second_stage_summary": None,
        "atom_cutoff_summary": None,
    }


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=False, allow_nan=False) + "\n",
        encoding="utf-8")


def load_baseline(path: Path) -> dict[str, Any]:
    baseline = read_json(path)
    require(type(baseline.get("schema_version")) is int and baseline["schema_version"] == SCHEMA_VERSION,
            f"Unsupported baseline schema version: {baseline.get('schema_version')!r}.")
    hashes = baseline.get("input_hashes")
    require(isinstance(hashes, dict) and set(hashes) == {"model", "map", "manifest"} and all(
        isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) for value in hashes.values()),
        "Baseline must pin model, map and manifest SHA-256 identities.")
    require(baseline.get("command_arguments") == COMMAND_ARGUMENT_TEMPLATE,
            "Baseline command arguments do not match the fixed benchmark command.")
    serial_ids = baseline.get("serial_ids")
    require(isinstance(serial_ids, list) and len(serial_ids) == EXPECTED_ATOM_COUNT and
            all(type(value) is int for value in serial_ids) and serial_ids == sorted(set(serial_ids)),
            "Baseline must contain 168 unique sorted serial IDs.")
    require_fields(baseline, {"expected_atom_count": int, "maximum_atoms_per_cluster": int,
                             "minimum_cluster_count": int, "generation_contract": dict}, "baseline")
    require(baseline["expected_atom_count"] == EXPECTED_ATOM_COUNT and
            baseline["maximum_atoms_per_cluster"] > 0 and baseline["minimum_cluster_count"] > 0 and
            baseline["maximum_atoms_per_cluster"] * baseline["minimum_cluster_count"] >= EXPECTED_ATOM_COUNT,
            "Invalid baseline atom cutoff expectations.")
    require_fields(baseline["generation_contract"], {
        "settings": dict, "kernel": dict, "atom_count": int, "fallback_charge_count": int,
        "charge_semantics": str,
    }, "generation_contract")
    settings = baseline["generation_contract"]["settings"]
    require(settings.get("blurring_width") == float(
        COMMAND_ARGUMENT_TEMPLATE[COMMAND_ARGUMENT_TEMPLATE.index("-r") + 1]) and
        settings.get("exclude_hydrogen") is True and settings.get("only_backbone") is False,
        "Generation settings do not match the fixed analysis command.")
    require(baseline.get("truth") == TRUTH_DEFINITION, "Unsupported baseline truth definition.")
    require(baseline.get("quality_gate") == {"status": "uncalibrated"} and
            "reference_quality_metrics" not in baseline, "Quality thresholds must remain uncalibrated.")
    return baseline


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the external fold-168 regression benchmark.")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--map", dest="map_path", type=Path, required=True)
    parser.add_argument("--simulation-manifest", type=Path,
                        help="Generation record; defaults to <map>.simulation.json")
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def run(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.output_dir / "run.log"
    actual_path = args.output_dir / "actual.json"
    report_path = args.output_dir / "report.json"
    input_paths = {
        "model": args.model.resolve(),
        "map": args.map_path.resolve(),
        "manifest": (args.simulation_manifest or Path(str(args.map_path) + ".simulation.json")).resolve(),
    }
    actual_hashes: dict[str, str] = {}
    actual = make_empty_actual(actual_hashes)
    command: list[str] = []
    errors: list[str] = []
    differences: list[str] = []
    wall_time_seconds: float | None = None
    log_text = ""
    scoring_status = "failed"
    gates = {
        "quality_gate": {"status": "uncalibrated", "passed": False},
        "iteration_gate": {"passed": None},
        "atom_cutoff_gate": {"passed": None},
    }

    try:
        baseline = load_baseline(args.baseline.resolve())
        executable = args.executable.resolve()
        if not executable.is_file():
            raise RegressionError(f"Benchmark executable does not exist: {executable}")
        actual_hashes = validate_input_hashes(input_paths, baseline["input_hashes"])
        actual["input_hashes"] = actual_hashes
        manifest = load_simulation_manifest(input_paths["manifest"], actual_hashes)
        validate_fixture(manifest, baseline)
        actual = make_empty_actual(actual_hashes)
        actual["generation_record"] = {name: value for name, value in manifest.items() if name != "atoms"}

        with tempfile.TemporaryDirectory(prefix="rhbm_fold_168_regression_") as temp_dir:
            temporary_database = Path(temp_dir) / "database.sqlite"
            if temporary_database.exists():
                raise RegressionError(
                    f"Temporary output database already exists: {temporary_database}")
            command = build_command(
                executable,
                temporary_database,
                input_paths["model"],
                input_paths["map"])
            start_time = time.perf_counter()
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
                cwd=temp_dir)
            wall_time_seconds = time.perf_counter() - start_time
            log_text = completed.stdout.decode("utf-8", errors="replace")
            log_path.write_text(log_text, encoding="utf-8")
            if completed.returncode != 0:
                raise RegressionError(
                    f"Benchmark command exited with status {completed.returncode}.")
            if not temporary_database.is_file():
                raise RegressionError(
                    "Benchmark command did not create the temporary output database.")

            # Do not attach truth to an execution whose inputs changed during fitting.
            validate_input_hashes(input_paths, actual_hashes)
            paired = pair_atom_truth(read_atom_results(temporary_database), manifest)
            metrics = calculate_quality_metrics(paired)
            actual["atoms"] = paired
            actual["quality_metrics"] = metrics
            actual["diagnostics"] = {
                "maximum_absolute_offset": max(abs(atom["intercept_mdpde"]) for atom in paired),
                "fallback_charge_count": manifest["fallback_charge_count"],
            }
            scoring_status = "complete"
            actual["second_stage_summary"] = parse_second_stage_summary(log_text)
            actual["atom_cutoff_summary"] = parse_atom_cutoff_summary(log_text)

        gates = evaluate_gates(baseline, actual)
        differences = [message for gate in gates.values() for message in gate["differences"]]
    except Exception as error:  # Preserve artifacts for all benchmark failures.
        errors.append(str(error))

    log_path.write_text(log_text, encoding="utf-8")
    write_json(actual_path, actual)
    passed = not errors and all(gate["passed"] is True for gate in gates.values())
    report = {
        "schema_version": SCHEMA_VERSION,
        "passed": passed,
        "truth_scoring": {"status": scoring_status},
        **gates,
        "stop_reason": (actual["second_stage_summary"] or {}).get("stop_reason"),
        "performance_gate": False,
        "wall_time_seconds": wall_time_seconds,
        "command": command,
        "errors": errors,
        "differences": differences,
        "artifacts": {
            "log": str(log_path),
            "actual": str(actual_path),
        },
    }
    write_json(report_path, report)

    if passed:
        print(
            "fold-168 regression passed; "
            f"wall time = {wall_time_seconds:.3f} s (observation only).")
        return 0
    print(f"fold-168 regression failed; see {report_path}.")
    for message in [*errors, *differences[:20]]:
        print(f"- {message}")
    if len(differences) > 20:
        print(f"- ... {len(differences) - 20} additional differences")
    return 1


def main() -> int:
    return run()


if __name__ == "__main__":
    raise SystemExit(main())
