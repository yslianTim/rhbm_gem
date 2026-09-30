"""Resolve current Joint simulation manifests and hash-pinned v1 records."""
from pathlib import Path
import json
import math
import re

from experiment_io import sha256_file

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = {"version": "sphere-fma-v1", "comparison": "squared_distance<=cutoff*cutoff",
           "coordinates": "fma(index,spacing,origin)", "squared_distance": "fma(dz,dz,fma(dy,dy,dx*dx))"}
POLICY = {"id": "element-scaled-v1", "oxygen": .8, "nitrogen": .9, "other": 1.}
_CHARGE_FAILURE_STATUSES = {"unsupported_residue", "unsupported_structure", "unsupported_spot",
                            "table_data_mismatch"}


def _require(condition, message):
    if not condition:
        raise RuntimeError(message)


def _require_fields(value, fields, context):
    _require(isinstance(value, dict), f"{context}: expected an object.")
    for name, field_type in fields.items():
        _require(type(value.get(name)) is field_type,
                 f"{context}.{name}: missing or invalid {field_type.__name__}.")


def _finite_number(value, context):
    _require(type(value) in (int, float), f"{context}: expected a finite number.")
    try:
        number = float(value)
    except OverflowError as error:
        raise RuntimeError(f"{context}: number outside float64 range.") from error
    _require(math.isfinite(number), f"{context}: non-finite number.")
    return number


def _require_vector(value, context, positive=False):
    _require(isinstance(value, list) and len(value) == 3,
             f"{context}: expected three coordinates.")
    for coordinate in value:
        number = _finite_number(coordinate, context)
        _require(not positive or number > 0.0, f"{context}: expected positive values.")


def _read_json(path):
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            _require(key not in result, f"Duplicate JSON key: {key}.")
            result[key] = value
        return result

    def check_finite(value):
        if isinstance(value, dict):
            for item in value.values():
                check_finite(item)
        elif isinstance(value, list):
            for item in value:
                check_finite(item)
        elif type(value) in (int, float):
            _finite_number(value, str(path))

    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique_object)
    except (OSError, ValueError) as error:
        raise RuntimeError(f"Failed to read JSON {path}: {error}") from error
    _require(isinstance(value, dict), f"{path}: expected a JSON object.")
    check_finite(value)
    return value


def _validate_atom_identity(atom):
    _require_fields(atom, {
        "serial_id": int, "sequence_id": int, "chain_id": str,
        "component_id": str, "atom_id": str, "alternate_indicator": str,
        "element": int, "structure": int, "position": list,
    }, "atom")
    _require(1 <= atom["element"] <= 118, "atom.element: invalid atomic number.")
    _require_vector(atom["position"], "atom.position")


def _load_v1_manifest(path, input_hashes):
    manifest = _read_json(path)
    _require_fields(manifest, {
        "schema_version": int, "generator": dict, "source": dict, "output": dict,
        "settings": dict, "kernel": dict, "execution": dict,
        "charge_semantics": str, "atom_count": int, "fallback_charge_count": int,
        "atoms": list,
    }, "manifest")
    _require(manifest["schema_version"] == 1, "Unsupported simulation manifest schema version.")
    _require(manifest["charge_semantics"] == "argument_passed_to_electric_potential",
             "Unsupported charge semantics.")
    _require_fields(manifest["generator"], {
        "version": str, "source_sha256": str, "configuration_sha256": str,
        "build_sha256": str,
    }, "generator")
    for name in ("source_sha256", "configuration_sha256", "build_sha256"):
        _require(re.fullmatch(r"[0-9a-f]{64}", manifest["generator"][name]) is not None,
                 f"generator.{name}: invalid SHA-256.")
    _require_fields(manifest["source"], {
        "model_path": str, "pdb_id": str, "model_sha256": str,
    }, "source")
    _require_fields(manifest["output"], {
        "map_file": str, "map_sha256": str, "format": str,
        "calculation_precision": str, "storage_precision": str,
    }, "output")
    _require(manifest["source"]["model_sha256"] == input_hashes["model"],
             "Manifest model SHA-256 mismatch.")
    _require(manifest["output"]["map_sha256"] == input_hashes["map"],
             "Manifest map SHA-256 mismatch.")
    for name, expected in (("format", "ccp4"), ("calculation_precision", "float64"),
                           ("storage_precision", "float32")):
        _require(manifest["output"][name] == expected, f"Unsupported output.{name}.")

    settings = manifest["settings"]
    _require_fields(settings, {
        "potential_model": str, "potential_model_code": int,
        "charge_mode": str, "charge_mode_code": int, "blurring_width_list": list,
        "grid_spacing": list, "grid_size": list, "origin": list,
        "coordinate_unit": str, "missing_charge_policy": str,
        "exclude_hydrogen": bool, "only_backbone": bool,
        "occupancy_applied": bool, "temperature_factor_applied": bool,
        "normalization_applied": bool,
    }, "settings")
    _require(settings["potential_model"] == "single_gaus" and settings["potential_model_code"] == 0,
             "Only single_gaus parameter truth is supported.")
    _require((settings["charge_mode"], settings["charge_mode_code"]) in
             (("neutral", 0), ("partial", 1), ("amber", 2)), "Invalid charge mode.")
    _require(settings["coordinate_unit"] == "angstrom" and
             settings["missing_charge_policy"] == "zero_with_status", "Unsupported generation policy.")
    _require(not any(settings[name] for name in
                     ("occupancy_applied", "temperature_factor_applied", "normalization_applied")),
             "Transformed simulation values do not support direct parameter truth.")
    width = _finite_number(settings.get("blurring_width"), "settings.blurring_width")
    _require(width > 0.0, "Blurring width must be positive.")
    _require(settings["blurring_width_list"] and all(
        _finite_number(w, "settings.blurring_width_list") > 0.0
        for w in settings["blurring_width_list"]), "Invalid blurring width list.")
    _require(width in settings["blurring_width_list"], "Blurring width absent from generation list.")
    _require(_finite_number(settings.get("cutoff_distance"), "settings.cutoff_distance") > 0.0,
             "Cutoff distance must be positive.")
    _require_vector(settings["grid_spacing"], "settings.grid_spacing", positive=True)
    _require_vector(settings["origin"], "settings.origin")
    _require(len(settings["grid_size"]) == 3 and all(type(n) is int and n > 0
             for n in settings["grid_size"]), "Invalid grid size.")
    kernel = manifest["kernel"]
    for name in ("charge_term_cutoff", "near_zero_distance", "effective_charge_width"):
        _require(_finite_number(kernel.get(name), f"kernel.{name}") > 0.0,
                 f"kernel.{name}: must be positive.")
    _require("minimum_charge_width" in kernel and kernel["minimum_charge_width"] is None and
             kernel["effective_charge_width"] == width, "Inconsistent single_gaus kernel width.")
    _require_fields(manifest["execution"], {
        "openmp_enabled": bool, "requested_job_count": int, "actual_job_count": int,
        "accumulation": str,
    }, "execution")
    _require(manifest["execution"]["requested_job_count"] > 0 and
             manifest["execution"]["actual_job_count"] > 0 and
             manifest["execution"]["accumulation"] == "z_planes_in_preparation_order",
             "Unsupported execution settings.")

    atoms = manifest["atoms"]
    _require(manifest["atom_count"] == len(atoms) and bool(atoms), "Invalid manifest atom count.")
    serial_ids = set()
    indices = set()
    fallback_count = 0
    for atom in atoms:
        _validate_atom_identity(atom)
        _require_fields(atom, {"preparation_index": int, "residue": int, "spot": int,
                               "lookup_status": str}, "atom")
        _require(atom["serial_id"] not in serial_ids, "Duplicate atom serial ID in manifest.")
        _require(atom["preparation_index"] not in indices, "Duplicate preparation index.")
        serial_ids.add(atom["serial_id"])
        indices.add(atom["preparation_index"])
        _require("table" in atom and "lookup_charge" in atom, "Missing charge lookup evidence.")
        table, status = atom["table"], atom["lookup_status"]
        charge = _finite_number(atom.get("charge_used"), "atom.charge_used")
        if settings["charge_mode"] == "neutral":
            _require(status == "neutral_mode" and table is None and
                     atom["lookup_charge"] is None and charge == 0.0, "Invalid neutral charge evidence.")
        else:
            allowed_tables = ("amber95",) if settings["charge_mode"] == "amber" else (
                "buried", "helix", "sheet")
            _require(table is None or table in allowed_tables, "Invalid charge table category.")
            if status == "found":
                _require(table is not None and
                         _finite_number(atom["lookup_charge"], "atom.lookup_charge") == charge,
                         "Successful lookup charge differs from charge_used.")
            else:
                _require(status in _CHARGE_FAILURE_STATUSES and atom["lookup_charge"] is None and
                         charge == 0.0, "Invalid zero fallback evidence.")
                fallback_count += 1
    _require(indices == set(range(len(atoms))), "Preparation indices must cover 0..atom_count-1.")
    _require(fallback_count == manifest["fallback_charge_count"], "Incorrect fallback charge count.")
    return manifest


def resolve(path, paths):
    value = _read_json(path)
    hashes = {key: sha256_file(source) for key, source in paths.items()}
    _require(value["source"]["model_sha256"] == hashes["model"] and
             value["output"]["map_sha256"] == hashes["map"], "Simulation input hash mismatch.")
    version = value.get("schema_version")
    if version == 1:
        fixture_path = "tests/fixtures/joint_component/simulation-contract.json"
        fixture = _read_json(ROOT / fixture_path)
        if hashes == fixture["input_hashes"]:
            _load_v1_manifest(path, hashes)
            widths = fixture.get("width_contract", {}).get("b",
                [value["settings"]["blurring_width"]] * len(value["atoms"]))
            return value, {"source": "hash-identified-v1", "fixture": fixture_path, "widths": widths}
        raise RuntimeError("Unknown v1 width contract; an explicit frozen input hash is required.")
    _require(type(version) is int and version == 2, "Unsupported simulation manifest version.")
    settings = value["settings"]; model = settings["potential_model"]; kernel = value["kernel"]
    _require(model in ("single_gaus", "five_gaus_charge", "single_gaus_user"), "Unknown kernel model.")
    _require(kernel["version"] == model+"-v1" and "effective_charge_width" not in kernel,
             "Ambiguous kernel version/global width.")
    _require(value["support"] == dict(SUPPORT, outer_cutoff=settings["cutoff_distance"]),
             "Wrong support/arithmetic contract.")
    _require(value["execution"]["accumulation"] == "z_planes_in_preparation_order",
             "Wrong accumulation order.")
    _require(kernel["width_policy"] == (POLICY if model == "single_gaus" else {"id": "model-specific-v1"}),
             "Wrong width policy.")
    atoms = value["atoms"]; base = settings["blurring_width"]
    _require(len(atoms) == value["atom_count"] and
             [atom["preparation_index"] for atom in atoms] == list(range(len(atoms))),
             "Wrong atom preparation order.")
    identities = []
    for atom in atoms:
        _validate_atom_identity(atom)
        identities.append(tuple(atom[key] for key in
            ("serial_id", "chain_id", "sequence_id", "component_id", "atom_id", "alternate_indicator")))
    _require(len(set(identities)) == len(atoms), "Duplicate atom identity.")
    for key in ("source_sha256", "configuration_sha256", "build_sha256"):
        _require(re.fullmatch(r"[0-9a-f]{64}", value["generator"][key]) is not None,
                 "Invalid generator hash.")
    widths = []
    for atom in atoms:
        width = base*(.8 if atom["element"] == 8 else .9 if atom["element"] == 7 else 1.)
        gaussian = width if model == "single_gaus" else None
        charge = width if model == "single_gaus" else max(base, 1e-5) if model == "five_gaus_charge" and base != 0 else None
        _require(atom["effective_gaussian_width"] == gaussian and
                 atom["effective_charge_width"] == charge, "Wrong effective atom widths.")
        _require(gaussian is None or math.isfinite(gaussian) and gaussian > 0,
                 "Invalid Gaussian width.")
        widths.append(gaussian)
    if model == "single_gaus":
        _require(kernel["near_zero_distance"] == 1e-5 and kernel["charge_term_cutoff"] == 2.5 and
                 kernel["minimum_charge_width"] is None, "Changed single_gaus kernel.")
    elif model == "five_gaus_charge" and base != 0:
        _require(kernel["near_zero_distance"] == 1e-5 and kernel["charge_term_cutoff"] == 3.0 and
                 kernel["minimum_charge_width"] == 1e-5, "Changed five_gaus kernel.")
    else:
        _require(all(kernel[key] is None for key in
                     ("near_zero_distance", "charge_term_cutoff", "minimum_charge_width")),
                 "Changed user kernel.")
    return value, {"source": "manifest-v2", "widths": widths, "support": value["support"]}
