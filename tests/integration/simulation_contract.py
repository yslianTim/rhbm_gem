"""Explicit v2 simulation semantics and hash-identified historical v1 records."""
from pathlib import Path
import math
import re
import fold_168_regression as fold

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = {"version": "sphere-fma-v1", "comparison": "squared_distance<=cutoff*cutoff",
           "coordinates": "fma(index,spacing,origin)", "squared_distance": "fma(dz,dz,fma(dy,dy,dx*dx))"}
POLICY = {"id": "element-scaled-v1", "oxygen": .8, "nitrogen": .9, "other": 1.}


def resolve(path, paths):
    value = fold.read_json(path); hashes = {k: fold.sha256_file(v) for k, v in paths.items()}
    fold.require(value["source"]["model_sha256"] == hashes["model"] and value["output"]["map_sha256"] == hashes["map"], "Simulation input hash mismatch.")
    version = value.get("schema_version")
    if version == 1:
        for fixture_path in ("joint_abc_coverage.json", "fold_168_simulation_baseline.json"):
            fixture = fold.read_json(ROOT/"tests/benchmarks"/fixture_path)
            if hashes == fixture["input_hashes"]:
                fold.load_simulation_manifest(path, hashes)
                widths = fixture.get("width_contract", {}).get("b", [value["settings"]["blurring_width"]]*len(value["atoms"]))
                return value, {"source": "hash-identified-v1", "fixture": fixture_path, "widths": widths}
        raise RuntimeError("Unknown v1 width contract; an explicit frozen input hash is required.")
    fold.require(type(version) is int and version == 2, "Unsupported simulation manifest version.")
    settings = value["settings"]; model = settings["potential_model"]; kernel = value["kernel"]
    fold.require(model in ("single_gaus", "five_gaus_charge", "single_gaus_user"), "Unknown kernel model.")
    fold.require(kernel["version"] == model+"-v1" and "effective_charge_width" not in kernel, "Ambiguous kernel version/global width.")
    fold.require(value["support"] == dict(SUPPORT, outer_cutoff=settings["cutoff_distance"]), "Wrong support/arithmetic contract.")
    fold.require(value["execution"]["accumulation"] == "z_planes_in_preparation_order", "Wrong accumulation order.")
    fold.require(kernel["width_policy"] == (POLICY if model == "single_gaus" else {"id": "model-specific-v1"}), "Wrong width policy.")
    atoms = value["atoms"]; base = settings["blurring_width"]
    fold.require(len(atoms) == value["atom_count"] and [a["preparation_index"] for a in atoms] == list(range(len(atoms))), "Wrong atom preparation order.")
    identities = []
    for atom in atoms:
        fold.validate_atom_identity(atom)
        identities.append(tuple(atom[k] for k in ("serial_id", "chain_id", "sequence_id", "component_id", "atom_id", "alternate_indicator")))
    fold.require(len(set(identities)) == len(atoms), "Duplicate atom identity.")
    for key in ("source_sha256", "configuration_sha256", "build_sha256"):
        fold.require(re.fullmatch(r"[0-9a-f]{64}", value["generator"][key]) is not None, "Invalid generator hash.")
    widths = []
    for a in atoms:
        width = base*(.8 if a["element"] == 8 else .9 if a["element"] == 7 else 1.)
        gaussian = width if model == "single_gaus" else None
        charge = width if model == "single_gaus" else max(base, 1e-5) if model == "five_gaus_charge" and base != 0 else None
        fold.require(a["effective_gaussian_width"] == gaussian and a["effective_charge_width"] == charge, "Wrong effective atom widths.")
        fold.require(gaussian is None or math.isfinite(gaussian) and gaussian > 0, "Invalid Gaussian width.")
        widths.append(gaussian)
    if model == "single_gaus":
        fold.require(kernel["near_zero_distance"] == 1e-5 and kernel["charge_term_cutoff"] == 2.5 and kernel["minimum_charge_width"] is None, "Changed single_gaus kernel.")
    elif model == "five_gaus_charge" and base != 0:
        fold.require(kernel["near_zero_distance"] == 1e-5 and kernel["charge_term_cutoff"] == 3.0 and kernel["minimum_charge_width"] == 1e-5, "Changed five_gaus kernel.")
    else:
        fold.require(all(kernel[k] is None for k in ("near_zero_distance", "charge_term_cutoff", "minimum_charge_width")), "Changed user kernel.")
    return value, {"source": "manifest-v2", "widths": widths, "support": value["support"]}
