"""Exercise the real CLI's joint analysis, persistence and export contract."""
from __future__ import annotations

import hashlib
import csv
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile


def main() -> int:
    executable = str(Path(sys.argv[1]).resolve())
    model = Path(__file__).resolve().parents[1] / "fixtures" / "test_model.cif"
    with tempfile.TemporaryDirectory(prefix="rhbm_joint_cli_") as directory:
        root = Path(directory)
        local_model = root / "input.cif"
        template_text = model.read_text()
        atom_marker = template_text.find("ATOM 1 C CA . ALA A 1")
        assert atom_marker >= 0
        atoms = []
        for index in range(25):
            serial = index + 1
            x, y = 3 + 7 * (index % 5), 3 + 7 * (index // 5)
            atoms.append(f"ATOM {serial} C CA . ALA A {serial} {x} {y} 3.0 1.0 0.0 1")
        local_model.write_text(template_text[:atom_marker] + "\n".join(atoms) + "\n#\n")
        model = local_model

        def run(*args: object, expected_exit: int | None = 0) -> None:
            result = subprocess.run([executable, *map(str, args)], capture_output=True, text=True)
            if expected_exit is None:
                assert result.returncode != 0, result.stdout + result.stderr
            else:
                assert result.returncode == expected_exit, result.stdout + result.stderr

        run("map_simulation", "-a", model, "-o", root, "--potential-model", "single",
            "--blurring-width", ".5", "-g", ".5", "-v", "0")
        map_path = next(root.glob("*.map"))
        generator = json.loads(Path(str(map_path) + ".simulation.json").read_text())["generator"]
        database = root / "joint.sqlite"
        run("potential_analysis", "--select-chain", "A", "-a", model,
            "-m", map_path, "-d", database, "-k", "example", "-v", "0", expected_exit=0)
        run("potential_analysis", "--estimator", "joint-components",
            "--select-chain", "A", "-a", model, "-m", map_path, "-d", database, "-k", "fixed", "-v", "0",
            expected_exit=0)
        model_hash = hashlib.sha256(model.read_bytes()).hexdigest()
        map_hash = hashlib.sha256(map_path.read_bytes()).hexdigest()
        # The same path with different bytes must produce a different fingerprint.
        model.write_text(model.read_text() + "\n# provenance test\n")
        run("potential_analysis", "--estimator", "joint-components", "--select-chain", "A", "-a", model,
            "-m", map_path, "-d", database, "-k", "changed", "-v", "0", expected_exit=0)
        changed_hash = hashlib.sha256(model.read_bytes()).hexdigest()
        run("potential_analysis", "--estimator", "joint-components", "--select-chain", "A", "-a", model,
            "-m", map_path, "-d", database, "-k", "example.contributions", "-v", "0", expected_exit=0)
        map_path.unlink()
        model.unlink()
        run("result_dump", "--printer", "joint", "-d", database, "-k", "example", "-o", root)
        collision = root / "collision"
        collision.mkdir()
        run("result_dump", "--printer", "joint", "-d", database, "-k", "example", "example.contributions",
            "-o", collision, expected_exit=None)
        assert not list(collision.glob("joint_*")), "Collision must fail before writing any result"
        saved = json.loads((root / "joint_result_example.json").read_text())
        with sqlite3.connect(database) as connection:
            assert connection.execute("PRAGMA user_version").fetchone()[0] == 20
            payload = connection.execute("SELECT result_json FROM model_joint_result WHERE key_tag='example'").fetchone()[0]
            fixed = json.loads(connection.execute("SELECT result_json FROM model_joint_result WHERE key_tag='fixed'").fetchone()[0])
            changed = json.loads(connection.execute("SELECT result_json FROM model_joint_result WHERE key_tag='changed'").fetchone()[0])
            stages = json.loads(connection.execute("SELECT result_json FROM model_stage_result WHERE key_tag='example'").fetchone()[0])
            assert stages["version"] == 4
            assert len(stages["atoms"]) == 25
            assert all(atom["peeled"] is None and atom["second_native"] is None for atom in stages["atoms"])
            for group in stages["groups"]:
                if group["summary"] and group["summary"]["inference"]:
                    assert "members" not in group["summary"]["inference"]
            assert all(atom["second"]["source"]["method"] == 3 for atom in stages["atoms"])
            for atom in stages["atoms"]:
                if atom["second"]["source"]["role"] == 1:  # target
                    assert atom["peeling"]["mode"] == "grid-consistent"
                else:
                    assert atom["peeling"] is None
                    assert atom["post_joint_local"] is None
            assert connection.execute("SELECT COUNT(*) FROM sqlite_master WHERE name IN ('model_atom_local_potential','model_atom_posterior','model_atom_group_potential')").fetchone()[0] == 0
        assert saved == json.loads(payload)
        assert saved["schema_version"] == 20
        assert saved["parameterization_contract"] == "target-only-pure-halo-support-exclusion-v1"
        assert saved["metadata"]["solver"]["search_method"] == "fixed-neighbor"
        default_solver = saved["metadata"]["solver"]
        assert {key: default_solver[key] for key in (
            "search_method", "fixed_neighbor_core_atoms",
            "fixed_neighbor_policy_version", "fixed_neighbor_maximum_sweeps",
            "fixed_neighbor_order", "fixed_neighbor_local_search",
        )} == {
            "search_method": "fixed-neighbor",
            "fixed_neighbor_core_atoms": 12,
            "fixed_neighbor_policy_version": "fixed-neighbor-production-v4",
            "fixed_neighbor_maximum_sweeps": 200,
            "fixed_neighbor_order": "forward",
            "fixed_neighbor_local_search": "active-set-lm",
        }
        fixed_solver = fixed["metadata"]["solver"]
        assert {key: fixed_solver[key] for key in (
            "search_method", "fixed_neighbor_core_atoms",
            "fixed_neighbor_policy_version",
            "fixed_neighbor_maximum_sweeps", "fixed_neighbor_order", "fixed_neighbor_local_search",
        )} == {
            "search_method": "fixed-neighbor",
            "fixed_neighbor_core_atoms": 12,
            "fixed_neighbor_policy_version": "fixed-neighbor-production-v4",
            "fixed_neighbor_maximum_sweeps": 200,
            "fixed_neighbor_order": "forward",
            "fixed_neighbor_local_search": "active-set-lm",
        }
        assert saved["selection_domain"]["target_indices"] == list(range(25))
        assert saved["selection_domain"]["contributor_policy"] == "pure-halo-support-exclusion-v1"
        assert saved["atom_ids"] == [str(index) for index in range(1, 26)]
        assert "fixed_halos" not in saved
        if saved["halo_median_reconstruction"] is not None:
            assert saved["halo_median_reconstruction"]["halos"] == []
        assert saved["initialization"]["data_scope"] == "contributor-local-sampling-may-read-outside-target-domain"
        csv_lines = (root / "joint_atoms_example.csv").read_text().splitlines()
        assert not (root / "joint_atoms_example.contributions.csv").exists()
        assert csv_lines[0].endswith(",SelectionRole,Parameterization,ContributionGroupID,TargetRuntimeConvergence,PointQuality")
        csv_rows = {row["SelectionRole"]: row for row in csv.DictReader(csv_lines)}
        assert csv_rows["target"]["Parameterization"] == "JointEstimated"
        assert csv_rows["target"]["StateAvailable"] == "1"
        assert csv_rows["target"]["TargetRuntimeConvergence"] == "passed"
        calibration = next(item for item in saved["evidence"] if item["name"] == "global-width-calibration")
        assert calibration["status"] == "passed"
        assert "global-first-stage-p95-k4" in calibration["reason"]
        assert abs(calibration["threshold"] - 4.0 * calibration["value"]) <= 1e-12
        assert 'target_runtime_convergence' in saved
        metadata = saved["metadata"]
        assert metadata["model_sha256"] == model_hash
        assert metadata["map_sha256"] == map_hash
        assert changed["metadata"]["model_sha256"] == changed_hash != model_hash
        assert metadata["map_normalization"]["requested"] is True
        assert metadata["map_normalization"]["applied"] is True
        assert metadata["map_normalization"]["divisor"] > 0
        assert metadata["units"]["A"] == "fit-map-unit*angstrom^3"
        assert metadata["units"]["C"] == "fit-map-unit*angstrom"
        for key in ("source_sha256", "configuration_sha256", "build_sha256"):
            assert len(metadata["software"][key]) == 64
        assert metadata["software"] == generator
        assert "regular_certificate" not in saved
        assert all("regular_certificate" not in component for component in saved["components"])
        assert saved["components"] and saved["assembled_state"] is not None
        for component in saved["components"]:
            assert component["search_completed"]
            assert component["state"] is not None
            assert component["stop_reason"] == "block-stationary"
            assert 0 < component["search_sweeps"] <= 200
            diagnostic = component["candidate_diagnostic"]
            assert diagnostic["convergence_criteria_id"] == "practical-v2"
            assert diagnostic["global_ac_kkt_threshold"] == 1e-7
            assert diagnostic["global_width_gradient_threshold"] == 1e-9
            assert "confirmed_log_width_change_threshold" not in diagnostic
            assert diagnostic["candidate_log_width_change_inf_norm"] >= 0
            point = component["target_evidence"]
            assert point["contract"] == "target-point-validation-v1"
            assert all(key not in point for key in ("ranks", "atoms", "column_scales"))
            assert point["point_atoms"]
        assert "prediction" not in saved and "observations" not in saved
        assert (root / "joint_atoms_example.csv").read_text().startswith("AtomID,ComponentID,A,B,C,")
        run("result_dump", "--printer", "gaus", "-d", database, "-k", "example", "-o", root)
        help_text = subprocess.run([executable, "--help"], capture_output=True, text=True, check=True).stdout
        if "umap_embedding" in help_text:
            rejected = subprocess.run([executable, "umap_embedding", "-d", str(database), "--model-key", "example"],
                                      capture_output=True, text=True)
            assert rejected.returncode != 0
            assert "All 3 selected UMAP feature columns are constant" in rejected.stdout + rejected.stderr, (
                rejected.stdout + rejected.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
