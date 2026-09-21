"""Exercise the real CLI's joint analysis, persistence and export contract."""
from __future__ import annotations

import hashlib
import shutil
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
        shutil.copyfile(model, local_model)
        local_model.write_text(local_model.read_text().rsplit("#", 1)[0] + "ATOM 2 C CB . ALA A 1 1.2 0.0 0.0 1.0 0.0 1\n#\n")
        model = local_model

        def run(*args: object, succeeds: bool = True) -> None:
            result = subprocess.run([executable, *map(str, args)], capture_output=True, text=True)
            assert (result.returncode == 0) == succeeds, result.stdout + result.stderr

        run("map_simulation", "-a", model, "-o", root, "--potential-model", "single",
            "--blurring-width", ".5", "-g", ".3", "-v", "0")
        map_path = next(root.glob("*.map"))
        generator = json.loads(Path(str(map_path) + ".simulation.json").read_text())["generator"]
        database = root / "joint.sqlite"
        run("potential_analysis", "--estimator", "joint-components", "--only-backbone", "true", "-a", model,
            "-m", map_path, "-d", database, "-k", "example", "--map-normalization", "false", "-v", "0")
        model_hash = hashlib.sha256(model.read_bytes()).hexdigest()
        map_hash = hashlib.sha256(map_path.read_bytes()).hexdigest()
        # The same path with different bytes must produce a different fingerprint.
        model.write_text(model.read_text() + "\n# provenance test\n")
        run("potential_analysis", "--estimator", "joint-components", "--only-backbone", "true", "-a", model,
            "-m", map_path, "-d", database, "-k", "changed", "--map-normalization", "false", "-v", "0")
        changed_hash = hashlib.sha256(model.read_bytes()).hexdigest()
        map_path.unlink()
        model.unlink()
        run("result_dump", "--printer", "joint", "-d", database, "-k", "example", "-o", root)
        saved = json.loads((root / "joint_result_example.json").read_text())
        with sqlite3.connect(database) as connection:
            assert connection.execute("PRAGMA user_version").fetchone()[0] == 17
            payload = connection.execute("SELECT result_json FROM model_joint_result WHERE key_tag='example'").fetchone()[0]
            changed = json.loads(connection.execute("SELECT result_json FROM model_joint_result WHERE key_tag='changed'").fetchone()[0])
        assert saved == json.loads(payload)
        assert saved["schema_version"] == 3
        assert saved["selection_domain"]["target_indices"] == [0]
        assert saved["atom_ids"] == ["1", "2"]
        assert saved["initialization"]["data_scope"] == "contributor-local-sampling-may-read-outside-target-domain"
        csv = (root / "joint_atoms_example.csv").read_text().splitlines()
        assert csv[0].endswith(",SelectionRole")
        assert csv[1].endswith(",target") and csv[2].endswith(",halo")
        metadata = saved["metadata"]
        assert metadata["model_sha256"] == model_hash
        assert metadata["map_sha256"] == map_hash
        assert changed["metadata"]["model_sha256"] == changed_hash != model_hash
        assert metadata["map_normalization"] == {"requested": False, "applied": False, "divisor": 1}
        assert metadata["units"]["A"] == "fit-map-unit*angstrom^3"
        assert metadata["units"]["C"] == "fit-map-unit*angstrom"
        for key in ("source_sha256", "configuration_sha256", "build_sha256"):
            assert len(metadata["software"][key]) == 64
        assert metadata["software"] == generator
        assert saved["regular_certificate"] == "not-run"
        assert saved["components"] and saved["assembled_state"] is not None
        assert "prediction" not in saved and "observations" not in saved
        assert (root / "joint_atoms_example.csv").read_text().startswith("AtomID,ComponentID,A,B,C,")
        run("result_dump", "--printer", "gaus", "-d", database, "-k", "example", succeeds=False)
        help_text = subprocess.run([executable, "--help"], capture_output=True, text=True, check=True).stdout
        if "umap_embedding" in help_text:
            rejected = subprocess.run([executable, "umap_embedding", "-d", str(database), "--model-key", "example"],
                                      capture_output=True, text=True)
            assert rejected.returncode != 0
            assert "Joint result UMAP is not supported" in rejected.stdout + rejected.stderr
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
