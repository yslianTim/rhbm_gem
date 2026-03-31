from __future__ import annotations

import subprocess
import shutil
import sys
import tempfile
from pathlib import Path

import rhbm_gem_module as rgm


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def run_script(script_name: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "resources" / "examples" / "python" / script_name), *args],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def ensure_execute(result: rgm.CommandResult, step_name: str) -> None:
    details = "\n".join(
        f"[{issue.phase}/{issue.level}] {issue.option_name}: {issue.message}"
        for issue in result.issues
    )
    if result.outcome == rgm.CommandOutcome.Succeeded:
        return
    raise AssertionError(
        f"{step_name} failed (outcome={result.outcome}).\n{details}"
    )


def stage_quickstart_inputs(workdir: Path) -> None:
    data_dir = workdir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    model_path = data_dir / "6Z6U.cif"
    shutil.copyfile(PROJECT_ROOT / "tests" / "fixtures" / "test_model.cif", model_path)

    simulation_request = rgm.MapSimulationRequest()
    simulation_request.model_file_path = str(model_path)
    simulation_request.output_dir = str(data_dir)
    simulation_request.map_file_name = "emd_11103_additional"
    simulation_request.blurring_width_list = [1.50]
    ensure_execute(rgm.RunMapSimulation(simulation_request), "RunMapSimulation")

    generated_maps = sorted(data_dir.glob("emd_11103_additional_*.map"))
    assert len(generated_maps) == 1, generated_maps
    generated_maps[0].replace(data_dir / "emd_11103_additional.map")


def assert_quickstart_smoke() -> None:
    with tempfile.TemporaryDirectory(prefix="rhbm_quickstart_smoke_") as temp_dir:
        workdir = Path(temp_dir) / "quickstart_output"
        stage_quickstart_inputs(workdir)

        result = run_script("00_quickstart.py", "--workdir", str(workdir), "--skip-download")
        output = f"{result.stdout}\n{result.stderr}"

        assert result.returncode == 0, output
        assert "Potential analysis" in result.stdout, output
        assert "Result dump" in result.stdout, output
        assert "Quickstart completed successfully." in result.stdout, output

        db_file = workdir / "data" / "database.sqlite"
        dump_outputs = sorted(path for path in workdir.iterdir() if path.is_file())

        assert db_file.exists(), output
        assert dump_outputs, output


def assert_end_to_end_smoke() -> None:
    with tempfile.TemporaryDirectory(prefix="rhbm_examples_smoke_") as temp_dir:
        workdir = Path(temp_dir) / "pipeline_output"
        data_dir = workdir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        staged_examples = (
            ("6Z6U.cif", "emd_11103_additional.map"),
            ("6DRV.cif", "emd_8908_additional.map"),
            ("9GXM.cif", "emd_51667_additional_1.map"),
        )

        for model_name, map_name in staged_examples:
            model_path = data_dir / model_name
            shutil.copyfile(PROJECT_ROOT / "tests" / "fixtures" / "test_model.cif", model_path)

            simulation_request = rgm.MapSimulationRequest()
            simulation_request.model_file_path = str(model_path)
            simulation_request.output_dir = str(data_dir)
            simulation_request.map_file_name = Path(map_name).stem
            simulation_request.blurring_width_list = [1.50]
            ensure_execute(rgm.RunMapSimulation(simulation_request), "RunMapSimulation")

            generated_maps = sorted(data_dir.glob(f"{Path(map_name).stem}_*.map"))
            assert len(generated_maps) == 1, generated_maps
            generated_maps[0].replace(data_dir / map_name)

        result = run_script("01_end_to_end_from_three_examples.py", "--workdir", str(workdir), "--skip-download")

        output = f"{result.stdout}\n{result.stderr}"
        assert result.returncode == 0, output
        assert "Three-example workflow completed successfully." in result.stdout, output

        db_file = workdir / "data" / "database.sqlite"
        dump_outputs = sorted(path for path in workdir.iterdir() if path.is_file())

        assert db_file.exists(), output
        assert dump_outputs, output


def main() -> int:
    assert_quickstart_smoke()
    assert_end_to_end_smoke()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
