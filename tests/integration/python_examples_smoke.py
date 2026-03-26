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


def ensure_execute(report: rgm.ExecutionReport, step_name: str) -> None:
    details = "\n".join(
        f"[{issue.phase}/{issue.level}] {issue.option_name}: {issue.message}"
        for issue in report.validation_issues
    )
    if report.prepared and report.executed:
        return
    raise AssertionError(
        f"{step_name} failed (prepared={report.prepared}, executed={report.executed}).\n{details}"
    )


def stage_quickstart_inputs(workdir: Path) -> None:
    data_dir = workdir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    model_path = data_dir / "6Z6U.cif"
    shutil.copyfile(PROJECT_ROOT / "tests" / "fixtures" / "test_model.cif", model_path)

    simulation_request = rgm.MapSimulationRequest()
    simulation_request.model_file_path = str(model_path)
    simulation_request.common.folder_path = str(data_dir)
    simulation_request.map_file_name = "emd_11103_additional"
    simulation_request.blurring_width_list = "1.50"
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
        result = run_script(
            "01_end_to_end_from_test_data.py",
            "--sampling-size",
            "10",
            "--workdir",
            str(workdir),
        )

        output = f"{result.stdout}\n{result.stderr}"
        assert result.returncode == 0, output
        assert "Pipeline completed successfully." in result.stdout, output

        maps_dir = workdir / "maps"
        dump_dir = workdir / "dump"
        db_file = workdir / "demo.sqlite"

        assert db_file.exists(), output
        assert sorted(maps_dir.glob("*.map")), output
        assert sorted(dump_dir.glob("*")), output


def main() -> int:
    assert_quickstart_smoke()
    assert_end_to_end_smoke()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
