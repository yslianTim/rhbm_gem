from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def run_script(script_name: str, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "examples" / "python" / script_name), *args],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def assert_quickstart_smoke() -> None:
    result = run_script("00_quickstart.py")
    output = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == 0, output
    assert "Imported module: rhbm_gem_module" in result.stdout, output
    assert "Available run functions:" in result.stdout, output


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
