from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    project_root = Path(__file__).resolve().parents[2]
    script = project_root / "examples" / "python" / "01_end_to_end_from_test_data.py"

    with tempfile.TemporaryDirectory(prefix="rhbm_examples_smoke_") as temp_dir:
        workdir = Path(temp_dir) / "pipeline_output"
        result = subprocess.run(
            [
                sys.executable,
                str(script),
                "--sampling-size",
                "10",
                "--workdir",
                str(workdir),
            ],
            cwd=project_root,
            capture_output=True,
            text=True,
            check=False,
        )

        output = f"{result.stdout}\n{result.stderr}"
        assert result.returncode == 0, output
        assert "Pipeline completed successfully." in result.stdout, output

        maps_dir = workdir / "maps"
        dump_dir = workdir / "dump"
        db_file = workdir / "demo.sqlite"
        map_outputs = sorted(maps_dir.glob("*.map"))
        dump_outputs = sorted(dump_dir.glob("*"))

        assert db_file.exists(), output
        assert map_outputs, output
        assert dump_outputs, output
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
