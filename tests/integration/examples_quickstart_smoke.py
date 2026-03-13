from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    project_root = Path(__file__).resolve().parents[2]
    script = project_root / "examples" / "python" / "00_quickstart.py"
    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=project_root,
        capture_output=True,
        text=True,
        check=False,
    )

    output = f"{result.stdout}\n{result.stderr}"
    assert result.returncode == 0, output
    assert "Imported module: rhbm_gem_module" in result.stdout, output
    assert "Available run functions:" in result.stdout, output
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
