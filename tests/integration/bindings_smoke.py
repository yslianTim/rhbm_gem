import rhbm_gem_module as m

from pathlib import Path
import re


def _manifest_entries() -> list[str]:
    manifest_path = Path(__file__).resolve().parents[2] / "include" / "rhbm_gem" / "core" / "command" / "CommandList.def"
    pattern = re.compile(
        r"RHBM_GEM_COMMAND\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*,"
    )
    return [match.group(1) for match in pattern.finditer(manifest_path.read_text(encoding="utf-8"))]


REQUEST_NAMES = ["CommonCommandRequest"] + [f"{stem}Request" for stem in _manifest_entries()]
RUN_FUNCTIONS = [f"Run{stem}" for stem in _manifest_entries()]


def main() -> int:
    assert len(REQUEST_NAMES) >= 2
    assert len(RUN_FUNCTIONS) >= 1
    assert hasattr(m, "LogLevel")
    assert hasattr(m, "ValidationPhase")
    assert hasattr(m, "ValidationIssue")
    assert hasattr(m, "ExecutionReport")
    assert hasattr(m, "PrinterType")
    assert hasattr(m.PrinterType, "ATOM_OUTLIER")
    assert hasattr(m, "TesterType")
    assert hasattr(m.TesterType, "BENCHMARK")
    assert all(hasattr(m, name) for name in REQUEST_NAMES)
    assert all(callable(getattr(m, name)) for name in RUN_FUNCTIONS)

    req = m.MapSimulationRequest()
    assert hasattr(req, "common")
    assert hasattr(req, "model_file_path")
    assert hasattr(req, "blurring_width_list")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
