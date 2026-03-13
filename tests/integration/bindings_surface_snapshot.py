import rhbm_gem_module as m

from pathlib import Path
import re


def _manifest_entries() -> list[str]:
    manifest_path = Path(__file__).resolve().parents[2] / "include" / "rhbm_gem" / "core" / "command" / "CommandList.def"
    pattern = re.compile(
        r"RHBM_GEM_COMMAND\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*,"
    )
    return [match.group(1) for match in pattern.finditer(manifest_path.read_text(encoding="utf-8"))]


EXPECTED_REQUEST_TYPES = {"CommonCommandRequest"} | {f"{stem}Request" for stem in _manifest_entries()}
EXPECTED_RUN_FUNCTIONS = {f"Run{stem}" for stem in _manifest_entries()}


def main() -> int:
    for name in EXPECTED_REQUEST_TYPES:
        assert hasattr(m, name), f"Missing request type: {name}"

    for name in EXPECTED_RUN_FUNCTIONS:
        assert hasattr(m, name), f"Missing run function: {name}"

    assert hasattr(m, "ExecutionReport")
    report = m.ExecutionReport()
    for field in ("prepared", "executed", "validation_issues"):
        assert hasattr(report, field), f"ExecutionReport missing field: {field}"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
