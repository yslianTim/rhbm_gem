import rhbm_gem_module as m

from command_manifest_expectations import expected_request_names, expected_run_functions

EXPECTED_REQUEST_TYPES = set(expected_request_names())
EXPECTED_RUN_FUNCTIONS = set(expected_run_functions())


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
