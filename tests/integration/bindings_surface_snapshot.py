import rhbm_gem_module as m


EXPECTED_REQUEST_TYPES = {
    "CommonCommandRequest",
    "PotentialAnalysisRequest",
    "PotentialDisplayRequest",
    "ResultDumpRequest",
    "MapSimulationRequest",
    "MapVisualizationRequest",
    "PositionEstimationRequest",
    "HRLModelTestRequest",
}

EXPECTED_RUN_FUNCTIONS = {
    "RunPotentialAnalysis",
    "RunPotentialDisplay",
    "RunResultDump",
    "RunMapSimulation",
    "RunMapVisualization",
    "RunPositionEstimation",
    "RunHRLModelTest",
}


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
