import rhbm_gem_module as m

REQUEST_NAMES = [
    "CommonCommandRequest",
    "PotentialAnalysisRequest",
    "PotentialDisplayRequest",
    "ResultDumpRequest",
    "MapSimulationRequest",
    "MapVisualizationRequest",
    "PositionEstimationRequest",
    "HRLModelTestRequest",
]

RUN_FUNCTIONS = [
    "RunPotentialAnalysis",
    "RunPotentialDisplay",
    "RunResultDump",
    "RunMapSimulation",
    "RunMapVisualization",
    "RunPositionEstimation",
    "RunHRLModelTest",
]


def main() -> int:
    assert len(REQUEST_NAMES) == 8
    assert len(RUN_FUNCTIONS) == 7
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
