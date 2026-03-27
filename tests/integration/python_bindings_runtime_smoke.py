from __future__ import annotations

import os
from pathlib import Path

import rhbm_gem_module as m


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENTAL_FEATURE_ENABLED = (
    os.environ.get("RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE", "OFF").upper() == "ON"
)

EXPECTED_COMMON_FIELDS = {
    "thread_size",
    "verbose_level",
    "database_path",
    "folder_path",
}


def assert_module_surface() -> None:
    assert hasattr(m, "LogLevel")
    assert hasattr(m, "ValidationPhase")
    assert hasattr(m, "ValidationIssue")
    assert hasattr(m, "ExecutionReport")
    assert hasattr(m, "PrinterType")
    assert hasattr(m.PrinterType, "ATOM_OUTLIER")
    assert hasattr(m, "TesterType")
    assert hasattr(m.TesterType, "BENCHMARK")

    if EXPERIMENTAL_FEATURE_ENABLED:
        assert hasattr(m, "MapVisualizationRequest")
        assert hasattr(m, "PositionEstimationRequest")
        assert hasattr(m, "RunMapVisualization")
        assert hasattr(m, "RunPositionEstimation")
    else:
        assert not hasattr(m, "MapVisualizationRequest")
        assert not hasattr(m, "PositionEstimationRequest")
        assert not hasattr(m, "RunMapVisualization")
        assert not hasattr(m, "RunPositionEstimation")


def assert_request_objects_are_usable() -> None:
    request_types = [
        m.PotentialAnalysisRequest,
        m.PotentialDisplayRequest,
        m.ResultDumpRequest,
        m.MapSimulationRequest,
        m.HRLModelTestRequest,
    ]
    if EXPERIMENTAL_FEATURE_ENABLED:
        request_types.extend(
            [
                m.MapVisualizationRequest,
                m.PositionEstimationRequest,
            ]
        )

    for request_type in request_types:
        request = request_type()
        common = request.common
        missing = [field for field in EXPECTED_COMMON_FIELDS if not hasattr(common, field)]
        assert not missing, f"{request_type.__name__}.common missing fields: {missing}"

    simulation = m.MapSimulationRequest()
    simulation.model_file_path = str(PROJECT_ROOT / "tests" / "fixtures" / "test_model.cif")
    simulation.common.folder_path = "runtime_smoke_output"
    simulation.blurring_width_list = "1.50"

    assert Path(simulation.model_file_path).name == "test_model.cif"
    assert Path(simulation.common.folder_path) == Path("runtime_smoke_output")
    assert simulation.blurring_width_list == "1.50"


def has_issue(report, option_name: str, phase) -> bool:
    return any(
        issue.option_name == option_name and issue.phase == phase
        for issue in report.validation_issues
    )


def assert_execution_report_runtime_behavior() -> None:
    report = m.RunMapSimulation(m.MapSimulationRequest())
    assert isinstance(report, m.ExecutionReport)
    assert not report.prepared
    assert not report.executed
    assert has_issue(report, "--model", m.ValidationPhase.Parse)

    analysis = m.PotentialAnalysisRequest()
    analysis.saved_key_tag = ""
    report = m.RunPotentialAnalysis(analysis)
    assert isinstance(report, m.ExecutionReport)
    assert not report.prepared
    assert has_issue(report, "--save-key", m.ValidationPhase.Parse)


def main() -> int:
    assert_module_surface()
    assert_request_objects_are_usable()
    assert_execution_report_runtime_behavior()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
