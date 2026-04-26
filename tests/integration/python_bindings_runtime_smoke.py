from __future__ import annotations

import os
from pathlib import Path

import rhbm_gem_module as m


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENTAL_FEATURE_ENABLED = (
    os.environ.get("RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE", "OFF").upper() == "ON"
)

EXPECTED_COMMON_FIELDS = {
    "job_count",
    "verbosity",
    "output_dir",
}

DATABASE_REQUEST_TYPES = {
    "PotentialAnalysisRequest",
    "PotentialDisplayRequest",
    "ResultDumpRequest",
}


def assert_module_surface() -> None:
    assert hasattr(m, "ValidationIssue")
    assert hasattr(m, "CommandResult")
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
        m.RHBMTestRequest,
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
        missing = [field for field in EXPECTED_COMMON_FIELDS if not hasattr(request, field)]
        assert not missing, f"{request_type.__name__} missing fields: {missing}"
        has_database_path = hasattr(request, "database_path")
        assert has_database_path == (request_type.__name__ in DATABASE_REQUEST_TYPES), (
            f"{request_type.__name__} database_path presence mismatch"
        )

    simulation = m.MapSimulationRequest()
    simulation.model_file_path = str(PROJECT_ROOT / "tests" / "fixtures" / "test_model.cif")
    simulation.output_dir = "runtime_smoke_output"
    simulation.blurring_width_list = [1.50]

    assert Path(simulation.model_file_path).name == "test_model.cif"
    assert Path(simulation.output_dir) == Path("runtime_smoke_output")
    assert simulation.blurring_width_list == [1.50]

    analysis = m.PotentialAnalysisRequest()
    for field_name in (
        "training_alpha_min",
        "training_alpha_max",
        "training_alpha_step",
    ):
        assert hasattr(analysis, field_name), f"PotentialAnalysisRequest missing {field_name}"
    analysis.training_alpha_min = 0.2
    analysis.training_alpha_max = 0.8
    analysis.training_alpha_step = 0.2
    assert analysis.training_alpha_min == 0.2
    assert analysis.training_alpha_max == 0.8
    assert analysis.training_alpha_step == 0.2


def has_issue(report, option_name: str) -> bool:
    return any(
        issue.option_name == option_name
        for issue in report.issues
    )


def assert_command_result_runtime_behavior() -> None:
    report = m.RunMapSimulation(m.MapSimulationRequest())
    assert isinstance(report, m.CommandResult)
    assert not report.succeeded
    assert has_issue(report, "--model")

    analysis = m.PotentialAnalysisRequest()
    analysis.saved_key_tag = ""
    report = m.RunPotentialAnalysis(analysis)
    assert isinstance(report, m.CommandResult)
    assert not report.succeeded
    assert has_issue(report, "--save-key")


def main() -> int:
    assert_module_surface()
    assert_request_objects_are_usable()
    assert_command_result_runtime_behavior()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
