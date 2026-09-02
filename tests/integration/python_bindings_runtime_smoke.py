from __future__ import annotations

import os
import tempfile
from pathlib import Path

import rhbm_gem_module as m


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENTAL_FEATURE_ENABLED = (
    os.environ.get("RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE", "OFF").upper() == "ON"
)
UMAP_ENABLED = os.environ.get("RHBM_GEM_ENABLE_UMAP", "ON").upper() == "ON"

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
    assert hasattr(m, "CommandDiagnostic")
    assert hasattr(m, "CommandResult")
    assert hasattr(m, "RunCommand")
    assert hasattr(m, "PainterType")
    assert not hasattr(m.PainterType, "QSCORE")
    assert hasattr(m, "PrinterType")
    assert hasattr(m.PrinterType, "ATOM_OUTLIER")
    assert hasattr(m, "TesterType")
    assert hasattr(m.TesterType, "BENCHMARK")
    assert hasattr(m.TesterType, "ATOMIC_MODEL")
    assert hasattr(m, "SphereSamplingMethod")
    assert hasattr(m.SphereSamplingMethod, "FIBONACCI_DETERMINISTIC")

    if EXPERIMENTAL_FEATURE_ENABLED:
        assert hasattr(m, "MapVisualizationRequest")
        assert hasattr(m, "PositionEstimationRequest")
    else:
        assert not hasattr(m, "MapVisualizationRequest")
        assert not hasattr(m, "PositionEstimationRequest")

    if UMAP_ENABLED:
        assert hasattr(m, "UmapEmbeddingRequest")
    else:
        assert not hasattr(m, "UmapEmbeddingRequest")

    for old_run_name in (
        "RunPotentialAnalysis",
        "RunPotentialDisplay",
        "RunResultDump",
        "RunMapSimulation",
        "RunRHBMTest",
        "RunMapVisualization",
        "RunPositionEstimation",
    ):
        assert not hasattr(m, old_run_name)


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
    if UMAP_ENABLED:
        request_types.append(m.UmapEmbeddingRequest)

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
    assert simulation.exclude_hydrogen is False
    assert simulation.only_backbone is False
    simulation.exclude_hydrogen = True
    simulation.only_backbone = True
    assert simulation.exclude_hydrogen is True
    assert simulation.only_backbone is True

    analysis = m.PotentialAnalysisRequest()
    assert analysis.map_normalization_flag is True
    assert analysis.exclude_hydrogen is False
    assert analysis.only_backbone is False
    analysis.map_normalization_flag = False
    analysis.exclude_hydrogen = True
    analysis.only_backbone = True
    assert analysis.map_normalization_flag is False
    assert analysis.exclude_hydrogen is True
    assert analysis.only_backbone is True
    for field_name in (
        "map_normalization_flag",
        "exclude_hydrogen",
        "only_backbone",
        "sampling_method",
    ):
        assert hasattr(analysis, field_name), f"PotentialAnalysisRequest missing {field_name}"
    for removed_field_name in (
        "sampling_profile_choice",
        "sampling_size",
        "sampling_range_min",
        "sampling_range_max",
        "sampling_height",
        "training_alpha_min",
        "training_alpha_max",
        "training_alpha_step",
        "training_report_dir",
        "training_alpha_flag",
        "alpha_r",
        "alpha_g",
    ):
        assert not hasattr(analysis, removed_field_name), (
            f"PotentialAnalysisRequest still exposes {removed_field_name}"
        )
    assert (
        analysis.sampling_method
        == m.SphereSamplingMethod.FIBONACCI_DETERMINISTIC
    )
    analysis.sampling_method = m.SphereSamplingMethod.VOLUME_UNIFORM_RANDOM
    assert (
        analysis.sampling_method
        == m.SphereSamplingMethod.VOLUME_UNIFORM_RANDOM
    )

    display = m.PotentialDisplayRequest()
    assert not hasattr(display, "map_file_path")
    display.painter_choice = m.PainterType.GAUS
    assert display.painter_choice == m.PainterType.GAUS

    if UMAP_ENABLED:
        umap = m.UmapEmbeddingRequest()
        assert umap.num_neighbors == 15
        assert umap.min_dist == 0.1
        assert umap.num_epochs == 0
        assert umap.random_seed == 42


def has_issue(report, option_name: str) -> bool:
    return any(
        issue.option_name == option_name
        for issue in report.issues
    )


def assert_command_result_runtime_behavior() -> None:
    report = m.RunCommand(m.MapSimulationRequest())
    assert isinstance(report, m.CommandResult)
    assert not report.succeeded
    assert has_issue(report, "-a,--model")

    analysis = m.PotentialAnalysisRequest()
    analysis.saved_key_tag = ""
    report = m.RunCommand(analysis)
    assert isinstance(report, m.CommandResult)
    assert not report.succeeded
    assert has_issue(report, "-k,--save-key")


def assert_umap_runtime_behavior() -> None:
    if not UMAP_ENABLED:
        return

    header = (
        "serial id,residue,spot,neighbor count,"
        "signal peeling ratio,tail peeling ratio,"
        "amplitude 1st,amplitude 2nd,amplitude 3rd,"
        "width 1st,width 2nd,width 3rd,"
        "offset 1st,offset 2nd,offset 3rd,"
        "amplitude rank 1st,amplitude rank 2nd,amplitude rank 3rd,"
        "width rank 1st,width rank 2nd,width rank 3rd,"
        "offset rank 1st,offset rank 2nd,offset rank 3rd"
    )
    with tempfile.TemporaryDirectory(prefix="rhbm_umap_python_") as temp_dir:
        workdir = Path(temp_dir)
        input_path = workdir / "local_fitting_result_python.csv"
        rows = [header]
        for observation in range(6):
            features = [
                (observation + 1) * (feature + 2)
                + (observation * observation + 3 * feature) % (feature + 3)
                for feature in range(21)
            ]
            rows.append(
                ",".join(
                    [str(observation + 1), "ALA", "CA"]
                    + [str(value) for value in features]
                )
            )
        input_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

        request = m.UmapEmbeddingRequest()
        request.input_csv_path = input_path
        request.output_dir = workdir / "output"
        request.num_neighbors = 3
        request.min_dist = 0.2
        request.num_epochs = 5
        request.random_seed = 99
        request.job_count = 1

        report = m.RunCommand(request)
        assert report.succeeded
        output_path = Path(request.output_dir) / "umap_embedding_python.csv"
        output_rows = output_path.read_text(encoding="utf-8").splitlines()
        assert len(output_rows) == 7
        assert all(len(row.split(",")) == 26 for row in output_rows)


def main() -> int:
    assert_module_surface()
    assert_request_objects_are_usable()
    assert_command_result_runtime_behavior()
    assert_umap_runtime_behavior()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
