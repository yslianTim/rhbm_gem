"""Check the current Joint experiment surface stays within its cleanup contract."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CURRENT_EXPERIMENT_SOURCES = {
    "joint_fixed_neighbor.cpp",
    "joint_offline_diagnostic.cpp",
    "joint_partial_selection.cpp",
    "joint_postprocessing_benchmark.cpp",
    "joint_statistical_experiment.cpp",
}

CURRENT_PYTHON_EXPERIMENTS = {
    "joint_benchmark.py",
    "joint_offline_diagnostic.py",
    "joint_partial_selection.py",
    "joint_statistical_experiment.py",
}

CURRENT_PYTHON_VALIDATION_RUNNERS = {
    "joint_component_audit.py",
    "joint_component_runtime.py",
    "joint_workflow_cli_smoke.py",
}

RETIRED_REPOSITORY_FILES = (
    "tests/experiments/joint_sparse_benchmark.cpp",
    "tests/core/joint_component/FixedBBlockCoordinate_test.cpp",
    "tests/core/joint_component/OperatorSearch_test.cpp",
    "tests/core/joint_component/ProfileOperator_test.cpp",
    "tests/core/joint_component/ProjectedTailQr_test.cpp",
    "src/core/detail/joint_component/SparseFactor.hpp",
    "src/core/detail/joint_component/SparseFactor.cpp",
)

RETIRED_PYTHON_CAMPAIGNS = (
    "tests/integration/joint_fixed_b_scaling.py",
    "tests/integration/joint_fixed_b_scaling_test.py",
    "tests/integration/joint_fixed_neighbor_outer_core_qualification.py",
    "tests/integration/joint_fixed_neighbor_outer_core_support.py",
    "tests/integration/joint_fixed_neighbor_outer_core_qualification_test.py",
    "tests/integration/joint_fixed_neighbor_order.py",
    "tests/integration/joint_fixed_neighbor_order_test.py",
    "tests/integration/joint_fixed_neighbor_stationarity.py",
    "tests/integration/joint_fixed_neighbor_stationarity_test.py",
    "tests/integration/joint_fixed_neighbor_prepared_block.py",
    "tests/integration/joint_fixed_neighbor_prepared_block_test.py",
)

ACTIVE_BACKEND_TERMS = (
    "SPQR",
    "SuiteSparse",
    "CHOLMOD",
    "RHBM_GEM_JOINT_SPQR",
    "RHBM_GEM_JOINT_SPARSE_BACKEND",
    "SparseBackend",
    "ActiveSparseBackend",
    "SparseBackendEnabled",
    "SPQR_ORDERING",
    "SpqrOrdering",
    "SparseFactor",
)

RETIRED_DERIVATIVE_CONTROLS = (
    "JacobianReductionKindForTesting",
    "JacobianReductionForTesting",
    "PrepareDerivativeForTesting",
    "ReduceDerivativeForTesting",
    "DerivativeTileRowsForTesting",
)

RETIRED_FIXED_B_SOLVER = (
    "FixedBBlockPolicy",
    "FixedBBlockOrder",
    "FixedBBlockResult",
    "SearchFixedBBlocks",
    "FixedBBlockRecord",
    "FixedBBlockSweep",
)

RETIRED_FIXTURES = (
    "joint_fixed_neighbor_inexact_baseline.json",
)

def main():
    integration = ROOT / "tests" / "integration"
    cmake = (ROOT / "tests" / "CMakeLists.txt").read_text()
    fixed_neighbor_source = (ROOT / "tests/experiments/joint_fixed_neighbor.cpp").read_text()
    benchmark_driver = (integration / "joint_benchmark.py").read_text()
    runtime_runner = (integration / "joint_component_runtime.py").read_text()
    audit_source = (ROOT / "tests/tools/joint_component_audit.cpp").read_text()
    audit_runner = (integration / "joint_component_audit.py").read_text()
    gitignore = (ROOT / ".gitignore").read_text()
    active_sources = "\n".join(path.read_text() for path in integration.glob("*.py"))
    active_paths = [path for path in (ROOT / "src/core").rglob("*") if path.is_file()]
    active_paths.extend(path for path in (ROOT / "include/rhbm_gem/core").rglob("*") if path.is_file())
    active_paths.extend(path for path in (ROOT / "tests" / "experiments").rglob("*")
                        if path.is_file())
    active_paths.extend(path for path in (ROOT / "tests" / "core" / "joint_component").rglob("*")
                        if path.is_file())
    active_paths.extend(path for path in integration.rglob("*")
                        if path.is_file() and path.name != Path(__file__).name)
    active_code = "\n".join(
        path.read_text(errors="ignore")
        for path in active_paths
        if path.suffix in {".c", ".cc", ".cpp", ".h", ".hh", ".hpp", ".hxx", ".py", ".sh"}
    )

    experiment_dir = ROOT / "tests" / "experiments"
    experiment_sources = {path.name for path in experiment_dir.iterdir()
                          if path.is_file() and path.suffix in {".c", ".cc", ".cpp"}}
    python_experiments = {path.name for path in integration.glob("joint_*.py")
                          if "__main__" in path.read_text(errors="ignore")
                          and not path.name.endswith("_test.py")
                          and path.name not in CURRENT_PYTHON_VALIDATION_RUNNERS}
    assert experiment_sources == CURRENT_EXPERIMENT_SOURCES, (
        f"tests/experiments executable sources must match the current whitelist: "
        f"unexpected={sorted(experiment_sources-CURRENT_EXPERIMENT_SOURCES)}, "
        f"missing={sorted(CURRENT_EXPERIMENT_SOURCES-experiment_sources)}")
    assert python_experiments == CURRENT_PYTHON_EXPERIMENTS, (
        f"current Joint Python experiment drivers must match the whitelist: "
        f"unexpected={sorted(python_experiments-CURRENT_PYTHON_EXPERIMENTS)}, "
        f"missing={sorted(CURRENT_PYTHON_EXPERIMENTS-python_experiments)}")
    missing_fixtures = [name for name in RETIRED_FIXTURES
                        if (ROOT / "tests" / "fixtures" / name).exists()]
    present_files = [name for name in RETIRED_REPOSITORY_FILES
                     if (ROOT / name).exists()]
    present_campaigns = [name for name in RETIRED_PYTHON_CAMPAIGNS
                         if (ROOT / name).exists()]
    test_paths = [path for path in (ROOT / "tests").rglob("*") if path.is_file()]
    retired_test_name_parts = (
        "spqrbenchmark", "benchmarkspqr", "spqrorder", "orderspqr",
        "sparsebackend", "backendselector",
    )
    present_retired_tests = [
        str(path.relative_to(ROOT))
        for path in test_paths
        if any(part in path.stem.lower().replace("_", "").replace("-", "")
               for part in retired_test_name_parts)
    ]
    assert not missing_fixtures, f"retired Joint fixtures returned: {missing_fixtures}"
    assert not present_files, f"retired Operator files returned: {present_files}"
    assert not present_campaigns, f"retired qualification drivers/tests returned: {present_campaigns}"
    assert not present_retired_tests, f"retired sparse tests returned: {present_retired_tests}"
    figure_root = ROOT / "docs/developer/figures"
    figure_dirs = {path.name for path in figure_root.iterdir() if path.is_dir()}
    unexpected_figure_dirs = sorted(figure_dirs - {"joint-validation"})
    assert figure_dirs == {"joint-validation"}, (
        f"closed Joint evidence directories returned: {unexpected_figure_dirs}")
    assert "joint_fixed_b_block_experiment" not in cmake
    assert "joint_fixed_b_scaling" not in cmake
    assert "joint_component_benchmark" not in cmake
    assert "commands.add_parser('benchmark')" not in runtime_runner
    for token in (
        "joint_fixed_b_observer_scalability",
        "joint_fixed_b_256_scalability",
        "joint_fixed_neighbor_outer_core_qualification_test",
        "joint_fixed_neighbor_stationarity_test",
        "joint_fixed_neighbor_order_test",
        "joint_fixed_neighbor_prepared_block_test",
        "joint_component_two_step_regression",
    ):
        assert token not in cmake, f"retired Joint CTest returned: {token}"
    for token in (
        "--neighbor-forward", "--neighbor-reverse",
        "--scaling-forward", "--scaling-reverse", "OUTER_CORE_ATOMS",
    ):
        assert token not in fixed_neighbor_source, f"retired benchmark mode returned: {token}"
    for token in ("--fixed-core-atoms", "fixed_neighbor_policy = custom"):
        assert token not in benchmark_driver, f"custom benchmark policy returned: {token}"
    for token in ("two-step-fixture", "--two-step-only"):
        assert token not in audit_source and token not in audit_runner, (
            f"historical two-step audit mode returned: {token}")
    for token in (
        "FixedNeighborLocalWork",
        "TwoAcceptedUpdates",
        "FullLocalSearch",
        "--fixed-local-work",
        "--inexact-one",
        "--inexact-two",
    ):
        assert token not in active_code, f"retired FixedNeighbor token remains: {token}"
    for token in ACTIVE_BACKEND_TERMS:
        assert token.lower() not in active_code.lower(), f"retired sparse backend token remains: {token}"
        assert token.lower() not in cmake.lower(), f"retired sparse backend token remains in tests CMake: {token}"
    src_cmake = (ROOT / "src" / "CMakeLists.txt").read_text()
    for token in ACTIVE_BACKEND_TERMS:
        assert token.lower() not in src_cmake.lower(), f"retired sparse backend token remains in src CMake: {token}"
    for token in RETIRED_DERIVATIVE_CONTROLS:
        assert token not in active_code, f"retired derivative control remains: {token}"
    for token in RETIRED_FIXED_B_SOLVER:
        assert token not in active_code, f"retired Fixed-B solver remains: {token}"
    historical_raw_path = "docs/developer/figures/" + "joint-fixed-neighbor-scaling-r1/individual-results"
    assert historical_raw_path not in active_sources
    assert "/docs/developer/figures/**/individual-results/" in gitignore
    assert "/docs/developer/figures/**/runs.json" in gitignore
    assert "/docs/developer/figures/**/*.progress.json" in gitignore
    assert (ROOT / "docs/developer/joint-experiments.md").is_file()
    assert (ROOT / "docs/developer/joint-fixed-neighbor-experimental.md").is_file()
    print("joint experiment cleanup contract passed")


if __name__ == "__main__":
    main()
