"""Check the current Joint experiment surface stays within its cleanup contract."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

RETIRED_DRIVERS = (
    "joint_fixed_neighbor_core_size.py",
    "joint_fixed_neighbor_outer_core_endpoint.py",
    "joint_fixed_neighbor_outer_core_frontier.py",
    "joint_optimized_frontier.py",
    "joint_fixed_neighbor_analysis.py",
    "joint_fixed_neighbor_attribution.py",
    "joint_fixed_neighbor_attribution_campaign.py",
    "joint_fixed_neighbor_endpoint_decomposition.py",
    "joint_fixed_neighbor_hybrid.py",
    "joint_fixed_neighbor_local_attribution.py",
    "joint_fixed_neighbor_local_certification.py",
    "joint_fixed_neighbor_local_schwarz.py",
    "joint_fixed_neighbor_local_schwarz_local_work.py",
    "joint_fixed_neighbor_numeric_reuse.py",
    "joint_fixed_neighbor_operator_diagnostic.py",
    "joint_fixed_neighbor_operator_route.py",
    "joint_fixed_neighbor_workspace_timing.py",
    "joint_fixed_neighbor_inexact.py",
    "joint_fixed_neighbor_inexact_campaign.py",
    "joint_fixed_neighbor_inexact_qualification.py",
)

RETIRED_TESTS = (
    "joint_route_frontier_test.py",
    "joint_optimized_frontier_test.py",
    "joint_factor_ownership_analysis_test.py",
    "joint_projected_tail_analysis_test.py",
    "joint_projected_width_analysis_test.py",
    "joint_fixed_neighbor_analysis_test.py",
    "joint_fixed_neighbor_endpoint_decomposition_test.py",
    "joint_fixed_neighbor_local_certification_test.py",
    "joint_fixed_neighbor_attribution_test.py",
    "joint_fixed_neighbor_local_attribution_test.py",
    "joint_fixed_neighbor_numeric_reuse_test.py",
    "joint_fixed_neighbor_workspace_timing_test.py",
    "joint_fixed_neighbor_operator_diagnostic_test.py",
    "joint_fixed_neighbor_local_schwarz_test.py",
    "joint_fixed_neighbor_local_schwarz_local_work_test.py",
    "joint_fixed_neighbor_hybrid_test.py",
    "joint_fixed_neighbor_operator_route_test.py",
    "joint_fixed_neighbor_inexact_test.py",
    "joint_fixed_neighbor_inexact_qualification_test.py",
)

RETIRED_REPOSITORY_FILES = (
    "tests/experiments/joint_sparse_benchmark.cpp",
    "tests/core/joint_component/OperatorSearch_test.cpp",
    "tests/core/joint_component/ProfileOperator_test.cpp",
)

RETIRED_FIXTURES = (
    "joint_fixed_neighbor_inexact_baseline.json",
)

CLOSED_MODES = (
    "--decompose", "--certified-local", "--qualification", "--attribution",
    "--matched-control", "--local-legacy", "--operator-diagnostic",
    "--local-operator-identity", "--local-operator-diagonal",
    "--local-operator-schwarz", "--local-operator-identity-search",
    "--local-operator-diagonal-search", "--local-operator-schwarz-search",
    "--local-operator-schwarz-two", "--local-operator-schwarz-full",
    "--local-operator-schwarz-polish-one", "--local-operator-schwarz-polish-two",
    "--tile-1024", "--tile-2048", "--tile-4096", "--tile-8192", "--tile-16384",
    "--inexact-one", "--inexact-one-search", "--inexact-one-endpoint",
    "--inexact-two", "--full-attribution",
    "--joint-search", "operator-pcg", "--operator-rank",
    "--operator-rank-seconds", "--operator-rank-work-entries",
    "--operator-rank-workspace-mib", "--schwarz-core-atoms",
    "--schwarz-overlap-hops", "--schwarz-max-block-atoms",
    "--schwarz-storage-mib", "--schwarz-scratch-mib",
)


def main():
    integration = ROOT / "tests" / "integration"
    cmake = (ROOT / "tests" / "CMakeLists.txt").read_text()
    experiment = (ROOT / "tests" / "experiments" / "joint_fixed_neighbor.cpp").read_text()
    gitignore = (ROOT / ".gitignore").read_text()
    active_sources = "\n".join(path.read_text() for path in integration.glob("*.py"))
    active_paths = [path for path in (ROOT / "src/core").rglob("*") if path.is_file()]
    active_paths.extend(path for path in (ROOT / "include/rhbm_gem/core").rglob("*") if path.is_file())
    active_paths.extend(path for path in (ROOT / "tests" / "experiments").rglob("*")
                        if path.is_file())
    active_paths.extend(path for path in integration.glob("*.py")
                        if path.name != Path(__file__).name)
    active_code = "\n".join(
        path.read_text(errors="ignore")
        for path in active_paths
        if path.suffix in {".cpp", ".hpp", ".py"}
    )

    missing = [name for name in RETIRED_DRIVERS + RETIRED_TESTS
               if (integration / name).exists()]
    missing_fixtures = [name for name in RETIRED_FIXTURES
                        if (ROOT / "tests" / "fixtures" / name).exists()]
    present_files = [name for name in RETIRED_REPOSITORY_FILES
                     if (ROOT / name).exists()]
    assert not missing, f"retired Joint files returned: {missing}"
    assert not missing_fixtures, f"retired Joint fixtures returned: {missing_fixtures}"
    assert not present_files, f"retired Operator files returned: {present_files}"
    assert not any(name.removesuffix(".py") in cmake for name in RETIRED_TESTS)
    assert not any(mode in experiment for mode in CLOSED_MODES)
    for token in (
        "FixedNeighborLocalWork",
        "TwoAcceptedUpdates",
        "FullLocalSearch",
        "--fixed-local-work",
        "--inexact-one",
        "--inexact-two",
    ):
        assert token not in active_code, f"retired FixedNeighbor token remains: {token}"
    historical_raw_path = "docs/developer/figures/" + "joint-fixed-neighbor-scaling-r1/individual-results"
    assert historical_raw_path not in active_sources
    assert "/docs/developer/figures/**/individual-results/" in gitignore
    assert "/docs/developer/figures/**/runs.json" in gitignore
    assert "/docs/developer/figures/**/*.progress.json" in gitignore
    assert (integration / "joint_fixed_neighbor_outer_core_qualification.py").is_file()
    assert (integration / "joint_fixed_neighbor_outer_core_support.py").is_file()
    assert (ROOT / "docs/developer/joint-experiments.md").is_file()
    assert (ROOT / "docs/developer/joint-fixed-neighbor-experimental.md").is_file()
    print("joint experiment cleanup contract passed")


if __name__ == "__main__":
    main()
