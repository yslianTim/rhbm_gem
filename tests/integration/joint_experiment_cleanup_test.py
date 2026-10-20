"""Check the current Joint experiment surface stays within its cleanup contract."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CURRENT_EXPERIMENT_SOURCES = {
    "joint_fullabc_width_study.cpp",
    "joint_fixed_width_oracle.cpp",
    "joint_fixed_neighbor.cpp",
    "joint_local_ac_kkt.cpp",
    "joint_local_replay_analysis.cpp",
    "joint_partial_selection.cpp",
    "joint_postprocessing_benchmark.cpp",
    "joint_width_trajectory.cpp",
    "joint_width_repro_no_instrumentation.cpp",
    "joint_production_bounded_synthetic.cpp",
    "joint_label_free_halo_policy.cpp",
    "joint_halo_geometry_census.cpp",
    "joint_halo_fullabc_validation.cpp",
}

CURRENT_PYTHON_EXPERIMENTS = {
    "joint_benchmark.py",
    "joint_partial_selection.py",
}

CURRENT_PYTHON_VALIDATION_RUNNERS = {
    "joint_component_runtime.py",
    "joint_workflow_cli_smoke.py",
}

PRODUCTION_REGRESSION_OWNERS = (
    "tests/core/joint_component/ActiveSetBoundedLM_test.cpp",
    "tests/core/joint_component/WidthCalibration_test.cpp",
    "tests/core/joint_component/Numerics_test.cpp",
    "tests/core/joint_component/FixedNeighborBlockCoordinate_test.cpp",
    "tests/core/joint_component/LocalConditioning_test.cpp",
    "tests/core/joint_component/ParameterStability_test.cpp",
    "tests/core/joint_component/PureHaloSupportExclusion_test.cpp",
    "tests/core/joint_component/Runtime_test.cpp",
    "tests/core/joint_component/Progress_test.cpp",
    "tests/core/joint_component/TargetEvidence_test.cpp",
    "tests/data/DataObjectPersistence_test.cpp",
    "tests/integration/CommandApiPipeline_test.cpp",
    "tests/experiments/joint_production_bounded_synthetic.cpp",
    "tests/experiments/joint_width_repro_no_instrumentation.cpp",
    "tests/experiments/joint_fixed_neighbor.cpp",
)

RETIRED_REPOSITORY_FILES = (
    "tests/experiments/joint_sparse_benchmark.cpp",
    "tests/core/joint_component/FixedBBlockCoordinate_test.cpp",
    "tests/core/joint_component/OperatorSearch_test.cpp",
    "tests/core/joint_component/ProfileOperator_test.cpp",
    "tests/core/joint_component/ProjectedTailQr_test.cpp",
    "src/core/detail/joint_component/SparseFactor.hpp",
    "src/core/detail/joint_component/SparseFactor.cpp",
    "tests/core/joint_component/BoundedRank_test.cpp",
    "tests/core/joint_component/ProjectedWidthResearch_test.cpp",
    "tests/core/joint_component/OfflineComponents_test.cpp",
    "tests/tools/joint_component_audit.cpp",
    "tests/integration/joint_component_audit.py",
    "tests/experiments/joint_offline_diagnostic.cpp",
    "tests/integration/joint_offline_diagnostic.py",
    "tests/integration/joint_offline_diagnostic_test.py",
    "tests/integration/joint_component_evidence_test.py",
    "tests/integration/joint_offline_support.py",
    "tests/support/BoundedRankResearch.cpp",
    "tests/support/BoundedRankResearch.hpp",
    "tests/support/ProjectedWidthResearch.cpp",
    "tests/support/ProjectedWidthResearch.hpp",
    "tests/support/MatrixFreeProjectedWidth.cpp",
    "tests/support/MatrixFreeProjectedWidth.hpp",
    "tests/support/StructuredProjectionResearch.cpp",
    "tests/support/StructuredProjectionResearch.hpp",
    "tests/support/JointOfflineTools.cpp",
    "tests/support/JointOfflineTools.hpp",
    "tests/support/JointLocalEvidence.cpp",
    "tests/support/JointLocalEvidence.hpp",
    "tests/support/JointDerivativeAudit.cpp",
    "tests/support/JointDerivativeAudit.hpp",
    "tests/experiments/joint_width_trust_experiment.cpp",
    "tests/experiments/run_joint_width_trust_paired.py",
    "tests/experiments/run_joint_global_width_bound_paired.py",
    "tests/integration/joint_global_width_bound_experiment_test.py",
    "tests/support/JointWidthTrustPolicy.cpp",
    "tests/support/JointWidthTrustPolicy.hpp",
    "tests/core/joint_component/WidthTrustPolicy_test.cpp",
    "tests/experiments/joint_dynamic_halo_mean.cpp",
    "tests/experiments/joint_dynamic_halo_joint.cpp",
    "tests/experiments/joint_dynamic_halo_synthetic_robustness.cpp",
    "tests/experiments/joint_halo_pruning_validation.cpp",
    "tests/experiments/joint_halo_observation_exclusion.cpp",
    "tests/experiments/joint_halo_observation_exclusion_fit.cpp",
    "tests/experiments/joint_halo_observation_exclusion_metrics.cpp",
    "tests/experiments/joint_halo_support_exclusion_synthetic.cpp",
    "tests/experiments/joint_halo_support_exclusion_validation.cpp",
    "tests/experiments/joint_halo_support_weighting_reference.cpp",
    "tests/experiments/fixed_halo_first_stage_diagnostic.cpp",
    "tests/experiments/joint_statistical_experiment.cpp",
    "tests/experiments/summarize_joint_halo_noise.py",
    "tests/core/joint_component/DynamicHaloMean_test.cpp",
    "tests/support/DynamicHaloMean.cpp",
    "tests/support/DynamicHaloMean.hpp",
    "tests/integration/joint_statistical_experiment.py",
    "tests/integration/joint_statistical_experiment_test.py",
    "tests/experiments/joint_fixed_halo_validation.cpp",
    "tests/core/joint_component/FixedHaloInitialization_test.cpp",
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

CURRENT_JOINT_DOCUMENTS = (
    "docs/developer/joint-benchmark.md",
    "docs/developer/joint-component-runtime.md",
    "docs/developer/joint-experiments.md",
    "docs/developer/joint-fixed-neighbor-experimental.md",
    "docs/developer/joint-operator-search.md",
    "docs/developer/joint-component-sparse-backend.md",
    "docs/developer/joint-fixed-actions.md",
    "docs/developer/joint-profile-operator.md",
    "docs/developer/joint-fixed-neighbor-workspace-attribution.md",
    "docs/developer/joint-component-evidence.md",
    "docs/developer/commands/potential_analysis.md",
)

RETIRED_LOCAL_SEARCH_DESCRIPTIONS = (
    "local legacycompact",
    "local `legacycompact`",
    "legacycompact profile solver",
    "legacycompact profile search",
    "legacycompact remains an active local",
    "legacycompact remains the local",
    "legacycompact is still active inside",
    "local search | legacycompact",
    "local search | `legacycompact`",
    "local search: legacycompact",
    "local search: `legacycompact`",
    "legacy-compact local search",
    "`legacy-compact` local search",
    "fixed_neighbor_local_search=legacy-compact",
    "fixed_neighbor_local_search: legacy-compact",
)

def main():
    integration = ROOT / "tests" / "integration"
    cmake = (ROOT / "tests" / "CMakeLists.txt").read_text()
    fixed_neighbor_source = (ROOT / "tests/experiments/joint_fixed_neighbor.cpp").read_text()
    benchmark_driver = (integration / "joint_benchmark.py").read_text()
    runtime_runner = (integration / "joint_component_runtime.py").read_text()
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

    # Check current surfaces only. The JSON decoder, persistence fixtures and
    # historical documentation intentionally retain old provenance values.
    core_source = "\n".join(
        path.read_text(errors="ignore")
        for path in (ROOT / "src/core").rglob("*")
        if path.is_file() and path.suffix in {".c", ".cc", ".cpp", ".h", ".hh", ".hpp", ".hxx"}
    )
    assert '"legacy-compact"' not in core_source, (
        "current Joint core must not emit the legacy local-search label")
    assert "'fixed_neighbor_local_search': 'profile-lm'" in benchmark_driver
    assert "'fixed_neighbor_local_search': 'LegacyCompact'" not in benchmark_driver
    assert "'fixed_neighbor_local_search': 'legacy-compact'" not in benchmark_driver
    current_joint_tests = (
        ROOT / "tests/core/joint_component/Progress_test.cpp",
        integration / "joint_benchmark_test.py",
        integration / "joint_workflow_cli_smoke.py",
    )
    current_test_text = "\n".join(path.read_text() for path in current_joint_tests)
    assert "local-search=profile-lm" in current_test_text
    assert "fixed_neighbor_local_search" in current_test_text
    assert "fixed_neighbor_local_search': 'LegacyCompact'" not in current_test_text
    assert 'fixed_neighbor_local_search": "legacy-compact"' not in current_test_text
    for relative_path in CURRENT_JOINT_DOCUMENTS:
        document = (ROOT / relative_path).read_text().lower()
        assert "profile lm" in document, (
            f"current Joint documentation omits the Profile LM name: {relative_path}")
        for phrase in RETIRED_LOCAL_SEARCH_DESCRIPTIONS:
            assert phrase not in document, (
                f"current Joint documentation describes the local solver as LegacyCompact: "
                f"{relative_path}: {phrase}")

    persistence_tests = (ROOT / "tests/data/DataObjectPersistence_test.cpp").read_text()
    assert 'old_current_decoded.metadata.solver->fixed_neighbor_local_search,"legacy-compact"' in persistence_tests
    assert 'legacy_provenance.search_method="legacy-compact"' in persistence_tests
    assert 'operator_provenance.search_method="operator-pcg"' in persistence_tests

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
    missing_regression_owners = [
        path for path in PRODUCTION_REGRESSION_OWNERS if not (ROOT / path).is_file()
    ]
    unregistered_regression_owners = [
        path for path in PRODUCTION_REGRESSION_OWNERS
        if Path(path).name not in cmake
    ]
    assert not missing_regression_owners, (
        f"Production regression owners were removed: {missing_regression_owners}")
    assert not unregistered_regression_owners, (
        f"Production regression owners lost their CMake registration: {unregistered_regression_owners}")
    assert "add_test(NAME joint_production_bounded_synthetic" in cmake, (
        "bound-active Production synthetic must remain a normal CTest")
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
    for token in (
        "joint_width_trust_experiment",
        "joint_global_width_bound_experiment",
        "JointWidthTrustPolicy",
        "width_trial_trust_for_testing",
        "width_trust_policy_rejected",
    ):
        assert token not in cmake, f"retired Width-trust CMake entry returned: {token}"
        assert token not in active_code, f"retired Width-trust tool or hook returned: {token}"
    for token in (
        "joint_dynamic_halo_",
        "joint_halo_pruning_validation",
        "joint_halo_observation_exclusion",
        "joint_halo_support_exclusion_synthetic",
        "joint_halo_support_exclusion_validation",
        "joint_halo_support_weighting_reference",
        "joint_statistical_experiment",
        "DynamicHaloMean",
    ):
        assert token not in cmake, f"retired Halo research CMake entry returned: {token}"
        assert token not in active_code, f"retired Halo research tool or helper returned: {token}"
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
    for token in (
        "BuildJointProblem",
        "BuildFixedHaloResearchProblem",
        "EstimateJointComponentsFixedHalo",
        "WithFixedHalos",
        "ReplayFixedHaloFullModel",
    ):
        assert token not in active_code, (
            f"retired FixedHalo estimation entry point returned: {token}")
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
