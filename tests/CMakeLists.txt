if(TARGET GTest::gtest_main)
    set(GTEST_MAIN_TARGET GTest::gtest_main)
else()
    set(GTEST_MAIN_TARGET gtest_main)
endif()

add_custom_target(tests_all)

function(rhbm_gem_collect_gtest_patterns out_var)
    set(patterns)

    foreach(source IN LISTS ARGN)
        set(source_path "${CMAKE_CURRENT_LIST_DIR}/${source}")
        if(NOT EXISTS "${source_path}")
            message(FATAL_ERROR "Missing test source for pattern extraction: ${source}")
        endif()

        file(STRINGS "${source_path}" test_macro_lines
            REGEX "^[ \t]*(TYPED_)?TEST(_F|_P)?\\("
        )

        foreach(test_line IN LISTS test_macro_lines)
            string(REGEX REPLACE
                "^[ \t]*(TYPED_)?TEST(_F|_P)?\\(([A-Za-z_][A-Za-z0-9_]*)[ \t]*,.*"
                "\\3"
                suite_name
                "${test_line}"
            )
            if(NOT suite_name STREQUAL test_line)
                list(APPEND patterns "*${suite_name}.*")
            endif()
        endforeach()
    endforeach()

    list(REMOVE_DUPLICATES patterns)
    set(${out_var} "${patterns}" PARENT_SCOPE)
endfunction()

function(add_rhbm_gtest_group test_name)
    set(options)
    set(oneValueArgs DOMAIN)
    set(multiValueArgs INTENTS SOURCES)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT ARG_DOMAIN)
        message(FATAL_ERROR "add_rhbm_gtest_group requires DOMAIN for ${test_name}")
    endif()
    if(NOT ARG_SOURCES)
        message(FATAL_ERROR "add_rhbm_gtest_group requires SOURCES for ${test_name}")
    endif()

    rhbm_gem_collect_gtest_patterns(test_patterns ${ARG_SOURCES})
    if(NOT test_patterns)
        message(FATAL_ERROR "No GoogleTest suites found for ${test_name}")
    endif()
    list(JOIN test_patterns ":" gtest_filter)

    add_test(
        NAME ${test_name}
        COMMAND $<TARGET_FILE:rhbm_tests>
            --gtest_color=yes
            --gtest_filter=${gtest_filter}
    )

    set(test_labels "domain:${ARG_DOMAIN}")
    list(APPEND test_labels "core")
    foreach(intent IN LISTS ARG_INTENTS)
        list(APPEND test_labels "intent:${intent}")
    endforeach()
    set_tests_properties(
        ${test_name}
        PROPERTIES LABELS "${test_labels}"
    )
endfunction()

set(RHBM_GEM_CORE_COMMAND_TEST_SOURCES
    core/command/CommandScenarios_test.cpp
    core/command/MapSimulation_test.cpp
)

if(RHBM_GEM_ENABLE_UMAP)
    list(APPEND RHBM_GEM_CORE_COMMAND_TEST_SOURCES
        core/command/UmapEmbeddingCommand_test.cpp
    )
endif()

set(CORE_COMMAND_TEST_SOURCES
    core/command/CommandRunnerLifecycle_test.cpp
    core/command/CommandValidationHelpers_test.cpp
    ${RHBM_GEM_CORE_COMMAND_TEST_SOURCES}
)

set(CORE_SAMPLER_TEST_SOURCES
    core/MapSampler_test.cpp
)

set(CORE_QSCORE_TEST_SOURCES
    core/QScoreHelper_test.cpp
)

set(CORE_ESTIMATOR_TEST_SOURCES
    core/second_stage/GraphAndBackground_test.cpp
    core/second_stage/SeedAndState_test.cpp
    core/second_stage/SolverAndProposal_test.cpp
    core/second_stage/CandidateAcceptance_test.cpp
    core/second_stage/Recovery_test.cpp
    core/second_stage/ProductionFitting_test.cpp
    core/second_stage/ConvergenceAndFinalization_test.cpp
    core/second_stage/PerformanceLogging_test.cpp
    core/EstimatorTester_test.cpp
    core/TestDataFactory_test.cpp
)

set(CORE_CONTRACT_TEST_SOURCES
    core/contract/CommandCatalog_test.cpp
    core/contract/CommandExecutionContract_test.cpp
    core/contract/DocumentationLayout_test.cpp
)

set(DATA_RUNTIME_TEST_SOURCES
    data/DataObjectAssemblySpatialQuery_test.cpp
    data/DataObjectFileIO_test.cpp
    data/DataObjectMapBehavior_test.cpp
    data/DataObjectModelAnalysis_test.cpp
    data/DataObjectPrecisionContract_test.cpp
    data/DataObjectImportRegression_test.cpp
    data/RHBMTypes_test.cpp
    data/LocalPotentialSeries_test.cpp
)

set(DATA_SCHEMA_TEST_SOURCES
    data/DataObjectPersistence_test.cpp
    data/DataObjectSchemaLifecycle_test.cpp
    data/DataObjectSchemaValidation_test.cpp
)

set(UTILS_ALGORITHM_TEST_SOURCES
    utils/algorithm/ConvergenceAlgorithm_test.cpp
    utils/algorithm/KDTreeAlgorithm_test.cpp
    utils/algorithm/RobustLoss_test.cpp
    utils/algorithm/ScaleReferenceTracker_test.cpp
    utils/algorithm/WeightedRidgeSolver_test.cpp
)

set(UTILS_MATH_TEST_SOURCES
    utils/math/ArrayHelper_test.cpp
    utils/math/EigenValidation_test.cpp
    utils/math/EigenHelper_test.cpp
    utils/math/ElectricPotential_test.cpp
    utils/math/GaussianModel3D_test.cpp
    utils/math/NumericValidation_test.cpp
    utils/math/GridSampler_test.cpp
    utils/math/SphereSampler_test.cpp
)

set(UTILS_DOMAIN_TEST_SOURCES
    utils/domain/AtomSelector_test.cpp
    utils/domain/ChemicalDataHelper_test.cpp
    utils/domain/ChimeraXHelper_test.cpp
    utils/domain/ComponentHelper_test.cpp
    utils/domain/FilePathHelper_test.cpp
    utils/domain/LoggerLevelBehavior_test.cpp
    utils/domain/LoggerOutput_test.cpp
    utils/domain/LoggerProgressTimer_test.cpp
    utils/domain/LocalPainter_test.cpp
    utils/domain/ROOTHelper_test.cpp
    utils/domain/SampleFilter_test.cpp
    utils/domain/StringHelperFormatting_test.cpp
    utils/domain/StringHelperParsing_test.cpp
    utils/domain/StringHelperTokenization_test.cpp
)

set(UTILS_HRL_TEST_SOURCES
    utils/hrl/GaussianEstimationTypes_test.cpp
    utils/hrl/LinearizationService_test.cpp
    utils/hrl/RHBMTrainer_test.cpp
    utils/hrl/RHBMHelper_test.cpp
    utils/hrl/MDPDERegression_test.cpp
)

set(JOINT_COMPONENT_CONTRACT_TEST_SOURCES
    core/joint_component/BlockCoordinate_test.cpp
    core/joint_component/Search_test.cpp
    core/joint_component/Observation_test.cpp
    core/joint_component/CompactAssessment_test.cpp
    core/joint_component/Components_test.cpp
)

set(JOINT_COMPONENT_WORKFLOW_TEST_SOURCES
    core/joint_component/Runtime_test.cpp
    core/joint_component/Progress_test.cpp
    core/joint_component/PartialSelection_test.cpp
    core/joint_component/ObservableProfile_test.cpp
)

set(JOINT_COMPONENT_NUMERICAL_TEST_SOURCES
    core/joint_component/FixedNeighborBlockCoordinate_test.cpp
    core/joint_component/Numerics_test.cpp
)

set(JOINT_COMPONENT_TEST_SOURCES
    ${JOINT_COMPONENT_CONTRACT_TEST_SOURCES}
    ${JOINT_COMPONENT_WORKFLOW_TEST_SOURCES}
    ${JOINT_COMPONENT_NUMERICAL_TEST_SOURCES}
)

set(INTEGRATION_COMMAND_TEST_SOURCES
    integration/CommandApiPipeline_test.cpp
)

set(RHBM_GEM_ALL_GTEST_SOURCES
    ${CORE_COMMAND_TEST_SOURCES}
    ${CORE_ESTIMATOR_TEST_SOURCES}
    ${CORE_QSCORE_TEST_SOURCES}
    ${CORE_SAMPLER_TEST_SOURCES}
    ${CORE_CONTRACT_TEST_SOURCES}
    ${DATA_RUNTIME_TEST_SOURCES}
    ${DATA_SCHEMA_TEST_SOURCES}
    ${UTILS_ALGORITHM_TEST_SOURCES}
    ${UTILS_MATH_TEST_SOURCES}
    ${UTILS_DOMAIN_TEST_SOURCES}
    ${UTILS_HRL_TEST_SOURCES}
    ${JOINT_COMPONENT_TEST_SOURCES}
    ${INTEGRATION_COMMAND_TEST_SOURCES}
)

get_target_property(_rhbm_gem_sources rhbm_gem SOURCES)
set(_rhbm_gem_test_sources)
foreach(_source IN LISTS _rhbm_gem_sources)
    if(IS_ABSOLUTE "${_source}")
        list(APPEND _rhbm_gem_test_sources "${_source}")
    else()
        list(APPEND _rhbm_gem_test_sources "${PROJECT_SOURCE_DIR}/src/${_source}")
    endif()
endforeach()

add_library(rhbm_gem_test STATIC
    ${_rhbm_gem_test_sources}
    support/SecondStageWorkCapture.cpp
    support/SecondStageBackgroundCapture.cpp
    support/SolverFailureCapture.cpp
    support/MDPDETestSupport.cpp
    support/EndpointRefinementTestSupport.cpp
)
foreach(_property IN ITEMS
    COMPILE_FEATURES INTERFACE_COMPILE_FEATURES
    COMPILE_DEFINITIONS INTERFACE_COMPILE_DEFINITIONS
    INCLUDE_DIRECTORIES INTERFACE_INCLUDE_DIRECTORIES
    SYSTEM_INCLUDE_DIRECTORIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
    COMPILE_OPTIONS INTERFACE_COMPILE_OPTIONS
    POSITION_INDEPENDENT_CODE)
    get_target_property(_value rhbm_gem ${_property})
    if(DEFINED _value AND NOT _value MATCHES "-NOTFOUND$")
        set_property(TARGET rhbm_gem_test PROPERTY ${_property} "${_value}")
    endif()
endforeach()
target_link_libraries(rhbm_gem_test PUBLIC rhbm_gem_dependencies PRIVATE CompilerFlags)
if(RHBM_GEM_ENABLE_UMAP)
    target_link_libraries(rhbm_gem_test PRIVATE "$<BUILD_INTERFACE:libscran::umappp>")
    target_compile_definitions(rhbm_gem_test PUBLIC RHBM_GEM_ENABLE_UMAP)
endif()
target_include_directories(rhbm_gem_test PRIVATE "${PROJECT_SOURCE_DIR}/tests")
target_compile_definitions(rhbm_gem_test PRIVATE RHBM_GEM_TEST_INSTRUMENTATION)

add_executable(rhbm_tests ${RHBM_GEM_ALL_GTEST_SOURCES})
target_sources(rhbm_tests PRIVATE support/SecondStageTestSupport.cpp)
target_compile_definitions(rhbm_tests PRIVATE RHBM_GEM_TEST_INSTRUMENTATION)
set_target_properties(rhbm_tests PROPERTIES
    OUTPUT_NAME "${PROJECT_DISPLAY_NAME}-TEST"
)
target_link_libraries(rhbm_tests PRIVATE
    rhbm_gem_test
    CompilerFlags
    ${GTEST_MAIN_TARGET}
)
target_include_directories(rhbm_tests PRIVATE
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src/core"
    "${PROJECT_SOURCE_DIR}/src/data"
    "${PROJECT_SOURCE_DIR}/tests"
)
add_dependencies(tests_all rhbm_tests)

# Test-only helpers contain no historical runner or multiprecision orchestration.
add_library(joint_test_support STATIC
    support/JointSyntheticWorkload.cpp
    support/JointBlockCoordinate.cpp
    support/JointTestGeometry.cpp support/JointTestNumerics.cpp support/JointDenseReference.cpp
    support/JointTestContext.cpp support/JointComponentChecks.cpp
    support/JointFixtureSupport.cpp)
target_compile_features(joint_test_support PUBLIC cxx_std_20)
target_link_libraries(joint_test_support PUBLIC rhbm_gem_dependencies CompilerFlags)
target_include_directories(joint_test_support PUBLIC
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/tests")
target_link_libraries(rhbm_tests PRIVATE joint_test_support)

add_executable(joint_component_runtime tools/joint_component_runtime.cpp support/JointComponentRuntime.cpp)
target_link_libraries(joint_component_runtime PRIVATE joint_test_support rhbm_gem)
add_dependencies(tests_all joint_component_runtime)

if(RHBM_GEM_BUILD_RESEARCH_TOOLS)
    add_custom_target(research_tools_all)
endif()

if(RHBM_GEM_BUILD_BENCHMARKS)
    add_custom_target(benchmarks_all)
    add_executable(joint_partial_selection experiments/joint_partial_selection.cpp)
    target_link_libraries(joint_partial_selection PRIVATE joint_test_support rhbm_gem)
    add_dependencies(benchmarks_all joint_partial_selection)

    add_executable(joint_postprocessing_benchmark experiments/joint_postprocessing_benchmark.cpp)
    target_link_libraries(joint_postprocessing_benchmark PRIVATE joint_test_support rhbm_gem)
    add_dependencies(benchmarks_all joint_postprocessing_benchmark)


    add_executable(joint_fixed_neighbor_experiment experiments/joint_fixed_neighbor.cpp)
    target_link_libraries(joint_fixed_neighbor_experiment PRIVATE joint_test_support rhbm_gem_test CompilerFlags)
    target_compile_definitions(joint_fixed_neighbor_experiment PRIVATE RHBM_GEM_TEST_INSTRUMENTATION)
    add_dependencies(benchmarks_all joint_fixed_neighbor_experiment)
endif()

if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
    add_library(joint_offline_support STATIC support/JointOfflineAudit.cpp
        support/JointPrecisionAudit.cpp support/JointOfflineTools.cpp
        support/JointLocalEvidence.cpp support/JointDerivativeAudit.cpp)
    target_link_libraries(joint_offline_support PUBLIC joint_test_support)
endif()

if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
    add_executable(joint_component_audit tools/joint_component_audit.cpp)
    target_link_libraries(joint_component_audit PRIVATE joint_offline_support rhbm_gem)
    add_executable(joint_offline_diagnostic experiments/joint_offline_diagnostic.cpp)
    target_link_libraries(joint_offline_diagnostic PRIVATE joint_offline_support rhbm_gem)
    add_executable(joint_offline_tests core/joint_component/Precision_test.cpp
        core/joint_component/OfflineComponents_test.cpp)
    target_link_libraries(joint_offline_tests PRIVATE joint_offline_support rhbm_gem ${GTEST_MAIN_TARGET})
    add_test(NAME joint_component_offline_tests COMMAND joint_offline_tests)
    set_tests_properties(joint_component_offline_tests PROPERTIES
        LABELS "offline;domain:core;intent:algorithm;joint:offline;long" TIMEOUT 3600
        ENVIRONMENT "OMP_NUM_THREADS=1")
    add_dependencies(tests_all joint_component_audit joint_offline_diagnostic joint_offline_tests)
endif()

if(RHBM_GEM_BUILD_RESEARCH_TOOLS)
    add_executable(joint_statistical_experiment experiments/joint_statistical_experiment.cpp)
    target_link_libraries(joint_statistical_experiment PRIVATE joint_test_support rhbm_gem)
    add_dependencies(research_tools_all joint_statistical_experiment)
endif()

add_rhbm_gtest_group(rhbm_tests_joint_component_contract
    DOMAIN core INTENTS algorithm contract SOURCES ${JOINT_COMPONENT_CONTRACT_TEST_SOURCES})
set_tests_properties(rhbm_tests_joint_component_contract PROPERTIES
    LABELS "core;domain:core;intent:algorithm;intent:contract;joint:contract" TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_rhbm_gtest_group(rhbm_tests_joint_component_workflow
    DOMAIN core INTENTS algorithm SOURCES ${JOINT_COMPONENT_WORKFLOW_TEST_SOURCES})
set_tests_properties(rhbm_tests_joint_component_workflow PROPERTIES
    LABELS "core;domain:core;intent:algorithm;joint:workflow" TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=1")

rhbm_gem_collect_gtest_patterns(joint_numerical_patterns ${JOINT_COMPONENT_NUMERICAL_TEST_SOURCES})
list(JOIN joint_numerical_patterns ":" joint_numerical_filter)
string(APPEND joint_numerical_filter
    "-JointComponentNumericsTest.CompactJacobianGramIdentityDiagnosticsAcrossLattices")
add_test(NAME rhbm_tests_joint_component_numerical
    COMMAND $<TARGET_FILE:rhbm_tests> --gtest_color=yes --gtest_filter=${joint_numerical_filter})
set_tests_properties(rhbm_tests_joint_component_numerical PROPERTIES
    LABELS "core;domain:core;intent:algorithm;joint:numerical" TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_test(NAME joint_component_gram_scalability
    COMMAND $<TARGET_FILE:rhbm_tests> --gtest_color=yes
        --gtest_filter=JointComponentNumericsTest.CompactJacobianGramIdentityDiagnosticsAcrossLattices)
set_tests_properties(joint_component_gram_scalability
    PROPERTIES LABELS "core;domain:core;intent:algorithm;joint:scalability;joint:offline;long" TIMEOUT 3600
        ENVIRONMENT "OMP_NUM_THREADS=1")

set(joint_persistence_test_patterns
    DataObjectPersistenceTest.JointResultsRoundTripCopyClearAndReplaceAtomically
    DataObjectPersistenceTest.JointSolverProvenanceRoundTripsAndMissingFieldRemainsReadable
    DataObjectPersistenceTest.FixedNeighborPolicyVersionsPreserveLegacyAndRejectRetiredV2Claims
    DataObjectPersistenceTest.CurrentFixedNeighborSolverProvenanceRoundTrips
    DataObjectPersistenceTest.HistoricalV2SolverProvenanceRemainsReadable
    DataObjectPersistenceTest.MalformedCurrentSolverProvenanceIsRejected
    DataObjectPersistenceTest.JointContributorSubsetRoundTripsAndRejectsForeignOrHydrogenIds
    DataObjectPersistenceTest.JointCodecRejectsMalformedMappingsAndPreservesMissingValues
    DataObjectPersistenceTest.JointMetadataRejectsInvalidValuesAndOldDocumentsWithoutWrites
    DataObjectPersistenceTest.JointSelectionMetadataValidationAndCsvRoles
    DataObjectPersistenceTest.JointNeutralRoundTripPreservesSourcesGeometryAndUnavailableStates
    DataObjectPersistenceTest.LegacyJointSnapshotMapsPointsWithoutInventingDerivedEvidence
    DataObjectPersistenceTest.CanonicalNativeDiagnosticsAndMissingStatesRoundTripExactly
)
list(JOIN joint_persistence_test_patterns ":" joint_persistence_test_filter)
add_test(NAME joint_persistence_contract_test
    COMMAND $<TARGET_FILE:rhbm_tests> --gtest_color=yes --gtest_filter=${joint_persistence_test_filter})
set_tests_properties(joint_persistence_contract_test PROPERTIES
    LABELS "core;domain:data;intent:contract;joint:workflow" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1")

add_rhbm_gtest_group(rhbm_tests_core_command
    DOMAIN core
    INTENTS command validation
    SOURCES ${CORE_COMMAND_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_core_sampler
    DOMAIN core
    INTENTS sampling
    SOURCES ${CORE_SAMPLER_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_core_qscore
    DOMAIN core
    INTENTS algorithm
    SOURCES ${CORE_QSCORE_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_core_estimator
    DOMAIN core
    INTENTS algorithm
    SOURCES ${CORE_ESTIMATOR_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_core_contract
    DOMAIN core
    INTENTS contract
    SOURCES ${CORE_CONTRACT_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_data_runtime
    DOMAIN data
    INTENTS io
    SOURCES ${DATA_RUNTIME_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_data_schema
    DOMAIN data
    INTENTS schema
    SOURCES ${DATA_SCHEMA_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_utils_algorithm
    DOMAIN utils
    INTENTS algorithm
    SOURCES ${UTILS_ALGORITHM_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_utils_math
    DOMAIN utils
    INTENTS algorithm
    SOURCES ${UTILS_MATH_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_utils_domain
    DOMAIN utils
    INTENTS validation
    SOURCES ${UTILS_DOMAIN_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_utils_hrl
    DOMAIN utils
    INTENTS algorithm
    SOURCES ${UTILS_HRL_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_integration_command
    DOMAIN integration
    INTENTS command
    SOURCES ${INTEGRATION_COMMAND_TEST_SOURCES}
)

if(NOT DEFINED Python_EXECUTABLE OR Python_EXECUTABLE STREQUAL "")
    find_package(Python QUIET COMPONENTS Interpreter)
endif()

if(DEFINED Python_EXECUTABLE AND NOT Python_EXECUTABLE STREQUAL "")
    add_test(NAME experiment_support_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/experiment_support_test.py")
    add_test(NAME experiment_runner_dependency_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/runner_dependency_test.py")
    add_test(NAME joint_benchmark_contract_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_benchmark_test.py")
    add_test(NAME joint_experiment_cleanup_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_experiment_cleanup_test.py")
    set_tests_properties(experiment_support_test experiment_runner_dependency_test
        PROPERTIES LABELS "core;domain:integration;intent:algorithm")
    set_tests_properties(joint_benchmark_contract_test
        PROPERTIES LABELS "core;domain:integration;intent:algorithm;joint:numerical" TIMEOUT 300
        ENVIRONMENT "OMP_NUM_THREADS=1")
    set_tests_properties(joint_experiment_cleanup_test
        PROPERTIES LABELS "core;domain:integration;intent:algorithm;joint:contract" TIMEOUT 180
        ENVIRONMENT "OMP_NUM_THREADS=1")
    add_test(NAME joint_component_evidence_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_evidence_test.py")
    add_test(NAME joint_component_runtime_runner_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime_test.py")
    add_test(NAME joint_component_regression
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime.py"
            regression --executable $<TARGET_FILE:joint_component_runtime> --work-dir "${PROJECT_BINARY_DIR}/joint-regression")
    add_test(NAME joint_component_physical_smoke
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime.py"
            physical-smoke --executable $<TARGET_FILE:joint_component_runtime>)
    set_tests_properties(joint_component_evidence_test
        PROPERTIES LABELS "core;domain:integration;intent:algorithm;joint:contract" TIMEOUT 180
        ENVIRONMENT "OMP_NUM_THREADS=1")
    set_tests_properties(
        joint_component_runtime_runner_test
        joint_component_regression joint_component_physical_smoke
        PROPERTIES LABELS "core;domain:integration;intent:algorithm;joint:workflow" TIMEOUT 600
        ENVIRONMENT "OMP_NUM_THREADS=1")

    if(RHBM_GEM_BUILD_BENCHMARKS)
        add_test(NAME joint_benchmark_smoke
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_benchmark_test.py"
                --smoke "${PROJECT_BINARY_DIR}")
        set_tests_properties(joint_benchmark_smoke PROPERTIES
            LABELS "benchmark;domain:integration;intent:benchmark;joint:numerical" TIMEOUT 300
            ENVIRONMENT "OMP_NUM_THREADS=1")
    endif()

    if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
        add_test(NAME joint_offline_diagnostic_smoke
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_offline_diagnostic.py"
                --work-dir "${PROJECT_BINARY_DIR}/joint-offline-diagnostic"
                --executable $<TARGET_FILE:joint_offline_diagnostic>
                --output "${PROJECT_BINARY_DIR}/joint-offline-diagnostic/weak-halo.json")
        set_tests_properties(joint_offline_diagnostic_smoke
            PROPERTIES LABELS "offline;domain:integration;intent:algorithm;joint:offline;long" TIMEOUT 3600
            ENVIRONMENT "OMP_NUM_THREADS=1")
    endif()

    if(RHBM_GEM_BUILD_RESEARCH_TOOLS)
        add_test(NAME joint_statistical_experiment_contract_test
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_statistical_experiment_test.py")
        add_test(NAME joint_statistical_experiment_smoke
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_statistical_experiment.py"
                --smoke --work-dir "${PROJECT_BINARY_DIR}/joint-statistical-smoke"
                --executable $<TARGET_FILE:joint_statistical_experiment>)
        set_tests_properties(joint_statistical_experiment_contract_test joint_statistical_experiment_smoke
            PROPERTIES LABELS "research;domain:integration;intent:algorithm;joint:research")
    endif()

    if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
        add_test(NAME joint_component_offline_regression
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_audit.py"
                --executable $<TARGET_FILE:joint_component_audit> --work-dir "${PROJECT_BINARY_DIR}/joint-offline")
        set_tests_properties(joint_component_offline_regression
            PROPERTIES LABELS "offline;domain:integration;intent:algorithm;joint:offline;long" TIMEOUT 3600
            ENVIRONMENT "OMP_NUM_THREADS=1")
    endif()
    if(RHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS)
        add_test(NAME joint_component_extended_regression
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime.py"
                regression --extended --executable $<TARGET_FILE:joint_component_runtime> --work-dir "${PROJECT_BINARY_DIR}/joint-extended")
        set_tests_properties(joint_component_extended_regression
            PROPERTIES LABELS "extended;domain:integration;intent:algorithm;joint:offline;long" TIMEOUT 3600
            ENVIRONMENT "OMP_NUM_THREADS=1")
    endif()
endif()

if(UNIX)
    find_program(BASH_EXECUTABLE NAMES bash)
    if(BASH_EXECUTABLE)
        add_test(
            NAME cli_quickstart_smoke_test
            COMMAND "${BASH_EXECUTABLE}"
                "${PROJECT_SOURCE_DIR}/tests/integration/cli_quickstart_smoke.sh"
                "${PROJECT_SOURCE_DIR}"
        )
        set_tests_properties(
            cli_quickstart_smoke_test
            PROPERTIES LABELS "core;domain:integration;intent:command"
        )
    endif()
endif()

if(DEFINED Python_EXECUTABLE)
    add_dependencies(tests_all rhbm_gem_cli)
    add_test(NAME joint_workflow_cli_smoke_test
        COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/integration/joint_workflow_cli_smoke.py
            $<TARGET_FILE:rhbm_gem_cli>)
    set_tests_properties(joint_workflow_cli_smoke_test
        PROPERTIES LABELS "core;domain:integration;intent:command;joint:workflow" TIMEOUT 300
        ENVIRONMENT "OMP_NUM_THREADS=1")
    add_test(
        NAME command_scaffold_smoke_test
        COMMAND ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/command_scaffold_smoke.py
            ${PROJECT_SOURCE_DIR}
    )
    set_tests_properties(
        command_scaffold_smoke_test
        PROPERTIES LABELS "core;domain:integration;intent:contract"
    )
endif()

if(BUILD_PYTHON_BINDINGS AND TARGET rhbm_gem_module AND DEFINED Python_EXECUTABLE)
    add_dependencies(tests_all rhbm_gem_module)

    function(add_python_integration_test test_name script_name labels)
        add_test(
        NAME ${test_name}
        COMMAND ${CMAKE_COMMAND} -E env
                "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
                "RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE=${RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE}"
                "RHBM_GEM_ENABLE_UMAP=${RHBM_GEM_ENABLE_UMAP}"
                ${Python_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/tests/integration/${script_name}
        )
        set_tests_properties(
            ${test_name}
            PROPERTIES LABELS "${labels}"
        )
    endfunction()

    add_python_integration_test(
        bindings_runtime_smoke_test
        python_bindings_runtime_smoke.py
        "domain:integration;intent:bindings"
    )
    add_python_integration_test(
        python_examples_smoke_test
        python_examples_smoke.py
        "domain:integration;intent:bindings"
    )
endif()

get_target_property(_production_sources rhbm_gem SOURCES)
foreach(_source IN LISTS _production_sources)
    if(_source MATCHES "(^|/)tests/support/")
        message(FATAL_ERROR "Test support source leaked into rhbm_gem: ${_source}")
    endif()
endforeach()
get_target_property(_production_definitions rhbm_gem COMPILE_DEFINITIONS)
if(_production_definitions MATCHES "RHBM_GEM_TEST_INSTRUMENTATION")
    message(FATAL_ERROR "Test instrumentation leaked into rhbm_gem")
endif()

get_target_property(_tests_all_dependencies tests_all MANUALLY_ADDED_DEPENDENCIES)
foreach(_optional_target IN ITEMS
    joint_partial_selection
    joint_postprocessing_benchmark
    joint_fixed_neighbor_experiment
    joint_statistical_experiment)
    if(_optional_target IN_LIST _tests_all_dependencies)
        message(FATAL_ERROR "Optional target ${_optional_target} must not be a tests_all dependency")
    endif()
endforeach()

if(RHBM_GEM_BUILD_BENCHMARKS)
    foreach(_target IN ITEMS joint_partial_selection joint_postprocessing_benchmark joint_fixed_neighbor_experiment)
        if(NOT TARGET ${_target})
            message(FATAL_ERROR "RHBM_GEM_BUILD_BENCHMARKS did not create ${_target}")
        endif()
    endforeach()
else()
    foreach(_target IN ITEMS joint_partial_selection joint_postprocessing_benchmark joint_fixed_neighbor_experiment)
        if(TARGET ${_target})
            message(FATAL_ERROR "${_target} exists while RHBM_GEM_BUILD_BENCHMARKS=OFF")
        endif()
    endforeach()
endif()

if(RHBM_GEM_BUILD_RESEARCH_TOOLS)
    foreach(_target IN ITEMS joint_statistical_experiment research_tools_all)
        if(NOT TARGET ${_target})
            message(FATAL_ERROR "RHBM_GEM_BUILD_RESEARCH_TOOLS did not create ${_target}")
        endif()
    endforeach()
elseif(TARGET joint_statistical_experiment OR TARGET research_tools_all)
    message(FATAL_ERROR "Research targets exist while RHBM_GEM_BUILD_RESEARCH_TOOLS=OFF")
endif()

if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
    if(NOT TARGET joint_component_audit OR NOT TARGET joint_offline_tests)
        message(FATAL_ERROR "RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS did not create its audit targets")
    endif()
elseif(TARGET joint_component_audit OR TARGET joint_offline_tests)
    message(FATAL_ERROR "Offline audit targets exist while RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=OFF")
endif()

if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
    if(NOT TARGET joint_offline_diagnostic OR NOT TARGET joint_offline_support)
        message(FATAL_ERROR "RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS did not create its diagnostic support")
    endif()
elseif(TARGET joint_offline_diagnostic OR TARGET joint_offline_support)
    message(FATAL_ERROR "Offline diagnostic targets exist while RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS=OFF")
endif()
