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
    core/second_stage/Observation_test.cpp
    core/second_stage/NumericalAuditNeutrality_test.cpp
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
    utils/hrl/MDPDEExperiment_test.cpp
)

set(JOINT_COMPONENT_TEST_SOURCES
    core/joint_component/FreeDesignRank_test.cpp
    core/joint_component/ProfileOperator_test.cpp
    core/joint_component/OperatorSearch_test.cpp
    core/joint_component/Numerics_test.cpp
    core/joint_component/Search_test.cpp
    core/joint_component/Components_test.cpp
    core/joint_component/Observation_test.cpp
    core/joint_component/Runtime_test.cpp
    core/joint_component/PartialSelection_test.cpp
    core/joint_component/ObservableProfile_test.cpp
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

add_executable(rhbm_tests ${RHBM_GEM_ALL_GTEST_SOURCES})
target_sources(rhbm_tests PRIVATE support/SecondStageTestSupport.cpp)
target_compile_definitions(rhbm_tests PRIVATE RHBM_GEM_TEST_INSTRUMENTATION)
set_target_properties(rhbm_tests PROPERTIES
    OUTPUT_NAME "${PROJECT_DISPLAY_NAME}-TEST"
)
target_link_libraries(rhbm_tests PRIVATE
    rhbm_gem
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
rhbm_gem_apply_observation_definitions(rhbm_tests)

add_dependencies(tests_all rhbm_tests)

# Test-only helpers contain no historical runner or multiprecision orchestration.
add_library(joint_test_support STATIC
    support/JointOperatorWorkload.cpp
    support/JointTestGeometry.cpp support/JointTestNumerics.cpp support/JointDenseReference.cpp
    support/JointTestContext.cpp support/JointComponentChecks.cpp
    support/JointFixtureSupport.cpp)
target_link_libraries(joint_test_support PUBLIC rhbm_gem CompilerFlags)
target_include_directories(joint_test_support PUBLIC "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/tests")
target_link_libraries(rhbm_tests PRIVATE joint_test_support)

add_executable(mdpde_experiment experiments/mdpde_experiment.cpp support/ForwardModelExperiment.cpp)
target_link_libraries(mdpde_experiment PRIVATE rhbm_gem CompilerFlags)
target_include_directories(mdpde_experiment PRIVATE "${PROJECT_SOURCE_DIR}/src" "${PROJECT_SOURCE_DIR}/tests")
target_compile_definitions(mdpde_experiment PRIVATE RHBM_GEM_TEST_INSTRUMENTATION)

add_executable(joint_component_runtime experiments/joint_component_runtime.cpp support/JointComponentRuntime.cpp)
target_link_libraries(joint_component_runtime PRIVATE joint_test_support)
add_dependencies(tests_all joint_component_runtime)

add_executable(joint_partial_selection experiments/joint_partial_selection.cpp)
target_link_libraries(joint_partial_selection PRIVATE joint_test_support)
add_dependencies(tests_all joint_partial_selection)

add_executable(joint_sparse_benchmark experiments/joint_sparse_benchmark.cpp)
target_link_libraries(joint_sparse_benchmark PRIVATE joint_test_support)
target_compile_definitions(joint_sparse_benchmark PRIVATE RHBM_GEM_TEST_INSTRUMENTATION)
add_dependencies(tests_all joint_sparse_benchmark)

add_executable(joint_postprocessing_benchmark experiments/joint_postprocessing_benchmark.cpp)
target_link_libraries(joint_postprocessing_benchmark PRIVATE joint_test_support)
add_dependencies(tests_all joint_postprocessing_benchmark)

add_executable(joint_component_benchmark experiments/joint_component_benchmark.cpp)
target_link_libraries(joint_component_benchmark PRIVATE rhbm_gem CompilerFlags)
add_dependencies(tests_all joint_component_benchmark)

if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
    add_library(joint_offline_support STATIC support/JointOfflineAudit.cpp
        support/JointPrecisionAudit.cpp support/JointOfflineTools.cpp
        support/JointLocalEvidence.cpp support/JointDerivativeAudit.cpp)
    target_link_libraries(joint_offline_support PUBLIC joint_test_support)
    add_executable(joint_component_audit experiments/joint_component_audit.cpp)
    target_link_libraries(joint_component_audit PRIVATE joint_offline_support)
    add_executable(joint_validation experiments/joint_validation.cpp)
    target_link_libraries(joint_validation PRIVATE joint_offline_support)
    add_executable(joint_offline_tests core/joint_component/Precision_test.cpp
        core/joint_component/OfflineComponents_test.cpp)
    target_link_libraries(joint_offline_tests PRIVATE joint_offline_support ${GTEST_MAIN_TARGET})
    add_test(NAME joint_component_offline_tests COMMAND joint_offline_tests)
    set_tests_properties(joint_component_offline_tests PROPERTIES LABELS "domain:core;intent:algorithm;joint:offline")
    add_dependencies(tests_all joint_component_audit joint_offline_tests)
endif()

add_rhbm_gtest_group(rhbm_tests_joint_component
    DOMAIN core INTENTS algorithm contract SOURCES ${JOINT_COMPONENT_TEST_SOURCES})
set_property(TEST rhbm_tests_joint_component APPEND PROPERTY LABELS "joint:runtime")

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
    add_test(NAME joint_component_evidence_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_evidence_test.py")
    add_test(NAME joint_component_runtime_runner_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime_test.py")
    add_test(NAME joint_bounded_validation_runner_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_bounded_validation_test.py")
    set_tests_properties(joint_bounded_validation_runner_test PROPERTIES ENVIRONMENT
        "JOINT_BOUNDED_TEST_EXECUTABLE=$<TARGET_FILE:joint_sparse_benchmark>;JOINT_BOUNDED_TEST_BACKEND=${RHBM_GEM_JOINT_SPARSE_BACKEND}")
    add_test(NAME joint_fixed_validation_runner_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_fixed_validation_test.py")
    add_test(NAME joint_search_validation_runner_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_search_validation_test.py")
    add_test(NAME joint_operator_validation_runner_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_operator_validation_test.py")
    add_test(NAME joint_component_regression
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime.py"
            regression --executable $<TARGET_FILE:joint_component_runtime> --work-dir "${PROJECT_BINARY_DIR}/joint-regression")
    add_test(NAME joint_component_physical_smoke
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime.py"
            physical-smoke --executable $<TARGET_FILE:joint_component_runtime>)
    set_tests_properties(joint_component_evidence_test joint_component_runtime_runner_test joint_component_regression joint_component_physical_smoke joint_operator_validation_runner_test joint_search_validation_runner_test joint_fixed_validation_runner_test joint_bounded_validation_runner_test
        PROPERTIES LABELS "domain:integration;intent:algorithm;joint:runtime")
    if(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS)
        add_test(NAME joint_validation_runner_test
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_validation_test.py")
        set_tests_properties(joint_validation_runner_test PROPERTIES LABELS "domain:integration;intent:algorithm;joint:offline")
        add_test(NAME joint_component_offline_regression
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_audit.py"
                --executable $<TARGET_FILE:joint_component_audit> --work-dir "${PROJECT_BINARY_DIR}/joint-offline")
        add_test(NAME joint_component_two_step_regression
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_audit.py"
                --two-step-only --executable $<TARGET_FILE:joint_component_audit> --work-dir "${PROJECT_BINARY_DIR}/joint-two-step")
        set_tests_properties(joint_component_offline_regression joint_component_two_step_regression PROPERTIES LABELS "domain:integration;intent:algorithm;joint:offline")
    endif()
    if(RHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS)
        add_test(NAME joint_component_extended_regression
            COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/joint_component_runtime.py"
                regression --extended --executable $<TARGET_FILE:joint_component_runtime> --work-dir "${PROJECT_BINARY_DIR}/joint-extended")
        set_tests_properties(joint_component_extended_regression PROPERTIES LABELS "domain:integration;intent:algorithm;joint:extended")
    endif()
    add_test(NAME endpoint_refinement_runner_test
        COMMAND "${Python_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/integration/endpoint_refinement_test.py")
    set_tests_properties(endpoint_refinement_runner_test
        PROPERTIES LABELS "domain:integration;intent:algorithm")
    add_test(
        NAME second_stage_audit_parser_test
        COMMAND "${Python_EXECUTABLE}"
            "${PROJECT_SOURCE_DIR}/tests/integration/second_stage_audit_test.py"
    )
    set_tests_properties(
        second_stage_audit_parser_test
        PROPERTIES LABELS "domain:integration;intent:algorithm"
    )


    add_test(
        NAME fold_168_regression_runner_test
        COMMAND "${Python_EXECUTABLE}"
            "${PROJECT_SOURCE_DIR}/tests/integration/fold_168_regression_test.py"
    )
    set_tests_properties(
        fold_168_regression_runner_test
        PROPERTIES LABELS "domain:integration;intent:algorithm"
    )
endif()

if(RHBM_GEM_ENABLE_FOLD_168_REGRESSION)
    if(NOT DEFINED Python_EXECUTABLE OR Python_EXECUTABLE STREQUAL "")
        find_package(Python REQUIRED COMPONENTS Interpreter)
    endif()

    foreach(input_path_variable IN ITEMS
        RHBM_GEM_FOLD_168_MODEL
        RHBM_GEM_FOLD_168_MAP)
        if("${${input_path_variable}}" STREQUAL "" OR
            NOT EXISTS "${${input_path_variable}}" OR
            IS_DIRECTORY "${${input_path_variable}}")
            message(FATAL_ERROR
                "${input_path_variable} must name an existing file when "
                "RHBM_GEM_ENABLE_FOLD_168_REGRESSION=ON")
        endif()
    endforeach()

    add_dependencies(tests_all rhbm_gem_cli)
    set(fold_168_test_name fold_168_simulation_regression)
    set(fold_168_test_script fold_168_regression.py)
    set(fold_168_output_directory fold_168)
    add_test(
        NAME ${fold_168_test_name}
        COMMAND "${Python_EXECUTABLE}"
            "${PROJECT_SOURCE_DIR}/tests/integration/${fold_168_test_script}"
            --executable $<TARGET_FILE:rhbm_gem_cli>
            --model "${RHBM_GEM_FOLD_168_MODEL}"
            --map "${RHBM_GEM_FOLD_168_MAP}"
            --baseline
                "${PROJECT_SOURCE_DIR}/tests/benchmarks/fold_168_simulation_baseline.json"
            --output-dir
                "${PROJECT_BINARY_DIR}/benchmark-results/${fold_168_output_directory}"
    )
    set_tests_properties(
        ${fold_168_test_name}
        PROPERTIES
            LABELS "domain:integration;intent:algorithm;benchmark:external"
            RUN_SERIAL TRUE
            # Full Debug measurements can exceed five minutes; this is not a quality gate.
            TIMEOUT 900
    )
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
            PROPERTIES LABELS "domain:integration;intent:command"
        )
    endif()
endif()

if(DEFINED Python_EXECUTABLE)
    add_dependencies(tests_all rhbm_gem_cli)
    add_test(NAME joint_workflow_cli_smoke_test
        COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/integration/joint_workflow_cli_smoke.py
            $<TARGET_FILE:rhbm_gem_cli>)
    set_tests_properties(joint_workflow_cli_smoke_test PROPERTIES LABELS "domain:integration;intent:command")
    add_test(
        NAME command_scaffold_smoke_test
        COMMAND ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/command_scaffold_smoke.py
            ${PROJECT_SOURCE_DIR}
    )
    set_tests_properties(
        command_scaffold_smoke_test
        PROPERTIES LABELS "domain:integration;intent:contract"
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
