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
    core/command/CommandValidationScenarios_test.cpp
    core/command/CommandWorkflowScenarios_test.cpp
)

set(CORE_COMMAND_TEST_SOURCES
    core/command/CommandBaseLifecycle_test.cpp
    core/command/CommandValidationHelpers_test.cpp
    core/command/MapSampling_test.cpp
    ${RHBM_GEM_CORE_COMMAND_TEST_SOURCES}
)

set(RHBM_GEM_EXPERIMENTAL_TEST_SOURCES)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)
    list(APPEND CORE_COMMAND_TEST_SOURCES
        ${RHBM_GEM_EXPERIMENTAL_TEST_SOURCES}
    )
endif()

set(CORE_CONTRACT_TEST_SOURCES
    core/contract/CommandCatalog_test.cpp
    core/contract/CommandExecutionContract_test.cpp
    core/contract/DocumentationLayout_test.cpp
    core/contract/PotentialPlotBuilder_test.cpp
    core/contract/PublicHeaderSurface_test.cpp
)

set(DATA_RUNTIME_TEST_SOURCES
    data/DataObjectRuntime_test.cpp
)

set(DATA_SCHEMA_TEST_SOURCES
    data/DataObjectSchema_test.cpp
)

set(UTILS_MATH_TEST_SOURCES
    utils/math/ArrayStats_test.cpp
    utils/math/EigenMatrixUtility_test.cpp
    utils/math/ElectricPotential_test.cpp
    utils/math/GausLinearTransformHelper_test.cpp
    utils/math/KDTreeAlgorithm_test.cpp
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
    utils/domain/ROOTHelper_test.cpp
    utils/domain/StringHelperFormatting_test.cpp
    utils/domain/StringHelperParsing_test.cpp
    utils/domain/StringHelperTokenization_test.cpp
)

set(UTILS_HRL_TEST_SOURCES
    utils/hrl/HRLAlphaTrainer_test.cpp
    utils/hrl/HRLDataTransform_test.cpp
    utils/hrl/HRLGroupEstimator_test.cpp
    utils/hrl/HRLModelAlgorithms_test.cpp
)

set(INTEGRATION_COMMAND_TEST_SOURCES
    integration/CommandApiPipeline_test.cpp
)

set(RHBM_GEM_ALL_GTEST_SOURCES
    ${CORE_COMMAND_TEST_SOURCES}
    ${CORE_CONTRACT_TEST_SOURCES}
    ${DATA_RUNTIME_TEST_SOURCES}
    ${DATA_SCHEMA_TEST_SOURCES}
    ${UTILS_MATH_TEST_SOURCES}
    ${UTILS_DOMAIN_TEST_SOURCES}
    ${UTILS_HRL_TEST_SOURCES}
    ${INTEGRATION_COMMAND_TEST_SOURCES}
)

add_executable(rhbm_tests ${RHBM_GEM_ALL_GTEST_SOURCES})
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
    "${PROJECT_SOURCE_DIR}/src/core"
    "${PROJECT_SOURCE_DIR}/src/data"
    "${PROJECT_SOURCE_DIR}/tests"
)

add_dependencies(tests_all rhbm_tests)

add_rhbm_gtest_group(rhbm_tests_core_command
    DOMAIN core
    INTENTS command validation
    SOURCES ${CORE_COMMAND_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_core_contract
    DOMAIN core
    INTENTS contract
    SOURCES ${CORE_CONTRACT_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_data_runtime
    DOMAIN data
    INTENTS io contract
    SOURCES ${DATA_RUNTIME_TEST_SOURCES}
)

add_rhbm_gtest_group(rhbm_tests_data_schema
    DOMAIN data
    INTENTS schema
    SOURCES ${DATA_SCHEMA_TEST_SOURCES}
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
