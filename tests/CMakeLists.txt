if(TARGET GTest::gtest_main)
    set(GTEST_MAIN_TARGET GTest::gtest_main)
else()
    set(GTEST_MAIN_TARGET gtest_main)
endif()

add_custom_target(tests_all)

function(add_rhbm_gtest_target target_name)
    set(options)
    set(oneValueArgs DOMAIN)
    set(multiValueArgs INTENTS SOURCES)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT ARG_DOMAIN)
        message(FATAL_ERROR "add_rhbm_gtest_target requires DOMAIN for target ${target_name}")
    endif()
    if(NOT ARG_SOURCES)
        message(FATAL_ERROR "add_rhbm_gtest_target requires SOURCES for target ${target_name}")
    endif()

    add_executable(${target_name} ${ARG_SOURCES})

    target_link_libraries(${target_name} PRIVATE
        rhbm_gem
        CompilerFlags
        ${GTEST_MAIN_TARGET}
    )

    target_include_directories(${target_name} PRIVATE
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/src/core"
        "${PROJECT_SOURCE_DIR}/src/core/internal"
        "${PROJECT_SOURCE_DIR}/src/gui/internal"
        "${PROJECT_SOURCE_DIR}/src/data"
        "${PROJECT_SOURCE_DIR}/src/data/internal"
        "${PROJECT_SOURCE_DIR}/tests"
    )

    set(test_labels "domain:${ARG_DOMAIN}")
    foreach(intent IN LISTS ARG_INTENTS)
        list(APPEND test_labels "intent:${intent}")
    endforeach()
    string(REPLACE ";" "\\;" escaped_test_labels "${test_labels}")

    gtest_discover_tests(${target_name}
        DISCOVERY_MODE PRE_TEST
        PROPERTIES LABELS "${escaped_test_labels}"
    )

    add_dependencies(tests_all ${target_name})
endfunction()

file(GLOB RHBM_GEM_CORE_COMMAND_TEST_SOURCES CONFIGURE_DEPENDS
    RELATIVE "${CMAKE_CURRENT_LIST_DIR}"
    "${CMAKE_CURRENT_LIST_DIR}/core/command/*Command_test.cpp"
)

# core
set(CORE_COMMAND_TEST_SOURCES
    core/command/CommandBase_test.cpp
    core/command/CommandOptionHelper_test.cpp
    core/command/CommandScalarValidationHelper_test.cpp
    core/command/CommandValidationIssue_test.cpp
    core/command/CommandApi_test.cpp
    ${RHBM_GEM_CORE_COMMAND_TEST_SOURCES}
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    list(APPEND CORE_COMMAND_TEST_SOURCES
        core/command/PotentialAnalysisExperimentalBondGate_test.cpp
    )
endif()

add_rhbm_gtest_target(rhbm_tests_core_command
    DOMAIN core
    INTENTS command validation
    SOURCES ${CORE_COMMAND_TEST_SOURCES}
)

add_rhbm_gtest_target(rhbm_tests_core_contract
    DOMAIN core
    INTENTS contract
    SOURCES
        core/contract/CommandCatalog_test.cpp
        core/contract/CommandExecutionContract_test.cpp
        core/contract/EnumOptionTraits_test.cpp
        core/contract/PotentialPlotBuilder_test.cpp
        core/contract/PublicHeaderSurface_test.cpp
)

# data
add_rhbm_gtest_target(rhbm_tests_data_runtime
    DOMAIN data
    INTENTS io contract
    SOURCES
        data/DataObjectRuntime_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_data_schema
    DOMAIN data
    INTENTS schema
    SOURCES
        data/DataObjectSchema_test.cpp
)

if(RHBM_GEM_LEGACY_V1_SUPPORT)
    target_compile_definitions(rhbm_tests_data_schema PRIVATE RHBM_GEM_LEGACY_V1_SUPPORT)
endif()

# utils
add_rhbm_gtest_target(rhbm_tests_utils_math
    DOMAIN utils
    INTENTS algorithm
    SOURCES
        utils/math/ArrayStats_test.cpp
        utils/math/EigenMatrixUtility_test.cpp
        utils/math/ElectricPotential_test.cpp
        utils/math/GausLinearTransformHelper_test.cpp
        utils/math/KDTreeAlgorithm_test.cpp
        utils/math/SphereSampler_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_utils_domain
    DOMAIN utils
    INTENTS validation
    SOURCES
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

add_rhbm_gtest_target(rhbm_tests_utils_hrl
    DOMAIN utils
    INTENTS algorithm
    SOURCES
        utils/hrl/HRLAlphaTrainer_test.cpp
        utils/hrl/HRLDataTransform_test.cpp
        utils/hrl/HRLGroupEstimator_test.cpp
        utils/hrl/HRLModelAlgorithms_test.cpp
)

# integration
add_rhbm_gtest_target(rhbm_tests_integration_command
    DOMAIN integration
    INTENTS command
    SOURCES
        integration/CommandApiPipeline_test.cpp
)

if(BUILD_PYTHON_BINDINGS AND TARGET rhbm_gem_module AND DEFINED Python_EXECUTABLE)
    add_dependencies(tests_all rhbm_gem_module)

    add_test(
        NAME bindings_smoke_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_smoke.py
    )
    add_test(
        NAME bindings_invalid_enum_validation_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_invalid_enum_validation.py
    )
    add_test(
        NAME examples_quickstart_smoke_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/examples_quickstart_smoke.py
    )
    add_test(
        NAME examples_end_to_end_smoke_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/examples_end_to_end_smoke.py
    )
    add_test(
        NAME bindings_common_setters_surface_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_common_setters_surface.py
    )
    add_test(
        NAME bindings_surface_snapshot_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_surface_snapshot.py
    )

    set_tests_properties(
        bindings_smoke_test
        bindings_invalid_enum_validation_test
        examples_quickstart_smoke_test
        examples_end_to_end_smoke_test
        bindings_common_setters_surface_test
        bindings_surface_snapshot_test
        PROPERTIES LABELS "domain:integration;intent:bindings"
    )
endif()
