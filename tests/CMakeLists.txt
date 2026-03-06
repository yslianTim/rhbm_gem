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
        core
        utils
        data
        CompilerFlags
        ${GTEST_MAIN_TARGET}
    )

    target_include_directories(${target_name} PRIVATE
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/src/core"
        "${PROJECT_SOURCE_DIR}/tests"
    )

    set_target_properties(${target_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
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

set(CORE_COMMAND_TEST_SOURCES
    core/command/Application_test.cpp
    core/command/CommandBase_test.cpp
    core/command/CommandCLIHelper_test.cpp
    core/command/CommandOptionHelper_test.cpp
    core/command/CommandScalarValidationHelper_test.cpp
    core/command/CommandValidationIssue_test.cpp
    core/command/HRLModelTestCommand_test.cpp
    core/command/MapSimulationCommand_test.cpp
    core/command/MapVisualizationCommand_test.cpp
    core/command/PositionEstimationCommand_test.cpp
    core/command/PotentialAnalysisCommand_test.cpp
    core/command/PotentialDisplayCommand_test.cpp
    core/command/ResultDumpCommand_test.cpp
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
        core/contract/BuiltInCommandBindingName_test.cpp
        core/contract/BuiltInCommandCatalog_test.cpp
        core/contract/CommandCommonOptions_test.cpp
        core/contract/CommandDescriptorShape_test.cpp
        core/contract/CommandExecutionContract_test.cpp
        core/contract/DocsSync_test.cpp
        core/contract/EnumOptionTraits_test.cpp
        core/contract/PublicHeaderSurface_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_data_io
    DOMAIN data
    INTENTS io contract
    SOURCES
        data/io/DataObjectDispatchIterationArchitecture_test.cpp
        data/io/DataObjectManager_test.cpp
        data/io/MapFileIO_test.cpp
        data/schema/DataObjectIOContract_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_data_schema
    DOMAIN data
    INTENTS schema
    SOURCES
        data/schema/DataObjectSchemaBootstrap_test.cpp
        data/schema/DataObjectSchemaValidation_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_data_migration
    DOMAIN data
    INTENTS migration schema
    SOURCES
        data/schema/DataObjectSchemaMigration_test.cpp
)

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

add_test(
    NAME logging_style_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/LoggingStyleGuard.cmake
)
set_tests_properties(logging_style_guard PROPERTIES
    LABELS "domain:integration;intent:style"
)

add_test(
    NAME parameterless_void_style_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/ParameterlessVoidStyleGuard.cmake
)
set_tests_properties(parameterless_void_style_guard PROPERTIES
    LABELS "domain:integration;intent:style"
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

    set_tests_properties(
        bindings_smoke_test
        bindings_invalid_enum_validation_test
        PROPERTIES LABELS "domain:integration;intent:bindings"
    )
endif()
