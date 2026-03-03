set(TEST_SOURCES
    Application_test.cpp
    ArrayStats_test.cpp
    AtomSelector_test.cpp
    BindingSurfaceParity_test.cpp
    BuiltInCommandCatalog_test.cpp
    ChemicalDataHelper_test.cpp
    ChimeraXHelper_test.cpp
    CommandBase_test.cpp
    CommandCLIHelper_test.cpp
    CommandExecutionContract_test.cpp
    CommandOptionHelper_test.cpp
    CommandRegistrationPolicy_test.cpp
    CommandRegistry_test.cpp
    CommandSurface_test.cpp
    CommandValidationIssue_test.cpp
    ComponentHelper_test.cpp
    ContributorGuidanceSync_test.cpp
    DataObjectManager_test.cpp
    DocsSync_test.cpp
    EigenMatrixUtility_test.cpp
    ElectricPotential_test.cpp
    EnumOptionTraits_test.cpp
    FilePathHelper_test.cpp
    GausLinearTransformHelper_test.cpp
    HRLAlphaTrainer_test.cpp
    HRLDataTransform_test.cpp
    HRLGroupEstimator_test.cpp
    HRLModelAlgorithms_test.cpp
    KDTreeAlgorithm_test.cpp
    Logger_test.cpp
    MapFileIO_test.cpp
    MapSimulationCommand_test.cpp
    MapVisualizationCommand_test.cpp
    PositionEstimationCommand_test.cpp
    PreparedStateInvalidation_test.cpp
    PotentialDisplayCommand_test.cpp
    ResultDumpCommand_test.cpp
    ROOTHelper_test.cpp
    SphereSampler_test.cpp
    StringHelper_test.cpp
)

set(TEST_EXEC_NAME unit_tests)
add_executable(${TEST_EXEC_NAME} ${TEST_SOURCES})

target_link_libraries(${TEST_EXEC_NAME} PRIVATE
    core
    utils
    data
    CompilerFlags
)

if(TARGET GTest::gtest_main)
    set(GTEST_MAIN_TARGET GTest::gtest_main)
else()
    set(GTEST_MAIN_TARGET gtest_main)
endif()
target_link_libraries(${TEST_EXEC_NAME} PRIVATE ${GTEST_MAIN_TARGET})

target_include_directories(${TEST_EXEC_NAME} PRIVATE "${PROJECT_BINARY_DIR}")

set_target_properties(${TEST_EXEC_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

gtest_discover_tests(${TEST_EXEC_NAME}
    DISCOVERY_MODE PRE_TEST
)

add_executable(command_catalog_binding_names_helper
    CommandCatalogBindingNamesHelper.cpp
)
target_link_libraries(command_catalog_binding_names_helper PRIVATE
    core
    CompilerFlags
)
add_dependencies(${TEST_EXEC_NAME} command_catalog_binding_names_helper)
set_target_properties(command_catalog_binding_names_helper PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

add_test(
    NAME logging_style_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/LoggingStyleGuard.cmake
)

add_test(
    NAME parameterless_void_style_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/ParameterlessVoidStyleGuard.cmake
)

if(BUILD_PYTHON_BINDINGS AND TARGET rhbm_gem_module AND DEFINED Python_EXECUTABLE)
    add_test(
        NAME bindings_smoke_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            -c
            "import subprocess, sys; import rhbm_gem_module as m; public_names = [line.strip() for line in subprocess.check_output([sys.argv[1]], text=True).splitlines() if line.strip()]; assert public_names; assert hasattr(m, 'LogLevel'); assert hasattr(m, 'ValidationPhase'); assert hasattr(m, 'ValidationIssue'); assert hasattr(m.PrinterType, 'ATOM_OUTLIER'); assert all(hasattr(m, name) for name in public_names); assert all(all(hasattr(getattr(m, name)(), method) for method in ('PrepareForExecution', 'HasValidationErrors', 'GetValidationIssues')) for name in public_names)"
            "$<TARGET_FILE:command_catalog_binding_names_helper>"
    )
endif()
