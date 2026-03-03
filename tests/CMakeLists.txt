set(TEST_SOURCES
    Application_test.cpp
    ArrayStats_test.cpp
    AtomSelector_test.cpp
    BindingSurfaceParity_test.cpp
    BuiltInCommandCatalog_test.cpp
    ChemicalDataHelper_test.cpp
    ChimeraXHelper_test.cpp
    CommandBase_test.cpp
    CommandExecutionContract_test.cpp
    CommandRegistry_test.cpp
    CommandSurface_test.cpp
    CommandValidationIssue_test.cpp
    ComponentHelper_test.cpp
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
            "import rhbm_gem_module as m; assert hasattr(m, 'LogLevel'); assert hasattr(m, 'ValidationPhase'); assert hasattr(m, 'ValidationIssue'); assert hasattr(m.PrinterType, 'ATOM_OUTLIER'); assert hasattr(m, 'PotentialAnalysisCommand'); assert hasattr(m, 'PotentialDisplayCommand'); assert hasattr(m, 'ResultDumpCommand'); assert hasattr(m, 'MapSimulationCommand'); assert not hasattr(m, 'MapVisualizationCommand'); assert not hasattr(m, 'PositionEstimationCommand'); assert not hasattr(m, 'HRLModelTestCommand'); cmd = m.MapSimulationCommand(); assert hasattr(cmd, 'PrepareForExecution'); assert hasattr(cmd, 'HasValidationErrors'); assert hasattr(cmd, 'GetValidationIssues')"
    )
endif()
