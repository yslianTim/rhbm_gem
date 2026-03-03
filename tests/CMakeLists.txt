set(TEST_SOURCES
    ArrayStats_test.cpp
    AtomSelector_test.cpp
    ChemicalDataHelper_test.cpp
    ChimeraXHelper_test.cpp
    ComponentHelper_test.cpp
    DataObjectManager_test.cpp
    EigenMatrixUtility_test.cpp
    ElectricPotential_test.cpp
    FilePathHelper_test.cpp
    GausLinearTransformHelper_test.cpp
    HRLAlphaTrainer_test.cpp
    HRLDataTransform_test.cpp
    HRLGroupEstimator_test.cpp
    HRLModelAlgorithms_test.cpp
    KDTreeAlgorithm_test.cpp
    Logger_test.cpp
    MapFileIO_test.cpp
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
