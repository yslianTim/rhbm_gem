## Top-level CMake configuration for the test suite
include(GoogleTest OPTIONAL RESULT_VARIABLE GTEST_CMAKE)
if(NOT GTEST_CMAKE OR GTEST_CMAKE STREQUAL "GTEST_CMAKE-NOTFOUND")
    message(FATAL_ERROR "GoogleTest.cmake not found. Ensure GTest is installed and GTest_DIR is set correctly.")
endif()

file(GLOB TEST_SOURCES *.cpp)

set(TEST_BINARY ${CMAKE_PROJECT_NAME}_test)

add_executable(${TEST_BINARY} ${TEST_SOURCES})

target_link_libraries(${TEST_BINARY} PRIVATE
    core
    utils
    data
    CompilerFlags
    GTest::gtest
    GTest::gtest_main
    GTest::gmock_main
)

set_target_properties(${TEST_BINARY} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

gtest_discover_tests(${TEST_BINARY})
