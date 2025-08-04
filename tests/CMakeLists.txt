find_package(GTest CONFIG REQUIRED)
include(GoogleTest)

file(GLOB TEST_SOURCES *.cpp)

set(TEST_BINARY ${CMAKE_PROJECT_NAME}_test)

add_executable(${TEST_BINARY} ${TEST_SOURCES})

target_include_directories(${TEST_BINARY} PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${TEST_BINARY} PRIVATE
    core
    utils
    data
    CompilerFlags
    GTest::gtest_main
)

set_target_properties(${TEST_BINARY} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

gtest_discover_tests(${TEST_BINARY})
