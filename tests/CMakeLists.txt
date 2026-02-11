file(GLOB TEST_SOURCES CONFIGURE_DEPENDS *.cpp)

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

gtest_discover_tests(${TEST_EXEC_NAME})
