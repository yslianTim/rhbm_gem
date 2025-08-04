## Top-level CMake configuration for the test suite
include(GoogleTest)

function(add_project_test name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE
        core
        utils
        data
        CompilerFlags
        GTest::gtest_main
        GTest::gmock
    )
    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    )
    gtest_discover_tests(${name})
endfunction()

# Automatically add each test subdirectory containing a CMakeLists.txt
# New test modules can simply create a subdirectory with its own CMakeLists.txt
file(GLOB CONFIGURE_DEPENDS CHILD_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} */CMakeLists.txt)
foreach(child ${CHILD_DIRS})
    get_filename_component(dir ${child} DIRECTORY)
    add_subdirectory(${dir})
endforeach()
