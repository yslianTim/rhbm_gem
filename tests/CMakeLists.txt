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
        gui_adapter
        CompilerFlags
        ${GTEST_MAIN_TARGET}
    )

    target_include_directories(${target_name} PRIVATE
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/src/core"
        "${PROJECT_SOURCE_DIR}/src/core/internal"
        "${PROJECT_SOURCE_DIR}/src/data"
        "${PROJECT_SOURCE_DIR}/src/data/internal"
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

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/core_tests.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/data_tests.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils_tests.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/integration_tests.cmake)
