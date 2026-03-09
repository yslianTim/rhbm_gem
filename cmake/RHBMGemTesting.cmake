include("${CMAKE_CURRENT_LIST_DIR}/RHBMGemThirdPartyGuard.cmake")

if(BUILD_TESTING)
    include(CTest)
    enable_testing()
    # Do not export test framework artifacts in normal project installs.
    set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest install targets" FORCE)
    set(INSTALL_GMOCK OFF CACHE BOOL "Disable GoogleMock install targets" FORCE)

    rhbm_gem_require_bundled_dependency("GoogleTest" "third_party/googletest/CMakeLists.txt")
    message(STATUS "Using bundled GoogleTest")
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    add_subdirectory(third_party/googletest)
    set(GTEST_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest/googletest")

    if(DEFINED GTEST_SOURCE_DIR)
        list(APPEND CMAKE_MODULE_PATH "${GTEST_SOURCE_DIR}/cmake")
    endif()

    include(GoogleTest)
    add_subdirectory(tests)

    if(ENABLE_COVERAGE)
        find_program(GCOV_EXECUTABLE NAMES gcov)
        if(NOT GCOV_EXECUTABLE)
            message(FATAL_ERROR "ENABLE_COVERAGE is ON but gcov was not found in PATH")
        endif()

        set(COVERAGE_OUTPUT_DIR "${PROJECT_BINARY_DIR}/coverage")
        add_custom_target(coverage
            COMMAND ${CMAKE_COMMAND}
                -DCOVERAGE_ACTION=clean
                -DCOVERAGE_BUILD_DIR:PATH=${PROJECT_BINARY_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GcovCoverageReport.cmake
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
            COMMAND ${CMAKE_COMMAND}
                -DCOVERAGE_ACTION=report
                -DCOVERAGE_BUILD_DIR:PATH=${PROJECT_BINARY_DIR}
                -DCOVERAGE_PROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
                -DCOVERAGE_OUTPUT_DIR:PATH=${COVERAGE_OUTPUT_DIR}
                -DCOVERAGE_INCLUDE_TESTS:BOOL=${COVERAGE_INCLUDE_TESTS}
                -DGCOV_EXECUTABLE:FILEPATH=${GCOV_EXECUTABLE}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GcovCoverageReport.cmake
            BYPRODUCTS
                "${COVERAGE_OUTPUT_DIR}/coverage_summary.txt"
                "${COVERAGE_OUTPUT_DIR}/coverage_detail.txt"
            DEPENDS tests_all
            USES_TERMINAL
            COMMENT "Running tests and generating gcov text coverage reports"
        )
    endif()
endif()

