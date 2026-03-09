include("${CMAKE_CURRENT_LIST_DIR}/RHBMGemFetchDeps.cmake")

if(BUILD_TESTING)
    include(CTest)
    enable_testing()
    # Do not export test framework artifacts in normal project installs.
    set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest install targets" FORCE)
    set(INSTALL_GMOCK OFF CACHE BOOL "Disable GoogleMock install targets" FORCE)

    if(USE_SYSTEM_LIBS)
        find_package(GTest QUIET)
        if(GTest_FOUND)
            message(STATUS "Using system GoogleTest")
            set(HAVE_GTEST TRUE)
        else()
            message(STATUS "System GoogleTest not found, will fetch via FetchContent")
            set(HAVE_GTEST FALSE)
        endif()
    else()
        message(STATUS "Forced to use FetchContent for GoogleTest")
        set(HAVE_GTEST FALSE)
    endif()

    if(NOT HAVE_GTEST)
        rhbm_gem_fetch_googletest()
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
