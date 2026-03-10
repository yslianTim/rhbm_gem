if(BUILD_TESTING)
    include(CTest)

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

add_custom_target(lint_repo
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -DPROJECT_BINARY_DIR:PATH=${PROJECT_BINARY_DIR}
        -DRHBM_GEM_INCLUDE_INSTALL_SMOKE:BOOL=OFF
        -P ${PROJECT_SOURCE_DIR}/cmake/RepoGuards.cmake
    USES_TERMINAL
    COMMENT "Running repository lint checks (style/structure/hygiene)"
)

add_custom_target(lint_install_smoke
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -DPROJECT_BINARY_DIR:PATH=${PROJECT_BINARY_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/InstallConsumerSmoke.cmake
    USES_TERMINAL
    COMMENT "Running install consumer smoke test"
)

function(rhbm_gem_add_lint_build_dependencies target_name)
    if(TARGET rhbm_gem)
        add_dependencies(${target_name} rhbm_gem)
    endif()
    if(TARGET rhbm_gem_cli)
        add_dependencies(${target_name} rhbm_gem_cli)
    endif()
    if(TARGET ${PROJECT_NAME}_gui)
        add_dependencies(${target_name} ${PROJECT_NAME}_gui)
    endif()
    if(TARGET rhbm_gem_module)
        add_dependencies(${target_name} rhbm_gem_module)
    endif()
endfunction()

rhbm_gem_add_lint_build_dependencies(lint_install_smoke)

add_custom_target(lint_all
    DEPENDS
        lint_repo
        lint_install_smoke
)
