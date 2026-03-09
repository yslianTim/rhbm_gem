add_rhbm_gtest_target(rhbm_tests_integration_command
    DOMAIN integration
    INTENTS command
    SOURCES
        integration/GuiCommandExecutorPipeline_test.cpp
)

add_test(
    NAME logging_style_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/LoggingStyleGuard.cmake
)
set_tests_properties(logging_style_guard PROPERTIES
    LABELS "domain:integration;intent:style"
)

add_test(
    NAME parameterless_void_style_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/ParameterlessVoidStyleGuard.cmake
)
set_tests_properties(parameterless_void_style_guard PROPERTIES
    LABELS "domain:integration;intent:style"
)

add_test(
    NAME repository_hygiene_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/RepositoryHygieneGuard.cmake
)
set_tests_properties(repository_hygiene_guard PROPERTIES
    LABELS "domain:integration;intent:style"
)

add_test(
    NAME fixture_tracking_guard
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/TestFixtureTrackingGuard.cmake
)
set_tests_properties(fixture_tracking_guard PROPERTIES
    LABELS "domain:integration;intent:style"
)

if(BUILD_PYTHON_BINDINGS AND TARGET rhbm_gem_module AND DEFINED Python_EXECUTABLE)
    add_dependencies(tests_all rhbm_gem_module)

    add_test(
        NAME bindings_smoke_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_smoke.py
    )
    add_test(
        NAME bindings_invalid_enum_validation_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_invalid_enum_validation.py
    )

    set_tests_properties(
        bindings_smoke_test
        bindings_invalid_enum_validation_test
        PROPERTIES LABELS "domain:integration;intent:bindings"
    )
endif()
