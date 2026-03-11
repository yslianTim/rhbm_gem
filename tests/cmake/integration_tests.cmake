add_rhbm_gtest_target(rhbm_tests_integration_command
    DOMAIN integration
    INTENTS command
    SOURCES
        integration/GuiCommandExecutorPipeline_test.cpp
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
    add_test(
        NAME examples_quickstart_smoke_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/examples_quickstart_smoke.py
    )
    add_test(
        NAME examples_end_to_end_smoke_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/examples_end_to_end_smoke.py
    )
    add_test(
        NAME bindings_common_setters_surface_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_common_setters_surface.py
    )
    add_test(
        NAME bindings_surface_snapshot_test
        COMMAND ${CMAKE_COMMAND} -E env
            "PYTHONPATH=$<TARGET_FILE_DIR:rhbm_gem_module>"
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/integration/bindings_surface_snapshot.py
    )

    set_tests_properties(
        bindings_smoke_test
        bindings_invalid_enum_validation_test
        examples_quickstart_smoke_test
        examples_end_to_end_smoke_test
        bindings_common_setters_surface_test
        bindings_surface_snapshot_test
        PROPERTIES LABELS "domain:integration;intent:bindings"
    )
endif()
