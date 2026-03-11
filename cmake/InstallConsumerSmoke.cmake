include("${PROJECT_SOURCE_DIR}/cmake/GuardCommon.cmake")
rhbm_guard_require_project_source_dir()
rhbm_guard_require_project_binary_dir()

set(SMOKE_ROOT "${PROJECT_BINARY_DIR}/consumer_smoke")
set(INSTALL_PREFIX "${SMOKE_ROOT}/install")
set(CONSUMER_BUILD_DIR "${SMOKE_ROOT}/build")
set(CONSUMER_SOURCE_DIR "${PROJECT_SOURCE_DIR}/tests/cmake/consumer_smoke")

file(REMOVE_RECURSE "${SMOKE_ROOT}")
file(MAKE_DIRECTORY "${SMOKE_ROOT}")

execute_process(
    COMMAND "${CMAKE_COMMAND}" --install "${PROJECT_BINARY_DIR}" --prefix "${INSTALL_PREFIX}"
    RESULT_VARIABLE INSTALL_RESULT
    OUTPUT_VARIABLE INSTALL_STDOUT
    ERROR_VARIABLE INSTALL_STDERR
)
if(NOT INSTALL_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Install consumer smoke test failed at install step.\n"
        "stdout:\n${INSTALL_STDOUT}\n"
        "stderr:\n${INSTALL_STDERR}\n")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}"
        -S "${CONSUMER_SOURCE_DIR}"
        -B "${CONSUMER_BUILD_DIR}"
        -DCMAKE_PREFIX_PATH:PATH=${INSTALL_PREFIX}
        -DRHBM_GEM_DIR:PATH=${INSTALL_PREFIX}/lib/cmake/RHBM_GEM
        -DCMAKE_BUILD_TYPE=Release
    RESULT_VARIABLE CONFIGURE_RESULT
    OUTPUT_VARIABLE CONFIGURE_STDOUT
    ERROR_VARIABLE CONFIGURE_STDERR
)
if(NOT CONFIGURE_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Install consumer smoke test failed at configure step.\n"
        "stdout:\n${CONFIGURE_STDOUT}\n"
        "stderr:\n${CONFIGURE_STDERR}\n")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${CONSUMER_BUILD_DIR}"
    RESULT_VARIABLE BUILD_RESULT
    OUTPUT_VARIABLE BUILD_STDOUT
    ERROR_VARIABLE BUILD_STDERR
)
if(NOT BUILD_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Install consumer smoke test failed at build step.\n"
        "stdout:\n${BUILD_STDOUT}\n"
        "stderr:\n${BUILD_STDERR}\n")
endif()
