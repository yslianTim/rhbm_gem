include("${CMAKE_CURRENT_LIST_DIR}/GuardCommon.cmake")
rhbm_guard_require_project_source_dir()

find_program(RHBM_GEM_PYTHON_EXECUTABLE NAMES python3 python)
if(NOT RHBM_GEM_PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Command manifest guard failed: python interpreter was not found.")
endif()

execute_process(
    COMMAND
        "${RHBM_GEM_PYTHON_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/scripts/developer/check_command_sync.py"
    RESULT_VARIABLE RHBM_GEM_COMMAND_MANIFEST_RESULT
    OUTPUT_VARIABLE RHBM_GEM_COMMAND_MANIFEST_STDOUT
    ERROR_VARIABLE RHBM_GEM_COMMAND_MANIFEST_STDERR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)

if(NOT RHBM_GEM_COMMAND_MANIFEST_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Command manifest guard failed.\n"
        "${RHBM_GEM_COMMAND_MANIFEST_STDOUT}\n"
        "${RHBM_GEM_COMMAND_MANIFEST_STDERR}")
endif()
