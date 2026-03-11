include(CMakeParseArguments)

function(rhbm_guard_require_project_source_dir)
    if(NOT DEFINED PROJECT_SOURCE_DIR)
        message(FATAL_ERROR "PROJECT_SOURCE_DIR is required")
    endif()
endfunction()

function(rhbm_guard_require_project_binary_dir)
    if(NOT DEFINED PROJECT_BINARY_DIR)
        message(FATAL_ERROR "PROJECT_BINARY_DIR is required")
    endif()
endfunction()

function(rhbm_guard_find_git out_git_executable)
    find_program(_rhbm_guard_git_executable NAMES git)
    if(NOT _rhbm_guard_git_executable)
        message(FATAL_ERROR "git executable was not found")
    endif()
    set(${out_git_executable} "${_rhbm_guard_git_executable}" PARENT_SCOPE)
endfunction()

function(rhbm_guard_get_tracked_files git_executable out_tracked_files)
    if("${git_executable}" STREQUAL "")
        message(FATAL_ERROR "git executable path is required")
    endif()

    execute_process(
        COMMAND "${git_executable}" -C "${PROJECT_SOURCE_DIR}" ls-files
        OUTPUT_VARIABLE _rhbm_guard_tracked_files_raw
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _rhbm_guard_git_ls_result
    )
    if(NOT _rhbm_guard_git_ls_result EQUAL 0)
        message(FATAL_ERROR "Failed to query tracked files via git ls-files")
    endif()

    set(_rhbm_guard_tracked_files)
    if(NOT _rhbm_guard_tracked_files_raw STREQUAL "")
        string(REPLACE "\n" ";" _rhbm_guard_tracked_files "${_rhbm_guard_tracked_files_raw}")
    endif()
    set(${out_tracked_files} "${_rhbm_guard_tracked_files}" PARENT_SCOPE)
endfunction()

function(rhbm_guard_collect_relative_files out_files)
    set(options)
    set(oneValueArgs BASE_DIR)
    set(multiValueArgs ROOTS GLOBS)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT ARG_BASE_DIR)
        message(FATAL_ERROR "rhbm_guard_collect_relative_files requires BASE_DIR")
    endif()
    if(NOT ARG_ROOTS)
        message(FATAL_ERROR "rhbm_guard_collect_relative_files requires ROOTS")
    endif()
    if(NOT ARG_GLOBS)
        message(FATAL_ERROR "rhbm_guard_collect_relative_files requires GLOBS")
    endif()

    set(_rhbm_guard_files)
    foreach(_rhbm_guard_root IN LISTS ARG_ROOTS)
        if(NOT EXISTS "${_rhbm_guard_root}")
            continue()
        endif()

        foreach(_rhbm_guard_glob IN LISTS ARG_GLOBS)
            file(GLOB_RECURSE _rhbm_guard_root_files
                RELATIVE "${ARG_BASE_DIR}"
                "${_rhbm_guard_root}/${_rhbm_guard_glob}")
            list(APPEND _rhbm_guard_files ${_rhbm_guard_root_files})
        endforeach()
    endforeach()

    if(_rhbm_guard_files)
        list(REMOVE_DUPLICATES _rhbm_guard_files)
    endif()
    set(${out_files} "${_rhbm_guard_files}" PARENT_SCOPE)
endfunction()
