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

function(rhbm_guard_check_logging_style)
    set(_logging_guard_sources
        "${PROJECT_SOURCE_DIR}/src"
        "${PROJECT_SOURCE_DIR}/include"
    )

    set(_stream_pattern "\\bstd::cout\\b|\\bstd::cerr\\b|\\bprintf\\(|\\bfprintf\\(|\\bputs\\(|\\bstd::clog\\b")
    set(_chatty_debug_pattern "Logger::Log\\(LogLevel::Debug,.*called")
    set(_allowed_stream_file "${PROJECT_SOURCE_DIR}/src/utils/domain/Logger.cpp")

    rhbm_guard_collect_relative_files(_logging_guard_files
        BASE_DIR "${PROJECT_SOURCE_DIR}"
        ROOTS ${_logging_guard_sources}
        GLOBS "*"
    )

    foreach(_relative_file IN LISTS _logging_guard_files)
        set(_file "${PROJECT_SOURCE_DIR}/${_relative_file}")
        if(IS_DIRECTORY "${_file}")
            continue()
        endif()

        file(READ "${_file}" _content)

        if(NOT _file STREQUAL _allowed_stream_file)
            string(REGEX MATCH "${_stream_pattern}" _stream_match "${_content}")
            if(_stream_match)
                message(FATAL_ERROR
                    "Logging style guard failed: direct stream output found in ${_relative_file}")
            endif()
        endif()

        string(REGEX MATCH "${_chatty_debug_pattern}" _chatty_match "${_content}")
        if(_chatty_match)
            message(FATAL_ERROR
                "Logging style guard failed: chatty debug trace found in ${_relative_file}")
        endif()
    endforeach()
endfunction()

function(rhbm_guard_check_parameterless_void)
    set(_parameterless_void_guard_sources
        "${PROJECT_SOURCE_DIR}/src"
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/tests"
        "${PROJECT_SOURCE_DIR}/bindings"
    )

    set(_parameterless_void_signature_pattern
        "\\([ \t]*void[ \t]*\\)[ \t]*(const|override|final|noexcept|requires|[:;{=]|$)")

    rhbm_guard_collect_relative_files(_parameterless_void_guard_files
        BASE_DIR "${PROJECT_SOURCE_DIR}"
        ROOTS ${_parameterless_void_guard_sources}
        GLOBS "*"
    )

    foreach(_relative_file IN LISTS _parameterless_void_guard_files)
        set(_file "${PROJECT_SOURCE_DIR}/${_relative_file}")
        if(IS_DIRECTORY "${_file}")
            continue()
        endif()

        file(STRINGS "${_file}" _parameterless_void_guard_lines)
        set(_line_number 0)
        foreach(_line IN LISTS _parameterless_void_guard_lines)
            math(EXPR _line_number "${_line_number} + 1")
            string(REGEX MATCH "${_parameterless_void_signature_pattern}" _void_match "${_line}")
            if(_void_match)
                message(FATAL_ERROR
                    "Parameterless void style guard failed: signature-style '(void)' found in "
                    "${_relative_file}:${_line_number}")
            endif()
        endforeach()
    endforeach()
endfunction()

function(rhbm_guard_check_hygiene)
    rhbm_guard_find_git(GIT_EXECUTABLE)
    rhbm_guard_get_tracked_files("${GIT_EXECUTABLE}" TRACKED_FILE_LIST)

    set(HYGIENE_VIOLATIONS)
    foreach(TRACKED_FILE IN LISTS TRACKED_FILE_LIST)
        if(TRACKED_FILE MATCHES "^CMakeFiles/")
            list(APPEND HYGIENE_VIOLATIONS "${TRACKED_FILE} (generated CMake metadata must not be tracked)")
        elseif(TRACKED_FILE MATCHES "^build/")
            list(APPEND HYGIENE_VIOLATIONS "${TRACKED_FILE} (build artifacts must not be tracked)")
        elseif(TRACKED_FILE MATCHES "(^|/)\\.DS_Store$")
            list(APPEND HYGIENE_VIOLATIONS "${TRACKED_FILE} (OS metadata must not be tracked)")
        elseif(TRACKED_FILE MATCHES "^[^/]+\\.sqlite$")
            list(APPEND HYGIENE_VIOLATIONS "${TRACKED_FILE} (root-level sqlite files must be moved to data/local/)")
        elseif(TRACKED_FILE MATCHES "^[^/]+\\.db$")
            list(APPEND HYGIENE_VIOLATIONS "${TRACKED_FILE} (root-level db files must be moved to data/local/)")
        elseif(TRACKED_FILE MATCHES "^[^/]+\\.sh$")
            list(APPEND HYGIENE_VIOLATIONS "${TRACKED_FILE} (root-level scripts must move under scripts/)")
        elseif(TRACKED_FILE MATCHES "^data/local/")
            list(APPEND HYGIENE_VIOLATIONS "${TRACKED_FILE} (data/local/ is reserved for ignored local datasets)")
        endif()
    endforeach()

    if(HYGIENE_VIOLATIONS)
        list(JOIN HYGIENE_VIOLATIONS "\n  - " HYGIENE_VIOLATIONS_TEXT)
        message(FATAL_ERROR
            "Repository hygiene guard failed.\n"
            "  - ${HYGIENE_VIOLATIONS_TEXT}\n"
            "Move generated/local artifacts to ignored locations before committing.")
    endif()
endfunction()

function(rhbm_guard_check_fixture_tracking)
    rhbm_guard_find_git(GIT_EXECUTABLE)

    rhbm_guard_collect_relative_files(TEST_SOURCE_FILES
        BASE_DIR "${PROJECT_SOURCE_DIR}"
        ROOTS "${PROJECT_SOURCE_DIR}/tests"
        GLOBS "*.cpp" "*.hpp" "*.py"
    )

    set(REFERENCED_FIXTURES)
    foreach(TEST_SOURCE_FILE IN LISTS TEST_SOURCE_FILES)
        file(READ "${PROJECT_SOURCE_DIR}/${TEST_SOURCE_FILE}" TEST_SOURCE_CONTENT)
        string(REGEX MATCHALL "TestDataPath\\(\"[^\"]+\"\\)" FIXTURE_CALLS "${TEST_SOURCE_CONTENT}")
        foreach(FIXTURE_CALL IN LISTS FIXTURE_CALLS)
            string(REGEX REPLACE "^TestDataPath\\(\"([^\"]+)\"\\)$" "\\1" FIXTURE_REL_PATH "${FIXTURE_CALL}")
            list(APPEND REFERENCED_FIXTURES "${FIXTURE_REL_PATH}")
        endforeach()
    endforeach()

    list(REMOVE_DUPLICATES REFERENCED_FIXTURES)

    set(FIXTURE_VIOLATIONS)
    foreach(FIXTURE_REL_PATH IN LISTS REFERENCED_FIXTURES)
        set(FIXTURE_ABS_PATH "${PROJECT_SOURCE_DIR}/tests/data/${FIXTURE_REL_PATH}")
        file(RELATIVE_PATH FIXTURE_REPO_PATH "${PROJECT_SOURCE_DIR}" "${FIXTURE_ABS_PATH}")

        if(NOT EXISTS "${FIXTURE_ABS_PATH}")
            list(APPEND FIXTURE_VIOLATIONS "${FIXTURE_REPO_PATH} (referenced but missing)")
            continue()
        endif()

        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" ls-files --error-unmatch "${FIXTURE_REPO_PATH}"
            RESULT_VARIABLE FIXTURE_TRACK_RESULT
            OUTPUT_QUIET
            ERROR_QUIET
        )

        if(NOT FIXTURE_TRACK_RESULT EQUAL 0)
            list(APPEND FIXTURE_VIOLATIONS "${FIXTURE_REPO_PATH} (referenced but not tracked)")
        endif()
    endforeach()

    if(FIXTURE_VIOLATIONS)
        list(JOIN FIXTURE_VIOLATIONS "\n  - " FIXTURE_VIOLATIONS_TEXT)
        message(FATAL_ERROR
            "Fixture tracking guard failed.\n"
            "  - ${FIXTURE_VIOLATIONS_TEXT}\n"
            "All files referenced by TestDataPath(...) must exist under tests/data/ and be tracked by git.")
    endif()
endfunction()

function(rhbm_guard_check_structure)
    rhbm_guard_find_git(GIT_EXECUTABLE)

    file(GLOB CORE_ROOT_HPP RELATIVE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/core/*.hpp")
    if(CORE_ROOT_HPP)
        list(JOIN CORE_ROOT_HPP "\n  - " CORE_ROOT_HPP_TEXT)
        message(FATAL_ERROR
            "Structure guard failed.\n"
            "  - src/core root must not contain forwarding headers.\n"
            "  - ${CORE_ROOT_HPP_TEXT}\n"
            "Move these headers into src/core/internal/, src/core/workflow/, or remove them.")
    endif()

    file(GLOB PUBLIC_CORE_FLAT_HPP
        RELATIVE "${PROJECT_SOURCE_DIR}"
        "${PROJECT_SOURCE_DIR}/include/rhbm_gem/core/*.hpp")
    if(PUBLIC_CORE_FLAT_HPP)
        list(JOIN PUBLIC_CORE_FLAT_HPP "\n  - " PUBLIC_CORE_FLAT_HPP_TEXT)
        message(FATAL_ERROR
            "Structure guard failed.\n"
            "  - include/rhbm_gem/core must use subfolders (command/ or painter/).\n"
            "  - ${PUBLIC_CORE_FLAT_HPP_TEXT}")
    endif()

    file(GLOB PUBLIC_DATA_FLAT_HPP
        RELATIVE "${PROJECT_SOURCE_DIR}"
        "${PROJECT_SOURCE_DIR}/include/rhbm_gem/data/*.hpp")
    if(PUBLIC_DATA_FLAT_HPP)
        list(JOIN PUBLIC_DATA_FLAT_HPP "\n  - " PUBLIC_DATA_FLAT_HPP_TEXT)
        message(FATAL_ERROR
            "Structure guard failed.\n"
            "  - include/rhbm_gem/data must use subfolders (object/, io/, dispatch/).\n"
            "  - ${PUBLIC_DATA_FLAT_HPP_TEXT}")
    endif()

    file(GLOB PUBLIC_UTILS_FLAT_HPP
        RELATIVE "${PROJECT_SOURCE_DIR}"
        "${PROJECT_SOURCE_DIR}/include/rhbm_gem/utils/*.hpp")
    if(PUBLIC_UTILS_FLAT_HPP)
        list(JOIN PUBLIC_UTILS_FLAT_HPP "\n  - " PUBLIC_UTILS_FLAT_HPP_TEXT)
        message(FATAL_ERROR
            "Structure guard failed.\n"
            "  - include/rhbm_gem/utils must use subfolders (domain/, math/, hrl/).\n"
            "  - ${PUBLIC_UTILS_FLAT_HPP_TEXT}")
    endif()

    file(GLOB CORE_ROOT_COMMAND_CPP RELATIVE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/core/*Command.cpp")
    if(CORE_ROOT_COMMAND_CPP)
        list(JOIN CORE_ROOT_COMMAND_CPP "\n  - " CORE_ROOT_COMMAND_CPP_TEXT)
        message(FATAL_ERROR
            "Structure guard failed.\n"
            "  - src/core root must not contain *Command.cpp files.\n"
            "  - ${CORE_ROOT_COMMAND_CPP_TEXT}\n"
            "Move these sources under src/core/command/.")
    endif()

    rhbm_guard_get_tracked_files("${GIT_EXECUTABLE}" TRACKED_FILE_LIST)

    file(GLOB_RECURSE PUBLIC_HEADERS RELATIVE "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/rhbm_gem/*.hpp")
    set(PUBLIC_INCLUDE_VIOLATIONS)

    foreach(TRACKED_FILE IN LISTS TRACKED_FILE_LIST)
        if(NOT TRACKED_FILE MATCHES "^(src|tests|bindings|examples)/.*\\.(hpp|cpp)$")
            continue()
        endif()

        set(FILE_PATH "${PROJECT_SOURCE_DIR}/${TRACKED_FILE}")
        if(NOT EXISTS "${FILE_PATH}")
            continue()
        endif()
        file(READ "${FILE_PATH}" FILE_CONTENT)

        foreach(PUBLIC_HEADER IN LISTS PUBLIC_HEADERS)
            get_filename_component(PUBLIC_BASENAME "${PUBLIC_HEADER}" NAME)
            string(REGEX REPLACE "([][+.*^$(){}|\\\\])" "\\\\\\1" ESCAPED_PUBLIC_BASENAME "${PUBLIC_BASENAME}")
            string(REGEX MATCH "#include[ \t]*[<\"]${ESCAPED_PUBLIC_BASENAME}[>\"]" MATCHED_INCLUDE "${FILE_CONTENT}")
            if(MATCHED_INCLUDE)
                list(APPEND PUBLIC_INCLUDE_VIOLATIONS
                    "${TRACKED_FILE} (must include <rhbm_gem/.../${PUBLIC_BASENAME}> instead of bare ${PUBLIC_BASENAME})")
                break()
            endif()
        endforeach()
    endforeach()

    if(PUBLIC_INCLUDE_VIOLATIONS)
        list(JOIN PUBLIC_INCLUDE_VIOLATIONS "\n  - " PUBLIC_INCLUDE_VIOLATIONS_TEXT)
        message(FATAL_ERROR
            "Structure guard failed.\n"
            "  - Found non-namespaced includes for public headers.\n"
            "  - ${PUBLIC_INCLUDE_VIOLATIONS_TEXT}\n"
            "Use <rhbm_gem/...> include paths for public headers.")
    endif()
endfunction()

function(rhbm_guard_check_absolute_paths)
    rhbm_guard_find_git(GIT_EXECUTABLE)
    rhbm_guard_get_tracked_files("${GIT_EXECUTABLE}" TRACKED_FILE_LIST)

    set(ABSOLUTE_PATH_VIOLATIONS)
    foreach(TRACKED_FILE IN LISTS TRACKED_FILE_LIST)
        if(NOT TRACKED_FILE MATCHES "^(src|include|bindings|examples|scripts|cmake)/.*\\.(hpp|cpp|h|py|cmake|sh)$")
            continue()
        endif()
        if(TRACKED_FILE STREQUAL "cmake/RepoGuards.cmake")
            continue()
        endif()

        set(FILE_PATH "${PROJECT_SOURCE_DIR}/${TRACKED_FILE}")
        if(NOT EXISTS "${FILE_PATH}")
            continue()
        endif()

        file(READ "${FILE_PATH}" FILE_CONTENT)

        string(REGEX MATCH "/Users/" HAS_USERS_PATH "${FILE_CONTENT}")
        if(HAS_USERS_PATH)
            list(APPEND ABSOLUTE_PATH_VIOLATIONS
                "${TRACKED_FILE} (contains hardcoded /Users/ absolute path)")
            continue()
        endif()

        string(REGEX MATCH "/home/" HAS_HOME_PATH "${FILE_CONTENT}")
        if(HAS_HOME_PATH)
            list(APPEND ABSOLUTE_PATH_VIOLATIONS
                "${TRACKED_FILE} (contains hardcoded /home/ absolute path)")
            continue()
        endif()

        string(REGEX MATCH "(^|[\"' \t(])[A-Za-z]:[/\\\\]" HAS_WINDOWS_DRIVE_PATH "${FILE_CONTENT}")
        if(HAS_WINDOWS_DRIVE_PATH)
            list(APPEND ABSOLUTE_PATH_VIOLATIONS
                "${TRACKED_FILE} (contains hardcoded Windows drive absolute path)")
        endif()
    endforeach()

    if(ABSOLUTE_PATH_VIOLATIONS)
        list(JOIN ABSOLUTE_PATH_VIOLATIONS "\n  - " ABSOLUTE_PATH_VIOLATIONS_TEXT)
        message(FATAL_ERROR
            "Absolute path guard failed.\n"
            "  - ${ABSOLUTE_PATH_VIOLATIONS_TEXT}\n"
            "Use command options, environment variables, or project-relative paths instead.")
    endif()
endfunction()

function(rhbm_guard_run_install_smoke)
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
endfunction()

rhbm_guard_require_project_source_dir()

if(NOT DEFINED RHBM_GEM_GUARD_MODE)
    set(RHBM_GEM_GUARD_MODE "repo")
endif()

if(RHBM_GEM_GUARD_MODE STREQUAL "repo")
    message(STATUS "Running guard: logging_style")
    rhbm_guard_check_logging_style()

    message(STATUS "Running guard: parameterless_void")
    rhbm_guard_check_parameterless_void()

    message(STATUS "Running guard: repository_hygiene")
    rhbm_guard_check_hygiene()

    message(STATUS "Running guard: fixture_tracking")
    rhbm_guard_check_fixture_tracking()

    message(STATUS "Running guard: structure")
    rhbm_guard_check_structure()

    message(STATUS "Running guard: absolute_paths")
    rhbm_guard_check_absolute_paths()
elseif(RHBM_GEM_GUARD_MODE STREQUAL "install_smoke")
    message(STATUS "Running guard: install_smoke")
    rhbm_guard_run_install_smoke()
else()
    message(FATAL_ERROR
        "Unsupported RHBM_GEM_GUARD_MODE='${RHBM_GEM_GUARD_MODE}'. "
        "Expected 'repo' or 'install_smoke'.")
endif()
