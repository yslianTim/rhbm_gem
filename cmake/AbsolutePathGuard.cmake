include("${CMAKE_CURRENT_LIST_DIR}/GuardCommon.cmake")
rhbm_guard_require_project_source_dir()
rhbm_guard_find_git(GIT_EXECUTABLE)
rhbm_guard_get_tracked_files("${GIT_EXECUTABLE}" TRACKED_FILE_LIST)

set(ABSOLUTE_PATH_VIOLATIONS)
foreach(TRACKED_FILE IN LISTS TRACKED_FILE_LIST)
    if(NOT TRACKED_FILE MATCHES "^(src|include|bindings|examples|scripts|cmake)/.*\\.(hpp|cpp|h|py|cmake|sh)$")
        continue()
    endif()
    if(TRACKED_FILE STREQUAL "cmake/AbsolutePathGuard.cmake")
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
