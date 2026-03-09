if(NOT DEFINED PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "PROJECT_SOURCE_DIR is required")
endif()

find_program(GIT_EXECUTABLE NAMES git)
if(NOT GIT_EXECUTABLE)
    message(FATAL_ERROR "git executable was not found")
endif()

execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" ls-files
    OUTPUT_VARIABLE TRACKED_FILES_RAW
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE GIT_LS_RESULT
)
if(NOT GIT_LS_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to query tracked files via git ls-files")
endif()

set(TRACKED_FILE_LIST)
if(NOT TRACKED_FILES_RAW STREQUAL "")
    string(REPLACE "\n" ";" TRACKED_FILE_LIST "${TRACKED_FILES_RAW}")
endif()

set(ABSOLUTE_PATH_VIOLATIONS)
foreach(TRACKED_FILE IN LISTS TRACKED_FILE_LIST)
    if(NOT TRACKED_FILE MATCHES "^(src|include|bindings|examples|scripts|cmake)/.*\\.(hpp|cpp|h|py|cmake|sh)$")
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
