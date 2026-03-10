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
