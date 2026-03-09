if(NOT DEFINED PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "PROJECT_SOURCE_DIR is required")
endif()

find_program(GIT_EXECUTABLE NAMES git)
if(NOT GIT_EXECUTABLE)
    message(FATAL_ERROR "git executable was not found")
endif()

file(GLOB_RECURSE TEST_SOURCE_FILES
    "${PROJECT_SOURCE_DIR}/tests/*.cpp"
    "${PROJECT_SOURCE_DIR}/tests/*.hpp"
    "${PROJECT_SOURCE_DIR}/tests/*.py"
)

set(REFERENCED_FIXTURES)
foreach(TEST_SOURCE_FILE IN LISTS TEST_SOURCE_FILES)
    file(READ "${TEST_SOURCE_FILE}" TEST_SOURCE_CONTENT)
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
