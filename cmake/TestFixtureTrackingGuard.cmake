include("${CMAKE_CURRENT_LIST_DIR}/GuardCommon.cmake")
rhbm_guard_require_project_source_dir()
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
