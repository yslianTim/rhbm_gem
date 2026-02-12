cmake_minimum_required(VERSION 3.15)

if(NOT DEFINED COVERAGE_ACTION)
    set(COVERAGE_ACTION "report")
endif()

if(NOT DEFINED COVERAGE_BUILD_DIR)
    message(FATAL_ERROR "COVERAGE_BUILD_DIR is required")
endif()

if(COVERAGE_ACTION STREQUAL "clean")
    file(GLOB_RECURSE GCDA_FILES "${COVERAGE_BUILD_DIR}/*.gcda")
    list(LENGTH GCDA_FILES GCDA_COUNT)
    if(GCDA_COUNT GREATER 0)
        file(REMOVE ${GCDA_FILES})
    endif()
    message(STATUS "Coverage cleanup complete: removed ${GCDA_COUNT} .gcda files")
    return()
endif()

if(NOT COVERAGE_ACTION STREQUAL "report")
    message(FATAL_ERROR "Unsupported COVERAGE_ACTION: ${COVERAGE_ACTION}")
endif()

if(NOT DEFINED COVERAGE_PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "COVERAGE_PROJECT_SOURCE_DIR is required for report action")
endif()

if(NOT DEFINED COVERAGE_OUTPUT_DIR)
    set(COVERAGE_OUTPUT_DIR "${COVERAGE_BUILD_DIR}/coverage")
endif()

if(NOT DEFINED COVERAGE_INCLUDE_TESTS)
    set(COVERAGE_INCLUDE_TESTS OFF)
endif()

if(NOT DEFINED GCOV_EXECUTABLE OR GCOV_EXECUTABLE STREQUAL "")
    find_program(GCOV_EXECUTABLE NAMES gcov)
endif()

if(NOT GCOV_EXECUTABLE)
    message(FATAL_ERROR "gcov executable was not found")
endif()

function(percent_to_basis_points INPUT_PERCENT OUTPUT_VAR)
    set(PERCENT_BASIS_POINTS 0)
    if(INPUT_PERCENT MATCHES "^([0-9]+)\\.([0-9]+)$")
        set(PERCENT_INT "${CMAKE_MATCH_1}")
        set(PERCENT_FRAC "${CMAKE_MATCH_2}")
        string(CONCAT PERCENT_FRAC_PADDED "${PERCENT_FRAC}" "00")
        string(SUBSTRING "${PERCENT_FRAC_PADDED}" 0 2 PERCENT_FRAC_2)
        math(EXPR PERCENT_BASIS_POINTS "${PERCENT_INT} * 100 + ${PERCENT_FRAC_2}")
    elseif(INPUT_PERCENT MATCHES "^[0-9]+$")
        math(EXPR PERCENT_BASIS_POINTS "${INPUT_PERCENT} * 100")
    endif()
    set(${OUTPUT_VAR} "${PERCENT_BASIS_POINTS}" PARENT_SCOPE)
endfunction()

function(format_basis_points INPUT_BASIS_POINTS OUTPUT_VAR)
    math(EXPR PERCENT_INT "${INPUT_BASIS_POINTS} / 100")
    math(EXPR PERCENT_FRAC "${INPUT_BASIS_POINTS} % 100")
    if(PERCENT_FRAC LESS 10)
        set(PERCENT_FRAC_STR "0${PERCENT_FRAC}")
    else()
        set(PERCENT_FRAC_STR "${PERCENT_FRAC}")
    endif()
    set(${OUTPUT_VAR} "${PERCENT_INT}.${PERCENT_FRAC_STR}" PARENT_SCOPE)
endfunction()

function(should_include_file INPUT_FILE INCLUDE_TESTS OUTPUT_VAR)
    set(INCLUDE_FILE FALSE)
    if(INPUT_FILE MATCHES "^third_party/")
        set(INCLUDE_FILE FALSE)
    elseif(INPUT_FILE MATCHES "^tests/")
        if(INCLUDE_TESTS)
            set(INCLUDE_FILE TRUE)
        endif()
    elseif(INPUT_FILE MATCHES "^src/")
        set(INCLUDE_FILE TRUE)
    endif()
    set(${OUTPUT_VAR} "${INCLUDE_FILE}" PARENT_SCOPE)
endfunction()

file(MAKE_DIRECTORY "${COVERAGE_OUTPUT_DIR}")
set(COVERAGE_DETAIL_FILE "${COVERAGE_OUTPUT_DIR}/coverage_detail.txt")
set(COVERAGE_SUMMARY_FILE "${COVERAGE_OUTPUT_DIR}/coverage_summary.txt")

file(GLOB_RECURSE GCNO_FILES "${COVERAGE_BUILD_DIR}/*.gcno")
list(SORT GCNO_FILES)

set(DETAIL_CONTENT "# file|line_coverage|lines\n")
set(SEEN_FILES)
set(SCANNED_GCNO_COUNT 0)
set(CONSIDERED_FILE_COUNT 0)
set(TOTAL_WEIGHTED_BASIS_X_LINES 0)
set(TOTAL_LINES 0)

foreach(GCNO_FILE IN LISTS GCNO_FILES)
    math(EXPR SCANNED_GCNO_COUNT "${SCANNED_GCNO_COUNT} + 1")
    execute_process(
        COMMAND "${GCOV_EXECUTABLE}" -n -r -s "${COVERAGE_PROJECT_SOURCE_DIR}" "${GCNO_FILE}"
        WORKING_DIRECTORY "${COVERAGE_OUTPUT_DIR}"
        RESULT_VARIABLE GCOV_RESULT
        OUTPUT_VARIABLE GCOV_OUTPUT
        ERROR_VARIABLE GCOV_ERROR
    )

    if(NOT GCOV_RESULT EQUAL 0)
        continue()
    endif()

    string(REPLACE "\r\n" "\n" GCOV_OUTPUT "${GCOV_OUTPUT}")
    string(REPLACE "\r" "\n" GCOV_OUTPUT "${GCOV_OUTPUT}")
    string(REPLACE "\n" ";" GCOV_LINES "${GCOV_OUTPUT}")

    set(CURRENT_FILE "")
    foreach(GCOV_LINE IN LISTS GCOV_LINES)
        if(GCOV_LINE MATCHES "^File '([^']+)'$")
            set(CURRENT_FILE "${CMAKE_MATCH_1}")
            if(IS_ABSOLUTE "${CURRENT_FILE}")
                file(RELATIVE_PATH CURRENT_FILE "${COVERAGE_PROJECT_SOURCE_DIR}" "${CURRENT_FILE}")
            endif()
        elseif(GCOV_LINE MATCHES "^Lines executed:([0-9]+\\.?[0-9]*)% of ([0-9]+)$")
            if(CURRENT_FILE STREQUAL "")
                continue()
            endif()

            should_include_file("${CURRENT_FILE}" "${COVERAGE_INCLUDE_TESTS}" INCLUDE_FILE)
            if(NOT INCLUDE_FILE)
                continue()
            endif()

            list(FIND SEEN_FILES "${CURRENT_FILE}" SEEN_INDEX)
            if(NOT SEEN_INDEX EQUAL -1)
                continue()
            endif()

            set(FILE_LINE_PERCENT "${CMAKE_MATCH_1}")
            set(FILE_LINE_COUNT "${CMAKE_MATCH_2}")

            percent_to_basis_points("${FILE_LINE_PERCENT}" FILE_LINE_PERCENT_BASIS)
            math(EXPR FILE_WEIGHTED_BASIS_X_LINES
                "${FILE_LINE_PERCENT_BASIS} * ${FILE_LINE_COUNT}")

            math(EXPR TOTAL_WEIGHTED_BASIS_X_LINES
                "${TOTAL_WEIGHTED_BASIS_X_LINES} + ${FILE_WEIGHTED_BASIS_X_LINES}")
            math(EXPR TOTAL_LINES "${TOTAL_LINES} + ${FILE_LINE_COUNT}")
            math(EXPR CONSIDERED_FILE_COUNT "${CONSIDERED_FILE_COUNT} + 1")
            list(APPEND SEEN_FILES "${CURRENT_FILE}")

            string(APPEND DETAIL_CONTENT
                "${CURRENT_FILE}|${FILE_LINE_PERCENT}%|${FILE_LINE_COUNT}\n")
        endif()
    endforeach()
endforeach()

if(TOTAL_LINES GREATER 0)
    math(EXPR WEIGHTED_COVERAGE_BASIS_POINTS
        "${TOTAL_WEIGHTED_BASIS_X_LINES} / ${TOTAL_LINES}")
else()
    set(WEIGHTED_COVERAGE_BASIS_POINTS 0)
endif()

format_basis_points("${WEIGHTED_COVERAGE_BASIS_POINTS}" WEIGHTED_COVERAGE_PERCENT)

if(COVERAGE_INCLUDE_TESTS)
    set(INCLUDE_TESTS_LABEL "ON")
else()
    set(INCLUDE_TESTS_LABEL "OFF")
endif()

string(CONCAT SUMMARY_CONTENT
    "Coverage Summary (gcov)\n"
    "Scanned .gcno files: ${SCANNED_GCNO_COUNT}\n"
    "Included files: ${CONSIDERED_FILE_COUNT}\n"
    "Total included lines: ${TOTAL_LINES}\n"
    "Weighted line coverage: ${WEIGHTED_COVERAGE_PERCENT}%\n"
    "Include tests: ${INCLUDE_TESTS_LABEL}\n"
    "Detail report: ${COVERAGE_DETAIL_FILE}\n"
)

file(WRITE "${COVERAGE_DETAIL_FILE}" "${DETAIL_CONTENT}")
file(WRITE "${COVERAGE_SUMMARY_FILE}" "${SUMMARY_CONTENT}")

message(STATUS "Coverage detail report generated: ${COVERAGE_DETAIL_FILE}")
message(STATUS "Coverage summary report generated: ${COVERAGE_SUMMARY_FILE}")
message(STATUS "Weighted line coverage: ${WEIGHTED_COVERAGE_PERCENT}%")
