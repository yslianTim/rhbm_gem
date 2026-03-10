if(NOT DEFINED PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "PROJECT_SOURCE_DIR is required")
endif()

find_program(GIT_EXECUTABLE NAMES git)
if(NOT GIT_EXECUTABLE)
    message(FATAL_ERROR "git executable was not found")
endif()

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
