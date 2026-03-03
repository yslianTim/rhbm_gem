set(_parameterless_void_guard_sources
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/include"
    "${PROJECT_SOURCE_DIR}/tests"
    "${PROJECT_SOURCE_DIR}/bindings"
)

set(_parameterless_void_signature_pattern
    "\\([ \t]*void[ \t]*\\)[ \t]*(const|override|final|noexcept|requires|[:;{=]|$)")

foreach(_root IN LISTS _parameterless_void_guard_sources)
    if(NOT EXISTS "${_root}")
        continue()
    endif()

    file(GLOB_RECURSE _parameterless_void_guard_files
        RELATIVE "${PROJECT_SOURCE_DIR}"
        "${_root}/*")

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
endforeach()
