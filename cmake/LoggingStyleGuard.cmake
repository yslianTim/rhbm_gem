set(_logging_guard_sources
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/include"
)

set(_stream_pattern "\\bstd::cout\\b|\\bstd::cerr\\b|\\bprintf\\(|\\bfprintf\\(|\\bputs\\(|\\bstd::clog\\b")
set(_chatty_debug_pattern "Logger::Log\\(LogLevel::Debug,.*called")
set(_allowed_stream_file "${PROJECT_SOURCE_DIR}/src/utils/Logger.cpp")

foreach(_root IN LISTS _logging_guard_sources)
    file(GLOB_RECURSE _logging_guard_files
        RELATIVE "${PROJECT_SOURCE_DIR}"
        "${_root}/*")

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
endforeach()
