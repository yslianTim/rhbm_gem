set(_rhbm_gem_source_dir "${CMAKE_CURRENT_LIST_DIR}/..")
get_filename_component(_rhbm_gem_source_dir "${_rhbm_gem_source_dir}" ABSOLUTE)

set(_rhbm_gem_presets
    debug
    release
)

foreach(_rhbm_gem_preset IN LISTS _rhbm_gem_presets)
    message(STATUS "[${_rhbm_gem_preset}] configuring")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" --preset "${_rhbm_gem_preset}"
        WORKING_DIRECTORY "${_rhbm_gem_source_dir}"
        RESULT_VARIABLE _rhbm_gem_configure_result
    )
    if(NOT _rhbm_gem_configure_result STREQUAL "0")
        message(FATAL_ERROR
            "Configure failed for preset '${_rhbm_gem_preset}' "
            "(result: ${_rhbm_gem_configure_result}).")
    endif()

    message(STATUS "[${_rhbm_gem_preset}] building")
    execute_process(
        COMMAND "${CMAKE_COMMAND}" --build --preset "${_rhbm_gem_preset}"
        WORKING_DIRECTORY "${_rhbm_gem_source_dir}"
        RESULT_VARIABLE _rhbm_gem_build_result
    )
    if(NOT _rhbm_gem_build_result STREQUAL "0")
        message(FATAL_ERROR
            "Build failed for preset '${_rhbm_gem_preset}' "
            "(result: ${_rhbm_gem_build_result}).")
    endif()
endforeach()

message(STATUS "All Debug and Release presets built successfully.")
