function(rhbm_gem_set_install_rpath target_name install_rpath)
    if("${install_rpath}" STREQUAL "")
        return()
    endif()

    set_target_properties(${target_name} PROPERTIES
        INSTALL_RPATH "${install_rpath}"
    )
endfunction()

function(rhbm_gem_apply_library_install_rpath target_name)
    if(DEFINED RHBM_GEM_LIBRARY_INSTALL_RPATH)
        rhbm_gem_set_install_rpath(${target_name} "${RHBM_GEM_LIBRARY_INSTALL_RPATH}")
    endif()
endfunction()

function(rhbm_gem_apply_executable_install_rpath target_name)
    if(DEFINED RHBM_GEM_EXECUTABLE_INSTALL_RPATH)
        rhbm_gem_set_install_rpath(${target_name} "${RHBM_GEM_EXECUTABLE_INSTALL_RPATH}")
    endif()
endfunction()
