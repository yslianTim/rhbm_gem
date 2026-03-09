function(rhbm_gem_require_bundled_dependency dependency_name check_path)
    set(_dependency_path "${PROJECT_SOURCE_DIR}/${check_path}")
    if(NOT EXISTS "${_dependency_path}")
        message(FATAL_ERROR
            "Bundled ${dependency_name} is required but '${check_path}' is missing.\n"
            "Restore the bundled dependency source under third_party/, or install\n"
            "the system package and configure with -DUSE_SYSTEM_LIBS=ON.")
    endif()
endfunction()
