function(rhbm_gem_require_bundled_dependency dependency_name check_path)
    set(_dependency_path "${PROJECT_SOURCE_DIR}/${check_path}")
    if(NOT EXISTS "${_dependency_path}")
        message(FATAL_ERROR
            "Bundled ${dependency_name} is required but '${check_path}' is missing.\n"
            "This project expects third-party dependencies via git submodules.\n"
            "Run: git submodule update --init --recursive")
    endif()
endfunction()

