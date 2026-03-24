cmake_minimum_required(VERSION 3.18)
project(RHBM_GEM
    VERSION 1.0
    DESCRIPTION RHBM-GEM
    HOMEPAGE_URL https://github.com/yslianTim/rhbm_gem
    LANGUAGES C CXX
)

set(PROJECT_DISPLAY_NAME "RHBM-GEM")
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Keep installed binaries relocatable with project-local shared libraries.
if(APPLE)
    set(RHBM_GEM_EXECUTABLE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
    set(RHBM_GEM_LIBRARY_INSTALL_RPATH "@loader_path")
elseif(UNIX)
    set(RHBM_GEM_EXECUTABLE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    set(RHBM_GEM_LIBRARY_INSTALL_RPATH "$ORIGIN")
endif()

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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CTest)
include(RHBMGemOptions)
include(RHBMGemDependencies)

# Keep all runtime targets (CLI, GUI, tests) in a single predictable folder.
set(RHBM_GEM_RUNTIME_OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${RHBM_GEM_RUNTIME_OUTPUT_DIR}")
foreach(config IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER "${config}" config_upper)
    set("CMAKE_RUNTIME_OUTPUT_DIRECTORY_${config_upper}" "${RHBM_GEM_RUNTIME_OUTPUT_DIR}")
endforeach()

add_subdirectory(src)

include(RHBMGemInstall)
include(RHBMGemDevTools)
