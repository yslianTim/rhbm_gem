cmake_minimum_required(VERSION 3.18)
message(STATUS "----------------------------")
message(STATUS "---   RHBM-GEM Project   ---")
message(STATUS "----------------------------")
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
    set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# Centralize all primary output artifacts under the top-level build directory.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
if(CMAKE_CONFIGURATION_TYPES)
    foreach(_config IN LISTS CMAKE_CONFIGURATION_TYPES)
        string(TOUPPER "${_config}" _config_upper)
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${_config_upper} "${PROJECT_BINARY_DIR}")
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${_config_upper} "${PROJECT_BINARY_DIR}")
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_config_upper} "${PROJECT_BINARY_DIR}")
    endforeach()
    unset(_config_upper)
    unset(_config)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(RHBMGemOptions)
include(RHBMGemDependencies)

add_subdirectory(src)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(bindings)
endif()

include(RHBMGemInstall)
include(RHBMGemTesting)
