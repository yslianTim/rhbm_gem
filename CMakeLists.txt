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

# Avoid writing entries into the user package registry from bundled deps.
set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON)

# Keep installed binaries relocatable with project-local shared libraries.
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(RHBMGemOptions)
include(RHBMGemDependencies)

add_subdirectory(src)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(bindings)
endif()

include(RHBMGemInstall)
include(RHBMGemTesting)
