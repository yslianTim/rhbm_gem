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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CTest)
include(RHBMGemOptions)
include(RHBMGemTargetHelpers)
include(RHBMGemDependencies)

add_subdirectory(src)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(bindings)
endif()

include(RHBMGemInstall)
include(RHBMGemTesting)
include(RHBMGemQuality)
