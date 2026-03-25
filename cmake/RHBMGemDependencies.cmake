set(RHBM_GEM_EIGEN3_URL "https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz")
set(RHBM_GEM_EIGEN3_URL_HASH "SHA256=e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec")

set(RHBM_GEM_CLI11_URL "https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.5.0.tar.gz")
set(RHBM_GEM_CLI11_URL_HASH "SHA256=17e02b4cddc2fa348e5dbdbb582c59a3486fa2b2433e70a0c3bacb871334fd55")

set(RHBM_GEM_GTEST_URL "https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz")
set(RHBM_GEM_GTEST_URL_HASH "SHA256=65fab701d9829d38cb77c14acdc431d2108bfdbf8979e40eb8ae567edf10b27c")

set(RHBM_GEM_PYBIND11_URL "https://github.com/pybind/pybind11/archive/refs/tags/v3.0.2.tar.gz")
set(RHBM_GEM_PYBIND11_URL_HASH "SHA256=2f20a0af0b921815e0e169ea7fec63909869323581b89d7de1553468553f6a2d")

set(RHBM_GEM_SQLITE3_URL "https://www.sqlite.org/2025/sqlite-amalgamation-3490100.zip")
set(RHBM_GEM_SQLITE3_URL_HASH "SHA256=6cebd1d8403fc58c30e93939b246f3e6e58d0765a5cd50546f16c00fd805d2c3")

set(RHBM_GEM_BOOST_URL "https://archives.boost.io/release/1.90.0/source/boost_1_90_0.tar.bz2")
set(RHBM_GEM_BOOST_URL_HASH "SHA256=49551aff3b22cbc5c5a9ed3dbc92f0e23ea50a0f7325b0d198b705e8ee3fc305")
set(RHBM_GEM_BOOST_FALLBACK_VERSION "1.90.0")

if(RHBM_GEM_DEP_PROVIDER STREQUAL "FETCH")
    include(FetchContent)
endif()

function(rhbm_gem_populate_content dep_name dep_url dep_hash out_source_dir)
    FetchContent_Declare(${dep_name}
        URL "${dep_url}"
        URL_HASH "${dep_hash}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(${dep_name})
    FetchContent_GetProperties(${dep_name})
    set(${out_source_dir} "${${dep_name}_SOURCE_DIR}" PARENT_SCOPE)
endfunction()

function(rhbm_gem_prepare_openmp_for_appleclang)
    if(NOT APPLE OR NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        return()
    endif()

    if(DEFINED OpenMP_ROOT)
        if(EXISTS "${OpenMP_ROOT}/include/omp.h")
            return()
        endif()
        return()
    endif()

    set(_rhbm_openmp_root "")
    set(_rhbm_openmp_candidates)
    if(DEFINED ENV{HOMEBREW_PREFIX})
        list(APPEND _rhbm_openmp_candidates "$ENV{HOMEBREW_PREFIX}/opt/libomp")
    endif()
    list(APPEND _rhbm_openmp_candidates
        "/opt/homebrew/opt/libomp"
        "/usr/local/opt/libomp"
    )

    foreach(_candidate_root IN LISTS _rhbm_openmp_candidates)
        if(EXISTS "${_candidate_root}/include/omp.h"
           AND (EXISTS "${_candidate_root}/lib/libomp.dylib"
                OR EXISTS "${_candidate_root}/lib/libomp.a"))
            set(_rhbm_openmp_root "${_candidate_root}")
            break()
        endif()
    endforeach()

    if(_rhbm_openmp_root)
        message(STATUS "AppleClang detected, probing OpenMP in: ${_rhbm_openmp_root}")
        set(OpenMP_ROOT "${_rhbm_openmp_root}" PARENT_SCOPE)
    endif()
endfunction()

function(rhbm_gem_link_boost_dependency target_name)
    if(RHBM_GEM_DEP_PROVIDER STREQUAL "SYSTEM")
        message(STATUS "Using system Boost package")
        if(TARGET Boost::headers)
            target_link_libraries(${target_name} INTERFACE Boost::headers)
        elseif(TARGET Boost::boost)
            target_link_libraries(${target_name} INTERFACE Boost::boost)
        elseif(Boost_INCLUDE_DIRS)
            target_include_directories(${target_name} SYSTEM INTERFACE ${Boost_INCLUDE_DIRS})
        else()
            message(FATAL_ERROR
                "System Boost package was configured but no Boost include target/directory was exported.")
        endif()
    else()
        message(STATUS "Using fetched Boost fallback (v${RHBM_GEM_BOOST_FALLBACK_VERSION})")
        target_include_directories(${target_name} SYSTEM INTERFACE
            "$<BUILD_INTERFACE:${RHBM_GEM_BOOST_INCLUDE_DIR}>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )
    endif()
endfunction()

add_library(rhbm_gem_dependencies INTERFACE)
set(RHBM_GEM_WITH_OPENMP FALSE)
set(RHBM_GEM_WITH_ROOT FALSE)
set(RHBM_GEM_OPENMP_ROOT "")

if(RHBM_GEM_DEP_PROVIDER STREQUAL "SYSTEM")
    message(STATUS "Dependency provider: SYSTEM")
    find_package(Eigen3 REQUIRED)
    find_package(CLI11 REQUIRED)
    find_package(SQLite3 REQUIRED)
    if(POLICY CMP0167)
        cmake_policy(PUSH)
        cmake_policy(SET CMP0167 OLD)
        find_package(Boost REQUIRED)
        cmake_policy(POP)
    else()
        find_package(Boost REQUIRED)
    endif()

    target_link_libraries(rhbm_gem_dependencies INTERFACE
        Eigen3::Eigen
        CLI11::CLI11
        SQLite::SQLite3
    )
else()
    message(STATUS "Dependency provider: FETCH")

    rhbm_gem_populate_content(
        rhbm_gem_eigen3
        "${RHBM_GEM_EIGEN3_URL}"
        "${RHBM_GEM_EIGEN3_URL_HASH}"
        RHBM_GEM_EIGEN3_SOURCE_DIR
    )

    rhbm_gem_populate_content(
        rhbm_gem_cli11
        "${RHBM_GEM_CLI11_URL}"
        "${RHBM_GEM_CLI11_URL_HASH}"
        RHBM_GEM_CLI11_SOURCE_DIR
    )

    rhbm_gem_populate_content(
        rhbm_gem_sqlite3
        "${RHBM_GEM_SQLITE3_URL}"
        "${RHBM_GEM_SQLITE3_URL_HASH}"
        RHBM_GEM_SQLITE3_SOURCE_DIR
    )
    if(NOT EXISTS "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.c"
        OR NOT EXISTS "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.h")
        message(FATAL_ERROR
            "Fetched SQLite3 source at '${RHBM_GEM_SQLITE3_SOURCE_DIR}' "
            "does not contain sqlite3.c/sqlite3.h.")
    endif()

    rhbm_gem_populate_content(
        rhbm_gem_boost
        "${RHBM_GEM_BOOST_URL}"
        "${RHBM_GEM_BOOST_URL_HASH}"
        RHBM_GEM_BOOST_SOURCE_DIR
    )
    if(NOT EXISTS "${RHBM_GEM_BOOST_SOURCE_DIR}/boost")
        message(FATAL_ERROR
            "Fetched Boost source at '${RHBM_GEM_BOOST_SOURCE_DIR}' "
            "does not contain the expected boost/ headers.")
    endif()

    set(RHBM_GEM_EIGEN_INCLUDE_DIR "${RHBM_GEM_EIGEN3_SOURCE_DIR}")
    set(RHBM_GEM_CLI11_INCLUDE_DIR "${RHBM_GEM_CLI11_SOURCE_DIR}/include")
    set(RHBM_GEM_BOOST_INCLUDE_DIR "${RHBM_GEM_BOOST_SOURCE_DIR}")

    target_include_directories(rhbm_gem_dependencies SYSTEM INTERFACE
        "$<BUILD_INTERFACE:${RHBM_GEM_EIGEN_INCLUDE_DIR}>"
        "$<BUILD_INTERFACE:${RHBM_GEM_CLI11_INCLUDE_DIR}>"
        "$<BUILD_INTERFACE:${RHBM_GEM_SQLITE3_SOURCE_DIR}>"
        "$<BUILD_INTERFACE:${RHBM_GEM_BOOST_INCLUDE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
endif()

rhbm_gem_link_boost_dependency(rhbm_gem_dependencies)

set(_rhbm_gem_min_openmp_version "4.5")
if(RHBM_GEM_OPENMP_MODE STREQUAL "OFF")
    message(STATUS "OpenMP features disabled by RHBM_GEM_OPENMP_MODE=OFF")
elseif(RHBM_GEM_OPENMP_MODE STREQUAL "AUTO")
    rhbm_gem_prepare_openmp_for_appleclang()
    find_package(OpenMP QUIET COMPONENTS CXX)
elseif(RHBM_GEM_OPENMP_MODE STREQUAL "ON")
    rhbm_gem_prepare_openmp_for_appleclang()
    find_package(OpenMP REQUIRED COMPONENTS CXX)
else()
    message(FATAL_ERROR "Unsupported RHBM_GEM_OPENMP_MODE='${RHBM_GEM_OPENMP_MODE}'")
endif()

if(OpenMP_CXX_FOUND)
    if(DEFINED OpenMP_CXX_VERSION AND OpenMP_CXX_VERSION VERSION_LESS _rhbm_gem_min_openmp_version)
        if(RHBM_GEM_OPENMP_MODE STREQUAL "ON")
            message(FATAL_ERROR
                "OpenMP version ${OpenMP_CXX_VERSION} is below required ${_rhbm_gem_min_openmp_version} "
                "while RHBM_GEM_OPENMP_MODE=ON.")
        endif()
        message(WARNING
            "Found OpenMP version ${OpenMP_CXX_VERSION} but need at least ${_rhbm_gem_min_openmp_version}. "
            "Disabling OpenMP features.")
    else()
        message(STATUS "OpenMP found (version ${OpenMP_CXX_VERSION}). Enabling OpenMP features.")
        target_compile_definitions(rhbm_gem_dependencies INTERFACE USE_OPENMP)
        target_link_libraries(rhbm_gem_dependencies INTERFACE OpenMP::OpenMP_CXX)
        set(RHBM_GEM_WITH_OPENMP TRUE)
        if(DEFINED OpenMP_ROOT)
            set(RHBM_GEM_OPENMP_ROOT "${OpenMP_ROOT}")
        endif()
    endif()
elseif(RHBM_GEM_OPENMP_MODE STREQUAL "AUTO")
    message(STATUS "OpenMP support not found, using serial version")
endif()

if(RHBM_GEM_ROOT_MODE STREQUAL "OFF")
    message(STATUS "ROOT features disabled by RHBM_GEM_ROOT_MODE=OFF")
elseif(RHBM_GEM_ROOT_MODE STREQUAL "AUTO")
    find_package(ROOT 6.28 QUIET COMPONENTS Core Hist Gpad RIO)
elseif(RHBM_GEM_ROOT_MODE STREQUAL "ON")
    find_package(ROOT 6.28 REQUIRED COMPONENTS Core Hist Gpad RIO)
else()
    message(FATAL_ERROR "Unsupported RHBM_GEM_ROOT_MODE='${RHBM_GEM_ROOT_MODE}'")
endif()

if(ROOT_FOUND)
    message(STATUS "ROOT library found, enabling ROOT features")
    target_compile_definitions(rhbm_gem_dependencies INTERFACE HAVE_ROOT)
    target_link_libraries(rhbm_gem_dependencies INTERFACE ROOT::Core ROOT::Hist ROOT::Gpad ROOT::RIO)
    set(RHBM_GEM_WITH_ROOT TRUE)
elseif(RHBM_GEM_ROOT_MODE STREQUAL "AUTO")
    message(STATUS "ROOT library not found, disabling ROOT features")
endif()

if(RHBM_GEM_BUILD_GUI)
    find_package(Qt6 QUIET COMPONENTS Core Widgets)
    if(Qt6_FOUND)
        message(STATUS "Qt6 found, enabling GUI executable build")
        add_library(rhbm_gem_gui_dependencies INTERFACE)
        target_link_libraries(rhbm_gem_gui_dependencies INTERFACE Qt6::Core Qt6::Widgets)
    else()
        message(STATUS "Qt6 not found, skipping GUI executable build")
    endif()
endif()

if(BUILD_PYTHON_BINDINGS)
    set(PYBIND11_FINDPYTHON ON)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

    if(RHBM_GEM_DEP_PROVIDER STREQUAL "SYSTEM")
        find_package(pybind11 CONFIG REQUIRED)
    else()
        set(PYBIND11_TEST OFF CACHE BOOL "Disable pybind11 tests" FORCE)
        FetchContent_Declare(rhbm_gem_pybind11
            URL "${RHBM_GEM_PYBIND11_URL}"
            URL_HASH "${RHBM_GEM_PYBIND11_URL_HASH}"
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        message(STATUS "Fetching pybind11 (v3.0.2) via FetchContent")
        FetchContent_MakeAvailable(rhbm_gem_pybind11)
    endif()
endif()

if(BUILD_TESTING)
    # Do not export test framework artifacts in normal project installs.
    set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest install targets" FORCE)
    set(INSTALL_GMOCK OFF CACHE BOOL "Disable GoogleMock install targets" FORCE)

    if(RHBM_GEM_DEP_PROVIDER STREQUAL "SYSTEM")
        find_package(GTest REQUIRED)
    else()
        if(MSVC)
            set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        endif()
        FetchContent_Declare(rhbm_gem_googletest
            URL "${RHBM_GEM_GTEST_URL}"
            URL_HASH "${RHBM_GEM_GTEST_URL_HASH}"
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        message(STATUS "Fetching GoogleTest (v1.17.0) via FetchContent")
        FetchContent_MakeAvailable(rhbm_gem_googletest)
    endif()
endif()
