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

    set(RHBM_GEM_EIGEN_INCLUDE_DIR "")
    set(RHBM_GEM_CLI11_INCLUDE_DIR "")
    set(RHBM_GEM_SQLITE3_SOURCE_DIR "")
    set(RHBM_GEM_BOOST_INCLUDE_DIR "")
else()
    message(STATUS "Dependency provider: FETCH")
    include(FetchContent)

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
endif()

if(BUILD_PYTHON_BINDINGS)
    set(PYBIND11_FINDPYTHON ON)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

    if(RHBM_GEM_DEP_PROVIDER STREQUAL "SYSTEM")
        find_package(pybind11 CONFIG REQUIRED)
    else()
        include(FetchContent)
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
        include(FetchContent)
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_Declare(rhbm_gem_googletest
            URL "${RHBM_GEM_GTEST_URL}"
            URL_HASH "${RHBM_GEM_GTEST_URL_HASH}"
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        message(STATUS "Fetching GoogleTest (v1.17.0) via FetchContent")
        FetchContent_MakeAvailable(rhbm_gem_googletest)
    endif()
endif()
