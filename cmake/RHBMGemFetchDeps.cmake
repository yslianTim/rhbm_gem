include(FetchContent)

set(RHBM_GEM_EIGEN3_URL "https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz")
set(RHBM_GEM_EIGEN3_URL_HASH "SHA256=e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec")

set(RHBM_GEM_GTEST_URL "https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz")
set(RHBM_GEM_GTEST_URL_HASH "SHA256=65fab701d9829d38cb77c14acdc431d2108bfdbf8979e40eb8ae567edf10b27c")

set(RHBM_GEM_PYBIND11_URL "https://github.com/pybind/pybind11/archive/refs/tags/v3.0.2.tar.gz")
set(RHBM_GEM_PYBIND11_URL_HASH "SHA256=2f20a0af0b921815e0e169ea7fec63909869323581b89d7de1553468553f6a2d")

function(rhbm_gem_fetch_eigen3 out_source_dir)
    FetchContent_Declare(rhbm_gem_eigen3
        URL "${RHBM_GEM_EIGEN3_URL}"
        URL_HASH "${RHBM_GEM_EIGEN3_URL_HASH}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    message(STATUS "Fetching Eigen3 (5.0.1) via FetchContent")
    FetchContent_MakeAvailable(rhbm_gem_eigen3)
    FetchContent_GetProperties(rhbm_gem_eigen3)
    set(${out_source_dir} "${rhbm_gem_eigen3_SOURCE_DIR}" PARENT_SCOPE)
endfunction()

function(rhbm_gem_fetch_googletest)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest install targets" FORCE)
    set(INSTALL_GMOCK OFF CACHE BOOL "Disable GoogleMock install targets" FORCE)

    FetchContent_Declare(rhbm_gem_googletest
        URL "${RHBM_GEM_GTEST_URL}"
        URL_HASH "${RHBM_GEM_GTEST_URL_HASH}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    message(STATUS "Fetching GoogleTest (v1.17.0) via FetchContent")
    FetchContent_MakeAvailable(rhbm_gem_googletest)
endfunction()

function(rhbm_gem_fetch_pybind11)
    set(PYBIND11_TEST OFF CACHE BOOL "Disable pybind11 tests" FORCE)

    FetchContent_Declare(rhbm_gem_pybind11
        URL "${RHBM_GEM_PYBIND11_URL}"
        URL_HASH "${RHBM_GEM_PYBIND11_URL_HASH}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    message(STATUS "Fetching pybind11 (v3.0.2) via FetchContent")
    FetchContent_MakeAvailable(rhbm_gem_pybind11)
endfunction()
