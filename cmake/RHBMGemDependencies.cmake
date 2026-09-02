set(RHBM_GEM_EIGEN3_MIN_VERSION "5.0.0")
set(RHBM_GEM_EIGEN3_MAX_EXCLUSIVE_VERSION "6.0.0")
set(RHBM_GEM_EIGEN3_VERSION_RANGE
    "${RHBM_GEM_EIGEN3_MIN_VERSION}...<${RHBM_GEM_EIGEN3_MAX_EXCLUSIVE_VERSION}")
set(RHBM_GEM_EIGEN3_FETCH_VERSION "5.0.0")
set(RHBM_GEM_EIGEN3_URL
    "https://gitlab.com/libeigen/eigen/-/archive/${RHBM_GEM_EIGEN3_FETCH_VERSION}/eigen-${RHBM_GEM_EIGEN3_FETCH_VERSION}.tar.gz")
set(RHBM_GEM_EIGEN3_URL_HASH "SHA256=315c881e19e17542a7d428c5aa37d113c89b9500d350c433797b730cd449c056")

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

set(RHBM_GEM_UMAPPP_VERSION "3.3.2")
set(RHBM_GEM_UMAPPP_URL "https://github.com/libscran/umappp/archive/refs/tags/v3.3.2.tar.gz")
set(RHBM_GEM_UMAPPP_URL_HASH "SHA256=43504eda3994f13d613d2597ca7a2ce9d0f0b1774e2aea3427f54d9fb4726e88")

set(RHBM_GEM_AARAND_URL "https://github.com/LTLA/aarand/archive/refs/tags/v1.1.0.tar.gz")
set(RHBM_GEM_AARAND_URL_HASH "SHA256=af0bc29e38a02a23a95e0ab988f42510c73fbeb89c6f22162fa6a98f1b863dbe")
set(RHBM_GEM_IRLBA_URL "https://github.com/libscran/irlba/archive/refs/tags/v3.1.0.tar.gz")
set(RHBM_GEM_IRLBA_URL_HASH "SHA256=2648a1be541963a5d3856ef932fc329ebdb61af1a386e3c548518c32bc1ab302")
set(RHBM_GEM_SUBPAR_URL "https://github.com/LTLA/subpar/archive/refs/tags/v0.5.0.tar.gz")
set(RHBM_GEM_SUBPAR_URL_HASH "SHA256=e86fc2af25625653cfaa0eca583f935ae3c5d1104868ff55df508472972cf8ac")
set(RHBM_GEM_SANISIZER_URL "https://github.com/LTLA/sanisizer/archive/refs/tags/v0.2.0.tar.gz")
set(RHBM_GEM_SANISIZER_URL_HASH "SHA256=2b5b5edd304d0c1453615cf490a5b5c16e451c10560876448cf24ed4ba6d0328")
set(RHBM_GEM_KNNCOLLE_URL "https://github.com/knncolle/knncolle/archive/refs/tags/v3.1.0.tar.gz")
set(RHBM_GEM_KNNCOLLE_URL_HASH "SHA256=4bd997de930fc34ac45c6c407ef874c97e1c1750235abb4d6c5e25eab7d28a26")

if(RHBM_GEM_DEP_PROVIDER STREQUAL "FETCH")
    include(FetchContent)
endif()

function(rhbm_gem_declare_content dep_name dep_url dep_hash)
    FetchContent_Declare(${dep_name}
        URL "${dep_url}"
        URL_HASH "${dep_hash}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
endfunction()

function(rhbm_gem_populate_content dep_name dep_url dep_hash out_source_dir)
    rhbm_gem_declare_content(${dep_name} "${dep_url}" "${dep_hash}")
    FetchContent_MakeAvailable(${dep_name})
    FetchContent_GetProperties(${dep_name})
    set(${out_source_dir} "${${dep_name}_SOURCE_DIR}" PARENT_SCOPE)
endfunction()

function(rhbm_gem_validate_eigen3_dependency)
    if(NOT TARGET Eigen3::Eigen)
        message(FATAL_ERROR
            "Eigen3 was found but did not export the Eigen3::Eigen target.")
    endif()
    if(NOT DEFINED Eigen3_VERSION OR "${Eigen3_VERSION}" STREQUAL "")
        message(FATAL_ERROR
            "Eigen3 did not report a version; RHBM-GEM requires "
            "${RHBM_GEM_EIGEN3_VERSION_RANGE}.")
    endif()
    if("${Eigen3_VERSION}" VERSION_LESS "${RHBM_GEM_EIGEN3_MIN_VERSION}")
        message(FATAL_ERROR
            "Eigen3 ${Eigen3_VERSION} is below the supported range "
            "${RHBM_GEM_EIGEN3_VERSION_RANGE}.")
    endif()
    if(NOT "${Eigen3_VERSION}" VERSION_LESS "${RHBM_GEM_EIGEN3_MAX_EXCLUSIVE_VERSION}")
        message(FATAL_ERROR
            "Eigen3 ${Eigen3_VERSION} is outside the supported range "
            "${RHBM_GEM_EIGEN3_VERSION_RANGE}.")
    endif()
endfunction()

function(rhbm_gem_prepare_eigen3_compat_redirect)
    file(WRITE "${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/Eigen3Config.cmake" [=[
if(NOT TARGET Eigen3::Eigen)
    set(Eigen3_FOUND FALSE)
    set(Eigen3_NOT_FOUND_MESSAGE
        "The RHBM-GEM Eigen3 compatibility redirect requires an existing Eigen3::Eigen target.")
endif()
]=])
    write_basic_package_version_file(
        "${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/Eigen3ConfigVersion.cmake"
        VERSION "${Eigen3_VERSION}"
        COMPATIBILITY SameMajorVersion
        ARCH_INDEPENDENT
    )
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
    find_package(Eigen3 ${RHBM_GEM_EIGEN3_VERSION_RANGE} CONFIG REQUIRED)
    rhbm_gem_validate_eigen3_dependency()
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

    if(TARGET SQLite3::SQLite3)
        set(_rhbm_gem_sqlite_target SQLite3::SQLite3)
    elseif(TARGET SQLite::SQLite3)
        set(_rhbm_gem_sqlite_target SQLite::SQLite3)
    else()
        message(FATAL_ERROR
            "System SQLite3 package was found but did not export a supported imported target.")
    endif()

    target_link_libraries(rhbm_gem_dependencies INTERFACE
        Eigen3::Eigen
        CLI11::CLI11
        ${_rhbm_gem_sqlite_target}
    )

    if(RHBM_GEM_ENABLE_UMAP)
        rhbm_gem_prepare_eigen3_compat_redirect()
        message(STATUS "Using system umappp package")
        find_package(libscran_umappp CONFIG REQUIRED)
    endif()
else()
    message(STATUS "Dependency provider: FETCH")

    if(RHBM_GEM_ENABLE_UMAP)
        # Let umappp populate Eigen inside its excluded subtree. This preserves
        # upstream's export ordering without adding Eigen's package rules to
        # the parent installation.
        if(NOT DEFINED CMAKE_EXPORT_PACKAGE_REGISTRY)
            set(CMAKE_EXPORT_PACKAGE_REGISTRY OFF)
        endif()
        set(EIGEN_BUILD_CMAKE_PACKAGE ON CACHE BOOL
            "Build the Eigen CMake package required by umappp dependencies" FORCE)
        rhbm_gem_declare_content(
            eigen "${RHBM_GEM_EIGEN3_URL}" "${RHBM_GEM_EIGEN3_URL_HASH}")
    else()
        set(EIGEN_BUILD_CMAKE_PACKAGE OFF CACHE BOOL
            "Build the Eigen CMake package required by umappp dependencies" FORCE)
        rhbm_gem_populate_content(
            eigen
            "${RHBM_GEM_EIGEN3_URL}"
            "${RHBM_GEM_EIGEN3_URL_HASH}"
            RHBM_GEM_EIGEN3_SOURCE_DIR
        )
    endif()

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

    if(RHBM_GEM_ENABLE_UMAP)
        # umappp's upstream FetchContent declarations use moving branches for
        # most transitive dependencies. Declare the supported releases first
        # so CMake's first-declaration-wins rule keeps this build reproducible.
        rhbm_gem_declare_content(
            aarand "${RHBM_GEM_AARAND_URL}" "${RHBM_GEM_AARAND_URL_HASH}")
        rhbm_gem_declare_content(
            irlba "${RHBM_GEM_IRLBA_URL}" "${RHBM_GEM_IRLBA_URL_HASH}")
        rhbm_gem_declare_content(
            subpar "${RHBM_GEM_SUBPAR_URL}" "${RHBM_GEM_SUBPAR_URL_HASH}")
        rhbm_gem_declare_content(
            sanisizer "${RHBM_GEM_SANISIZER_URL}" "${RHBM_GEM_SANISIZER_URL_HASH}")
        rhbm_gem_declare_content(
            knncolle "${RHBM_GEM_KNNCOLLE_URL}" "${RHBM_GEM_KNNCOLLE_URL_HASH}")

        message(STATUS "Fetching umappp (v${RHBM_GEM_UMAPPP_VERSION}) via FetchContent")
        FetchContent_Populate(umappp
            URL "${RHBM_GEM_UMAPPP_URL}"
            URL_HASH "${RHBM_GEM_UMAPPP_URL_HASH}"
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            SOURCE_DIR "${FETCHCONTENT_BASE_DIR}/umappp-src"
            BINARY_DIR "${FETCHCONTENT_BASE_DIR}/umappp-build"
        )
        add_subdirectory(
            "${umappp_SOURCE_DIR}"
            "${umappp_BINARY_DIR}"
            EXCLUDE_FROM_ALL
        )

        FetchContent_GetProperties(eigen)
        set(RHBM_GEM_EIGEN3_SOURCE_DIR "${eigen_SOURCE_DIR}")
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

if(RHBM_GEM_ENABLE_UMAP)
    if(NOT TARGET libscran::umappp)
        message(FATAL_ERROR
            "UMAP support is enabled, but the libscran::umappp target is unavailable.")
    endif()
    message(STATUS "UMAP support enabled through libscran::umappp")
else()
    message(STATUS "UMAP support disabled by RHBM_GEM_ENABLE_UMAP=OFF")
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
