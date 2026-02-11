cmake_minimum_required(VERSION 3.15)
message(STATUS "----------------------------")
message(STATUS "---   RHBM-GEM Project   ---")
message(STATUS "----------------------------")
project(RHBM_GEM
    VERSION 1.0
    DESCRIPTION RHBM-GEM
    HOMEPAGE_URL xxx
    LANGUAGES C CXX
)

set(PROJECT_DISPLAY_NAME "RHBM-GEM")
include(GNUInstallDirs)

# Keep installed binaries relocatable with project-local shared libraries.
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Unit‑test build option
option(BUILD_TESTING "Build unit tests" ON)
option(ENABLE_COVERAGE "Enable gcov coverage instrumentation" OFF)
option(COVERAGE_INCLUDE_TESTS "Include tests in coverage summary" OFF)

if(ENABLE_COVERAGE AND NOT BUILD_TESTING)
    message(FATAL_ERROR "ENABLE_COVERAGE requires BUILD_TESTING=ON")
endif()

# Specify the C++ standard (Use interface libraries)
add_library(CompilerFlags INTERFACE)
target_compile_features(CompilerFlags INTERFACE cxx_std_17)
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")

if(ENABLE_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        message(FATAL_ERROR
            "ENABLE_COVERAGE supports GNU/Clang/AppleClang only. Detected: ${CMAKE_CXX_COMPILER_ID}")
    endif()
    message(STATUS "Coverage instrumentation enabled")
endif()

target_compile_options(CompilerFlags INTERFACE
    "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wpedantic;-Wconversion;-Wsign-conversion;-Wunused;-Wcast-align;-Wold-style-cast;-Wstrict-overflow;-Wno-psabi>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3;/Zc:__cplusplus>>"
)
# Constrain Eigen headers in Release builds to MPL-2.0-only subsets.
target_compile_definitions(CompilerFlags INTERFACE
    "$<$<CONFIG:Release>:EIGEN_MPL2_ONLY>"
)

if(ENABLE_COVERAGE)
    target_compile_options(CompilerFlags INTERFACE
        -O0
        -g
        --coverage
        -fprofile-arcs
        -ftest-coverage
    )
    target_link_options(CompilerFlags INTERFACE
        --coverage
    )
endif()

option(USE_SYSTEM_LIBS "Use system-installed libraries" ON)

# Locate Eigen package
if(USE_SYSTEM_LIBS)
    find_package(Eigen3 QUIET)
    if(Eigen3_FOUND)
        message(STATUS "Using system Eigen3")
        set(HAVE_EIGEN3 TRUE)
    else()
        message(STATUS "System Eigen3 not found, will use third_party version")
        set(HAVE_EIGEN3 FALSE)
    endif()
else()
    message(STATUS "Forced to use third_party libraries")
    set(HAVE_EIGEN3 FALSE)
endif()

if(NOT HAVE_EIGEN3)
    # Disable Eigen's own tests to avoid policy warnings
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "Disable Eigen tests" FORCE)
    # Temporarily store the original BUILD_TESTING value and turn it off to
    # avoid Eigen's test configuration interfering with ours.
    if(DEFINED BUILD_TESTING)
        set(_ORIG_BUILD_TESTING "${BUILD_TESTING}")
    else()
        unset(_ORIG_BUILD_TESTING)
    endif()
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    add_subdirectory(third_party/eigen)
    # Restore BUILD_TESTING to its original state.  If it was originally
    # undefined, remove it from the cache so that later calls such as
    # include(CTest) or option(BUILD_TESTING ...) can define it without being
    # forced to our temporary value.
    if(DEFINED _ORIG_BUILD_TESTING)
        set(BUILD_TESTING "${_ORIG_BUILD_TESTING}" CACHE BOOL "" FORCE)
    else()
        unset(BUILD_TESTING CACHE)
    endif()
endif()

# Locate CLI11 package
if(USE_SYSTEM_LIBS)
    find_package(CLI11 QUIET)
    if(CLI11_FOUND)
        message(STATUS "Using system CLI11")
        set(HAVE_CLI11 TRUE)
    else()
        message(STATUS "System CLI11 not found, will use third_party version")
        set(HAVE_CLI11 FALSE)
    endif()
else()
    set(HAVE_CLI11 FALSE)
endif()

if(NOT HAVE_CLI11)
    add_subdirectory(third_party/CLI11)
endif()

# Locate Pybind11 package
if(USE_SYSTEM_LIBS)
    set(PYBIND11_FINDPYTHON ON)
    find_package(PYBIND11 QUIET)
    if(PYBIND11_FOUND)
        message(STATUS "Using system Pybind11")
        set(HAVE_PYBIND11 TRUE)
    else()
        message(STATUS "System Pybind11 not found, will use third_party version")
        set(HAVE_PYBIND11 FALSE)
    endif()
else()
    set(HAVE_PYBIND11 FALSE)
endif()

if(NOT HAVE_PYBIND11)
    add_subdirectory(third_party/pybind11)
endif()

# Control where the static and shared libraries are built so that on windows
# we don't need to tinker with the path to run the executable
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# 
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_subdirectory(src)
add_subdirectory(bindings)

# Install public headers and legal notices.
install(DIRECTORY "${CMAKE_SOURCE_DIR}/include/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(
    FILES
        "${CMAKE_SOURCE_DIR}/LICENSE"
        "${CMAKE_SOURCE_DIR}/THIRD_PARTY_NOTICES.md"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}"
)

# -------------------------------------------------
# Unit tests (GoogleTest via FetchContent)
# -------------------------------------------------
if(BUILD_TESTING)
    include(CTest)
    enable_testing()
    # Do not export test framework artifacts in normal project installs.
    set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest install targets" FORCE)
    set(INSTALL_GMOCK OFF CACHE BOOL "Disable GoogleMock install targets" FORCE)

    if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest/CMakeLists.txt")
        message(STATUS "Using bundled GoogleTest")
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        add_subdirectory(third_party/googletest)
        set(GTEST_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/googletest/googletest")
    else()
        message(STATUS "Fetching GoogleTest")
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        v1.17.0
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
        set(GTEST_SOURCE_DIR "${googletest_SOURCE_DIR}/googletest")
    endif()

    if(DEFINED GTEST_SOURCE_DIR)
        list(APPEND CMAKE_MODULE_PATH "${GTEST_SOURCE_DIR}/cmake")
    endif()

    include(GoogleTest)
    add_subdirectory(tests)    # all test sources live here

    if(ENABLE_COVERAGE)
        find_program(GCOV_EXECUTABLE NAMES gcov)
        if(NOT GCOV_EXECUTABLE)
            message(FATAL_ERROR "ENABLE_COVERAGE is ON but gcov was not found in PATH")
        endif()

        set(COVERAGE_OUTPUT_DIR "${CMAKE_BINARY_DIR}/coverage")
        add_custom_target(coverage
            COMMAND ${CMAKE_COMMAND}
                -DCOVERAGE_ACTION=clean
                -DCOVERAGE_BUILD_DIR:PATH=${CMAKE_BINARY_DIR}
                -P ${CMAKE_SOURCE_DIR}/cmake/GcovCoverageReport.cmake
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
            COMMAND ${CMAKE_COMMAND}
                -DCOVERAGE_ACTION=report
                -DCOVERAGE_BUILD_DIR:PATH=${CMAKE_BINARY_DIR}
                -DCOVERAGE_PROJECT_SOURCE_DIR:PATH=${CMAKE_SOURCE_DIR}
                -DCOVERAGE_OUTPUT_DIR:PATH=${COVERAGE_OUTPUT_DIR}
                -DCOVERAGE_INCLUDE_TESTS:BOOL=${COVERAGE_INCLUDE_TESTS}
                -DGCOV_EXECUTABLE:FILEPATH=${GCOV_EXECUTABLE}
                -P ${CMAKE_SOURCE_DIR}/cmake/GcovCoverageReport.cmake
            BYPRODUCTS
                "${COVERAGE_OUTPUT_DIR}/coverage_summary.txt"
                "${COVERAGE_OUTPUT_DIR}/coverage_detail.txt"
            DEPENDS unit_tests
            USES_TERMINAL
            COMMENT "Running tests and generating gcov text coverage reports"
        )
    endif()
endif()
