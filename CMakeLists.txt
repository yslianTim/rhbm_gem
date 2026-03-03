cmake_minimum_required(VERSION 3.18)
message(STATUS "----------------------------")
message(STATUS "---   RHBM-GEM Project   ---")
message(STATUS "----------------------------")
project(RHBM_GEM
    VERSION 1.0
    DESCRIPTION RHBM-GEM
    HOMEPAGE_URL https://github.com/yslianTim/em_map_analysis
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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Unit‑test build option
option(BUILD_TESTING "Build unit tests" ON)
option(ENABLE_COVERAGE "Enable gcov coverage instrumentation" OFF)
option(COVERAGE_INCLUDE_TESTS "Include tests in coverage summary" OFF)
option(BUILD_PYTHON_BINDINGS "Build pybind11 Python extension module" ON)
option(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS
    "Enable experimental bond-analysis workflow in PotentialAnalysisCommand" OFF)
set(RHBM_GEM_OPENMP_MODE "AUTO" CACHE STRING "OpenMP feature mode: AUTO, ON, or OFF")
set_property(CACHE RHBM_GEM_OPENMP_MODE PROPERTY STRINGS AUTO ON OFF)
set(RHBM_GEM_BOOST_MODE "AUTO" CACHE STRING "Boost feature mode: AUTO, ON, or OFF")
set_property(CACHE RHBM_GEM_BOOST_MODE PROPERTY STRINGS AUTO ON OFF)
set(RHBM_GEM_BOOST_COMPONENTS "" CACHE STRING
    "Semicolon-separated Boost components required when Boost support is enabled (for example filesystem;system)")
set(RHBM_GEM_ROOT_MODE "AUTO" CACHE STRING "ROOT feature mode: AUTO, ON, or OFF")
set_property(CACHE RHBM_GEM_ROOT_MODE PROPERTY STRINGS AUTO ON OFF)
set(RHBM_GEM_PYTHON_INSTALL_LAYOUT "SITE_PREFIX" CACHE STRING
    "Python extension install layout: SITE_PREFIX or LIBDIR")
set_property(CACHE RHBM_GEM_PYTHON_INSTALL_LAYOUT PROPERTY STRINGS SITE_PREFIX LIBDIR)
set(RHBM_GEM_PYTHON_INSTALL_DIR "" CACHE PATH
    "Override installation directory for Python extension module (relative or absolute path)")

set(_RHBM_GEM_VALID_FEATURE_MODES AUTO ON OFF)
string(TOUPPER "${RHBM_GEM_OPENMP_MODE}" RHBM_GEM_OPENMP_MODE)
if(NOT RHBM_GEM_OPENMP_MODE IN_LIST _RHBM_GEM_VALID_FEATURE_MODES)
    message(FATAL_ERROR
        "Invalid RHBM_GEM_OPENMP_MODE='${RHBM_GEM_OPENMP_MODE}'. Valid values: AUTO, ON, OFF.")
endif()
set(RHBM_GEM_OPENMP_MODE "${RHBM_GEM_OPENMP_MODE}" CACHE STRING
    "OpenMP feature mode: AUTO, ON, or OFF" FORCE)

string(TOUPPER "${RHBM_GEM_BOOST_MODE}" RHBM_GEM_BOOST_MODE)
if(NOT RHBM_GEM_BOOST_MODE IN_LIST _RHBM_GEM_VALID_FEATURE_MODES)
    message(FATAL_ERROR
        "Invalid RHBM_GEM_BOOST_MODE='${RHBM_GEM_BOOST_MODE}'. Valid values: AUTO, ON, OFF.")
endif()
set(RHBM_GEM_BOOST_MODE "${RHBM_GEM_BOOST_MODE}" CACHE STRING
    "Boost feature mode: AUTO, ON, or OFF" FORCE)

string(TOUPPER "${RHBM_GEM_ROOT_MODE}" RHBM_GEM_ROOT_MODE)
if(NOT RHBM_GEM_ROOT_MODE IN_LIST _RHBM_GEM_VALID_FEATURE_MODES)
    message(FATAL_ERROR
        "Invalid RHBM_GEM_ROOT_MODE='${RHBM_GEM_ROOT_MODE}'. Valid values: AUTO, ON, OFF.")
endif()
set(RHBM_GEM_ROOT_MODE "${RHBM_GEM_ROOT_MODE}" CACHE STRING
    "ROOT feature mode: AUTO, ON, or OFF" FORCE)
unset(_RHBM_GEM_VALID_FEATURE_MODES)

set(_RHBM_GEM_VALID_PY_LAYOUTS SITE_PREFIX LIBDIR)
string(TOUPPER "${RHBM_GEM_PYTHON_INSTALL_LAYOUT}" RHBM_GEM_PYTHON_INSTALL_LAYOUT)
if(NOT RHBM_GEM_PYTHON_INSTALL_LAYOUT IN_LIST _RHBM_GEM_VALID_PY_LAYOUTS)
    message(FATAL_ERROR
        "Invalid RHBM_GEM_PYTHON_INSTALL_LAYOUT='${RHBM_GEM_PYTHON_INSTALL_LAYOUT}'. "
        "Valid values: SITE_PREFIX, LIBDIR.")
endif()
set(RHBM_GEM_PYTHON_INSTALL_LAYOUT "${RHBM_GEM_PYTHON_INSTALL_LAYOUT}" CACHE STRING
    "Python extension install layout: SITE_PREFIX or LIBDIR" FORCE)
unset(_RHBM_GEM_VALID_PY_LAYOUTS)

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
        set(RHBM_GEM_USE_SYSTEM_CLI11 TRUE)
    else()
        message(STATUS "System CLI11 not found, will use third_party version")
        set(HAVE_CLI11 FALSE)
        set(RHBM_GEM_USE_SYSTEM_CLI11 FALSE)
    endif()
else()
    set(HAVE_CLI11 FALSE)
    set(RHBM_GEM_USE_SYSTEM_CLI11 FALSE)
endif()

if(NOT HAVE_CLI11)
    add_subdirectory(third_party/CLI11)
endif()

add_library(rhbm_cli11 INTERFACE)
if(HAVE_CLI11)
    target_link_libraries(rhbm_cli11 INTERFACE CLI11::CLI11)
else()
    target_include_directories(rhbm_cli11 INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/CLI11/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/third_party/CLI11/include/CLI" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()

install(TARGETS rhbm_cli11 EXPORT RHBM_GEMTargets)

if(BUILD_PYTHON_BINDINGS)
    set(PYBIND11_FINDPYTHON ON)
    if(CMAKE_VERSION VERSION_LESS "3.18")
        find_package(Python REQUIRED COMPONENTS Interpreter Development)
    else()
        find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    endif()

    # Locate Pybind11 package
    if(USE_SYSTEM_LIBS)
        find_package(pybind11 CONFIG QUIET)
        if(pybind11_FOUND)
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
endif()

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# 
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_subdirectory(src)
if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(bindings)
endif()

if(NOT DEFINED RHBM_GEM_WITH_OPENMP)
    set(RHBM_GEM_WITH_OPENMP FALSE)
endif()
if(NOT DEFINED RHBM_GEM_WITH_BOOST)
    set(RHBM_GEM_WITH_BOOST FALSE)
endif()
if(NOT DEFINED RHBM_GEM_WITH_ROOT)
    set(RHBM_GEM_WITH_ROOT FALSE)
endif()
if(NOT DEFINED RHBM_GEM_USE_SYSTEM_SQLITE3)
    set(RHBM_GEM_USE_SYSTEM_SQLITE3 FALSE)
endif()
if(NOT DEFINED RHBM_GEM_USE_SYSTEM_CLI11)
    set(RHBM_GEM_USE_SYSTEM_CLI11 FALSE)
endif()

# Install public headers and legal notices.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(
    FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
        "${CMAKE_CURRENT_SOURCE_DIR}/THIRD_PARTY_NOTICES.md"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}"
)
install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/examples/python/"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/examples/python"
    FILES_MATCHING
    PATTERN "*.py"
)
install(
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test_model.cif"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/examples/python/data"
)

set(RHBM_GEM_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RHBM_GEMConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfig.cmake"
    INSTALL_DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(
    EXPORT RHBM_GEMTargets
    FILE RHBM_GEMTargets.cmake
    NAMESPACE RHBM_GEM::
    DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfigVersion.cmake"
    DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
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

        set(COVERAGE_OUTPUT_DIR "${PROJECT_BINARY_DIR}/coverage")
        add_custom_target(coverage
            COMMAND ${CMAKE_COMMAND}
                -DCOVERAGE_ACTION=clean
                -DCOVERAGE_BUILD_DIR:PATH=${PROJECT_BINARY_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GcovCoverageReport.cmake
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
            COMMAND ${CMAKE_COMMAND}
                -DCOVERAGE_ACTION=report
                -DCOVERAGE_BUILD_DIR:PATH=${PROJECT_BINARY_DIR}
                -DCOVERAGE_PROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
                -DCOVERAGE_OUTPUT_DIR:PATH=${COVERAGE_OUTPUT_DIR}
                -DCOVERAGE_INCLUDE_TESTS:BOOL=${COVERAGE_INCLUDE_TESTS}
                -DGCOV_EXECUTABLE:FILEPATH=${GCOV_EXECUTABLE}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GcovCoverageReport.cmake
            BYPRODUCTS
                "${COVERAGE_OUTPUT_DIR}/coverage_summary.txt"
                "${COVERAGE_OUTPUT_DIR}/coverage_detail.txt"
            DEPENDS unit_tests
            USES_TERMINAL
            COMMENT "Running tests and generating gcov text coverage reports"
        )
    endif()
endif()
