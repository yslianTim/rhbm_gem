# Core build options
option(BUILD_TESTING "Build unit tests" ON)
option(ENABLE_COVERAGE "Enable gcov coverage instrumentation" OFF)
option(COVERAGE_INCLUDE_TESTS "Include tests in coverage summary" OFF)
option(BUILD_PYTHON_BINDINGS "Build pybind11 Python extension module" ON)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(RHBM_GEM_BUILD_GUI
    "Attempt to build the Qt6 GUI executable when Qt6 is available" ON)
option(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS
    "Enable experimental bond-analysis workflow in PotentialAnalysisCommand" OFF)
option(RHBM_GEM_LEGACY_V1_SUPPORT
    "Enable migration support for legacy v1 SQLite schema" ON)

set(RHBM_GEM_DEP_PROVIDER "SYSTEM" CACHE STRING
    "Dependency provider mode: SYSTEM or FETCH")
set_property(CACHE RHBM_GEM_DEP_PROVIDER PROPERTY STRINGS SYSTEM FETCH)

set(RHBM_GEM_OPENMP_MODE "AUTO" CACHE STRING "OpenMP feature mode: AUTO, ON, or OFF")
set_property(CACHE RHBM_GEM_OPENMP_MODE PROPERTY STRINGS AUTO ON OFF)
set(RHBM_GEM_ROOT_MODE "AUTO" CACHE STRING "ROOT feature mode: AUTO, ON, or OFF")
set_property(CACHE RHBM_GEM_ROOT_MODE PROPERTY STRINGS AUTO ON OFF)
set(RHBM_GEM_PYTHON_INSTALL_LAYOUT "SITE_PREFIX" CACHE STRING
    "Python extension install layout: SITE_PREFIX or LIBDIR")
set_property(CACHE RHBM_GEM_PYTHON_INSTALL_LAYOUT PROPERTY STRINGS SITE_PREFIX LIBDIR)
set(RHBM_GEM_PYTHON_INSTALL_DIR "" CACHE PATH
    "Override installation directory for Python extension module (relative or absolute path)")

set(_RHBM_GEM_VALID_DEP_PROVIDERS SYSTEM FETCH)
string(TOUPPER "${RHBM_GEM_DEP_PROVIDER}" RHBM_GEM_DEP_PROVIDER)
if(NOT RHBM_GEM_DEP_PROVIDER IN_LIST _RHBM_GEM_VALID_DEP_PROVIDERS)
    message(FATAL_ERROR
        "Invalid RHBM_GEM_DEP_PROVIDER='${RHBM_GEM_DEP_PROVIDER}'. Valid values: SYSTEM, FETCH.")
endif()
set(RHBM_GEM_DEP_PROVIDER "${RHBM_GEM_DEP_PROVIDER}" CACHE STRING
    "Dependency provider mode: SYSTEM or FETCH" FORCE)
unset(_RHBM_GEM_VALID_DEP_PROVIDERS)

set(_RHBM_GEM_VALID_FEATURE_MODES AUTO ON OFF)
string(TOUPPER "${RHBM_GEM_OPENMP_MODE}" RHBM_GEM_OPENMP_MODE)
if(NOT RHBM_GEM_OPENMP_MODE IN_LIST _RHBM_GEM_VALID_FEATURE_MODES)
    message(FATAL_ERROR
        "Invalid RHBM_GEM_OPENMP_MODE='${RHBM_GEM_OPENMP_MODE}'. Valid values: AUTO, ON, OFF.")
endif()
set(RHBM_GEM_OPENMP_MODE "${RHBM_GEM_OPENMP_MODE}" CACHE STRING
    "OpenMP feature mode: AUTO, ON, or OFF" FORCE)

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

add_library(CompilerFlags INTERFACE)
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
