function(rhbm_gem_normalize_cache_enum cache_var description)
    set(_rhbm_gem_valid_values ${ARGN})
    set_property(CACHE ${cache_var} PROPERTY STRINGS ${_rhbm_gem_valid_values})

    string(TOUPPER "${${cache_var}}" _rhbm_gem_normalized_value)
    if(NOT _rhbm_gem_normalized_value IN_LIST _rhbm_gem_valid_values)
        string(REPLACE ";" ", " _rhbm_gem_valid_values_text "${_rhbm_gem_valid_values}")
        message(FATAL_ERROR
            "Invalid ${cache_var}='${${cache_var}}'. Valid values: ${_rhbm_gem_valid_values_text}.")
    endif()

    set(${cache_var} "${_rhbm_gem_normalized_value}" CACHE STRING "${description}" FORCE)
endfunction()

# Core build options
option(ENABLE_COVERAGE "Enable gcov coverage instrumentation" OFF)
option(COVERAGE_INCLUDE_TESTS "Include tests in coverage summary" OFF)
option(BUILD_PYTHON_BINDINGS "Build pybind11 Python extension module" ON)
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(RHBM_GEM_BUILD_GUI
    "Attempt to build the Qt6 GUI executable when Qt6 is available" ON)
option(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS
    "Enable experimental bond-analysis workflow in PotentialAnalysisCommand" OFF)

set(RHBM_GEM_DEP_PROVIDER "SYSTEM" CACHE STRING
    "Dependency provider mode: SYSTEM or FETCH")
set(RHBM_GEM_OPENMP_MODE "AUTO" CACHE STRING "OpenMP feature mode: AUTO, ON, or OFF")
set(RHBM_GEM_ROOT_MODE "AUTO" CACHE STRING "ROOT feature mode: AUTO, ON, or OFF")
set(RHBM_GEM_PYTHON_INSTALL_LAYOUT "SITE_PREFIX" CACHE STRING
    "Python extension install layout: SITE_PREFIX or LIBDIR")
set(RHBM_GEM_PYTHON_INSTALL_DIR "" CACHE PATH
    "Override installation directory for Python extension module (relative or absolute path)")

rhbm_gem_normalize_cache_enum(
    RHBM_GEM_DEP_PROVIDER
    "Dependency provider mode: SYSTEM or FETCH"
    SYSTEM FETCH
)
rhbm_gem_normalize_cache_enum(
    RHBM_GEM_OPENMP_MODE
    "OpenMP feature mode: AUTO, ON, or OFF"
    AUTO ON OFF
)
rhbm_gem_normalize_cache_enum(
    RHBM_GEM_ROOT_MODE
    "ROOT feature mode: AUTO, ON, or OFF"
    AUTO ON OFF
)
rhbm_gem_normalize_cache_enum(
    RHBM_GEM_PYTHON_INSTALL_LAYOUT
    "Python extension install layout: SITE_PREFIX or LIBDIR"
    SITE_PREFIX LIBDIR
)

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
