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
option(RHBM_GEM_ENABLE_UMAP
    "Enable the umap_embedding command through umappp" ON)
option(RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    "Enable experimental features across the project" OFF)
option(RHBM_GEM_ENABLE_FOLD_168_REGRESSION
    "Enable the external 168-atom simulation regression benchmark" OFF)
option(RHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS
    "Enable frozen 168-atom and additional joint component regressions" OFF)
option(RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS
    "Build independent joint component derivative and certification tools" OFF)
foreach(retired_option IN ITEMS RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT)
    if(DEFINED ${retired_option})
        message(FATAL_ERROR "${retired_option} is retired. Remove both old options from your command and cache (cmake -U RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE -U RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT), then use RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT.")
    endif()
endforeach()
option(RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT
    "Record existing second-stage decision evidence without additional solving" OFF)

set(RHBM_GEM_FOLD_168_MODEL "" CACHE FILEPATH
    "Path to the fixed fold-168 regression CIF model")
set(RHBM_GEM_FOLD_168_MAP "" CACHE FILEPATH
    "Path to the fixed fold-168 regression map")

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

if(RHBM_GEM_ENABLE_FOLD_168_REGRESSION AND NOT BUILD_TESTING)
    message(FATAL_ERROR
        "RHBM_GEM_ENABLE_FOLD_168_REGRESSION requires BUILD_TESTING=ON")
endif()

if((RHBM_GEM_ENABLE_JOINT_EXTENDED_TESTS OR RHBM_GEM_ENABLE_JOINT_OFFLINE_AUDITS) AND NOT BUILD_TESTING)
    message(FATAL_ERROR "Joint component test options require BUILD_TESTING=ON")
endif()

function(rhbm_gem_apply_observation_definitions target)
    if(RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT)
        target_compile_definitions(${target} PRIVATE RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT)
    endif()
endfunction()

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
