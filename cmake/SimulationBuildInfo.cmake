# Build-time fingerprints include uncommitted source edits, without depending on a Git checkout.
file(GLOB_RECURSE _simulation_sources CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/*.cpp" "${PROJECT_SOURCE_DIR}/src/*.hpp"
    "${PROJECT_SOURCE_DIR}/include/*.hpp" "${PROJECT_SOURCE_DIR}/cmake/*.cmake"
    "${PROJECT_SOURCE_DIR}/src/CMakeLists.txt")
list(APPEND _simulation_sources "${PROJECT_SOURCE_DIR}/CMakeLists.txt")
list(SORT _simulation_sources)
set(_simulation_build_dir "${PROJECT_BINARY_DIR}/generated/$<CONFIG>")
file(GENERATE OUTPUT "${_simulation_build_dir}/SimulationSources.txt"
    CONTENT "$<JOIN:${_simulation_sources},\n>\n")
file(GENERATE OUTPUT "${_simulation_build_dir}/SimulationConfiguration-$<COMPILE_LANGUAGE>.txt" CONTENT
"version=${PROJECT_VERSION}
compiler=${CMAKE_CXX_COMPILER};${CMAKE_CXX_COMPILER_ID};${CMAKE_CXX_COMPILER_VERSION}
system=${CMAKE_SYSTEM_NAME};${CMAKE_SYSTEM_PROCESSOR}
configuration=$<CONFIG>
flags=${CMAKE_CXX_FLAGS}
debug_flags=${CMAKE_CXX_FLAGS_DEBUG}
release_flags=${CMAKE_CXX_FLAGS_RELEASE}
relwithdebinfo_flags=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
minsizerel_flags=${CMAKE_CXX_FLAGS_MINSIZEREL}
boost=${Boost_VERSION_STRING}
eigen=${Eigen3_VERSION}
joint_sparse_backend=${RHBM_GEM_JOINT_SPARSE_BACKEND}
spqr=${SPQR_VERSION}
cholmod=${CHOLMOD_VERSION}
suitesparse_config=${SUITESPARSE_CONFIG_VERSION}
spqr_transitive_link=${RHBM_GEM_SPQR_BLAS_LINK}
spqr_blas_version=${RHBM_GEM_SPQR_BLAS_VERSION}
openmp=${RHBM_GEM_WITH_OPENMP};${OpenMP_CXX_VERSION};${OpenMP_CXX_FLAGS}
root=${RHBM_GEM_WITH_ROOT}
dependency_provider=${RHBM_GEM_DEP_PROVIDER}
compile_features=$<TARGET_PROPERTY:rhbm_gem,COMPILE_FEATURES>
compile_options=$<TARGET_GENEX_EVAL:rhbm_gem,$<TARGET_PROPERTY:rhbm_gem,COMPILE_OPTIONS>>
compile_definitions=$<TARGET_GENEX_EVAL:rhbm_gem,$<TARGET_PROPERTY:rhbm_gem,COMPILE_DEFINITIONS>>
include_directories=$<TARGET_GENEX_EVAL:rhbm_gem,$<TARGET_PROPERTY:rhbm_gem,INCLUDE_DIRECTORIES>>
link_options=$<TARGET_GENEX_EVAL:rhbm_gem,$<TARGET_PROPERTY:rhbm_gem,LINK_OPTIONS>>
" TARGET rhbm_gem)
add_custom_command(OUTPUT "${_simulation_build_dir}/SimulationBuildInfo.hpp"
    COMMAND ${CMAKE_COMMAND}
        "-DSOURCE_ROOT=${PROJECT_SOURCE_DIR}"
        "-DINPUT_LIST=${_simulation_build_dir}/SimulationSources.txt"
        "-DCONFIG_FILE=${_simulation_build_dir}/SimulationConfiguration-CXX.txt"
        "-DOUTPUT_FILE=${_simulation_build_dir}/SimulationBuildInfo.hpp"
        "-DPROJECT_VERSION=${PROJECT_VERSION}"
        -P "${PROJECT_SOURCE_DIR}/cmake/GenerateSimulationBuildInfo.cmake"
    DEPENDS ${_simulation_sources}
        "${_simulation_build_dir}/SimulationSources.txt"
        "${_simulation_build_dir}/SimulationConfiguration-CXX.txt"
    VERBATIM)
target_sources(rhbm_gem PRIVATE "${_simulation_build_dir}/SimulationBuildInfo.hpp")
target_include_directories(rhbm_gem PRIVATE "${_simulation_build_dir}")
