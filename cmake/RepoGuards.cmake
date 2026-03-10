if(NOT DEFINED PROJECT_SOURCE_DIR)
    message(FATAL_ERROR "PROJECT_SOURCE_DIR is required")
endif()
if(NOT DEFINED PROJECT_BINARY_DIR)
    message(FATAL_ERROR "PROJECT_BINARY_DIR is required")
endif()
if(NOT DEFINED RHBM_GEM_INCLUDE_INSTALL_SMOKE)
    set(RHBM_GEM_INCLUDE_INSTALL_SMOKE ON)
endif()

set(_RHBM_GEM_GUARD_SCRIPTS
    LoggingStyleGuard.cmake
    ParameterlessVoidStyleGuard.cmake
    RepositoryHygieneGuard.cmake
    TestFixtureTrackingGuard.cmake
    StructureGuard.cmake
    AbsolutePathGuard.cmake
)

foreach(_guard_script IN LISTS _RHBM_GEM_GUARD_SCRIPTS)
    message(STATUS "Running guard: ${_guard_script}")
    include("${PROJECT_SOURCE_DIR}/cmake/${_guard_script}")
endforeach()

if(RHBM_GEM_INCLUDE_INSTALL_SMOKE)
    message(STATUS "Running guard: InstallConsumerSmoke.cmake")
    include("${PROJECT_SOURCE_DIR}/cmake/InstallConsumerSmoke.cmake")
endif()

unset(_guard_script)
unset(_RHBM_GEM_GUARD_SCRIPTS)
