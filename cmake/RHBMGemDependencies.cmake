include("${CMAKE_CURRENT_LIST_DIR}/RHBMGemThirdPartyGuard.cmake")

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
    rhbm_gem_require_bundled_dependency("Eigen3" "third_party/eigen/CMakeLists.txt")
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
    # Restore BUILD_TESTING to its original state. If it was originally
    # undefined, remove it from the cache so later option() calls can define it.
    if(DEFINED _ORIG_BUILD_TESTING)
        set(BUILD_TESTING "${_ORIG_BUILD_TESTING}" CACHE BOOL "" FORCE)
    else()
        unset(BUILD_TESTING CACHE)
    endif()
endif()

add_library(rhbm_eigen INTERFACE)
if(HAVE_EIGEN3)
    target_link_libraries(rhbm_eigen INTERFACE Eigen3::Eigen)
else()
    target_include_directories(rhbm_eigen SYSTEM INTERFACE
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party/eigen>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    install(
        DIRECTORY "${PROJECT_SOURCE_DIR}/third_party/eigen/Eigen"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
    if(EXISTS "${PROJECT_SOURCE_DIR}/third_party/eigen/unsupported")
        install(
            DIRECTORY "${PROJECT_SOURCE_DIR}/third_party/eigen/unsupported"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        )
    endif()
endif()
install(TARGETS rhbm_eigen EXPORT RHBM_GEMTargets)

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
    rhbm_gem_require_bundled_dependency("CLI11" "third_party/CLI11/CMakeLists.txt")
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
    install(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/third_party/CLI11/include/CLI"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
endif()

install(TARGETS rhbm_cli11 EXPORT RHBM_GEMTargets)

if(BUILD_PYTHON_BINDINGS)
    set(PYBIND11_FINDPYTHON ON)
    if(CMAKE_VERSION VERSION_LESS "3.18")
        find_package(Python REQUIRED COMPONENTS Interpreter Development)
    else()
        find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    endif()

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
        rhbm_gem_require_bundled_dependency("pybind11" "third_party/pybind11/CMakeLists.txt")
        add_subdirectory(third_party/pybind11)
    endif()
endif()
