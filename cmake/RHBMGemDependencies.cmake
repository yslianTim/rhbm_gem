include("${CMAKE_CURRENT_LIST_DIR}/RHBMGemFetchDeps.cmake")

# Locate Eigen package
if(USE_SYSTEM_LIBS)
    find_package(Eigen3 QUIET)
    if(Eigen3_FOUND)
        message(STATUS "Using system Eigen3")
        set(HAVE_EIGEN3 TRUE)
    else()
        message(STATUS "System Eigen3 not found, will fetch via FetchContent")
        set(HAVE_EIGEN3 FALSE)
    endif()
else()
    message(STATUS "Forced to use FetchContent for Eigen3")
    set(HAVE_EIGEN3 FALSE)
endif()

if(NOT HAVE_EIGEN3)
    rhbm_gem_fetch_eigen3(RHBM_GEM_EIGEN3_SOURCE_DIR)
endif()

add_library(rhbm_eigen INTERFACE)
if(HAVE_EIGEN3)
    target_link_libraries(rhbm_eigen INTERFACE Eigen3::Eigen)
else()
    target_include_directories(rhbm_eigen SYSTEM INTERFACE
        "$<BUILD_INTERFACE:${RHBM_GEM_EIGEN3_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    install(
        DIRECTORY "${RHBM_GEM_EIGEN3_SOURCE_DIR}/Eigen"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
    if(EXISTS "${RHBM_GEM_EIGEN3_SOURCE_DIR}/unsupported")
        install(
            DIRECTORY "${RHBM_GEM_EIGEN3_SOURCE_DIR}/unsupported"
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
        message(STATUS "System CLI11 not found, will fetch via FetchContent")
        set(HAVE_CLI11 FALSE)
        set(RHBM_GEM_USE_SYSTEM_CLI11 FALSE)
    endif()
else()
    message(STATUS "Forced to use FetchContent for CLI11")
    set(HAVE_CLI11 FALSE)
    set(RHBM_GEM_USE_SYSTEM_CLI11 FALSE)
endif()

if(NOT HAVE_CLI11)
    rhbm_gem_fetch_cli11(RHBM_GEM_CLI11_SOURCE_DIR)
endif()

add_library(rhbm_cli11 INTERFACE)
if(HAVE_CLI11)
    target_link_libraries(rhbm_cli11 INTERFACE CLI11::CLI11)
else()
    target_include_directories(rhbm_cli11 INTERFACE
        "$<BUILD_INTERFACE:${RHBM_GEM_CLI11_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    install(
        DIRECTORY "${RHBM_GEM_CLI11_SOURCE_DIR}/include/CLI"
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
            message(STATUS "System Pybind11 not found, will fetch via FetchContent")
            set(HAVE_PYBIND11 FALSE)
        endif()
    else()
        message(STATUS "Forced to use FetchContent for Pybind11")
        set(HAVE_PYBIND11 FALSE)
    endif()

    if(NOT HAVE_PYBIND11)
        rhbm_gem_fetch_pybind11()
    endif()
endif()
