cmake_minimum_required(VERSION 3.15)
message(STATUS "-----------------------------")
message(STATUS "-- EP Map Analysis Project --")
message(STATUS "-----------------------------")
project(PotentialMap
    VERSION 1.0
    DESCRIPTION TestProject
    HOMEPAGE_URL xxx
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify the C++ standard (Use interface libraries)
add_library(CompilerFlags INTERFACE)
target_compile_features(CompilerFlags INTERFACE cxx_std_17)
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(CompilerFlags INTERFACE
    "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wpedantic;-Wconversion;-Wsign-conversion;-Wunused;-Wcast-align;-Wold-style-cast;-Wstrict-overflow;-Wno-psabi>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3;/Zc:__cplusplus>>"
)

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
    # Temporarily disable global testing to avoid Eigen's test configuration
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    add_subdirectory(third_party/eigen)
    set(BUILD_TESTING ON CACHE BOOL "" FORCE)
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

# -------------------------------------------------
# Unit‑test build option
# -------------------------------------------------
option(BUILD_TESTING "Build unit tests" ON)

add_subdirectory(src)
add_subdirectory(bindings)

# -------------------------------------------------
# Unit tests (GoogleTest via FetchContent)
# -------------------------------------------------
if(BUILD_TESTING)
    include(CTest)
    enable_testing()

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
endif()
