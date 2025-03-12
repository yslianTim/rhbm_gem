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

# Specify the C++ standard (Use interface libraries)
add_library(CompilerFlags INTERFACE)
target_compile_features(CompilerFlags INTERFACE cxx_std_17)
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(CompilerFlags INTERFACE
    "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused;-Wno-psabi>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
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
    add_subdirectory(third_party/Eigen)
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

# Control where the static and shared libraries are built so that on windows
# we don't need to tinker with the path to run the executable
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

# 
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(src)