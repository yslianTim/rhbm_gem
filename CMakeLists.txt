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

# Locate OpenMP package
find_package(OpenMP REQUIRED)

# Locate SQLite package
find_package(SQLite3 REQUIRED)

# Locate ROOT package
find_package(ROOT 6.28 REQUIRED COMPONENTS Core Hist Gpad RIO Tree)
list(APPEND ROOT_LIBS ROOT::Core ROOT::Hist ROOT::Gpad ROOT::RIO ROOT::Tree)

# Control where the static and shared libraries are built so that on windows
# we don't need to tinker with the path to run the executable
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

# 
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(third_party/CLI11)
add_subdirectory(third_party/eigen)

add_subdirectory(src)