## Top-level CMake configuration for the test suite
include(GoogleTest)

add_subdirectory(core)
add_subdirectory(data)
add_subdirectory(utils)
