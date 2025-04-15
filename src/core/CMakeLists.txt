# File dependencies
file(GLOB SOURCE *.cpp)


# 建立 core 模組
add_library(core SHARED ${SOURCE})
target_include_directories(core PUBLIC ${CMAKE_SOURCE_DIR}/include/core)
target_link_libraries(core PUBLIC
    CLI11::CLI11
    data
    CompilerFlags
)

set_target_properties(core PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)