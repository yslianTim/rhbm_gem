file(GLOB SOURCE *.cpp)

pybind11_add_module(cpp_module ${SOURCE})


target_link_libraries(cpp_module PRIVATE core data utils CompilerFlags)

set_target_properties(cpp_module PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)