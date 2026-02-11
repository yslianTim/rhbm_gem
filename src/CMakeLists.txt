add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(data)

set(MAIN_EXEC_TARGET ${PROJECT_NAME})
add_executable(${MAIN_EXEC_TARGET} main.cpp)
set_target_properties(${MAIN_EXEC_TARGET} PROPERTIES
    OUTPUT_NAME "${PROJECT_DISPLAY_NAME}"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

target_link_libraries(${MAIN_EXEC_TARGET} PUBLIC
    core
    utils
    data
)
target_link_libraries(${MAIN_EXEC_TARGET} PRIVATE
    CLI11::CLI11
    CompilerFlags
)

target_include_directories(${MAIN_EXEC_TARGET} PUBLIC "${PROJECT_BINARY_DIR}")

# Install rules
install(
    TARGETS ${MAIN_EXEC_TARGET}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
