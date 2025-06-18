add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(data)

set(MAIN_EXEC_NAME ${CMAKE_PROJECT_NAME}_run)
add_executable(${MAIN_EXEC_NAME} main.cpp)

target_link_libraries(${MAIN_EXEC_NAME} PUBLIC
    core
    utils
    data
    CompilerFlags
)

target_include_directories(${MAIN_EXEC_NAME} PUBLIC "${PROJECT_BINARY_DIR}")

# Install rules
install(TARGETS ${MAIN_EXEC_NAME} DESTINATION bin)
