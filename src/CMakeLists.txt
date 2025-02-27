# 加入各模組子目錄
add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(data)

# 建立主執行檔
set(MAIN_EXEC_NAME ${CMAKE_PROJECT_NAME}_run)
add_executable(${MAIN_EXEC_NAME} main.cpp)

# 將主執行檔連結至核心與工具模組
target_link_libraries(${MAIN_EXEC_NAME} PUBLIC
    core
    ${EIGEN_LIBS}
    ${ROOT_LIBS}
    ${OPENMP_LIBS}
    ${CLI_LIBS}
    ${SQLITE3_LIBS}
    CompilerFlags
)

target_include_directories(${MAIN_EXEC_NAME} PUBLIC "${PROJECT_BINARY_DIR}")

# Install rules
install(TARGETS ${MAIN_EXEC_NAME} DESTINATION bin)