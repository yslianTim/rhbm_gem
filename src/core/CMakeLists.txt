# File dependencies
file(GLOB SOURCE *.cpp)


# 建立 core 模組
add_library(core STATIC ${SOURCE})
target_include_directories(core PUBLIC ${CMAKE_SOURCE_DIR}/include/core)
target_link_libraries(core PUBLIC
    CLI11::CLI11
    data
)