# File dependencies
file(GLOB SOURCE *.cpp)


# 建立 data 模組
add_library(data SHARED ${SOURCE})
target_include_directories(data PUBLIC ${CMAKE_SOURCE_DIR}/include/data)
target_link_libraries(data PUBLIC
    SQLite::SQLite3
    utils
)