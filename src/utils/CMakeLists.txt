# File dependencies
file(GLOB SOURCE *.cpp)


# 建立 utils 模組
add_library(utils SHARED ${SOURCE})
target_include_directories(utils PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_link_libraries(utils PUBLIC ${EIGEN_LIBS})