# File dependencies
file(GLOB SOURCE *.cpp)


# 建立 utils 模組
add_library(utils SHARED ${SOURCE})
target_include_directories(utils PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_link_libraries(utils PUBLIC
    Eigen3::Eigen
    ${ROOT_LIBS}
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP support found")
    target_compile_definitions(utils PUBLIC USE_OPENMP)
    target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP support not found, using serial version")
endif()