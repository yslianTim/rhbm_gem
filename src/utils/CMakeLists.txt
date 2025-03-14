# File dependencies
file(GLOB SOURCE *.cpp)


# Build utils module
add_library(utils SHARED ${SOURCE})
target_include_directories(utils PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_link_libraries(utils PUBLIC Eigen3::Eigen)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP support found")
    target_compile_definitions(utils PUBLIC USE_OPENMP)
    target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP support not found, using serial version")
endif()

find_package(ROOT 6.28 QUIET COMPONENTS Core Hist Gpad RIO)
if(ROOT_FOUND)
    message(STATUS "ROOT library found, enabling ROOT features")
    target_compile_definitions(utils PUBLIC HAVE_ROOT)
    target_link_libraries(utils PUBLIC ROOT::Core ROOT::Hist ROOT::Gpad ROOT::RIO)
else()
    message(STATUS "ROOT library not found, disabling ROOT features")
endif()