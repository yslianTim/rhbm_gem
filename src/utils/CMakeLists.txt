# File dependencies
file(GLOB SOURCE *.cpp)


# Build utils module
add_library(utils SHARED ${SOURCE})
target_include_directories(utils PUBLIC ${CMAKE_SOURCE_DIR}/include/utils)
target_link_libraries(utils PUBLIC Eigen3::Eigen)

set(MIN_OPENMP_VERSION "4.5")
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    if(DEFINED OpenMP_CXX_VERSION AND OpenMP_CXX_VERSION VERSION_LESS MIN_OPENMP_VERSION)
        message(WARNING "Found OpenMP version ${OpenMP_CXX_VERSION} but need at least ${MIN_OPENMP_VERSION}. Disabling OpenMP features.")
    else()
        message(STATUS "OpenMP found (version ${OpenMP_CXX_VERSION}). Enabling OpenMP features.")
        target_compile_definitions(utils PUBLIC USE_OPENMP)
        target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
    endif()
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

set_target_properties(utils PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)