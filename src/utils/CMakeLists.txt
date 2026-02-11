set(UTILS_SOURCES
    AtomKeySystem.cpp
    AtomSelector.cpp
    BondKeySystem.cpp
    ChemicalDataHelper.cpp
    ComponentHelper.cpp
    ComponentKeySystem.cpp
    CylinderSampler.cpp
    ElectricPotential.cpp
    FilePathHelper.cpp
    GausLinearTransformHelper.cpp
    GridSampler.cpp
    HRLModelHelper.cpp
    HRLModelTester.cpp
    LocalPainter.cpp
    Logger.cpp
    ROOTHelper.cpp
    SphereSampler.cpp
)

# Build utils module
add_library(utils ${UTILS_SOURCES})
target_include_directories(utils PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/utils>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/utils>"
)
target_link_libraries(utils PUBLIC Eigen3::Eigen)
target_link_libraries(utils PRIVATE $<BUILD_INTERFACE:CompilerFlags>)

set(MIN_OPENMP_VERSION "4.5")
set(RHBM_GEM_ENABLE_OPENMP FALSE)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    if(DEFINED OpenMP_CXX_VERSION AND OpenMP_CXX_VERSION VERSION_LESS MIN_OPENMP_VERSION)
        message(WARNING "Found OpenMP version ${OpenMP_CXX_VERSION} but need at least ${MIN_OPENMP_VERSION}. Disabling OpenMP features.")
    else()
        message(STATUS "OpenMP found (version ${OpenMP_CXX_VERSION}). Enabling OpenMP features.")
        target_compile_definitions(utils PUBLIC USE_OPENMP)
        target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
        set(RHBM_GEM_ENABLE_OPENMP TRUE)
    endif()
else()
    message(STATUS "OpenMP support not found, using serial version")
endif()

set(RHBM_GEM_ENABLE_ROOT FALSE)
find_package(ROOT 6.28 QUIET COMPONENTS Core Hist Gpad RIO)
if(ROOT_FOUND)
    message(STATUS "ROOT library found, enabling ROOT features")
    target_compile_definitions(utils PUBLIC HAVE_ROOT)
    target_link_libraries(utils PUBLIC ROOT::Core ROOT::Hist ROOT::Gpad ROOT::RIO)
    set(RHBM_GEM_ENABLE_ROOT TRUE)
else()
    message(STATUS "ROOT library not found, disabling ROOT features")
endif()

set_target_properties(utils PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

install(
    TARGETS utils
    EXPORT RHBM_GEMTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

set(RHBM_GEM_WITH_OPENMP "${RHBM_GEM_ENABLE_OPENMP}" CACHE INTERNAL "RHBM_GEM was built with OpenMP support" FORCE)
set(RHBM_GEM_WITH_ROOT "${RHBM_GEM_ENABLE_ROOT}" CACHE INTERNAL "RHBM_GEM was built with ROOT support" FORCE)
