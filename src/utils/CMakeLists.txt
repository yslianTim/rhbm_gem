set(UTILS_DOMAIN_SOURCES
    domain/AtomKeySystem.cpp
    domain/AtomSelector.cpp
    domain/BondKeySystem.cpp
    domain/ChemicalDataHelper.cpp
    domain/ComponentHelper.cpp
    domain/ComponentKeySystem.cpp
    domain/FilePathHelper.cpp
    domain/LocalPainter.cpp
    domain/Logger.cpp
    domain/ROOTHelper.cpp
)

set(UTILS_MATH_SOURCES
    math/CylinderSampler.cpp
    math/ElectricPotential.cpp
    math/GausLinearTransformHelper.cpp
    math/GridSampler.cpp
    math/SphereSampler.cpp
)

set(UTILS_HRL_SOURCES
    hrl/HRLAlphaTrainer.cpp
    hrl/HRLDataTransform.cpp
    hrl/HRLGroupEstimator.cpp
    hrl/HRLModelAlgorithms.cpp
    hrl/HRLModelTester.cpp
)

set(UTILS_INTERNAL_OBJECT_TARGETS
    utils_domain
    utils_math
    utils_hrl
)

function(rhbm_configure_utils_object_target target_name)
    target_include_directories(${target_name} PUBLIC
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    target_link_libraries(${target_name} PUBLIC rhbm_eigen)
    target_link_libraries(${target_name} PRIVATE $<BUILD_INTERFACE:CompilerFlags>)
endfunction()

function(rhbm_utils_object_compile_definitions)
    foreach(target_name IN LISTS UTILS_INTERNAL_OBJECT_TARGETS)
        target_compile_definitions(${target_name} PRIVATE ${ARGN})
    endforeach()
endfunction()

function(rhbm_utils_object_link_libraries)
    foreach(target_name IN LISTS UTILS_INTERNAL_OBJECT_TARGETS)
        target_link_libraries(${target_name} PRIVATE ${ARGN})
    endforeach()
endfunction()

function(rhbm_utils_object_include_directories)
    foreach(target_name IN LISTS UTILS_INTERNAL_OBJECT_TARGETS)
        target_include_directories(${target_name} PRIVATE ${ARGN})
    endforeach()
endfunction()

add_library(utils_domain OBJECT ${UTILS_DOMAIN_SOURCES})
rhbm_configure_utils_object_target(utils_domain)

add_library(utils_math OBJECT ${UTILS_MATH_SOURCES})
rhbm_configure_utils_object_target(utils_math)

add_library(utils_hrl OBJECT ${UTILS_HRL_SOURCES})
rhbm_configure_utils_object_target(utils_hrl)

add_library(utils
    $<TARGET_OBJECTS:utils_domain>
    $<TARGET_OBJECTS:utils_math>
    $<TARGET_OBJECTS:utils_hrl>
)
target_include_directories(utils PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(utils PUBLIC rhbm_eigen)
target_link_libraries(utils PRIVATE
    utils_domain
    utils_math
    utils_hrl
    $<BUILD_INTERFACE:CompilerFlags>
)

set(MIN_OPENMP_VERSION "4.5")
set(RHBM_GEM_ENABLE_OPENMP FALSE)
if(RHBM_GEM_OPENMP_MODE STREQUAL "OFF")
    message(STATUS "OpenMP features disabled by RHBM_GEM_OPENMP_MODE=OFF")
elseif(RHBM_GEM_OPENMP_MODE STREQUAL "AUTO")
    find_package(OpenMP QUIET COMPONENTS CXX)
elseif(RHBM_GEM_OPENMP_MODE STREQUAL "ON")
    find_package(OpenMP REQUIRED COMPONENTS CXX)
else()
    message(FATAL_ERROR "Unsupported RHBM_GEM_OPENMP_MODE='${RHBM_GEM_OPENMP_MODE}'")
endif()

if(OpenMP_CXX_FOUND)
    if(DEFINED OpenMP_CXX_VERSION AND OpenMP_CXX_VERSION VERSION_LESS MIN_OPENMP_VERSION)
        if(RHBM_GEM_OPENMP_MODE STREQUAL "ON")
            message(FATAL_ERROR
                "OpenMP version ${OpenMP_CXX_VERSION} is below required ${MIN_OPENMP_VERSION} "
                "while RHBM_GEM_OPENMP_MODE=ON.")
        endif()
        message(WARNING
            "Found OpenMP version ${OpenMP_CXX_VERSION} but need at least ${MIN_OPENMP_VERSION}. "
            "Disabling OpenMP features.")
    else()
        message(STATUS "OpenMP found (version ${OpenMP_CXX_VERSION}). Enabling OpenMP features.")
        target_compile_definitions(utils PUBLIC USE_OPENMP)
        target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
        rhbm_utils_object_compile_definitions(USE_OPENMP)
        rhbm_utils_object_link_libraries(OpenMP::OpenMP_CXX)
        set(RHBM_GEM_ENABLE_OPENMP TRUE)
    endif()
else()
    if(RHBM_GEM_OPENMP_MODE STREQUAL "AUTO")
        message(STATUS "OpenMP support not found, using serial version")
    endif()
endif()

set(RHBM_GEM_ENABLE_BOOST FALSE)
set(_RHBM_GEM_BOOST_COMPONENTS)
set(_RHBM_GEM_BOOST_FIND_ARGS)
if(RHBM_GEM_BOOST_COMPONENTS)
    set(_RHBM_GEM_BOOST_COMPONENTS ${RHBM_GEM_BOOST_COMPONENTS})
    list(APPEND _RHBM_GEM_BOOST_FIND_ARGS COMPONENTS ${_RHBM_GEM_BOOST_COMPONENTS})
endif()

if(RHBM_GEM_BOOST_MODE STREQUAL "OFF")
    message(STATUS "Boost features disabled by RHBM_GEM_BOOST_MODE=OFF")
elseif(RHBM_GEM_BOOST_MODE STREQUAL "AUTO")
    find_package(Boost CONFIG QUIET ${_RHBM_GEM_BOOST_FIND_ARGS})
    if(NOT Boost_FOUND)
        if(POLICY CMP0167)
            cmake_policy(PUSH)
            cmake_policy(SET CMP0167 OLD)
            find_package(Boost QUIET ${_RHBM_GEM_BOOST_FIND_ARGS})
            cmake_policy(POP)
        else()
            find_package(Boost QUIET ${_RHBM_GEM_BOOST_FIND_ARGS})
        endif()
    endif()
elseif(RHBM_GEM_BOOST_MODE STREQUAL "ON")
    find_package(Boost CONFIG QUIET ${_RHBM_GEM_BOOST_FIND_ARGS})
    if(NOT Boost_FOUND)
        if(POLICY CMP0167)
            cmake_policy(PUSH)
            cmake_policy(SET CMP0167 OLD)
            find_package(Boost REQUIRED ${_RHBM_GEM_BOOST_FIND_ARGS})
            cmake_policy(POP)
        else()
            find_package(Boost REQUIRED ${_RHBM_GEM_BOOST_FIND_ARGS})
        endif()
    endif()
else()
    message(FATAL_ERROR "Unsupported RHBM_GEM_BOOST_MODE='${RHBM_GEM_BOOST_MODE}'")
endif()

if(Boost_FOUND)
    if(_RHBM_GEM_BOOST_COMPONENTS)
        message(STATUS "Boost found, enabling Boost features (components: ${_RHBM_GEM_BOOST_COMPONENTS})")
    else()
        message(STATUS "Boost found, enabling Boost features")
    endif()

    target_compile_definitions(utils PUBLIC HAVE_BOOST)
    rhbm_utils_object_compile_definitions(HAVE_BOOST)

    if(TARGET Boost::headers)
        target_link_libraries(utils PUBLIC Boost::headers)
        rhbm_utils_object_link_libraries(Boost::headers)
    elseif(TARGET Boost::boost)
        target_link_libraries(utils PUBLIC Boost::boost)
        rhbm_utils_object_link_libraries(Boost::boost)
    elseif(Boost_INCLUDE_DIRS)
        target_include_directories(utils PUBLIC ${Boost_INCLUDE_DIRS})
        rhbm_utils_object_include_directories(${Boost_INCLUDE_DIRS})
    endif()

    set(_RHBM_GEM_BOOST_COMPONENT_TARGETS)
    foreach(_boost_component IN LISTS _RHBM_GEM_BOOST_COMPONENTS)
        if(TARGET "Boost::${_boost_component}")
            list(APPEND _RHBM_GEM_BOOST_COMPONENT_TARGETS "Boost::${_boost_component}")
        endif()
    endforeach()
    if(_RHBM_GEM_BOOST_COMPONENT_TARGETS)
        target_link_libraries(utils PUBLIC ${_RHBM_GEM_BOOST_COMPONENT_TARGETS})
        rhbm_utils_object_link_libraries(${_RHBM_GEM_BOOST_COMPONENT_TARGETS})
    elseif(Boost_LIBRARIES)
        target_link_libraries(utils PUBLIC ${Boost_LIBRARIES})
        rhbm_utils_object_link_libraries(${Boost_LIBRARIES})
    endif()

    set(RHBM_GEM_ENABLE_BOOST TRUE)
else()
    if(RHBM_GEM_BOOST_MODE STREQUAL "AUTO")
        message(STATUS "Boost library not found, disabling Boost features")
    endif()
endif()

set(RHBM_GEM_ENABLE_ROOT FALSE)
if(RHBM_GEM_ROOT_MODE STREQUAL "OFF")
    message(STATUS "ROOT features disabled by RHBM_GEM_ROOT_MODE=OFF")
elseif(RHBM_GEM_ROOT_MODE STREQUAL "AUTO")
    find_package(ROOT 6.28 QUIET COMPONENTS Core Hist Gpad RIO)
elseif(RHBM_GEM_ROOT_MODE STREQUAL "ON")
    find_package(ROOT 6.28 REQUIRED COMPONENTS Core Hist Gpad RIO)
else()
    message(FATAL_ERROR "Unsupported RHBM_GEM_ROOT_MODE='${RHBM_GEM_ROOT_MODE}'")
endif()

if(ROOT_FOUND)
    message(STATUS "ROOT library found, enabling ROOT features")
    target_compile_definitions(utils PUBLIC HAVE_ROOT)
    target_link_libraries(utils PUBLIC ROOT::Core ROOT::Hist ROOT::Gpad ROOT::RIO)
    rhbm_utils_object_compile_definitions(HAVE_ROOT)
    rhbm_utils_object_link_libraries(ROOT::Core ROOT::Hist ROOT::Gpad ROOT::RIO)
    set(RHBM_GEM_ENABLE_ROOT TRUE)
else()
    if(RHBM_GEM_ROOT_MODE STREQUAL "AUTO")
        message(STATUS "ROOT library not found, disabling ROOT features")
    endif()
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
set(RHBM_GEM_WITH_BOOST "${RHBM_GEM_ENABLE_BOOST}" CACHE INTERNAL "RHBM_GEM was built with Boost support" FORCE)
set(RHBM_GEM_WITH_ROOT "${RHBM_GEM_ENABLE_ROOT}" CACHE INTERNAL "RHBM_GEM was built with ROOT support" FORCE)

unset(_RHBM_GEM_BOOST_COMPONENT_TARGETS)
unset(_RHBM_GEM_BOOST_FIND_ARGS)
unset(_RHBM_GEM_BOOST_COMPONENTS)
