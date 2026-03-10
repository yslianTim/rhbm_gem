function(rhbm_utils_prepare_openmp_for_appleclang)
    if(NOT APPLE OR NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        return()
    endif()

    if(DEFINED OpenMP_ROOT AND EXISTS "${OpenMP_ROOT}/include/omp.h")
        set(_rhbm_openmp_root "${OpenMP_ROOT}")
    else()
        set(_rhbm_openmp_candidates)
        if(DEFINED ENV{HOMEBREW_PREFIX})
            list(APPEND _rhbm_openmp_candidates "$ENV{HOMEBREW_PREFIX}/opt/libomp")
        endif()
        list(APPEND _rhbm_openmp_candidates
            "/opt/homebrew/opt/libomp"
            "/usr/local/opt/libomp"
        )

        foreach(_candidate_root IN LISTS _rhbm_openmp_candidates)
            if(EXISTS "${_candidate_root}/include/omp.h"
               AND (EXISTS "${_candidate_root}/lib/libomp.dylib"
                    OR EXISTS "${_candidate_root}/lib/libomp.a"))
                set(_rhbm_openmp_root "${_candidate_root}")
                break()
            endif()
        endforeach()
    endif()

    if(_rhbm_openmp_root)
        message(STATUS "AppleClang detected, probing OpenMP in: ${_rhbm_openmp_root}")
        set(OpenMP_ROOT "${_rhbm_openmp_root}" CACHE PATH
            "Root directory for OpenMP installation"
            FORCE
        )
        list(PREPEND CMAKE_PREFIX_PATH "${_rhbm_openmp_root}")
        list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    endif()
endfunction()

set(RHBM_GEM_LIBRARY_SOURCES
    core/command/Application.cpp
    core/painter/AtomPainter.cpp
    core/command/BuiltInCommandCatalog.cpp
    core/command/CommandBase.cpp
    core/painter/ComparisonPainter.cpp
    core/workflow/DataObjectWorkflowOps.cpp
    core/painter/DemoPainterRendering.cpp
    core/painter/DemoPainterWorkflow.cpp
    core/painter/GausPainterRendering.cpp
    core/painter/GausPainterFitting.cpp
    core/painter/PotentialPlotBuilder.cpp
    core/command/HRLModelTestCommand.cpp
    core/workflow/HRLModelTestCommandWorkflows.cpp
    core/command/MapSampling.cpp
    core/command/MapSimulationCommand.cpp
    core/command/MapVisualizationCommand.cpp
    core/painter/ModelPainterRender.cpp
    core/painter/ModelPainterAnnotation.cpp
    core/painter/ModelPainterStyle.cpp
    core/command/PositionEstimationCommand.cpp
    core/command/PotentialAnalysisCommand.cpp
    core/workflow/PotentialAnalysisAnalysisWorkflow.cpp
    core/workflow/PotentialAnalysisTrainingWorkflow.cpp
    core/workflow/PotentialAnalysisReportWorkflow.cpp
    core/command/PotentialDisplayCommand.cpp
    core/command/ResultDumpCommand.cpp
    core/workflow/PotentialAnalysisTrainingSupport.cpp
    core/workflow/ResultDumpCommandWorkflow.cpp

    data/dispatch/DataObjectDispatch.cpp
    data/object/AtomClassifier.cpp
    data/object/AtomObject.cpp
    data/object/AtomicModelDataBlock.cpp
    data/object/BondClassifier.cpp
    data/object/BondObject.cpp
    data/object/ChemicalComponentEntry.cpp
    data/object/GroupPotentialEntry.cpp
    data/object/LocalPotentialEntry.cpp
    data/object/MapObject.cpp
    data/object/ModelObject.cpp
    data/object/PotentialEntryQuery.cpp
    data/io/file/CCP4Format.cpp
    data/io/file/CifFormat.cpp
    data/io/file/FileFormatBackendFactory.cpp
    data/io/file/FileFormatRegistry.cpp
    data/io/file/FileProcessFactoryResolver.cpp
    data/io/file/MapAxisOrderHelper.cpp
    data/io/file/MapFileReader.cpp
    data/io/file/MapFileWriter.cpp
    data/io/file/MapObjectFactory.cpp
    data/io/file/MmCifLoopParser.cpp
    data/io/file/MmCifTokenizer.cpp
    data/io/file/ModelFileReader.cpp
    data/io/file/ModelFileWriter.cpp
    data/io/file/ModelObjectFactory.cpp
    data/io/file/MrcFormat.cpp
    data/io/file/PdbFormat.cpp
    data/io/sqlite/BuiltInDataObjectDaos.cpp
    data/io/sqlite/DataObjectDAOFactoryRegistry.cpp
    data/io/sqlite/DatabaseManager.cpp
    data/io/sqlite/ManagedStoreRegistry.cpp
    data/io/sqlite/MapObjectDAO.cpp
    data/io/sqlite/ModelAnalysisPersistence.cpp
    data/io/sqlite/ModelObjectDAOSqlite.cpp
    data/io/sqlite/ModelStructurePersistence.cpp
    data/schema/DatabaseSchemaManager.cpp
    data/io/DataIoServices.cpp
    data/io/DataObjectManager.cpp

    utils/domain/AtomKeySystem.cpp
    utils/domain/AtomSelector.cpp
    utils/domain/BondKeySystem.cpp
    utils/domain/ChemicalDataHelper.cpp
    utils/domain/ComponentHelper.cpp
    utils/domain/ComponentKeySystem.cpp
    utils/domain/FilePathHelper.cpp
    utils/domain/LocalPainter.cpp
    utils/domain/Logger.cpp
    utils/domain/ROOTHelper.cpp
    utils/math/CylinderSampler.cpp
    utils/math/ElectricPotential.cpp
    utils/math/GausLinearTransformHelper.cpp
    utils/math/GridSampler.cpp
    utils/math/SphereSampler.cpp
    utils/hrl/HRLAlphaTrainer.cpp
    utils/hrl/HRLDataTransform.cpp
    utils/hrl/HRLGroupEstimator.cpp
    utils/hrl/HRLModelAlgorithms.cpp
    utils/hrl/HRLModelTester.cpp

    gui/GuiCommandExecutor.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    list(APPEND RHBM_GEM_LIBRARY_SOURCES
        core/workflow/PotentialAnalysisBondWorkflow.cpp
    )
endif()

if(RHBM_GEM_LEGACY_V1_SUPPORT)
    list(APPEND RHBM_GEM_LIBRARY_SOURCES
        data/migration/legacy_v1/LegacyModelObjectReader.cpp
    )
endif()

if(RHBM_GEM_DEP_PROVIDER STREQUAL "FETCH")
    list(APPEND RHBM_GEM_LIBRARY_SOURCES
        "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.c"
    )
endif()

add_library(rhbm_gem ${RHBM_GEM_LIBRARY_SOURCES})

target_include_directories(rhbm_gem PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(rhbm_gem PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/internal"
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/internal"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils"
)

if(RHBM_GEM_DEP_PROVIDER STREQUAL "SYSTEM")
    target_link_libraries(rhbm_gem PUBLIC
        Eigen3::Eigen
        CLI11::CLI11
        SQLite::SQLite3
    )
else()
    target_include_directories(rhbm_gem SYSTEM PUBLIC
        "$<BUILD_INTERFACE:${RHBM_GEM_EIGEN_INCLUDE_DIR}>"
        "$<BUILD_INTERFACE:${RHBM_GEM_CLI11_INCLUDE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    )
    target_include_directories(rhbm_gem PRIVATE
        "${RHBM_GEM_SQLITE3_SOURCE_DIR}"
    )
endif()

target_link_libraries(rhbm_gem PRIVATE
    $<BUILD_INTERFACE:CompilerFlags>
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
endif()
if(RHBM_GEM_LEGACY_V1_SUPPORT)
    target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_LEGACY_V1_SUPPORT)
endif()

set(MIN_OPENMP_VERSION "4.5")
set(RHBM_GEM_ENABLE_OPENMP FALSE)
if(RHBM_GEM_OPENMP_MODE STREQUAL "OFF")
    message(STATUS "OpenMP features disabled by RHBM_GEM_OPENMP_MODE=OFF")
elseif(RHBM_GEM_OPENMP_MODE STREQUAL "AUTO")
    rhbm_utils_prepare_openmp_for_appleclang()
    find_package(OpenMP QUIET COMPONENTS CXX)
elseif(RHBM_GEM_OPENMP_MODE STREQUAL "ON")
    rhbm_utils_prepare_openmp_for_appleclang()
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
        target_compile_definitions(rhbm_gem PUBLIC USE_OPENMP)
        target_link_libraries(rhbm_gem PUBLIC OpenMP::OpenMP_CXX)
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

    target_compile_definitions(rhbm_gem PUBLIC HAVE_BOOST)

    if(TARGET Boost::headers)
        target_link_libraries(rhbm_gem PUBLIC Boost::headers)
    elseif(TARGET Boost::boost)
        target_link_libraries(rhbm_gem PUBLIC Boost::boost)
    elseif(Boost_INCLUDE_DIRS)
        target_include_directories(rhbm_gem PUBLIC ${Boost_INCLUDE_DIRS})
    endif()

    set(_RHBM_GEM_BOOST_COMPONENT_TARGETS)
    foreach(_boost_component IN LISTS _RHBM_GEM_BOOST_COMPONENTS)
        if(TARGET "Boost::${_boost_component}")
            list(APPEND _RHBM_GEM_BOOST_COMPONENT_TARGETS "Boost::${_boost_component}")
        endif()
    endforeach()
    if(_RHBM_GEM_BOOST_COMPONENT_TARGETS)
        target_link_libraries(rhbm_gem PUBLIC ${_RHBM_GEM_BOOST_COMPONENT_TARGETS})
    elseif(Boost_LIBRARIES)
        target_link_libraries(rhbm_gem PUBLIC ${Boost_LIBRARIES})
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
    target_compile_definitions(rhbm_gem PUBLIC HAVE_ROOT)
    target_link_libraries(rhbm_gem PUBLIC ROOT::Core ROOT::Hist ROOT::Gpad ROOT::RIO)
    set(RHBM_GEM_ENABLE_ROOT TRUE)
else()
    if(RHBM_GEM_ROOT_MODE STREQUAL "AUTO")
        message(STATUS "ROOT library not found, disabling ROOT features")
    endif()
endif()

install(
    TARGETS rhbm_gem
    EXPORT RHBM_GEMTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

set(RHBM_GEM_WITH_OPENMP "${RHBM_GEM_ENABLE_OPENMP}" CACHE INTERNAL "RHBM_GEM was built with OpenMP support" FORCE)
set(RHBM_GEM_WITH_BOOST "${RHBM_GEM_ENABLE_BOOST}" CACHE INTERNAL "RHBM_GEM was built with Boost support" FORCE)
set(RHBM_GEM_WITH_ROOT "${RHBM_GEM_ENABLE_ROOT}" CACHE INTERNAL "RHBM_GEM was built with ROOT support" FORCE)
if(DEFINED OpenMP_ROOT)
    set(RHBM_GEM_OPENMP_ROOT "${OpenMP_ROOT}" CACHE INTERNAL "OpenMP root used by RHBM_GEM build" FORCE)
endif()

set(MAIN_EXEC_TARGET ${PROJECT_NAME})
add_executable(${MAIN_EXEC_TARGET} main.cpp)
set_target_properties(${MAIN_EXEC_TARGET} PROPERTIES
    OUTPUT_NAME "${PROJECT_DISPLAY_NAME}"
)
target_link_libraries(${MAIN_EXEC_TARGET} PRIVATE
    rhbm_gem
    CompilerFlags
)
target_include_directories(${MAIN_EXEC_TARGET} PUBLIC "${PROJECT_BINARY_DIR}")

install(
    TARGETS ${MAIN_EXEC_TARGET}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

if(RHBM_GEM_BUILD_GUI)
    find_package(Qt6 QUIET COMPONENTS Core Widgets)
    if(Qt6_FOUND)
        message(STATUS "Qt6 found, enabling GUI executable build")

        set(GUI_EXEC_TARGET ${PROJECT_NAME}_gui)
        add_executable(${GUI_EXEC_TARGET}
            gui/main.cpp
            gui/MainWindow.cpp
        )
        set_target_properties(${GUI_EXEC_TARGET} PROPERTIES
            OUTPUT_NAME "${PROJECT_DISPLAY_NAME}-GUI"
        )
        target_link_libraries(${GUI_EXEC_TARGET} PRIVATE
            rhbm_gem
            Qt6::Core
            Qt6::Widgets
            CompilerFlags
        )
        target_include_directories(${GUI_EXEC_TARGET} PRIVATE
            "${PROJECT_SOURCE_DIR}/src/gui"
        )

        install(
            TARGETS ${GUI_EXEC_TARGET}
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        )
    else()
        message(STATUS "Qt6 not found, skipping GUI executable build")
    endif()
endif()

unset(_RHBM_GEM_BOOST_COMPONENT_TARGETS)
unset(_RHBM_GEM_BOOST_FIND_ARGS)
unset(_RHBM_GEM_BOOST_COMPONENTS)
