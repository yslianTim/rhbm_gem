add_library(rhbm_gem)
target_compile_features(rhbm_gem PUBLIC cxx_std_17)

set(RHBM_GEM_COMMAND_SOURCES
    core/command/HRLModelTestCommand.cpp
    core/command/MapSimulationCommand.cpp
    core/command/PotentialAnalysisCommand.cpp
    core/command/PotentialDisplayCommand.cpp
    core/command/ResultDumpCommand.cpp
)

set(RHBM_GEM_EXPERIMENTAL_SOURCES
    core/command/MapVisualizationCommand.cpp
    core/command/PositionEstimationCommand.cpp
    core/experimental/PotentialAnalysisBondWorkflow.cpp
)

target_sources(rhbm_gem PRIVATE
    # core
    core/command/CommandApi.cpp
    core/command/detail/CommandCli.cpp
    core/painter/AtomPainter.cpp
    core/command/detail/CommandBase.cpp
    core/painter/ComparisonPainter.cpp
    core/painter/DemoPainter.cpp
    core/painter/GausPainter.cpp
    core/painter/ModelPainter.cpp
    core/painter/PotentialPlotBuilder.cpp
    core/painter/detail/PainterModelAccess.cpp
    ${RHBM_GEM_COMMAND_SOURCES}

    # data
    data/object/AtomClassifier.cpp
    core/painter/AtomStyleCatalog.cpp
    data/object/AtomObject.cpp
    data/io/DataRepository.cpp
    data/io/file/ModelImportState.cpp
    data/object/BondClassifier.cpp
    core/painter/BondStyleCatalog.cpp
    data/object/BondObject.cpp
    data/object/ChemicalComponentEntry.cpp
    data/object/ModelAnalysisEditor.cpp
    data/object/ModelAnalysisView.cpp
    data/detail/GroupPotentialEntry.cpp
    data/detail/LocalPotentialEntry.cpp
    data/detail/ModelDerivedState.cpp
    data/detail/ModelObjectParts.cpp
    data/object/MapObject.cpp
    data/detail/MapSpatialIndex.cpp
    data/detail/ModelAnalysisData.cpp
    data/object/ModelObject.cpp
    data/io/file/CCP4Format.cpp
    data/io/file/CifFormat.cpp
    data/io/file/ModelMapFileIO.cpp
    data/io/file/MapAxisOrderHelper.cpp
    data/io/file/MrcFormat.cpp
    data/io/file/PdbFormat.cpp
    data/io/sqlite/ModelObjectStorage.cpp
    data/io/sqlite/SQLitePersistence.cpp

    # utils
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
    utils/math/GaussianPotentialSampler.cpp
    utils/math/GaussianResponseMath.cpp
    utils/math/GridSampler.cpp
    utils/math/SphereSampler.cpp
    utils/hrl/RHBMTrainer.cpp
    utils/hrl/RHBMHelper.cpp
    utils/hrl/GaussianLinearizationService.cpp
    utils/hrl/HRLModelTestDataFactory.cpp
    utils/hrl/HRLModelTester.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)
    target_sources(rhbm_gem PRIVATE
        ${RHBM_GEM_EXPERIMENTAL_SOURCES}
    )
endif()

if(RHBM_GEM_DEP_PROVIDER STREQUAL "FETCH")
    target_sources(rhbm_gem PRIVATE
        "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.c"
    )
endif()

target_include_directories(rhbm_gem PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(rhbm_gem PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core"
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils"
)

target_link_libraries(rhbm_gem PUBLIC rhbm_gem_dependencies)
target_link_libraries(rhbm_gem PRIVATE CompilerFlags)

if(WIN32 AND BUILD_SHARED_LIBS)
    set_target_properties(rhbm_gem PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
rhbm_gem_apply_library_install_rpath(rhbm_gem)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)
    target_compile_definitions(rhbm_gem PUBLIC RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)
endif()
install(
    TARGETS rhbm_gem rhbm_gem_dependencies
    EXPORT RHBM_GEMTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

# Avoid target-directory collisions on case-insensitive filesystems
# (e.g. `rhbm_gem` library vs `RHBM_GEM` executable).
set(MAIN_EXEC_TARGET rhbm_gem_cli)
add_executable(${MAIN_EXEC_TARGET} main.cpp)
set_target_properties(${MAIN_EXEC_TARGET} PROPERTIES
    OUTPUT_NAME "${PROJECT_DISPLAY_NAME}"
)
target_link_libraries(${MAIN_EXEC_TARGET} PRIVATE
    rhbm_gem
    CompilerFlags
)
target_include_directories(${MAIN_EXEC_TARGET} PRIVATE "${PROJECT_BINARY_DIR}")
rhbm_gem_apply_executable_install_rpath(${MAIN_EXEC_TARGET})

install(
    TARGETS ${MAIN_EXEC_TARGET}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()
