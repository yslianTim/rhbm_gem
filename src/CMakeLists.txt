add_library(rhbm_gem)
target_compile_features(rhbm_gem PUBLIC cxx_std_17)

set(RHBM_GEM_COMMAND_SOURCES
    core/command/HRLModelTestCommand.cpp
    core/command/MapSimulationCommand.cpp
    core/command/MapVisualizationCommand.cpp
    core/command/PositionEstimationCommand.cpp
    core/command/PotentialAnalysisCommand.cpp
    core/command/PotentialDisplayCommand.cpp
    core/command/ResultDumpCommand.cpp
)

target_sources(rhbm_gem PRIVATE
    # core
    core/command/CommandApi.cpp
    core/command/CommandDataSupport.cpp
    core/painter/AtomPainter.cpp
    core/command/CommandCatalog.cpp
    core/command/CommandBase.cpp
    core/painter/ComparisonPainter.cpp
    core/painter/DemoPainterAdvancedRender.cpp
    core/painter/DemoPainterRendering.cpp
    core/painter/DemoPainterWorkflow.cpp
    core/painter/GausPainterComponentRender.cpp
    core/painter/GausPainterRendering.cpp
    core/painter/GausPainterFitting.cpp
    core/painter/PotentialPlotBuilder.cpp
    core/painter/PotentialPlotBuilderSequenceGraphs.cpp
    ${RHBM_GEM_COMMAND_SOURCES}
    core/command/MapSampling.cpp
    core/painter/ModelPainterMapValueRender.cpp
    core/painter/ModelPainterRender.cpp
    core/painter/ModelPainterAnnotation.cpp
    core/painter/ModelPainterStyle.cpp

    # data
    data/object/AtomClassifier.cpp
    data/object/AtomObject.cpp
    data/io/file/ModelImportState.cpp
    data/object/BondClassifier.cpp
    data/object/BondObject.cpp
    data/object/ChemicalComponentEntry.cpp
    data/object/DataObjectDispatch.cpp
    data/object/GroupPotentialEntry.cpp
    data/object/LocalPotentialEntry.cpp
    data/object/MapObject.cpp
    data/object/ModelObject.cpp
    data/object/PotentialEntryQuery.cpp
    data/io/file/CCP4Format.cpp
    data/io/file/CifFormat.cpp
    data/io/file/FileIO.cpp
    data/io/file/MapAxisOrderHelper.cpp
    data/io/file/MrcFormat.cpp
    data/io/file/PdbFormat.cpp
    data/io/sqlite/ModelObjectStorage.cpp
    data/io/sqlite/SQLitePersistence.cpp
    data/io/DataObjectManager.cpp

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
    utils/math/GausLinearTransformHelper.cpp
    utils/math/GridSampler.cpp
    utils/math/SphereSampler.cpp
    utils/hrl/HRLAlphaTrainer.cpp
    utils/hrl/HRLDataTransform.cpp
    utils/hrl/HRLGroupEstimator.cpp
    utils/hrl/HRLModelAlgorithms.cpp
    utils/hrl/HRLModelTester.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    target_sources(rhbm_gem PRIVATE
        core/workflow/PotentialAnalysisBondWorkflow.cpp
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
    "${CMAKE_CURRENT_SOURCE_DIR}/core/internal"
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/internal"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils"
)

target_link_libraries(rhbm_gem PUBLIC rhbm_gem_dependencies)
target_link_libraries(rhbm_gem PRIVATE CompilerFlags)

if(WIN32 AND BUILD_SHARED_LIBS)
    set_target_properties(rhbm_gem PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
rhbm_gem_apply_library_install_rpath(rhbm_gem)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
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

if(TARGET rhbm_gem_gui_dependencies)
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
        rhbm_gem_gui_dependencies
        CompilerFlags
    )
    target_include_directories(${GUI_EXEC_TARGET} PRIVATE
        "${PROJECT_SOURCE_DIR}/src/gui"
        "${PROJECT_SOURCE_DIR}/src/core"
        "${PROJECT_SOURCE_DIR}/src/core/internal"
    )
    rhbm_gem_apply_executable_install_rpath(${GUI_EXEC_TARGET})

    install(
        TARGETS ${GUI_EXEC_TARGET}
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
endif()
