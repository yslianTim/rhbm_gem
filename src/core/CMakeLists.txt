set(CORE_SOURCES
    Application.cpp
    painter/AtomPainter.cpp
    BuiltInCommandCatalog.cpp
    CommandBase.cpp
    painter/ComparisonPainter.cpp
    workflow/DataObjectWorkflowOps.cpp
    DataObjectManager.cpp
    painter/DemoPainter.cpp
    painter/GausPainter.cpp
    HRLModelTestCommand.cpp
    workflow/HRLModelTestCommandWorkflows.cpp
    MapSampling.cpp
    MapSimulationCommand.cpp
    MapVisualizationCommand.cpp
    painter/ModelPainter.cpp
    PositionEstimationCommand.cpp
    PotentialAnalysisCommand.cpp
    workflow/PotentialAnalysisCommandWorkflow.cpp
    PotentialDisplayCommand.cpp
    ResultDumpCommand.cpp
    workflow/ResultDumpCommandWorkflow.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    list(APPEND CORE_SOURCES
        workflow/PotentialAnalysisBondWorkflow.cpp
    )
endif()

add_library(core ${CORE_SOURCES})
target_include_directories(core PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/core>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/core>"
)
target_include_directories(core PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(core PUBLIC
    utils
    data
    rhbm_cli11
)
target_link_libraries(core PRIVATE
    $<BUILD_INTERFACE:CompilerFlags>
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    target_compile_definitions(core PRIVATE RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
endif()

set_target_properties(core PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

install(
    TARGETS core
    EXPORT RHBM_GEMTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
