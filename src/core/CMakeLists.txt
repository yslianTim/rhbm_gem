set(CORE_SOURCES
    command/Application.cpp
    painter/AtomPainter.cpp
    command/BuiltInCommandCatalog.cpp
    command/CommandBase.cpp
    painter/ComparisonPainter.cpp
    workflow/DataObjectWorkflowOps.cpp
    command/DataObjectManager.cpp
    painter/DemoPainter.cpp
    painter/GausPainter.cpp
    command/HRLModelTestCommand.cpp
    workflow/HRLModelTestCommandWorkflows.cpp
    command/MapSampling.cpp
    command/MapSimulationCommand.cpp
    command/MapVisualizationCommand.cpp
    painter/ModelPainter.cpp
    command/PositionEstimationCommand.cpp
    command/PotentialAnalysisCommand.cpp
    workflow/PotentialAnalysisCommandWorkflow.cpp
    command/PotentialDisplayCommand.cpp
    command/ResultDumpCommand.cpp
    workflow/ResultDumpCommandWorkflow.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    list(APPEND CORE_SOURCES
        workflow/PotentialAnalysisBondWorkflow.cpp
    )
endif()

add_library(core ${CORE_SOURCES})
target_include_directories(core PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(core PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/internal"
    "${PROJECT_SOURCE_DIR}/src/data"
    "${PROJECT_SOURCE_DIR}/src/data/internal"
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
