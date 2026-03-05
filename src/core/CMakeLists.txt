set(CORE_SOURCES
    Application.cpp
    AtomPainter.cpp
    BuiltInCommandCatalog.cpp
    CommandBase.cpp
    ComparisonPainter.cpp
    DataObjectWorkflowVisitors.cpp
    DataObjectManager.cpp
    DemoPainter.cpp
    GausPainter.cpp
    HRLModelTestCommand.cpp
    MapInterpolationVisitor.cpp
    MapSimulationCommand.cpp
    MapVisualizationCommand.cpp
    ModelPainter.cpp
    PositionEstimationCommand.cpp
    PotentialAnalysisCommand.cpp
    PotentialDisplayCommand.cpp
    ResultDumpCommand.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    list(APPEND CORE_SOURCES
        PotentialAnalysisBondWorkflow.cpp
    )
endif()

add_library(core ${CORE_SOURCES})
target_include_directories(core PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/core>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/core>"
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
