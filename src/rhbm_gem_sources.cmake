include(${CMAKE_CURRENT_LIST_DIR}/CommandSources.generated.cmake)

set(RHBM_GEM_LIBRARY_SOURCES
    core/command/Application.cpp
    core/command/CommandApi.cpp
    core/command/CommandRuntimeRegistry.cpp
    core/painter/AtomPainter.cpp
    core/painter/RootPlotRenderBackend.cpp
    core/command/CommandCatalog.cpp
    core/command/CommandBase.cpp
    core/painter/ComparisonPainter.cpp
    core/workflow/DataObjectWorkflowOps.cpp
    core/painter/DemoPainterAdvancedRender.cpp
    core/painter/DemoPainterRendering.cpp
    core/painter/DemoPainterWorkflow.cpp
    core/painter/GausPainterComponentRender.cpp
    core/painter/GausPainterRendering.cpp
    core/painter/GausPainterFitting.cpp
    core/painter/PotentialPlotBuilder.cpp
    core/painter/PotentialPlotBuilderSequenceGraphs.cpp
    core/painter/PotentialPlotBuilderPlotFunctions.cpp
    ${RHBM_GEM_COMMAND_SOURCES}
    core/workflow/HRLModelTestCommandWorkflows.cpp
    core/command/MapSampling.cpp
    core/painter/ModelPainterMapValueRender.cpp
    core/painter/ModelPainterRender.cpp
    core/painter/ModelPainterAnnotation.cpp
    core/painter/ModelPainterStyle.cpp
    core/workflow/PotentialAnalysisAnalysisWorkflow.cpp
    core/workflow/PotentialAnalysisTrainingWorkflow.cpp
    core/workflow/PotentialAnalysisReportWorkflow.cpp
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
    data/io/file/CifFormatBondConstruction.cpp
    data/io/file/CifFormatCategoryLoaders.cpp
    data/io/file/CifFormatOrchestration.cpp
    data/io/file/CifFormatWriteBlock.cpp
    data/io/file/ICifCategoryParser.cpp
    data/io/file/FileFormatBackendFactory.cpp
    data/io/file/FileFormatRegistry.cpp
    data/io/file/FileIO.cpp
    data/io/file/MapAxisOrderHelper.cpp
    data/io/file/MmCifLoopParser.cpp
    data/io/file/MmCifTokenizer.cpp
    data/io/file/MrcFormat.cpp
    data/io/file/PdbFormat.cpp
    data/io/sqlite/DatabaseManager.cpp
    data/io/sqlite/MapObjectDAO.cpp
    data/io/sqlite/ModelAnalysisPersistence.cpp
    data/io/sqlite/ModelObjectDAOSqlite.cpp
    data/io/sqlite/ModelStructurePersistence.cpp
    data/schema/DatabaseSchemaManager.cpp
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
