add_library(rhbm_gem)
target_compile_features(rhbm_gem PUBLIC cxx_std_20)
target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_AUDIT_VERSION="${PROJECT_VERSION}")

set(RHBM_GEM_COMMAND_SOURCES
    core/command/detail/LocalFittingFeatures.cpp
    core/command/RHBMTestCommand.cpp
    core/command/detail/RHBMTestPlotting.cpp
    core/command/MapSimulationCommand.cpp
    core/command/detail/MapSimulation.cpp
    core/command/detail/SimulationManifest.cpp
    core/command/PotentialAnalysisCommand.cpp
    core/command/PotentialDisplayCommand.cpp
    core/command/ResultDumpCommand.cpp
)

set(RHBM_GEM_EXPERIMENTAL_SOURCES
    core/command/MapVisualizationCommand.cpp
    core/command/PositionEstimationCommand.cpp
)

target_sources(rhbm_gem PRIVATE
    # core
    core/EstimatorTester.cpp
    core/detail/joint_component/Numerics.cpp
    core/detail/joint_component/TiledDerivative.cpp
    core/detail/joint_component/LinearSolve.cpp
    core/detail/joint_component/CompactSvd.cpp
    core/detail/joint_component/SparseFactor.cpp
    core/detail/joint_component/Search.cpp
    core/detail/joint_component/Components.cpp
    core/JointComponentEstimator.cpp
    core/GaussianEstimator.cpp
    core/detail/FirstStageInitialization.cpp
    core/detail/PostFitPeeling.cpp
    core/detail/gaussian_fit/GaussianModelOperations.cpp
    core/detail/gaussian_fit/PreparedLocalGaussianFit.cpp
    core/detail/second_stage/SecondStageState.cpp
    core/detail/second_stage/observation/SecondStageObservation.cpp
    core/detail/second_stage/CouplingGraph.cpp
    core/detail/second_stage/JointFitting.cpp
    core/detail/second_stage/CandidateEvaluation.cpp
    core/detail/second_stage/ComponentAssembly.cpp
    core/detail/second_stage/CandidateTransaction.cpp
    core/detail/second_stage/CandidateTransactionLocal.cpp
    core/detail/second_stage/ObjectiveEvaluation.cpp
    core/detail/second_stage/SuspiciousUpdate.cpp
    core/detail/second_stage/FixedPointRecovery.cpp
    core/detail/second_stage/IterationProposal.cpp
    core/detail/second_stage/Quarantine.cpp
    core/detail/second_stage/CandidateTransactionBoundary.cpp
    core/detail/second_stage/DependencyPolish.cpp
    core/detail/second_stage/observation/SecondStageLogging.cpp
    core/detail/second_stage/observation/PerformanceCounters.cpp
    core/detail/second_stage/ConvergenceCertificate.cpp
    core/detail/second_stage/IterationProcess.cpp
    core/MapSampler.cpp
    core/QScoreHelper.cpp
    core/TestDataFactory.cpp
    core/command/CommandSystem.cpp
    core/painter/AtomPainter.cpp
    core/painter/ComparisonPainter.cpp
    core/painter/DemoPainter.cpp
    core/painter/GausPainter.cpp
    core/painter/detail/PotentialPlotBuilder.cpp
    ${RHBM_GEM_COMMAND_SOURCES}

    # data
    data/detail/AtomClassifier.cpp
    data/detail/JointStageAdapter.cpp
    data/object/AtomLocalPotentialView.cpp
    data/object/AtomObject.cpp
    data/io/JointAnalysisFileIO.cpp
    data/io/detail/JointResultJson.cpp
    data/io/DataRepository.cpp
    data/io/ModelMapFileIO.cpp
    data/io/file/ModelImportState.cpp
    data/object/BondObject.cpp
    data/object/ChemicalComponentEntry.cpp
    data/object/ModelAnalysisEditor.cpp
    data/object/ModelAnalysisView.cpp
    data/detail/ModelDerivedState.cpp
    data/detail/ModelObjectParts.cpp
    data/object/MapObject.cpp
    data/detail/ModelAnalysisData.cpp
    data/object/ModelObject.cpp
    data/io/file/CCP4Format.cpp
    data/io/file/CifFormat.cpp
    data/io/file/MapHelper.cpp
    data/io/file/MrcFormat.cpp
    data/io/file/PdbFormat.cpp
    data/io/sqlite/ModelObjectStorage.cpp

    # utils
    utils/domain/AtomKeySystem.cpp
    utils/domain/AtomSelector.cpp
    utils/domain/BondKeySystem.cpp
    utils/domain/ChemicalDataHelper.cpp
    utils/domain/ComponentHelper.cpp
    utils/domain/ComponentKeySystem.cpp
    utils/domain/FilePathHelper.cpp
    utils/domain/FileFingerprint.cpp
    utils/domain/LocalPainter.cpp
    utils/domain/Logger.cpp
    utils/domain/ROOTHelper.cpp
    utils/domain/SampleFilter.cpp
    utils/math/ElectricPotential.cpp
    utils/math/GaussianModel3D.cpp
    utils/math/GridSampler.cpp
    utils/math/SphereSampler.cpp
    utils/hrl/LocalPotentialSeries.cpp
    utils/hrl/RHBMTrainer.cpp
    utils/hrl/RHBMHelper.cpp
    utils/hrl/MDPDEEndpointRefinement.cpp
    utils/hrl/LinearizationService.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)
    target_sources(rhbm_gem PRIVATE
        ${RHBM_GEM_EXPERIMENTAL_SOURCES}
    )
endif()

if(RHBM_GEM_ENABLE_UMAP)
    target_sources(rhbm_gem PRIVATE
        core/command/UmapEmbeddingCommand.cpp
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
if(RHBM_GEM_JOINT_SPARSE_BACKEND STREQUAL "SPQR")
    rhbm_gem_prepare_openmp_for_appleclang()
    find_package(SPQR 4 CONFIG REQUIRED)
    if(NOT SPQR_VERSION VERSION_LESS "5")
        message(FATAL_ERROR "Joint SPQR backend requires SPQR 4.x")
    endif()
    target_link_libraries(rhbm_gem PRIVATE SuiteSparse::SPQR SuiteSparse::CHOLMOD)
    target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_JOINT_SPQR)
    set(RHBM_GEM_SPQR_BLAS_LINK "unreported-by-package")
    if(TARGET SuiteSparse::SPQR_static)
        get_target_property(RHBM_GEM_SPQR_BLAS_LINK SuiteSparse::SPQR_static INTERFACE_LINK_LIBRARIES)
        string(GENEX_STRIP "${RHBM_GEM_SPQR_BLAS_LINK}" RHBM_GEM_SPQR_BLAS_LINK)
    endif()
    set(RHBM_GEM_SPQR_BLAS_VERSION "${BLAS_VERSION}")
    if(APPLE)
        get_target_property(_spqr_library SuiteSparse::SPQR IMPORTED_LOCATION_RELEASE)
        if(NOT _spqr_library)
            get_target_property(_spqr_library SuiteSparse::SPQR IMPORTED_LOCATION)
        endif()
        execute_process(COMMAND /usr/bin/otool -L "${_spqr_library}"
            OUTPUT_VARIABLE _spqr_links ERROR_QUIET)
        string(REGEX MATCH "[^\n]*(Accelerate|libblas|libopenblas)[^\n]*" RHBM_GEM_SPQR_BLAS_VERSION "${_spqr_links}")
        string(STRIP "${RHBM_GEM_SPQR_BLAS_VERSION}" RHBM_GEM_SPQR_BLAS_VERSION)
        unset(_spqr_library)
        unset(_spqr_links)
    endif()
    if(NOT RHBM_GEM_SPQR_BLAS_VERSION)
        set(RHBM_GEM_SPQR_BLAS_VERSION "unreported-by-package")
    endif()
endif()
target_link_libraries(rhbm_gem PRIVATE "$<BUILD_INTERFACE:CompilerFlags>")
include(SimulationBuildInfo)
if(RHBM_GEM_ENABLE_UMAP)
    target_link_libraries(rhbm_gem PRIVATE
        "$<BUILD_INTERFACE:libscran::umappp>"
    )
    target_compile_definitions(rhbm_gem PUBLIC RHBM_GEM_ENABLE_UMAP)
endif()

if(WIN32 AND BUILD_SHARED_LIBS)
    set_target_properties(rhbm_gem PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
rhbm_gem_apply_library_install_rpath(rhbm_gem)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)
    target_compile_definitions(rhbm_gem PUBLIC RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE)
endif()
rhbm_gem_apply_observation_definitions(rhbm_gem)

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

if(BUILD_TESTING)
    target_sources(rhbm_gem PRIVATE "${PROJECT_SOURCE_DIR}/tests/support/SecondStageNumericalProbe.cpp"
        "${PROJECT_SOURCE_DIR}/tests/support/SolverFailureCapture.cpp"
        "${PROJECT_SOURCE_DIR}/tests/support/MDPDEExperiment.cpp"
        "${PROJECT_SOURCE_DIR}/tests/support/EndpointRefinementExperiment.cpp")
    target_include_directories(rhbm_gem PRIVATE "${PROJECT_SOURCE_DIR}/tests")
    target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_TEST_INSTRUMENTATION)
endif()
