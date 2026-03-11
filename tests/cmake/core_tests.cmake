set(CORE_COMMAND_TEST_SOURCES
    core/command/Application_test.cpp
    core/command/CommandBase_test.cpp
    core/command/CommandCLIHelper_test.cpp
    core/command/CommandOptionHelper_test.cpp
    core/command/CommandScalarValidationHelper_test.cpp
    core/command/CommandValidationIssue_test.cpp
    core/command/GuiCommandExecutor_test.cpp
    core/command/PotentialAnalysisCommand_test.cpp
    core/command/PotentialDisplayCommand_test.cpp
    core/command/ResultDumpCommand_test.cpp
    core/command/MapSimulationCommand_test.cpp
    core/command/MapVisualizationCommand_test.cpp
    core/command/PositionEstimationCommand_test.cpp
    core/command/HRLModelTestCommand_test.cpp
)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    list(APPEND CORE_COMMAND_TEST_SOURCES
        core/command/PotentialAnalysisExperimentalBondGate_test.cpp
    )
endif()

add_rhbm_gtest_target(rhbm_tests_core_command
    DOMAIN core
    INTENTS command validation
    SOURCES ${CORE_COMMAND_TEST_SOURCES}
)

add_rhbm_gtest_target(rhbm_tests_core_contract
    DOMAIN core
    INTENTS contract
    SOURCES
        core/contract/BuiltInCommandBindingName_test.cpp
        core/contract/BuiltInCommandCatalog_test.cpp
        core/contract/CommandScaffoldScript_test.cpp
        core/contract/CommandCommonOptions_test.cpp
        core/contract/CommandDescriptorShape_test.cpp
        core/contract/CommandExecutionContract_test.cpp
        core/contract/DocsSync_test.cpp
        core/contract/EnumOptionTraits_test.cpp
        core/contract/PublicHeaderSurface_test.cpp
)
