include(${CMAKE_CURRENT_LIST_DIR}/CoreCommandTests.generated.cmake)

set(CORE_COMMAND_TEST_SOURCES
    core/command/Application_test.cpp
    core/command/CommandBase_test.cpp
    core/command/CommandCLIHelper_test.cpp
    core/command/CommandOptionHelper_test.cpp
    core/command/CommandScalarValidationHelper_test.cpp
    core/command/CommandValidationIssue_test.cpp
    core/command/GuiCommandExecutor_test.cpp
    ${RHBM_GEM_GENERATED_CORE_COMMAND_TEST_SOURCES}
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
        core/contract/CommandCatalog_test.cpp
        core/contract/CommandExecutionContract_test.cpp
        core/contract/DocsSync_test.cpp
        core/contract/EnumOptionTraits_test.cpp
        core/contract/PublicHeaderSurface_test.cpp
)
