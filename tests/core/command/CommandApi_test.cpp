#include <gtest/gtest.h>

#include "support/CommandValidationAssertions.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>

namespace rg = rhbm_gem;

TEST(CommandApiTest, MapSimulationRequiresModelFilePath)
{
    const auto result{ rg::RunMapSimulation(rg::MapSimulationRequest{}) };

    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--model",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}
