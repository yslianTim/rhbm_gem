#include <gtest/gtest.h>

#include <filesystem>

#include "CommandTestHelpers.hpp"
#include "CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/PotentialAnalysisCommand.hpp"

namespace rg = rhbm_gem;

TEST(PotentialAnalysisCommandTest, NegativeAlphaValuesBecomeValidationErrors)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.alpha_r = -0.1;
    request.alpha_g = 0.0;
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(command, "--alpha-r", rg::ValidationPhase::Parse),
        nullptr);
    ASSERT_NE(
        command_test::FindValidationIssue(command, "--alpha-g", rg::ValidationPhase::Parse),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, SimulationRequiresPositiveResolutionAtPrepare)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.simulation_flag = true;
    request.simulated_map_resolution = 0.0;
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--sim-resolution",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, NonPositiveSamplingHeightBecomesValidationError)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.sampling_height = 0.0;
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--sampling-height",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, EmptySavedKeyBecomesValidationError)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.saved_key_tag = "";
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--save-key",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}
