#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include <rhbm_gem/core/PotentialAnalysisCommand.hpp>

namespace rg = rhbm_gem;

TEST(PotentialAnalysisCommandTest, NegativeAlphaValuesBecomeValidationErrors)
{
    rg::PotentialAnalysisCommand command;
    command.SetAlphaR(-0.1);
    command.SetAlphaG(0.0);

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
    rg::PotentialAnalysisCommand command;
    command.SetSimulationFlag(true);
    command.SetSimulatedMapResolution(0.0);

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
    rg::PotentialAnalysisCommand command;
    command.SetSamplingHeight(0.0);

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
    rg::PotentialAnalysisCommand command;
    command.SetSavedKeyTag("");

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--save-key",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}
