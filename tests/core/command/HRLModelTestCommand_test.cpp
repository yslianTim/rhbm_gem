#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>

namespace rg = rhbm_gem;

TEST(HRLModelTestCommandTest, InvalidTesterChoiceBecomesValidationError)
{
    rg::HRLModelTestCommand command;
    command.SetTesterChoice(static_cast<rg::TesterType>(999));

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--tester",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(HRLModelTestCommandTest, NegativeAlphaValuesBecomeValidationErrors)
{
    rg::HRLModelTestCommand command;
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

TEST(HRLModelTestCommandTest, FitRangeOrderingBecomesPrepareValidationError)
{
    rg::HRLModelTestCommand command;
    command.SetFitRangeMinimum(2.0);
    command.SetFitRangeMaximum(1.0);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--fit-range",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}
