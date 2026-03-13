#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/HRLModelTestCommand.hpp"

namespace rg = rhbm_gem;

TEST(HRLModelTestCommandTest, InvalidTesterChoiceBecomesValidationError)
{
    rg::HRLModelTestCommand command{};
    rg::HRLModelTestRequest request{};
    request.tester_choice = static_cast<rg::TesterType>(999);
    command.ApplyRequest(request);

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
    rg::HRLModelTestCommand command{};
    rg::HRLModelTestRequest request{};
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

TEST(HRLModelTestCommandTest, FitRangeOrderingBecomesPrepareValidationError)
{
    rg::HRLModelTestCommand command{};
    rg::HRLModelTestRequest request{};
    request.fit_range_min = 2.0;
    request.fit_range_max = 1.0;
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--fit-range",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}
