#include <gtest/gtest.h>

#include "support/CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/command/HRLModelTestCommand.hpp"

namespace rg = rhbm_gem;

TEST(HRLModelTestCommandTest, FitRangeOrderingBecomesPrepareValidationError)
{
    rg::HRLModelTestCommand command{rg::CommonOptionProfile::FileWorkflow};
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

TEST(HRLModelTestCommandTest, ReapplyingInvalidRequestAfterSuccessfulPrepareRefreshesState)
{
    rg::HRLModelTestCommand command{rg::CommonOptionProfile::FileWorkflow};
    rg::HRLModelTestRequest request{};
    request.fit_range_min = 0.1;
    request.fit_range_max = 1.0;
    command.ApplyRequest(request);

    ASSERT_TRUE(command.PrepareForExecution());

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
