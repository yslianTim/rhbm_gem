#include <gtest/gtest.h>

#include "support/CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/command/HRLModelTestCommand.hpp"

namespace rg = rhbm_gem;

TEST(HRLModelTestCommandTest, FitRangeOrderingBecomesPrepareValidationError)
{
    rg::HRLModelTestCommand command{};
    rg::HRLModelTestRequest request{};
    request.fit_range_min = 2.0;
    request.fit_range_max = 1.0;
    command.ApplyRequest(request);

    EXPECT_FALSE(command.Run());
    EXPECT_FALSE(command.WasPrepared());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--fit-range",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}
