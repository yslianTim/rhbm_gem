#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/PotentialAnalysisCommand.hpp"

namespace rg = rhbm_gem;

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

TEST(PotentialAnalysisCommandTest, InvertedSamplingRangeBecomesPrepareValidationError)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.sampling_range_min = 2.0;
    request.sampling_range_max = 1.0;
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--sampling-range",
            rg::ValidationPhase::Prepare,
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
