#include <gtest/gtest.h>

#include <algorithm>

#include "PotentialAnalysisCommand.hpp"

namespace rg = rhbm_gem;

TEST(PotentialAnalysisCommandTest, NegativeAlphaValuesBecomeValidationErrors)
{
    rg::PotentialAnalysisCommand command;
    command.SetAlphaR(-0.1);
    command.SetAlphaG(0.0);

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto alpha_r_issue{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--alpha-r";
            })
    };
    ASSERT_NE(alpha_r_issue, issues.end());
    EXPECT_EQ(alpha_r_issue->phase, rg::ValidationPhase::Parse);

    const auto alpha_g_issue{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--alpha-g";
            })
    };
    ASSERT_NE(alpha_g_issue, issues.end());
    EXPECT_EQ(alpha_g_issue->phase, rg::ValidationPhase::Parse);
}

TEST(PotentialAnalysisCommandTest, SimulationRequiresPositiveResolutionAtPrepare)
{
    rg::PotentialAnalysisCommand command;
    command.SetSimulationFlag(true);
    command.SetSimulatedMapResolution(0.0);

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--sim-resolution";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Prepare);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(PotentialAnalysisCommandTest, NonPositiveSamplingHeightBecomesValidationError)
{
    rg::PotentialAnalysisCommand command;
    command.SetSamplingHeight(0.0);

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--sampling-height";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(PotentialAnalysisCommandTest, EmptySavedKeyBecomesValidationError)
{
    rg::PotentialAnalysisCommand command;
    command.SetSavedKeyTag("");

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--save-key";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}
