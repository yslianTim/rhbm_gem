#include <gtest/gtest.h>

#include <algorithm>

#include "HRLModelTestCommand.hpp"

namespace rg = rhbm_gem;

TEST(HRLModelTestCommandTest, InvalidTesterChoiceBecomesValidationError)
{
    rg::HRLModelTestCommand command;
    command.SetTesterChoice(static_cast<rg::TesterType>(999));

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--tester";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(HRLModelTestCommandTest, NegativeAlphaValuesBecomeValidationErrors)
{
    rg::HRLModelTestCommand command;
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

TEST(HRLModelTestCommandTest, FitRangeOrderingBecomesPrepareValidationError)
{
    rg::HRLModelTestCommand command;
    command.SetFitRangeMinimum(2.0);
    command.SetFitRangeMaximum(1.0);

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--fit-range";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Prepare);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}
