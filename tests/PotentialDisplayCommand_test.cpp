#include <gtest/gtest.h>

#include <algorithm>

#include "PotentialDisplayCommand.hpp"

namespace rg = rhbm_gem;

TEST(PotentialDisplayCommandTest, MalformedReferenceModelKeyListBecomesValidationError)
{
    rg::PotentialDisplayCommand command;
    command.SetPainterChoice(rg::PainterType::MODEL);
    command.SetModelKeyTagList("model_key");

    EXPECT_NO_THROW(command.SetRefModelKeyTagListMap("invalid"));
    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--ref-model-keylist";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(PotentialDisplayCommandTest, InvalidPainterChoiceBecomesValidationError)
{
    rg::PotentialDisplayCommand command;
    command.SetPainterChoice(static_cast<rg::PainterType>(999));
    command.SetModelKeyTagList("model_key");

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--painter";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}
