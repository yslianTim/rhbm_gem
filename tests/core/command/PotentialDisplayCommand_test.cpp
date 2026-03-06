#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include "PotentialDisplayCommand.hpp"

namespace rg = rhbm_gem;

TEST(PotentialDisplayCommandTest, MalformedReferenceModelKeyListBecomesValidationError)
{
    rg::PotentialDisplayCommand command;
    command.SetPainterChoice(rg::PainterType::MODEL);
    command.SetModelKeyTagList("model_key");

    EXPECT_NO_THROW(command.SetRefModelKeyTagListMap("invalid"));
    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--ref-model-keylist",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialDisplayCommandTest, InvalidPainterChoiceBecomesValidationError)
{
    rg::PotentialDisplayCommand command;
    command.SetPainterChoice(static_cast<rg::PainterType>(999));
    command.SetModelKeyTagList("model_key");

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--painter",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialDisplayCommandTest, WellFormedReferenceModelKeyListPassesPrepareValidation)
{
    rg::PotentialDisplayCommand command;
    command.SetPainterChoice(rg::PainterType::MODEL);
    command.SetModelKeyTagList("model_key");
    command.SetRefModelKeyTagListMap("[with_charge]ref_a,ref_b;[no_charge]ref_c");

    EXPECT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(
        command_test::FindValidationIssue(
            command,
            "--ref-model-keylist",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}
