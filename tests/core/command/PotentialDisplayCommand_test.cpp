#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/PotentialDisplayCommand.hpp>

namespace rg = rhbm_gem;

TEST(PotentialDisplayCommandTest, MalformedReferenceModelKeyListBecomesValidationError) {
    rg::PotentialDisplayCommand command{};
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = "model_key";
    request.ref_model_key_tag_list = "invalid";
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--ref-model-keylist",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialDisplayCommandTest, InvalidPainterChoiceBecomesValidationError) {
    rg::PotentialDisplayCommand command{};
    rg::PotentialDisplayRequest request{};
    request.painter_choice = static_cast<rg::PainterType>(999);
    request.model_key_tag_list = "model_key";
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--painter",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialDisplayCommandTest, WellFormedReferenceModelKeyListPassesPrepareValidation) {
    rg::PotentialDisplayCommand command{};
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = "model_key";
    request.ref_model_key_tag_list = "[with_charge]ref_a,ref_b;[no_charge]ref_c";
    command.ApplyRequest(request);

    EXPECT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(
        command_test::FindValidationIssue(
            command,
            "--ref-model-keylist",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}
