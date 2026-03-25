#include <gtest/gtest.h>

#include "support/CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/PotentialDisplayCommand.hpp"

namespace rg = rhbm_gem;

TEST(PotentialDisplayCommandTest, MalformedReferenceModelKeyListBecomesValidationError) {
    rg::PotentialDisplayCommand command{rg::CommonOptionProfile::DatabaseWorkflow};
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

TEST(PotentialDisplayCommandTest, WellFormedReferenceModelKeyListPassesPrepareValidation) {
    rg::PotentialDisplayCommand command{rg::CommonOptionProfile::DatabaseWorkflow};
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
