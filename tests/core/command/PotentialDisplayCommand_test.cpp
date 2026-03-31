#include <gtest/gtest.h>

#include "support/CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/command/PotentialDisplayCommand.hpp"

namespace rg = rhbm_gem;

TEST(PotentialDisplayCommandTest, MalformedReferenceModelKeyListBecomesValidationError) {
    rg::PotentialDisplayCommand command{};
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups[""] = { "invalid" };
    command.ApplyRequest(request);

    EXPECT_FALSE(command.Run());
    EXPECT_FALSE(command.WasPrepared());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--ref-group",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialDisplayCommandTest, WellFormedReferenceModelKeyListPassesPrepareValidation) {
    rg::PotentialDisplayCommand command{};
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups = {
        { "with_charge", { "ref_a", "ref_b" } },
        { "no_charge", { "ref_c" } },
    };
    command.ApplyRequest(request);

    EXPECT_FALSE(command.Run());
    EXPECT_TRUE(command.WasPrepared());
    EXPECT_EQ(
        command_test::FindValidationIssue(
            command,
            "--ref-group",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}
