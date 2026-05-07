#include <gtest/gtest.h>

#include "support/CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandSystem.hpp>

namespace rg = rhbm_gem;

TEST(CommandValidationScenariosTest, MapSimulationRejectsAllInvalidBlurringWidthsAtPrepare)
{
    rg::MapSimulationRequest request{};
    request.blurring_width_list = { -1.0, 0.0 };

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--blurring-width"),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRequiresPositiveResolutionWhenSimulationEnabled)
{
    rg::PotentialAnalysisRequest request{};
    request.simulation_flag = true;
    request.simulated_map_resolution = 0.0;

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--sim-resolution"),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRejectsInvertedSamplingRangeAtPrepare)
{
    rg::PotentialAnalysisRequest request{};
    request.sampling_range_min = 2.0;
    request.sampling_range_max = 1.0;

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--sampling-range"),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRejectsInvertedTrainingAlphaRangeAtPrepare)
{
    rg::PotentialAnalysisRequest request{};
    request.training_alpha_min = 1.0;
    request.training_alpha_max = 0.5;

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--training-alpha-range"),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisCoercesInvalidTrainingAlphaStepAtParse)
{
    rg::PotentialAnalysisRequest request{};
    request.training_alpha_step = 0.0;

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--training-alpha-step"),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRejectsEmptySavedKeyAtParse)
{
    rg::PotentialAnalysisRequest request{};
    request.saved_key_tag = "";

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--save-key"),
        nullptr);
}

TEST(CommandValidationScenariosTest, RHBMTestRejectsInvertedFitRangeAtPrepare)
{
    rg::RHBMTestRequest request{};
    request.fit_range_min = 2.0;
    request.fit_range_max = 1.0;

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--fit-range"),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialDisplayRejectsMalformedReferenceGroups)
{
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups[""] = { "invalid" };

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--ref-group"),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialDisplayAllowsWellFormedReferenceGroupsPastPrepare)
{
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups = {
        { "with_charge", { "ref_a", "ref_b" } },
        { "no_charge", { "ref_c" } },
    };

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(
        command_test::FindValidationIssue(
            result.issues,
            "--ref-group"),
        nullptr);
}

TEST(CommandValidationScenariosTest, ResultDumpRequiresMapFileForMapPrinterAtPrepare)
{
    rg::ResultDumpRequest request{};
    request.printer_choice = rg::PrinterType::MAP_VALUE;
    request.model_key_tag_list = { "model" };

    const auto result{ rg::RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--map"),
        nullptr);
}
