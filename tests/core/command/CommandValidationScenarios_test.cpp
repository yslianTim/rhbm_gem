#include <gtest/gtest.h>

#include "support/CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/RHBMTestCommand.hpp"
#include "command/MapSimulationCommand.hpp"
#include "command/PotentialAnalysisCommand.hpp"
#include "command/PotentialDisplayCommand.hpp"
#include "command/ResultDumpCommand.hpp"

namespace rg = rhbm_gem;

namespace {

template <typename Command, typename Request>
bool ApplyAndRun(Command & command, const Request & request)
{
    command.ApplyRequest(request);
    return command.Run();
}

} // namespace

TEST(CommandValidationScenariosTest, MapSimulationRejectsAllInvalidBlurringWidthsAtPrepare)
{
    rg::MapSimulationCommand command{};
    rg::MapSimulationRequest request{};
    request.blurring_width_list = { -1.0, 0.0 };

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--blurring-width",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRequiresPositiveResolutionWhenSimulationEnabled)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.simulation_flag = true;
    request.simulated_map_resolution = 0.0;

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--sim-resolution",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRejectsInvertedSamplingRangeAtPrepare)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.sampling_range_min = 2.0;
    request.sampling_range_max = 1.0;

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--sampling-range",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRejectsInvertedTrainingAlphaRangeAtPrepare)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.training_alpha_min = 1.0;
    request.training_alpha_max = 0.5;

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--training-alpha-range",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisCoercesInvalidTrainingAlphaStepAtParse)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.training_alpha_step = 0.0;

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--training-alpha-step",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialAnalysisRejectsEmptySavedKeyAtParse)
{
    rg::PotentialAnalysisCommand command{};
    rg::PotentialAnalysisRequest request{};
    request.saved_key_tag = "";

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--save-key",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, RHBMTestRejectsInvertedFitRangeAtPrepare)
{
    rg::RHBMTestCommand command{};
    rg::RHBMTestRequest request{};
    request.fit_range_min = 2.0;
    request.fit_range_max = 1.0;

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--fit-range",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialDisplayRejectsMalformedReferenceGroups)
{
    rg::PotentialDisplayCommand command{};
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups[""] = { "invalid" };

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--ref-group",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, PotentialDisplayAllowsWellFormedReferenceGroupsPastPrepare)
{
    rg::PotentialDisplayCommand command{};
    rg::PotentialDisplayRequest request{};
    request.painter_choice = rg::PainterType::MODEL;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups = {
        { "with_charge", { "ref_a", "ref_b" } },
        { "no_charge", { "ref_c" } },
    };

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_TRUE(command.WasPrepared());
    EXPECT_EQ(
        command_test::FindValidationIssue(
            command,
            "--ref-group",
            std::nullopt,
            LogLevel::Error),
        nullptr);
}

TEST(CommandValidationScenariosTest, ResultDumpRequiresMapFileForMapPrinterAtPrepare)
{
    rg::ResultDumpCommand command{};
    rg::ResultDumpRequest request{};
    request.printer_choice = rg::PrinterType::MAP_VALUE;
    request.model_key_tag_list = { "model" };

    EXPECT_FALSE(ApplyAndRun(command, request));
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--map",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}
