#include <gtest/gtest.h>

#include <fstream>

#include "CommandTestHelpers.hpp"
#include "CommandValidationAssertions.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>

namespace rg = rhbm_gem;

TEST(CommandApiTest, MapSimulationRequiresModelFilePath)
{
    const auto result{ rg::RunMapSimulation(rg::MapSimulationRequest{}) };

    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--model",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(CommandApiTest, MapSimulationRejectsInvalidBlurringWidthList)
{
    rg::MapSimulationRequest request;
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = "-1.0,0.0";

    const auto result{ rg::RunMapSimulation(request) };
    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--blurring-width",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandApiTest, PotentialAnalysisValidatesSimulationResolution)
{
    command_test::ScopedTempDir temp_dir{ "command_api_potential_analysis" };
    const auto dummy_map_path{ temp_dir.path() / "dummy.map" };
    std::ofstream(dummy_map_path) << "dummy";

    rg::PotentialAnalysisRequest request;
    request.common.database_path = temp_dir.path() / "demo.sqlite";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.map_file_path = dummy_map_path;
    request.simulation_flag = true;
    request.simulated_map_resolution = 0.0;

    const auto result{ rg::RunPotentialAnalysis(request) };
    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--sim-resolution",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandApiTest, ResultDumpMapPathIsOnlyRequiredForMapPrinter)
{
    command_test::ScopedTempDir temp_dir{ "command_api_result_dump" };

    rg::ResultDumpRequest map_request;
    map_request.common.database_path = temp_dir.path() / "demo.sqlite";
    map_request.model_key_tag_list = "demo_model";
    map_request.printer_choice = rg::PrinterType::MAP_VALUE;

    const auto map_result{ rg::RunResultDump(map_request) };
    EXPECT_FALSE(map_result.prepared);
    EXPECT_FALSE(map_result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            map_result.validation_issues,
            "--map",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);

    rg::ResultDumpRequest gaus_request;
    gaus_request.common.database_path = temp_dir.path() / "demo.sqlite";
    gaus_request.model_key_tag_list = "demo_model";
    gaus_request.printer_choice = rg::PrinterType::GAUS_ESTIMATES;

    const auto gaus_result{ rg::RunResultDump(gaus_request) };
    EXPECT_TRUE(gaus_result.prepared);
    EXPECT_FALSE(gaus_result.executed);
    EXPECT_EQ(
        command_test::FindValidationIssue(
            gaus_result.validation_issues,
            "--map",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}
