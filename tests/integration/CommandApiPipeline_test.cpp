#include <gtest/gtest.h>

#include <filesystem>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>

namespace rg = rhbm_gem;

namespace {

std::filesystem::path FindGeneratedMap(const std::filesystem::path & directory)
{
    for (const auto & entry : std::filesystem::directory_iterator(directory))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".map")
        {
            return entry.path();
        }
    }
    return {};
}

std::size_t CountRegularFiles(const std::filesystem::path & directory)
{
    if (!std::filesystem::exists(directory)) return 0;

    std::size_t count{ 0 };
    for (const auto & entry : std::filesystem::directory_iterator(directory))
    {
        if (entry.is_regular_file())
        {
            ++count;
        }
    }
    return count;
}

} // namespace

TEST(CommandApiPipelineTest, ExecutesSimulationAnalysisAndDumpPipeline)
{
    command_test::ScopedTempDir temp_dir{ "command_executor_pipeline" };
    const auto maps_dir{ temp_dir.path() / "maps" };
    const auto analysis_output_dir{ temp_dir.path() / "analysis_output" };
    const auto dump_output_dir{ temp_dir.path() / "dump_output" };
    const auto database_path{ temp_dir.path() / "pipeline.sqlite" };

    std::filesystem::create_directories(maps_dir);
    std::filesystem::create_directories(analysis_output_dir);
    std::filesystem::create_directories(dump_output_dir);

    rg::MapSimulationRequest simulation_request;
    simulation_request.output_dir = maps_dir;
    simulation_request.model_file_path = command_test::TestDataPath("test_model.cif");
    simulation_request.map_file_name = "sim_map";
    simulation_request.blurring_width_list = { 1.50 };

    const auto simulation_result{
        rg::RunMapSimulation(simulation_request)
    };
    ASSERT_TRUE(simulation_result.succeeded);

    const auto generated_map_file{ FindGeneratedMap(maps_dir) };
    ASSERT_FALSE(generated_map_file.empty());
    ASSERT_TRUE(std::filesystem::exists(generated_map_file));

    rg::PotentialAnalysisRequest analysis_request;
    analysis_request.database_path = database_path;
    analysis_request.output_dir = analysis_output_dir;
    analysis_request.model_file_path = command_test::TestDataPath("test_model.cif");
    analysis_request.map_file_path = generated_map_file;
    analysis_request.saved_key_tag = "pipeline_model";
    analysis_request.sampling_size = 200;

    const auto analysis_result{
        rg::RunPotentialAnalysis(analysis_request)
    };
    ASSERT_TRUE(analysis_result.succeeded);

    rg::ResultDumpRequest dump_request;
    dump_request.database_path = database_path;
    dump_request.output_dir = dump_output_dir;
    dump_request.printer_choice = rg::PrinterType::GAUS_ESTIMATES;
    dump_request.model_key_tag_list = { "pipeline_model" };

    const auto dump_result{
        rg::RunResultDump(dump_request)
    };
    ASSERT_TRUE(dump_result.succeeded);
    EXPECT_GT(CountRegularFiles(dump_output_dir), 0u);
}

TEST(CommandApiPipelineTest, PotentialAnalysisTrainingEmitsRequestedAlphaRReport)
{
    command_test::ScopedTempDir temp_dir{ "command_executor_training_reports" };
    const auto maps_dir{ temp_dir.path() / "maps" };
    const auto analysis_output_dir{ temp_dir.path() / "analysis_output" };
    const auto training_report_dir{ temp_dir.path() / "training_reports" };
    const auto dump_output_dir{ temp_dir.path() / "dump_output" };
    const auto database_path{ temp_dir.path() / "training.sqlite" };

    std::filesystem::create_directories(maps_dir);
    std::filesystem::create_directories(analysis_output_dir);
    std::filesystem::create_directories(dump_output_dir);

    rg::MapSimulationRequest simulation_request;
    simulation_request.output_dir = maps_dir;
    simulation_request.model_file_path = command_test::TestDataPath("test_model.cif");
    simulation_request.map_file_name = "sim_map";
    simulation_request.blurring_width_list = { 1.50 };

    const auto simulation_result{
        rg::RunMapSimulation(simulation_request)
    };
    ASSERT_TRUE(simulation_result.succeeded);

    const auto generated_map_file{ FindGeneratedMap(maps_dir) };
    ASSERT_FALSE(generated_map_file.empty());
    ASSERT_TRUE(std::filesystem::exists(generated_map_file));

    rg::PotentialAnalysisRequest analysis_request;
    analysis_request.database_path = database_path;
    analysis_request.output_dir = analysis_output_dir;
    analysis_request.model_file_path = command_test::TestDataPath("test_model.cif");
    analysis_request.map_file_path = generated_map_file;
    analysis_request.saved_key_tag = "trained_model";
    analysis_request.training_alpha_flag = true;
    analysis_request.training_report_dir = training_report_dir;
    analysis_request.sampling_size = 600;

    const auto analysis_result{
        rg::RunPotentialAnalysis(analysis_request)
    };
    ASSERT_TRUE(analysis_result.succeeded);
    EXPECT_TRUE(std::filesystem::exists(training_report_dir / "alpha_r_bias.pdf"));

    rg::ResultDumpRequest dump_request;
    dump_request.database_path = database_path;
    dump_request.output_dir = dump_output_dir;
    dump_request.printer_choice = rg::PrinterType::GAUS_ESTIMATES;
    dump_request.model_key_tag_list = { "trained_model" };

    const auto dump_result{
        rg::RunResultDump(dump_request)
    };
    ASSERT_TRUE(dump_result.succeeded);
    EXPECT_GT(CountRegularFiles(dump_output_dir), 0u);
}
