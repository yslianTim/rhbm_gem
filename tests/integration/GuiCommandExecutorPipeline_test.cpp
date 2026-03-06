#include <gtest/gtest.h>

#include <filesystem>

#include "CommandTestHelpers.hpp"
#include "GuiCommandExecutor.hpp"

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

TEST(GuiCommandExecutorPipelineTest, ExecutesSimulationAnalysisAndDumpPipeline)
{
    command_test::ScopedTempDir temp_dir{ "gui_executor_pipeline" };
    const auto maps_dir{ temp_dir.path() / "maps" };
    const auto analysis_output_dir{ temp_dir.path() / "analysis_output" };
    const auto dump_output_dir{ temp_dir.path() / "dump_output" };
    const auto database_path{ temp_dir.path() / "pipeline.sqlite" };

    std::filesystem::create_directories(maps_dir);
    std::filesystem::create_directories(analysis_output_dir);
    std::filesystem::create_directories(dump_output_dir);

    rg::gui::MapSimulationRequest simulation_request;
    simulation_request.common.folder_path = maps_dir;
    simulation_request.model_file_path = command_test::TestDataPath("test_model.cif");
    simulation_request.map_file_name = "sim_map";
    simulation_request.blurring_width_list = "1.50";

    const auto simulation_result{
        rg::gui::GuiCommandExecutor::ExecuteMapSimulation(simulation_request)
    };
    ASSERT_TRUE(simulation_result.prepared);
    ASSERT_TRUE(simulation_result.executed);

    const auto generated_map_file{ FindGeneratedMap(maps_dir) };
    ASSERT_FALSE(generated_map_file.empty());
    ASSERT_TRUE(std::filesystem::exists(generated_map_file));

    rg::gui::PotentialAnalysisRequest analysis_request;
    analysis_request.common.database_path = database_path;
    analysis_request.common.folder_path = analysis_output_dir;
    analysis_request.model_file_path = command_test::TestDataPath("test_model.cif");
    analysis_request.map_file_path = generated_map_file;
    analysis_request.saved_key_tag = "gui_pipeline_model";
    analysis_request.sampling_size = 200;

    const auto analysis_result{
        rg::gui::GuiCommandExecutor::ExecutePotentialAnalysis(analysis_request)
    };
    ASSERT_TRUE(analysis_result.prepared);
    ASSERT_TRUE(analysis_result.executed);

    rg::gui::ResultDumpRequest dump_request;
    dump_request.common.database_path = database_path;
    dump_request.common.folder_path = dump_output_dir;
    dump_request.printer_choice = rg::PrinterType::GAUS_ESTIMATES;
    dump_request.model_key_tag_list = "gui_pipeline_model";

    const auto dump_result{
        rg::gui::GuiCommandExecutor::ExecuteResultDump(dump_request)
    };
    ASSERT_TRUE(dump_result.prepared);
    ASSERT_TRUE(dump_result.executed);
    EXPECT_GT(CountRegularFiles(dump_output_dir), 0u);
}
