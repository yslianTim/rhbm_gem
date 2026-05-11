#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandSystem.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

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

void ExpectSelectedAtomsHaveAlphaR(const rg::ModelObject & model, double alpha_r)
{
    const auto & atom_list{ model.GetSelectedAtoms() };
    ASSERT_GT(atom_list.size(), 0u);
    for (const auto * atom : atom_list)
    {
        EXPECT_DOUBLE_EQ(alpha_r, rg::LocalPotentialView::RequireFor(*atom).GetAlphaR());
    }
}

void ExpectSelectedAtomsShareAlphaR(const rg::ModelObject & model)
{
    const auto & atom_list{ model.GetSelectedAtoms() };
    ASSERT_GT(atom_list.size(), 0u);
    const auto alpha_r{ rg::LocalPotentialView::RequireFor(*atom_list.front()).GetAlphaR() };
    ASSERT_TRUE(std::isfinite(alpha_r));
    for (const auto * atom : atom_list)
    {
        EXPECT_DOUBLE_EQ(alpha_r, rg::LocalPotentialView::RequireFor(*atom).GetAlphaR());
    }
}

} // namespace

TEST(CommandApiPipelineTest, PotentialAnalysisRequestExposesTrainingAlphaGridControls)
{
    rg::PotentialAnalysisRequest request;
    EXPECT_DOUBLE_EQ(request.training_alpha_min, 0.0);
    EXPECT_DOUBLE_EQ(request.training_alpha_max, 1.0);
    EXPECT_DOUBLE_EQ(request.training_alpha_step, 0.1);

    request.training_alpha_min = 0.2;
    request.training_alpha_max = 0.8;
    request.training_alpha_step = 0.2;
    EXPECT_DOUBLE_EQ(request.training_alpha_min, 0.2);
    EXPECT_DOUBLE_EQ(request.training_alpha_max, 0.8);
    EXPECT_DOUBLE_EQ(request.training_alpha_step, 0.2);
}

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
        rg::RunCommand(simulation_request)
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
    analysis_request.alpha_r = 0.37;

    const auto analysis_result{
        rg::RunCommand(analysis_request)
    };
    ASSERT_TRUE(analysis_result.succeeded);

    rg::DataRepository repository{ database_path };
    auto model{ repository.LoadModel("pipeline_model") };
    ASSERT_NE(model, nullptr);
    ExpectSelectedAtomsHaveAlphaR(*model, analysis_request.alpha_r);

    rg::ResultDumpRequest dump_request;
    dump_request.database_path = database_path;
    dump_request.output_dir = dump_output_dir;
    dump_request.printer_choice = rg::PrinterType::GAUS_ESTIMATES;
    dump_request.model_key_tag_list = { "pipeline_model" };

    const auto dump_result{
        rg::RunCommand(dump_request)
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
        rg::RunCommand(simulation_request)
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

    const auto analysis_result{
        rg::RunCommand(analysis_request)
    };
    ASSERT_TRUE(analysis_result.succeeded);
    EXPECT_TRUE(std::filesystem::exists(training_report_dir / "alpha_r_bias.pdf"));

    rg::DataRepository repository{ database_path };
    auto model{ repository.LoadModel("trained_model") };
    ASSERT_NE(model, nullptr);
    ExpectSelectedAtomsShareAlphaR(*model);

    rg::ResultDumpRequest dump_request;
    dump_request.database_path = database_path;
    dump_request.output_dir = dump_output_dir;
    dump_request.printer_choice = rg::PrinterType::GAUS_ESTIMATES;
    dump_request.model_key_tag_list = { "trained_model" };

    const auto dump_result{
        rg::RunCommand(dump_request)
    };
    ASSERT_TRUE(dump_result.succeeded);
    EXPECT_GT(CountRegularFiles(dump_output_dir), 0u);
}
