#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <stdexcept>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rg = rhbm_gem;
namespace rgc = rhbm_gem::core;

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

void ExpectSelectedAtomsHaveFiniteNonNegativeAlphaR(const rg::ModelObject & model)
{
    const auto & atom_list{ model.GetSelectedAtoms() };
    ASSERT_GT(atom_list.size(), 0u);
    for (const auto * atom : atom_list)
    {
        const auto alpha_r{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetAlphaR(
                rg::LocalFittingStage::Third)
        };
        EXPECT_TRUE(std::isfinite(alpha_r));
        EXPECT_GE(alpha_r, 0.0);
    }
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

    rgc::MapSimulationRequest simulation_request;
    simulation_request.output_dir = maps_dir;
    simulation_request.model_file_path = command_test::TestDataPath("test_model.cif");
    simulation_request.map_file_name = "sim_map";
    simulation_request.blurring_width_list = { 1.50 };

    const auto simulation_result{
        rgc::RunCommand(simulation_request)
    };
    ASSERT_TRUE(simulation_result.succeeded);

    const auto generated_map_file{ FindGeneratedMap(maps_dir) };
    ASSERT_FALSE(generated_map_file.empty());
    ASSERT_TRUE(std::filesystem::exists(generated_map_file));

    rgc::PotentialAnalysisRequest analysis_request;
    analysis_request.database_path = database_path;
    analysis_request.output_dir = analysis_output_dir;
    analysis_request.model_file_path = command_test::TestDataPath("test_model.cif");
    analysis_request.map_file_path = generated_map_file;
    analysis_request.saved_key_tag = "pipeline/model test";

    const auto analysis_result{
        rgc::RunCommand(analysis_request)
    };
    ASSERT_TRUE(analysis_result.succeeded);
    EXPECT_TRUE(std::filesystem::exists(
        analysis_output_dir / "local_fitting_result_pipeline_model_test.csv"));

    rg::DataRepository repository{ database_path };
    auto model{ repository.LoadModel("pipeline/model test") };
    ASSERT_NE(model, nullptr);
    ExpectSelectedAtomsHaveFiniteNonNegativeAlphaR(*model);
    ASSERT_EQ(model->GetAtomList().size(), 1u);
    EXPECT_TRUE(std::isfinite(model->GetAtomList().front()->GetStandardQScore()));
    EXPECT_DOUBLE_EQ(
        model->GetStandardAverageQScore(),
        model->GetAtomList().front()->GetStandardQScore());

    const auto failure_output_dir{ temp_dir.path() / "analysis_failure_output" };
    std::filesystem::create_directories(
        failure_output_dir / "local_fitting_result_write_failure.csv");
    auto failure_request{ analysis_request };
    failure_request.output_dir = failure_output_dir;
    failure_request.saved_key_tag = "write/failure";
    const auto failure_result{ rgc::RunCommand(failure_request) };
    EXPECT_FALSE(failure_result.succeeded);
    EXPECT_THROW(
        (void)repository.LoadModel("write/failure"),
        std::runtime_error);

    rgc::ResultDumpRequest dump_request;
    dump_request.database_path = database_path;
    dump_request.output_dir = dump_output_dir;
    dump_request.printer_choice = rgc::PrinterType::GAUS_ESTIMATES;
    dump_request.model_key_tag_list = { "pipeline/model test" };

    const auto dump_result{
        rgc::RunCommand(dump_request)
    };
    ASSERT_TRUE(dump_result.succeeded);
    EXPECT_GT(CountRegularFiles(dump_output_dir), 0u);
}
