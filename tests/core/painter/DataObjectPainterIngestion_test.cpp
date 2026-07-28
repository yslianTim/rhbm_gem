#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <stdexcept>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/core/PainterFunctions.hpp>
#include "core/painter/detail/PainterModelValidation.hpp"
#include <rhbm_gem/data/io/DataRepository.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;
namespace rgc = rhbm_gem::core;

namespace {

std::shared_ptr<rg::ModelObject> LoadAnalyzedModelFixture(
    const std::filesystem::path & temp_root,
    const std::filesystem::path & map_path)
{
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto database_path{ temp_root / "db" / "database.sqlite" };

    rgc::PotentialAnalysisRequest request{};
    request.output_dir = temp_root / "analysis_out";
    request.database_path = database_path;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.saved_key_tag = "analyzed_model";

    const auto result{ rgc::RunCommand(request) };
    if (!result.succeeded)
    {
        throw std::runtime_error("Failed to build analysis-ready model fixture.");
    }

    rg::DataRepository repository{ database_path };
    auto model{ repository.LoadModel("analyzed_model") };
    return std::shared_ptr<rg::ModelObject>{ std::move(model) };
}

} // namespace

TEST(DataObjectPainterIngestionTest, PainterValidationAcceptsAnalysisReadyModels)
{
    command_test::ScopedTempDir temp_dir{ "painter_ingestion_analyzed" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_path{
        command_test::GenerateMapFile(temp_dir.path() / "map", model_path, "fixture_map")
    };
    auto model{ LoadAnalyzedModelFixture(temp_dir.path(), map_path) };
    ASSERT_NE(model, nullptr);

    EXPECT_NO_THROW(rg::painter_internal::RequireLocalAnalyzedModel(*model, "AtomPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "GausPainter"));
    EXPECT_NO_THROW(rgc::PaintQScore(
        rgc::ModelObjectList{ model.get() },
        temp_dir.path().string()));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "ComparisonPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "DemoPainter"));
}

TEST(DataObjectPainterIngestionTest, PotentialDisplayQScoreRunsWithoutMap)
{
    command_test::ScopedTempDir temp_dir{ "potential_display_qscore_without_map" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_path{
        command_test::GenerateMapFile(temp_dir.path() / "map", model_path, "fixture_map")
    };
    const auto model{ LoadAnalyzedModelFixture(temp_dir.path(), map_path) };
    ASSERT_NE(model, nullptr);

    const auto output_dir{ temp_dir.path() / "display_out" };
    rgc::PotentialDisplayRequest request{};
    request.database_path = temp_dir.path() / "db" / "database.sqlite";
    request.output_dir = output_dir;
    request.painter_choice = rgc::PainterType::QSCORE;
    request.model_key_tag_list = { "analyzed_model" };

    const auto result{ rgc::RunCommand(request) };

    ASSERT_TRUE(result.succeeded);
    bool found_qscore_output{ false };
    for (const auto & entry : std::filesystem::directory_iterator(output_dir))
    {
        const auto file_name{ entry.path().filename().string() };
        if (entry.is_regular_file()
            && file_name.find("average_qscore_to_sequence_summary_") == 0
            && entry.path().extension() == ".pdf")
        {
            found_qscore_output = true;
            break;
        }
    }
    EXPECT_TRUE(found_qscore_output);
}

TEST(DataObjectPainterIngestionTest, PainterFunctionEntrypointsRejectModelsWithoutAnalysisData)
{
    command_test::ScopedTempDir temp_dir{ "painter_ingestion_invalid" };
    auto model{ data_test::MakeModelWithBond() };
    ASSERT_NE(model, nullptr);
    model->SetKeyTag("raw_model");
    model->SelectAllAtoms();
    model->SelectAllBonds();

    const rgc::ModelObjectList model_objects{ model.get() };
    const rgc::ReferenceModelGroupMap reference_model_groups{
        { "ref", { model.get() } }
    };
    const auto output_folder{ temp_dir.path().string() };

    EXPECT_THROW(
        rgc::PaintAtom(model_objects, output_folder),
        std::runtime_error);
    EXPECT_THROW(
        rgc::PaintGaus(model_objects, output_folder),
        std::runtime_error);
    EXPECT_THROW(
        rgc::PaintQScore(model_objects, output_folder),
        std::runtime_error);
    EXPECT_THROW(
        rgc::PaintComparison(
            model_objects,
            reference_model_groups,
            output_folder),
        std::runtime_error);
    EXPECT_THROW(
        rgc::PaintDemo(model_objects, reference_model_groups, output_folder),
        std::runtime_error);
}
