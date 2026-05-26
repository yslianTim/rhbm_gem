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
    const std::filesystem::path & temp_root)
{
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_path{
        command_test::GenerateMapFile(temp_root / "map", model_path, "fixture_map")
    };
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
    auto model{ LoadAnalyzedModelFixture(temp_dir.path()) };
    ASSERT_NE(model, nullptr);

    EXPECT_NO_THROW(rg::painter_internal::RequireLocalAnalyzedModel(*model, "AtomPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "ModelPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "GausPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "ComparisonPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "DemoPainter"));
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
        rgc::PaintModel(model_objects, output_folder),
        std::runtime_error);
    EXPECT_THROW(
        rgc::PaintGaus(model_objects, output_folder),
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
