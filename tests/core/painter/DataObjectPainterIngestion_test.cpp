#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <memory>
#include <stdexcept>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/core/PainterFunctions.hpp>
#include "core/painter/detail/PainterModelValidation.hpp"
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;
namespace rgc = rhbm_gem::core;

namespace {

rg::MapObject MakeConstantMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (std::size_t i = 0; i < 8; ++i)
    {
        values[i] = 1.0f;
    }
    return rg::MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

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
    auto map{ rg::ReadMap(map_path) };
    ASSERT_NE(model, nullptr);
    ASSERT_NE(map, nullptr);

    EXPECT_NO_THROW(rg::painter_internal::RequireLocalAnalyzedModel(*model, "AtomPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "GausPainter"));
    EXPECT_NO_THROW(rgc::PaintQScore(
        rgc::ModelObjectList{ model.get() },
        *map,
        temp_dir.path().string()));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "ComparisonPainter"));
    EXPECT_NO_THROW(rg::painter_internal::RequireGroupedAnalyzedModel(*model, "DemoPainter"));
}

TEST(DataObjectPainterIngestionTest, PainterFunctionEntrypointsRejectModelsWithoutAnalysisData)
{
    command_test::ScopedTempDir temp_dir{ "painter_ingestion_invalid" };
    auto model{ data_test::MakeModelWithBond() };
    auto map{ MakeConstantMapObject() };
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
        rgc::PaintQScore(model_objects, map, output_folder),
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
