#include <gtest/gtest.h>

#include <filesystem>
#include <memory>
#include <stdexcept>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

std::shared_ptr<rg::ModelObject> LoadAnalyzedModelFixture(
    const std::filesystem::path & temp_root)
{
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_path{
        command_test::GenerateMapFile(temp_root / "map", model_path, "fixture_map")
    };
    const auto database_path{ temp_root / "db" / "database.sqlite" };

    rg::PotentialAnalysisRequest request{};
    request.output_dir = temp_root / "analysis_out";
    request.database_path = database_path;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.saved_key_tag = "analyzed_model";

    const auto result{ rg::RunPotentialAnalysis(request) };
    if (!result.succeeded)
    {
        throw std::runtime_error("Failed to build analysis-ready model fixture.");
    }

    rg::DataRepository repository{ database_path };
    auto model{ repository.LoadModel("analyzed_model") };
    return std::shared_ptr<rg::ModelObject>{ std::move(model) };
}

} // namespace

TEST(DataObjectPainterIngestionTest, PainterTypedIngestionUsesAnalysisReadyModelApis)
{
    command_test::ScopedTempDir temp_dir{ "painter_ingestion_analyzed" };
    auto model{ LoadAnalyzedModelFixture(temp_dir.path()) };
    ASSERT_NE(model, nullptr);

    rg::AtomPainter atom_painter;
    EXPECT_NO_THROW(atom_painter.AddModel(*model));

    rg::ModelPainter model_painter;
    EXPECT_NO_THROW(model_painter.AddModel(*model));

    rg::GausPainter gaus_painter;
    EXPECT_NO_THROW(gaus_painter.AddModel(*model));

    rg::ComparisonPainter comparison_painter;
    EXPECT_NO_THROW(comparison_painter.AddModel(*model));
    EXPECT_NO_THROW(comparison_painter.AddReferenceModel(*model, "ref"));

    rg::DemoPainter demo_painter;
    EXPECT_NO_THROW(demo_painter.AddModel(*model));
    EXPECT_NO_THROW(demo_painter.AddReferenceModel(*model, "ref"));
}

TEST(DataObjectPainterIngestionTest, PainterContractsRejectModelsWithoutAnalysisData)
{
    auto model{ data_test::MakeModelWithBond() };
    ASSERT_NE(model, nullptr);
    model->SetKeyTag("raw_model");
    model->SelectAllAtoms();
    model->SelectAllBonds();

    rg::AtomPainter atom_painter;
    EXPECT_THROW(atom_painter.AddModel(*model), std::runtime_error);

    rg::ModelPainter model_painter;
    EXPECT_THROW(model_painter.AddModel(*model), std::runtime_error);

    rg::GausPainter gaus_painter;
    EXPECT_THROW(gaus_painter.AddModel(*model), std::runtime_error);

    rg::ComparisonPainter comparison_painter;
    EXPECT_THROW(comparison_painter.AddModel(*model), std::runtime_error);
    EXPECT_THROW(comparison_painter.AddReferenceModel(*model, "ref"), std::runtime_error);

    rg::DemoPainter demo_painter;
    EXPECT_THROW(demo_painter.AddModel(*model), std::runtime_error);
    EXPECT_THROW(demo_painter.AddReferenceModel(*model, "ref"), std::runtime_error);
}
