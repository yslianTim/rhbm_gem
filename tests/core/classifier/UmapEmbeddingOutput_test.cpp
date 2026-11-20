#include <gtest/gtest.h>

#include "core/classifier/UmapEmbeddingAlgorithm.hpp"
#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/classifier/UmapEmbeddingOutput.hpp"
#include "core/classifier/UmapEmbeddingPlot.hpp"
#include "support/CommandTestHelpers.hpp"

#include <filesystem>
#include <fstream>
#include <string>

namespace rhbm_gem::core::classifier {

TEST(UmapEmbeddingOutputTest, MissingSecondEmbeddingLeavesEmptyCoordinatesAndHeaderOnlySummary)
{
    const command_test::ScopedTempDir dir{ "umap_missing_second_embedding" };
    const auto output_path{ GetEmbeddingOutputPath(dir.path(), "model") };
    const auto first_summary{ dir.path() / "umap_embedding_model_mean_shift_bins.csv" };
    const auto second_summary{ dir.path() / "umap_embedding_model_second_umap_mean_shift_bins.csv" };

    PreparedUmapInput prepared;
    prepared.feature_rows.resize(3);
    prepared.configured_spot_indices = { 1, 1, 1 };
    prepared.elements = { Element::CARBON, Element::CARBON, Element::CARBON };
    for (std::size_t point = 0; point < prepared.feature_rows.size(); ++point)
    {
        auto & row{ prepared.feature_rows[point] };
        row.serial_id = static_cast<int>(point + 1);
        row.residue = "ALA";
        row.spot = "CA";
    }
    UmapEmbeddingResult first;
    first.coordinates = { -1.0, 1.0, 0.0, 0.0, 1.0, -1.0 };
    first.mean_shift_coordinates = first.coordinates;
    first.binning = BinMeanShiftEmbedding(
        prepared.feature_rows, first.mean_shift_coordinates, GetUmapSpatialScale(first.coordinates));
    const UmapEmbeddingResult second;

    std::ofstream(second_summary) << "previous summary\n";
    std::string error;
    ASSERT_TRUE(WriteEmbeddingOutputs(output_path, prepared, first, second, error)) << error;
    std::ifstream points{ output_path };
    std::string line;
    ASSERT_TRUE(std::getline(points, line));
    for (std::size_t point = 0; point < prepared.feature_rows.size(); ++point)
    {
        ASSERT_TRUE(std::getline(points, line));
        EXPECT_TRUE(line.ends_with(",,"));
    }
    EXPECT_FALSE(std::getline(points, line));

    std::ifstream first_bins{ first_summary };
    std::ifstream second_bins{ second_summary };
    std::string first_header;
    ASSERT_TRUE(std::getline(first_bins, first_header));
    ASSERT_TRUE(std::getline(second_bins, line));
    EXPECT_EQ(line, first_header);
    EXPECT_FALSE(std::getline(second_bins, line));

    const auto second_plot{ dir.path() / "umap_embedding_model_second_umap.pdf" };
    const auto second_shift_plot{ dir.path() / "umap_embedding_model_second_umap_mean_shift.pdf" };
#ifdef HAVE_ROOT
    std::ofstream(second_plot) << "previous second plot";
    std::ofstream(second_shift_plot) << "previous second mean shift plot";
#endif
    ASSERT_TRUE(WriteEmbeddingPlots(output_path, prepared, first, second, error)) << error;
    EXPECT_FALSE(std::filesystem::exists(second_plot));
    EXPECT_FALSE(std::filesystem::exists(second_shift_plot));
#ifdef HAVE_ROOT
    EXPECT_TRUE(std::filesystem::is_regular_file(dir.path() / "umap_embedding_model.pdf"));
    EXPECT_TRUE(std::filesystem::is_regular_file(dir.path() / "umap_embedding_model_mean_shift.pdf"));
#else
    EXPECT_FALSE(std::filesystem::exists(dir.path() / "umap_embedding_model.pdf"));
    EXPECT_FALSE(std::filesystem::exists(dir.path() / "umap_embedding_model_mean_shift.pdf"));
#endif
}

} // namespace rhbm_gem::core::classifier
