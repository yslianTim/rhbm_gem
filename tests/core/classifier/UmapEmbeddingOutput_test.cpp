#include <gtest/gtest.h>
#include <boost/json.hpp>

#include "core/classifier/UmapEmbeddingAlgorithm.hpp"
#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/classifier/UmapEmbeddingOutput.hpp"
#include "core/classifier/UmapEmbeddingPlot.hpp"
#include "support/CommandTestHelpers.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
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
        row.features[4] = std::numeric_limits<double>::quiet_NaN();
        row.features[5] = 4.0;
        row.features[6] = -2.0;
    }
    UmapEmbeddingResult first;
    first.coordinates = { -1.0, 1.0, 0.0, 0.0, 1.0, -1.0 };
    first.mean_shift_coordinates = first.coordinates;
    first.binning = BinMeanShiftEmbedding(
        prepared.feature_rows, first.mean_shift_coordinates, GetUmapSpatialScale(first.coordinates));
    PopulateUmapBinFeatureMedians(prepared.feature_rows, first.binning);
    ASSERT_TRUE(BuildSecondUmapFeatures(prepared.feature_rows, first.binning, "model").values.empty());
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
    EXPECT_EQ(line, first_header + ",local carbon ratio,local nitrogen ratio,local oxygen ratio"
        ",global carbon ratio,global nitrogen ratio,global oxygen ratio");
    EXPECT_FALSE(std::getline(second_bins, line));
    for (std::size_t bin = 0; bin < prepared.feature_rows.size(); ++bin)
    {
        ASSERT_TRUE(std::getline(first_bins, line));
        EXPECT_TRUE(line.ends_with(",1,,4,-2,1,CA,1,1,,0,0,,0,0")) << line;
    }
    EXPECT_FALSE(std::getline(first_bins, line));

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

TEST(UmapEmbeddingOutputTest, CountsElementRatiosInSecondMeanShiftBins)
{
    const command_test::ScopedTempDir dir{ "umap_second_bin_element_ratios" };
    const auto output_path{ GetEmbeddingOutputPath(dir.path(), "model") };
    const auto second_summary{ dir.path() / "umap_embedding_model_second_umap_mean_shift_bins.csv" };

    PreparedUmapInput prepared;
    prepared.elements = { Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
        Element::SULFUR, Element::UNK, Element::OXYGEN };
    prepared.feature_rows.resize(prepared.elements.size());
    // Use the same configured spot for every element to exercise the Element-based count.
    prepared.configured_spot_indices.resize(prepared.elements.size(), 1);
    for (std::size_t point = 0; point < prepared.feature_rows.size(); ++point)
    {
        auto & row{ prepared.feature_rows[point] };
        row.serial_id = static_cast<int>(point + 1);
        row.residue = "ALA";
        row.spot = "CA";
    }
    prepared.feature_rows.back().features[4] = 1.0;

    UmapEmbeddingResult first;
    first.coordinates = { 0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0, 0.0, 6.0, 0.0 };
    first.mean_shift_coordinates = first.coordinates;
    first.binning = BinMeanShiftEmbedding(
        prepared.feature_rows, first.mean_shift_coordinates, GetUmapSpatialScale(first.coordinates));
    PopulateUmapBinFeatureMedians(prepared.feature_rows, first.binning);

    UmapEmbeddingResult second;
    second.zero_mad_fallback_features = BuildSecondUmapFeatures(
        prepared.feature_rows, first.binning, "model").zero_mad_fallback_features;
    second.coordinates = first.coordinates;
    second.mean_shift_coordinates = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0 };
    second.binning = BinMeanShiftEmbedding(
        prepared.feature_rows, second.mean_shift_coordinates, GetUmapSpatialScale(second.coordinates));

    std::string error;
    ASSERT_TRUE(WriteEmbeddingOutputs(output_path, prepared, first, second, error)) << error;
    std::ifstream metadata_file{ output_path.string() + ".metadata.json" };
    std::string metadata_line;
    ASSERT_TRUE(std::getline(metadata_file, metadata_line));
    EXPECT_EQ(boost::json::parse(metadata_line).at("second_umap_zero_mad_fallback_features").as_array(),
        (boost::json::array{ "amplitude joint median of bin", "amplitude joint" }));
    std::ifstream bins{ second_summary };
    std::string line;
    ASSERT_TRUE(std::getline(bins, line));
    ASSERT_TRUE(std::getline(bins, line));
    std::istringstream mixed_bin{ line };
    std::string field;
    for (std::size_t column = 0; column < 17; ++column)
    {
        ASSERT_TRUE(std::getline(mixed_bin, field, ','));
    }
    for (const double ratio : std::array{ 2.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0, 1.0, 0.5 })
    {
        ASSERT_TRUE(std::getline(mixed_bin, field, ','));
        EXPECT_DOUBLE_EQ(std::stod(field), ratio);
    }
    EXPECT_FALSE(std::getline(mixed_bin, field, ','));

    ASSERT_TRUE(std::getline(bins, line));
    EXPECT_TRUE(line.ends_with(",0,0,1,0,0,0.5")) << line;
    EXPECT_FALSE(std::getline(bins, line));
}

} // namespace rhbm_gem::core::classifier
