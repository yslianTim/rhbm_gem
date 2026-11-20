#include <gtest/gtest.h>

#include "core/classifier/UmapEmbeddingFeatures.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rhbm_gem::core::classifier {

TEST(UmapEmbeddingFeaturesTest, FiltersRequiredFeaturesAndUsesPopulationZScoreInRowOrder)
{
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    const std::array feature_values{
        std::array{ 1.0, 8.0, missing, 9.0, 10.0, 2.0, 3.0 },
        std::array{ 50.0, 5.0, 0.0, missing, 10.0, 2.0, 3.0 },
        std::array{ 3.0, 8.0, 1.0, 3.0, 10.0, 2.0, 3.0 },
        std::array{ 5.0, 8.0, 2.0, 6.0, 10.0, 2.0, 3.0 },
    };
    std::vector<std::unique_ptr<AtomObject>> atoms;
    std::vector<LocalFittingFeatureRow> rows;
    for (std::size_t index = 0; index < feature_values.size(); ++index)
    {
        const auto serial_id{ static_cast<int>(index + 1) };
        auto atom{ std::make_unique<AtomObject>() };
        atom->SetSerialID(serial_id);
        atom->SetElement(Element::CARBON);
        atoms.push_back(std::move(atom));
        LocalFittingFeatureRow row;
        row.serial_id = serial_id;
        row.residue = "ALA";
        row.spot = "CA";
        row.features = feature_values[index];
        rows.push_back(std::move(row));
    }
    rows[1].tail_peeling_ratio_assessment.reason = "fixture-tail-missing";
    const ModelObject model{ std::move(atoms) };

    std::string error;
    const auto prepared{ BuildAndStandardizeInput(model, std::move(rows), "fixture", 15, error) };
    ASSERT_TRUE(prepared) << error;
    EXPECT_EQ(prepared->effective_neighbors, 2);
    EXPECT_EQ(prepared->excluded_rows, 1u);
    ASSERT_EQ(prepared->exclusions.size(), 1u);
    EXPECT_EQ(prepared->exclusions.at("unavailable: tail peeling ratio (fixture-tail-missing)"), 1u);
    EXPECT_EQ(prepared->constant_features, (std::vector<std::string_view>{ "distance to closest neighbor" }));
    EXPECT_EQ(GetUmapFeatureNames(), (std::vector<std::string_view>{
        "neighbor distance sum in 2A", "distance to closest neighbor", "tail peeling ratio" }));

    ASSERT_EQ(prepared->feature_rows.size(), 3u);
    EXPECT_EQ(prepared->feature_rows[0].serial_id, 1);
    EXPECT_EQ(prepared->feature_rows[1].serial_id, 3);
    EXPECT_EQ(prepared->feature_rows[2].serial_id, 4);
    EXPECT_TRUE(std::isnan(prepared->feature_rows.front().features[2]));

    const auto z{ std::sqrt(1.5) };
    const std::array expected{ -z, 0.0, z, 0.0, 0.0, -z, z, 0.0, 0.0 };
    EXPECT_EQ(prepared->features.column_count, 3u);
    ASSERT_EQ(prepared->features.values.size(), expected.size());
    for (std::size_t index = 0; index < expected.size(); ++index)
    {
        EXPECT_NEAR(prepared->features.values[index], expected[index], 1e-15);
    }
}

TEST(UmapEmbeddingFeaturesTest, ImputesSecondFeaturesBeforePopulationZScore)
{
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    const std::array signal{ 0.0, 2.0, 6.0, 8.0, missing };
    const std::array offsets{ 1.0, 3.0, 1.0, 3.0, missing };
    std::vector<LocalFittingFeatureRow> rows(signal.size());
    for (std::size_t point = 0; point < rows.size(); ++point)
    {
        const double index{ static_cast<double>(point) };
        rows[point].features = { 2.0 * index, 100.0, signal[point], 1.0, missing,
            9.0 - 2.0 * index, offsets[point] };
    }
    EXPECT_EQ(GetSecondUmapFeatureNames(), (std::vector<std::string_view>{
        "first bin center x (diagonal-normalized)", "first bin center y (diagonal-normalized)",
        "neighbor distance sum in 2A", "signal peeling ratio",
        "width joint", "offset joint" }));

    UmapBinning bins;
    bins.scale = { {}, 10.0, 1.0 };
    bins.bins[{ 1, 2 }] = { { 0, 1, 2, 3, 4 }, 1 };
    bins.point_bin_numbers = { 1, 1, 1, 1, 1 };
    const auto features{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(features.column_count, 6u);
    ASSERT_EQ(features.values.size(), rows.size() * features.column_count);
    const std::array filled_signal{ 0.0, 2.0, 6.0, 8.0, 4.0 };
    const std::array filled_offsets{ 1.0, 3.0, 1.0, 3.0, 2.0 };
    for (std::size_t point = 0; point < rows.size(); ++point)
    {
        const auto offset{ point * features.column_count };
        const double centered{ 2.0 * static_cast<double>(point) - 4.0 };
        EXPECT_DOUBLE_EQ(features.values[offset], 0.075);
        EXPECT_DOUBLE_EQ(features.values[offset + 1], 0.125);
        EXPECT_NEAR(features.values[offset + 2], centered / std::sqrt(8.0), 1e-15);
        EXPECT_NEAR(features.values[offset + 3], (filled_signal[point] - 4.0) / std::sqrt(8.0), 1e-15);
        EXPECT_NEAR(features.values[offset + 4], -centered / std::sqrt(8.0), 1e-15);
        EXPECT_NEAR(features.values[offset + 5], (filled_offsets[point] - 2.0) / std::sqrt(0.8), 1e-15);
    }
}

TEST(UmapEmbeddingFeaturesTest, PreservesDistancesBetweenGeometricCentersIndependentlyOfBinNumbers)
{
    const std::vector<LocalFittingFeatureRow> rows(6);
    UmapBinning bins;
    bins.scale = { { -2.0, 1.0 }, 10.0, 1.0 };
    bins.bins = {
        { { 0, 0 }, { { 0, 1 }, 1 } },
        { { 1, 0 }, { { 2 }, 3 } },
        { { 0, 1 }, { { 3 }, 2 } },
        { { 1, 1 }, { { 4 }, 4 } },
        { { 7, 2 }, { { 5 }, 5 } },
    };
    bins.point_bin_numbers = { 1, 1, 3, 2, 4, 5 };
    const auto features{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(features.column_count, 6u);
    ASSERT_EQ(features.values.size(), rows.size() * 6);

    // Centers of the 0.5-wide squares with origin (-2, 1).
    const std::array centers{
        std::array{ -1.75, 1.25 }, std::array{ -1.75, 1.25 }, std::array{ -1.25, 1.25 },
        std::array{ -1.75, 1.75 }, std::array{ -1.25, 1.75 }, std::array{ 1.75, 2.25 },
    };
    for (std::size_t a = 0; a < rows.size(); ++a)
    {
        EXPECT_NEAR(features.values[6 * a], (centers[a][0] + 2.0) / 10.0, 1e-15);
        EXPECT_NEAR(features.values[6 * a + 1], (centers[a][1] - 1.0) / 10.0, 1e-15);
        for (std::size_t b = 0; b < a; ++b)
        {
            const double expected{
                std::hypot(centers[a][0] - centers[b][0], centers[a][1] - centers[b][1]) / 10.0
            };
            const double actual{ std::hypot(features.values[6 * a] - features.values[6 * b],
                features.values[6 * a + 1] - features.values[6 * b + 1]) };
            EXPECT_NEAR(actual, expected, 1e-15);
        }
    }
    for (auto & [index, bin] : bins.bins)
    {
        bin.number = 6 - bin.number;
        for (const auto point : bin.points) bins.point_bin_numbers[point] = bin.number;
    }
    EXPECT_EQ(BuildSecondUmapFeatures(rows, bins).values, features.values);
}

TEST(UmapEmbeddingFeaturesTest, SkipsSecondEmbeddingWhenAllSelectedColumnsAreConstant)
{
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    std::vector<LocalFittingFeatureRow> rows(3);
    for (std::size_t point = 0; point < rows.size(); ++point)
    {
        const auto unselected{ static_cast<double>(point) };
        rows[point].features = { 2.0, unselected, missing, unselected, 7.0, 9.0, 3.0 };
    }
    UmapBinning bins;
    bins.bins[{ 0, 0 }] = { { 0, 1, 2 }, 1 };
    bins.point_bin_numbers = { 1, 1, 1 };
    EXPECT_TRUE(BuildSecondUmapFeatures(rows, bins).values.empty());
    bins.scale = { {}, 10.0, 1.0 };
    EXPECT_TRUE(BuildSecondUmapFeatures(rows, bins).values.empty());

    rows[1].features[0] = 3.0;
    bins.scale = {};
    const auto features{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(features.values.size(), 18u);
    EXPECT_DOUBLE_EQ(features.values[6], 0.0);
    EXPECT_DOUBLE_EQ(features.values[7], 0.0);
    EXPECT_NEAR(features.values[8], std::sqrt(2.0), 1e-15);
    EXPECT_DOUBLE_EQ(features.values[9], 0.0);
}

} // namespace rhbm_gem::core::classifier
