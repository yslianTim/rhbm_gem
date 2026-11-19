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
    std::vector<LocalFittingFeatureRow> rows(signal.size());
    for (std::size_t point = 0; point < rows.size(); ++point)
    {
        rows[point].features = { 2.0, 100.0, signal[point], 1.0, missing, 9.0, 3.0 };
    }
    EXPECT_EQ(GetSecondUmapFeatureNames(), (std::vector<std::string_view>{
        "bin number", "neighbor distance sum in 2A", "signal peeling ratio",
        "amplitude joint", "width joint", "offset joint" }));

    const std::vector<std::size_t> bins{ 1, 4, 1, 1, 4 };
    const auto features{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(features.column_count, 6u);
    ASSERT_EQ(features.values.size(), rows.size() * features.column_count);
    const std::array filled_signal{ 0.0, 2.0, 6.0, 8.0, 4.0 };
    for (std::size_t point = 0; point < rows.size(); ++point)
    {
        const auto offset{ point * features.column_count };
        EXPECT_NEAR(features.values[offset], bins[point] == 1 ? -std::sqrt(2.0 / 3.0)
            : std::sqrt(3.0 / 2.0), 1e-15);
        EXPECT_DOUBLE_EQ(features.values[offset + 1], 0.0);
        EXPECT_NEAR(features.values[offset + 2], (filled_signal[point] - 4.0) / std::sqrt(8.0), 1e-15);
        EXPECT_DOUBLE_EQ(features.values[offset + 3], 0.0);
        EXPECT_DOUBLE_EQ(features.values[offset + 4], 0.0);
        EXPECT_DOUBLE_EQ(features.values[offset + 5], 0.0);
    }
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
    EXPECT_TRUE(BuildSecondUmapFeatures(rows, { 1, 1, 1 }).values.empty());
}

} // namespace rhbm_gem::core::classifier
