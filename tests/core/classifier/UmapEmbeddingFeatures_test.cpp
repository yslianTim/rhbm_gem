#include <gtest/gtest.h>

#include "core/classifier/UmapEmbeddingFeatures.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
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
        atom->SetElement(index == 1 ? Element::NITROGEN : Element::CARBON);
        atoms.push_back(std::move(atom));
        LocalFittingFeatureRow row;
        row.serial_id = serial_id;
        row.residue = "ALA";
        row.spot = "CA";
        row.features = feature_values[index];
        rows.push_back(std::move(row));
    }
    rows[1].tail_peeling_ratio_assessment.reason = "fixture-tail-missing";
    auto atom_without_features{ std::make_unique<AtomObject>() };
    atom_without_features->SetSerialID(5);
    atom_without_features->SetElement(Element::OXYGEN);
    atoms.push_back(std::move(atom_without_features));
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
    EXPECT_EQ(prepared->elements, (std::vector<Element>{ Element::CARBON, Element::CARBON, Element::CARBON }));
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

TEST(UmapEmbeddingFeaturesTest, ComputesIndependentFiniteFeatureMediansPerBin)
{
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    constexpr auto inf{ std::numeric_limits<double>::infinity() };
    constexpr auto large{ std::numeric_limits<double>::max() };
    const std::array source_values{
        std::array{ 3.0, 0.0, -3.0, inf, 5.0, missing, -6.0 },
        std::array{ 1.0, 0.0, missing, -2.0, -2.0, 0.0, missing },
        std::array{ 3.0, 0.0, inf, -4.0, -2.0, -inf, inf },
        std::array{ 3.0, 0.0, 6.0, 0.0, 0.0, 5.0, -inf },
        std::array{ 9.0, missing, -6.0, 2.0, inf, 1.0, missing },
        std::array{ missing, -inf, -inf, 4.0, missing, missing, missing },
        std::array{ -large, -2.0, 0.0, missing, large, -4.0, -2.0 },
        std::array{ large, -4.0, -1.0, 4.0, large / 2.0, 0.0, 2.0 },
        std::array{ inf, 7.0, 5.0, 0.0, -inf, 7.0, 2.0 },
    };
    std::vector<LocalFittingFeatureRow> rows(source_values.size());
    for (std::size_t point = 0; point < rows.size(); ++point)
        rows[point].features = source_values[point];
    UmapBinning bins;
    bins.bins = {
        { { 0, 0 }, { { 0 }, 1 } },
        { { 0, 1 }, { { 1, 2, 3, 4, 5 }, 2 } },
        { { 1, 0 }, { { 6, 7 }, 3 } },
        { { 1, 1 }, { { 8 }, 4 } },
    };
    PopulateUmapBinFeatureMedians(rows, bins);

    using Medians = std::array<std::optional<double>, kLocalFittingFeatureCount>;
    EXPECT_EQ(bins.bins.at({ 0, 0 }).feature_medians,
        (Medians{ 3.0, 0.0, -3.0, std::nullopt, 5.0, std::nullopt, -6.0 }));
    EXPECT_EQ(bins.bins.at({ 0, 1 }).feature_medians,
        (Medians{ 3.0, 0.0, 0.0, 0.0, -2.0, 1.0, std::nullopt }));
    EXPECT_EQ(bins.bins.at({ 1, 0 }).feature_medians,
        (Medians{ 0.0, -3.0, -0.5, 4.0, 0.75 * large, -2.0, 0.0 }));
    EXPECT_EQ(bins.bins.at({ 1, 1 }).feature_medians,
        (Medians{ std::nullopt, 7.0, 5.0, 0.0, std::nullopt, 7.0, 2.0 }));

    rows[0].features[0] = missing;
    rows[0].features[3] = -1.0;
    PopulateUmapBinFeatureMedians(rows, bins);
    EXPECT_EQ(bins.bins.at({ 0, 0 }).feature_medians,
        (Medians{ std::nullopt, 0.0, -3.0, -1.0, 5.0, std::nullopt, -6.0 }));
}

TEST(UmapEmbeddingFeaturesTest, ImputesAndStandardizesElevenColumnsInAtomRowOrder)
{
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    constexpr auto inf{ std::numeric_limits<double>::infinity() };
    const std::array source_values{
        std::array{ 0.0, 99.0, 0.0, 11.0, 0.0, 2.0, -4.0 },
        std::array{ 2.0, 98.0, 2.0, 13.0, 2.0, 4.0, -2.0 },
        std::array{ 4.0, 97.0, 6.0, 15.0, 4.0, 6.0, 0.0 },
        std::array{ 6.0, 96.0, missing, 17.0, 8.0, 10.0, 8.0 },
        std::array{ 8.0, 95.0, 8.0, 20.0, missing, 4.0, 2.0 },
        std::array{ 10.0, 94.0, inf, 22.0, inf, 4.0, 4.0 },
    };
    std::vector<LocalFittingFeatureRow> rows(source_values.size());
    for (std::size_t point = 0; point < rows.size(); ++point)
        rows[point].features = source_values[point];

    EXPECT_EQ(GetSecondUmapFeatureNames(), (std::vector<std::string_view>{
        "amplitude joint median of bin", "width joint median of bin", "offset joint median of bin",
        "neighbor distance sum in 2A median of bin", "distance to closest neighbor median of bin",
        "signal peeling ratio median of bin", "tail peeling ratio median of bin",
        "neighbor distance sum in 2A", "signal peeling ratio", "width joint", "offset joint" }));

    UmapBinning bins;
    bins.scale = { { -2.0, 1.0 }, 10.0, 1.0 };
    bins.bins = {
        { { 0, 0 }, { { 0, 1, 2 }, 1 } },
        { { 0, 1 }, { { 4, 5 }, 2 } },
        { { 1, 0 }, { { 3 }, 3 } },
    };
    bins.point_bin_numbers = { 1, 1, 1, 3, 2, 2 };
    PopulateUmapBinFeatureMedians(rows, bins);
    const auto features{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(features.column_count, 11u);
    ASSERT_EQ(features.values.size(), rows.size() * features.column_count);

    // The amplitude fill is 2 because the first bin contributes three rows.
    const std::array imputed_columns{
        std::array{ 2.0, 2.0, 2.0, 8.0, 2.0, 2.0 },
        std::array{ 4.0, 4.0, 4.0, 10.0, 4.0, 4.0 },
        std::array{ -2.0, -2.0, -2.0, 8.0, 3.0, 3.0 },
        std::array{ 2.0, 2.0, 2.0, 6.0, 9.0, 9.0 },
        std::array{ 98.0, 98.0, 98.0, 96.0, 94.5, 94.5 },
        std::array{ 2.0, 2.0, 2.0, 2.0, 8.0, 8.0 },
        std::array{ 13.0, 13.0, 13.0, 17.0, 21.0, 21.0 },
        std::array{ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 },
        std::array{ 0.0, 2.0, 6.0, 4.0, 8.0, 4.0 },
        std::array{ 2.0, 4.0, 6.0, 10.0, 4.0, 4.0 },
        std::array{ -4.0, -2.0, 0.0, 8.0, 2.0, 4.0 },
    };
    const std::array means{ 3.0, 5.0, 4.0 / 3.0, 5.0, 96.5, 4.0, 49.0 / 3.0,
        5.0, 4.0, 5.0, 4.0 / 3.0 };
    const std::array variances{ 5.0, 5.0, 125.0 / 9.0, 10.0, 2.5, 8.0, 116.0 / 9.0,
        35.0 / 3.0, 20.0 / 3.0, 19.0 / 3.0, 140.0 / 9.0 };
    for (std::size_t column = 0; column < features.column_count; ++column)
    {
        for (std::size_t point = 0; point < rows.size(); ++point)
        {
            SCOPED_TRACE(testing::Message() << "column " << column << ", point " << point);
            const auto value{ features.values[point * features.column_count + column] };
            EXPECT_TRUE(std::isfinite(value));
            EXPECT_NEAR(value, (imputed_columns[column][point] - means[column])
                / std::sqrt(variances[column]), 1e-14);
        }
    }

    for (auto & [index, bin] : bins.bins)
    {
        bin.number = 20 - 2 * bin.number;
        for (const auto point : bin.points) bins.point_bin_numbers[point] = bin.number;
    }
    EXPECT_EQ(BuildSecondUmapFeatures(rows, bins).values, features.values);

    UmapBinning relocated;
    relocated.point_bin_numbers = bins.point_bin_numbers;
    for (const auto & [index, bin] : bins.bins)
        relocated.bins[{ 9 - index[0], 3 - index[1] }] = bin;
    EXPECT_EQ(BuildSecondUmapFeatures(rows, relocated).values, features.values);
}

TEST(UmapEmbeddingFeaturesTest, RequiresVariationInAtLeastOneMedianOrAtomColumn)
{
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    std::vector<LocalFittingFeatureRow> rows(3);
    for (auto & row : rows) row.features = { 2.0, 3.0, missing, 4.0, 7.0, 9.0, 0.0 };
    UmapBinning bins;
    bins.bins = {
        { { 0, 0 }, { { 0, 2 }, 1 } },
        { { 1, 0 }, { { 1 }, 2 } },
    };
    PopulateUmapBinFeatureMedians(rows, bins);
    const auto constant{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(constant.column_count, 11u);
    EXPECT_TRUE(constant.values.empty());

    rows[1].features[1] = 8.0;
    PopulateUmapBinFeatureMedians(rows, bins);
    const auto median_variation{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(median_variation.values.size(), 33u);
    EXPECT_NEAR(median_variation.values[11 + 4], std::sqrt(2.0), 1e-14);

    rows[1].features[1] = 3.0;
    rows[0].features[6] = -1.0;
    rows[2].features[6] = 1.0;
    PopulateUmapBinFeatureMedians(rows, bins);
    const auto atom_variation{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(atom_variation.values.size(), 33u);
    EXPECT_NEAR(atom_variation.values[10], -std::sqrt(1.5), 1e-14);

    for (std::size_t index = 0; index < 33; ++index)
    {
        if (index % 11 != 4) EXPECT_DOUBLE_EQ(median_variation.values[index], 0.0);
        if (index % 11 != 10) EXPECT_DOUBLE_EQ(atom_variation.values[index], 0.0);
    }

    for (auto & row : rows) row.features.fill(missing);
    PopulateUmapBinFeatureMedians(rows, bins);
    EXPECT_TRUE(BuildSecondUmapFeatures(rows, bins).values.empty());
}

TEST(UmapEmbeddingFeaturesTest, StandardizesAdjacentFloatingPointValues)
{
    std::vector<LocalFittingFeatureRow> rows(8);
    for (auto & row : rows) row.features[6] = 1.0;
    rows.back().features[6] = std::nextafter(1.0, 2.0);
    UmapBinning bins;
    bins.bins = {
        { { 0, 0 }, { { 0, 1, 2, 3, 4, 5, 6 }, 1 } },
        { { 1, 0 }, { { 7 }, 2 } },
    };
    PopulateUmapBinFeatureMedians(rows, bins);
    const auto features{ BuildSecondUmapFeatures(rows, bins) };
    ASSERT_EQ(features.column_count, 11u);
    ASSERT_EQ(features.values.size(), rows.size() * features.column_count);
    for (std::size_t row = 0; row < rows.size(); ++row)
    {
        const double expected{ row == 7 ? std::sqrt(7.0) : -1.0 / std::sqrt(7.0) };
        EXPECT_NEAR(features.values[row * features.column_count + 2], expected, 1e-14);
        EXPECT_NEAR(features.values[row * features.column_count + 10], expected, 1e-14);
    }
}

} // namespace rhbm_gem::core::classifier
