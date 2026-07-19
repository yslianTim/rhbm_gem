#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>

namespace {

rhbm_gem::MapObject MakeMapObject(const std::vector<float> & values)
{
    const std::array<int, 3> grid_size{
        static_cast<int>(values.size()),
        1,
        1
    };
    const std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    const std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto map_values{ std::make_unique<float[]>(values.size()) };
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        map_values[i] = values[i];
    }
    return rhbm_gem::MapObject{
        grid_size,
        grid_spacing,
        origin,
        std::move(map_values)
    };
}

} // namespace

TEST(QScoreHelperTest, ReturnsAmplitudeAndOffsetForTypicalMap)
{
    const auto map{ MakeMapObject({ 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f }) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    constexpr float expected_standard_deviation{ 2.44948974f };
    constexpr float expected_offset{ 4.5f - expected_standard_deviation };
    constexpr float expected_height{ 8.0f - expected_offset };
    EXPECT_NEAR(height, expected_height, 1.0e-5f);
    EXPECT_NEAR(offset, expected_offset, 1.0e-5f);
}

TEST(QScoreHelperTest, CapsReferenceHighAtMeanPlusTenStandardDeviations)
{
    std::vector<float> values(121, 0.0f);
    values.back() = 1.0f;
    const auto map{ MakeMapObject(values) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    constexpr float expected_height{ 1.0f / 121.0f + 10.0f / 11.0f };
    EXPECT_NEAR(height, expected_height, 1.0e-5f);
    EXPECT_FLOAT_EQ(offset, 0.0f);
    EXPECT_LT(height + offset, map.GetMapValueMax());
}

TEST(QScoreHelperTest, ClampsOffsetToObservedMinimum)
{
    const auto map{ MakeMapObject({ 0.0f, 0.0f, 0.0f, 10.0f }) };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    EXPECT_FLOAT_EQ(height, 10.0f);
    EXPECT_FLOAT_EQ(offset, 0.0f);
}

TEST(QScoreHelperTest, ReturnsZeroAmplitudeForConstantMap)
{
    const auto map{
        MakeMapObject({ 3.5f, 3.5f, 3.5f, 3.5f, 3.5f, 3.5f, 3.5f, 3.5f })
    };

    const auto [height, offset]{
        rhbm_gem::core::GetReferenceGaussianParameters(map)
    };

    EXPECT_FLOAT_EQ(height, 0.0f);
    EXPECT_FLOAT_EQ(offset, 3.5f);
}
