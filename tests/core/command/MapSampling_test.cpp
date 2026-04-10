#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <tuple>
#include <vector>

#include "command/MapSampling.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rg = rhbm_gem;

namespace {

class SinglePointSampler
{
public:
    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & position) const
    {
        return { std::make_tuple(0.0f, position) };
    }
};

class ShiftedPointSampler
{
public:
    SamplingPointList GenerateSamplingPoints(
        const std::array<float, 3> & position,
        const std::array<float, 3> & direction_like_input) const
    {
        return { std::make_tuple(
            1.0f,
            std::array<float, 3>{
                position[0] + direction_like_input[0],
                position[1] + direction_like_input[1],
                position[2] + direction_like_input[2]
            }) };
    }
};

rg::MapObject MakeMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        values[i] = static_cast<float>(i + 1);
    }
    return rg::MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

} // namespace

TEST(MapSamplingTest, PositionOnlySamplerReturnsExpectedPointValueAndIsDeterministic)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    const auto first{
        rg::SampleMapValues(map, sampler, { 0.0f, 0.0f, 0.0f }) };
    const auto second{
        rg::SampleMapValues(map, sampler, { 1.0f, 1.0f, 1.0f }) };
    const auto first_again{
        rg::SampleMapValues(map, sampler, { 0.0f, 0.0f, 0.0f }) };

    ASSERT_EQ(first.size(), 1u);
    ASSERT_EQ(second.size(), 1u);
    ASSERT_EQ(first_again.size(), 1u);
    EXPECT_FLOAT_EQ(std::get<0>(first.front()), 0.0f);
    EXPECT_FLOAT_EQ(std::get<0>(second.front()), 0.0f);
    EXPECT_FLOAT_EQ(std::get<1>(first.front()), map.GetMapValue(0, 0, 0));
    EXPECT_FLOAT_EQ(std::get<1>(second.front()), map.GetMapValue(1, 1, 1));
    EXPECT_FLOAT_EQ(std::get<1>(first_again.front()), std::get<1>(first.front()));
}

TEST(MapSamplingTest, OrientedSamplerUsesDirectionLikeInput)
{
    auto map{ MakeMapObject() };
    ShiftedPointSampler sampler;

    const auto sampling_data{
        rg::SampleMapValues(
            map,
            sampler,
            { 0.0f, 0.0f, 0.0f },
            { 1.0f, 1.0f, 1.0f })
    };

    ASSERT_EQ(1u, sampling_data.size());
    EXPECT_FLOAT_EQ(1.0f, std::get<0>(sampling_data.front()));
    EXPECT_FLOAT_EQ(map.GetMapValue(1, 1, 1), std::get<1>(sampling_data.front()));
}
