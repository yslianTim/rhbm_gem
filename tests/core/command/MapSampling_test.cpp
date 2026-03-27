#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <tuple>
#include <vector>

#include "internal/command/MapSampling.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/SamplerBase.hpp>

namespace rg = rhbm_gem;

namespace {

class SinglePointSampler : public ::SamplerBase
{
public:
    void Print() const override {}

    std::vector<std::tuple<float, std::array<float, 3>>> GenerateSamplingPoints(
        const std::array<float, 3> & position,
        const std::array<float, 3> & axis_vector) const override
    {
        (void)axis_vector;
        return { std::make_tuple(0.0f, position) };
    }

    unsigned int GetSamplingSize() const override { return m_sampling_size; }
    void SetSamplingSize(unsigned int value) override { m_sampling_size = value; }

private:
    unsigned int m_sampling_size{ 1 };
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

TEST(MapSamplingTest, SampleMapValuesReturnsExpectedPointValueAndIsDeterministic)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    const auto first{
        rg::SampleMapValues(map, sampler, { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }) };
    const auto second{
        rg::SampleMapValues(map, sampler, { 1.0f, 1.0f, 1.0f }, { 0.0f, 1.0f, 0.0f }) };
    const auto first_again{
        rg::SampleMapValues(map, sampler, { 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }) };

    ASSERT_EQ(first.size(), 1u);
    ASSERT_EQ(second.size(), 1u);
    ASSERT_EQ(first_again.size(), 1u);
    EXPECT_FLOAT_EQ(std::get<0>(first.front()), 0.0f);
    EXPECT_FLOAT_EQ(std::get<0>(second.front()), 0.0f);
    EXPECT_FLOAT_EQ(std::get<1>(first.front()), map.GetMapValue(0, 0, 0));
    EXPECT_FLOAT_EQ(std::get<1>(second.front()), map.GetMapValue(1, 1, 1));
    EXPECT_FLOAT_EQ(std::get<1>(first_again.front()), std::get<1>(first.front()));
}
