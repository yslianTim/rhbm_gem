#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#include "command/MapSampling.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rg = rhbm_gem;

namespace {

class SinglePointSampler
{
public:
    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & position) const
    {
        return { SamplingPoint{ 0.0f, position } };
    }
};

class ShiftedPointSampler
{
public:
    SamplingPointList GenerateSamplingPoints(
        const std::array<float, 3> & position,
        const std::array<float, 3> & direction_like_input) const
    {
        return { SamplingPoint{
            1.0f,
            std::array<float, 3>{
                position[0] + direction_like_input[0],
                position[1] + direction_like_input[1],
                position[2] + direction_like_input[2]
            }
        } };
    }
};

class FixedPointSampler
{
public:
    SamplingPointList GenerateSamplingPoints(const std::array<float, 3> & position) const
    {
        return {
            SamplingPoint{ 0.0f, position },
            SamplingPoint{ 1.0f, { position[0] + 1.0f, position[1], position[2] } },
            SamplingPoint{ 1.0f, { position[0], position[1] + 1.0f, position[2] } },
            SamplingPoint{ 1.0f, { position[0] - 1.0f, position[1], position[2] } }
        };
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

std::unique_ptr<rg::ModelObject> MakeLinearNeighborModel()
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(2);

    auto center_atom{ std::make_unique<rg::AtomObject>() };
    center_atom->SetSerialID(1);
    center_atom->SetPosition(0.0f, 0.0f, 0.0f);

    auto neighbor_atom{ std::make_unique<rg::AtomObject>() };
    neighbor_atom->SetSerialID(2);
    neighbor_atom->SetPosition(1.0f, 0.0f, 0.0f);

    atom_list.emplace_back(std::move(center_atom));
    atom_list.emplace_back(std::move(neighbor_atom));
    return std::make_unique<rg::ModelObject>(std::move(atom_list));
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
    EXPECT_FLOAT_EQ(first.front().distance, 0.0f);
    EXPECT_FLOAT_EQ(second.front().distance, 0.0f);
    EXPECT_FLOAT_EQ(first.front().response, map.GetMapValue(0, 0, 0));
    EXPECT_FLOAT_EQ(second.front().response, map.GetMapValue(1, 1, 1));
    EXPECT_TRUE(first.front().HasPosition());
    EXPECT_FLOAT_EQ(first.front().position->at(0), 0.0f);
    EXPECT_FLOAT_EQ(first_again.front().response, first.front().response);
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
    EXPECT_FLOAT_EQ(1.0f, sampling_data.front().distance);
    EXPECT_FLOAT_EQ(map.GetMapValue(1, 1, 1), sampling_data.front().response);
    ASSERT_TRUE(sampling_data.front().HasPosition());
    EXPECT_FLOAT_EQ(1.0f, sampling_data.front().position->at(0));
}

TEST(MapSamplingTest, AtomSamplerWithZeroAngleMatchesPositionOnlySampling)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    const auto from_position{
        rg::SampleMapValues(map, sampler, atom->GetPosition()) };
    const auto from_atom{
        rg::SampleMapValues(map, sampler, *atom, 1.1, 0.0) };

    ASSERT_EQ(from_position.size(), from_atom.size());
    ASSERT_EQ(from_atom.size(), 1u);
    EXPECT_FLOAT_EQ(from_position.front().distance, from_atom.front().distance);
    EXPECT_FLOAT_EQ(from_position.front().response, from_atom.front().response);
    ASSERT_TRUE(from_atom.front().HasPosition());
    EXPECT_FLOAT_EQ(from_atom.front().position->at(0), 0.0f);
}

TEST(MapSamplingTest, AtomSamplerRejectsPointsWithinAngleThresholdOfNeighborDirections)
{
    auto map{ MakeMapObject() };
    FixedPointSampler sampler;
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    const auto sampling_data{
        rg::SampleMapValues(map, sampler, *atom, 1.1, 30.0) };

    ASSERT_EQ(sampling_data.size(), 3u);
    EXPECT_FLOAT_EQ(sampling_data.at(0).distance, 0.0f);
    ASSERT_TRUE(sampling_data.at(1).HasPosition());
    EXPECT_FLOAT_EQ(sampling_data.at(1).position->at(0), 0.0f);
    EXPECT_FLOAT_EQ(sampling_data.at(1).position->at(1), 1.0f);
    ASSERT_TRUE(sampling_data.at(2).HasPosition());
    EXPECT_FLOAT_EQ(sampling_data.at(2).position->at(0), -1.0f);
    EXPECT_FLOAT_EQ(sampling_data.at(2).position->at(1), 0.0f);
}

TEST(MapSamplingTest, AtomSamplerKeepsAllPointsWhenNeighborRadiusFindsNoNeighbors)
{
    auto map{ MakeMapObject() };
    FixedPointSampler sampler;
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    const auto sampling_data{
        rg::SampleMapValues(map, sampler, *atom, 0.5, 30.0) };

    ASSERT_EQ(sampling_data.size(), 4u);
}

TEST(MapSamplingTest, AtomSamplerAllowsDetachedAtomWhenAngleIsZeroButRejectsAngleFiltering)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    rg::AtomObject detached_atom;
    detached_atom.SetPosition(0.0f, 0.0f, 0.0f);

    const auto sampling_data{
        rg::SampleMapValues(map, sampler, detached_atom, 1.0, 0.0) };

    ASSERT_EQ(sampling_data.size(), 1u);
    EXPECT_THROW(
        (void)rg::SampleMapValues(map, sampler, detached_atom, 1.0, 30.0),
        std::runtime_error);
}

TEST(MapSamplingTest, AtomSamplerRejectsInvalidFilterParameters)
{
    auto map{ MakeMapObject() };
    SinglePointSampler sampler;
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    EXPECT_THROW(
        (void)rg::SampleMapValues(map, sampler, *atom, -0.1, 0.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::SampleMapValues(
            map,
            sampler,
            *atom,
            1.0,
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::SampleMapValues(map, sampler, *atom, 1.0, -1.0),
        std::invalid_argument);
    EXPECT_THROW(
        (void)rg::SampleMapValues(map, sampler, *atom, 1.0, 181.0),
        std::invalid_argument);
}
