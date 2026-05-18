#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <vector>

#include "command/detail/MapSampling.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

using namespace rhbm_gem;

namespace {

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

MapObject MakeMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        values[i] = static_cast<float>(i + 1);
    }
    return MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

std::unique_ptr<ModelObject> MakeLinearNeighborModel()
{
    std::vector<std::unique_ptr<AtomObject>> atom_list;
    atom_list.reserve(2);

    auto center_atom{ std::make_unique<AtomObject>() };
    center_atom->SetSerialID(1);
    center_atom->SetPosition(0.0f, 0.0f, 0.0f);

    auto neighbor_atom{ std::make_unique<AtomObject>() };
    neighbor_atom->SetSerialID(2);
    neighbor_atom->SetPosition(1.0f, 0.0f, 0.0f);

    atom_list.emplace_back(std::move(center_atom));
    atom_list.emplace_back(std::move(neighbor_atom));
    return std::make_unique<ModelObject>(std::move(atom_list));
}

} // namespace

TEST(MapSamplingTest, OrientedSamplerUsesDirectionLikeInput)
{
    auto map{ MakeMapObject() };
    ShiftedPointSampler sampler;

    const auto sampling_data{
        SampleMapValues(
            map,
            sampler,
            { 0.0f, 0.0f, 0.0f },
            { 1.0f, 1.0f, 1.0f })
    };

    ASSERT_EQ(1u, sampling_data.size());
    EXPECT_FLOAT_EQ(1.0f, sampling_data.front().point.distance);
    EXPECT_FLOAT_EQ(map.GetMapValue(1, 1, 1), sampling_data.front().response);
    EXPECT_FLOAT_EQ(1.0f, sampling_data.front().point.position.at(0));
}

TEST(MapSamplingTest, CommandAtomSamplerBuildsSphereSamplerFromProfileChoice)
{
    auto map{ MakeMapObject() };
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    const auto sampling_data{
        SampleAtomMapValues(
            map,
            *atom,
            SphereSamplingProfileChoice::FIBONACCI_DETERMINISTIC)
    };

    EXPECT_FALSE(sampling_data.empty());
}

TEST(MapSamplingTest, AtomSamplerRequiresAttachedAtomBeforeSampling)
{
    auto map{ MakeMapObject() };
    AtomObject detached_atom;
    detached_atom.SetPosition(0.0f, 0.0f, 0.0f);

    EXPECT_THROW(
        (void)SampleAtomMapValues(
            map,
            detached_atom,
            SphereSamplingProfileChoice::FIBONACCI_DETERMINISTIC),
        std::runtime_error);
}
