#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/math/GridSampler.hpp>

using namespace rhbm_gem;
using namespace rhbm_gem::core;

namespace {

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

std::size_t CountSelectedSamples(const LocalPotentialSampleList & sampling_data)
{
    std::size_t count{ 0 };
    for (const auto & sample : sampling_data)
    {
        if (sample.point.is_selected)
        {
            count++;
        }
    }
    return count;
}

} // namespace

TEST(MapSamplerTest, GridSamplerProducesMapSamples)
{
    auto map{ MakeMapObject() };
    GridSampler sampler;
    sampler.SetGridResolution(2);
    sampler.SetWindowSize(1.0f);
    sampler.SetReferenceUVector({ 1.0f, 0.0f, 0.0f });

    const auto sampling_data{
        SampleMapValues(
            map,
            sampler,
            { 0.5f, 0.5f, 0.5f },
            { 0.0f, 0.0f, 1.0f })
    };

    ASSERT_EQ(4u, sampling_data.size());
    EXPECT_NEAR(0.70710677f, sampling_data.front().point.distance, 1.0e-6f);
    EXPECT_FLOAT_EQ(0.0f, sampling_data.front().point.position.at(0));
    EXPECT_FLOAT_EQ(0.0f, sampling_data.front().point.position.at(1));
    EXPECT_FLOAT_EQ(0.5f, sampling_data.front().point.position.at(2));
}

TEST(MapSamplerTest, AtomSamplerUsesSamplingMethod)
{
    auto map{ MakeMapObject() };
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    const auto sampling_data{
        SampleAtomMapValues(
            map,
            *atom,
            SphereSamplingMethod::FibonacciDeterministic)
    };

    EXPECT_FALSE(sampling_data.empty());
}

TEST(MapSamplerTest, AtomSamplerKeepsRejectedSamplesAsUnselected)
{
    auto map{ MakeMapObject() };
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    const auto sampling_data{
        SampleAtomMapValues(
            map,
            *atom,
            SphereSamplingMethod::FibonacciDeterministic)
    };

    const auto selected_count{ CountSelectedSamples(sampling_data) };
    ASSERT_FALSE(sampling_data.empty());
    EXPECT_GT(selected_count, 0u);
    EXPECT_LT(selected_count, sampling_data.size());
}

TEST(MapSamplerTest, AtomSamplerRequiresAttachedAtomBeforeSampling)
{
    auto map{ MakeMapObject() };
    AtomObject detached_atom;
    detached_atom.SetPosition(0.0f, 0.0f, 0.0f);

    EXPECT_THROW(
        (void)SampleAtomMapValues(
            map,
            detached_atom,
            SphereSamplingMethod::FibonacciDeterministic),
        std::runtime_error);
}
