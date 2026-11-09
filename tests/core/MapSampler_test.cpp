#include <gtest/gtest.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <rhbm_gem/core/MapSampler.hpp>
#include "core/detail/MapInterpolation.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/math/GridSampler.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

using namespace rhbm_gem;
using namespace rhbm_gem::core;

namespace {

MapObject MakeMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<double, 3> grid_spacing{ 1.0, 1.0, 1.0 };
    std::array<double, 3> origin{ 0.0, 0.0, 0.0 };
    auto values{ std::make_unique<double[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        values[i] = static_cast<double>(i + 1);
    }
    return MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

MapObject MakeWorkflowMapObject()
{
    std::array<int, 3> grid_size{ 8, 8, 8 };
    std::array<double, 3> grid_spacing{ 1.0, 1.0, 1.0 };
    std::array<double, 3> origin{ -3.0, -3.0, -3.0 };
    const size_t voxel_size{
        static_cast<size_t>(grid_size.at(0) * grid_size.at(1) * grid_size.at(2))
    };
    auto values{ std::make_unique<double[]>(voxel_size) };
    for (size_t i = 0; i < voxel_size; ++i)
    {
        values[i] = 1.0;
    }
    return MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

std::unique_ptr<ModelObject> MakeLinearNeighborModel()
{
    std::vector<std::unique_ptr<AtomObject>> atom_list;
    atom_list.reserve(2);

    auto center_atom{ std::make_unique<AtomObject>() };
    center_atom->SetSerialID(1);
    center_atom->SetElement(Element::CARBON);
    center_atom->SetPosition(0.0, 0.0, 0.0);

    auto neighbor_atom{ std::make_unique<AtomObject>() };
    neighbor_atom->SetSerialID(2);
    neighbor_atom->SetElement(Element::CARBON);
    neighbor_atom->SetPosition(1.0, 0.0, 0.0);

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
    sampler.SetWindowSize(1.0);
    sampler.SetReferenceUVector({ 1.0, 0.0, 0.0 });

    const auto sampling_data{
        SampleMapValues(
            map,
            sampler,
            { 0.5, 0.5, 0.5 },
            { 0.0, 0.0, 1.0 })
    };

    ASSERT_EQ(4u, sampling_data.size());
    EXPECT_NEAR(0.70710677, sampling_data.front().point.distance, 1.0e-6);
    EXPECT_DOUBLE_EQ(0.0, sampling_data.front().point.position.at(0));
    EXPECT_DOUBLE_EQ(0.0, sampling_data.front().point.position.at(1));
    EXPECT_DOUBLE_EQ(0.5, sampling_data.front().point.position.at(2));
}

TEST(MapSamplerTest, SpiralAtomSamplesUseFirstStageMaskAndTricubicResponses)
{
    auto map{ MakeMapObject() };
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };
    constexpr auto method{ SphereSamplingMethod::SpiralDeterministic };
    const auto points{ sphere_sampler::GenerateSamplingPointList(atom->GetPositionRef(), method) };
    const auto samples{ SampleAtomMapValues(map, *atom, method) };

    ASSERT_EQ(samples.size(), 200u);
    ASSERT_EQ(samples.size(), points.size());
    bool saw_cone_only_rejection{ false };
    bool saw_ownership_only_rejection{ false };
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        const auto & point{ points[i] };
        const double x{ point.position[0] };
        // The sole neighbor is at (1, 0, 0), so ownership changes at x = 0.5.
        const bool owned_by_neighbor{ x > 0.5 };
        const bool inside_cone{ x / point.distance > std::sqrt(3.0) / 2.0 };
        saw_cone_only_rejection |= inside_cone && !owned_by_neighbor;
        saw_ownership_only_rejection |= owned_by_neighbor && !inside_cone;
        EXPECT_EQ(samples[i].point.position, point.position);
        EXPECT_EQ(samples[i].point.distance, point.distance);
        EXPECT_EQ(samples[i].point.is_selected, !owned_by_neighbor && !inside_cone);
        const auto stencil{ detail::MakeTricubicStencil(map, point.position) };
        const auto response{ detail::InterpolateTricubic(stencil, [&](const auto & node)
        {
            return map.GetMapValue(node[0], node[1], node[2]);
        }) };
        EXPECT_EQ(samples[i].response, response);
    }
    EXPECT_TRUE(saw_cone_only_rejection);
    EXPECT_TRUE(saw_ownership_only_rejection);
}

TEST(MapSamplerTest, AtomSamplerDefaultMatchesExplicitFibonacci)
{
    auto map{ MakeMapObject() };
    auto model{ MakeLinearNeighborModel() };
    const auto * atom{ model->GetAtomList().at(0).get() };

    const auto default_samples{ SampleAtomMapValues(map, *atom) };
    const auto explicit_samples{ SampleAtomMapValues(
        map, *atom, SphereSamplingMethod::FibonacciDeterministic) };

    ASSERT_EQ(default_samples.size(), explicit_samples.size());
    for (std::size_t i = 0; i < default_samples.size(); ++i)
    {
        EXPECT_EQ(default_samples[i].response, explicit_samples[i].response);
        EXPECT_EQ(default_samples[i].point.distance, explicit_samples[i].point.distance);
        EXPECT_EQ(default_samples[i].point.position, explicit_samples[i].point.position);
        EXPECT_EQ(default_samples[i].point.is_selected, explicit_samples[i].point.is_selected);
    }
}

TEST(MapSamplerTest, HydrogenNeighborDoesNotChangeSamplingMask)
{
    auto map{ MakeMapObject() };
    std::vector<std::unique_ptr<AtomObject>> baseline_atoms;
    auto baseline_center{ std::make_unique<AtomObject>() };
    baseline_center->SetSerialID(1);
    baseline_center->SetElement(Element::CARBON);
    baseline_center->SetPosition(0.0, 0.0, 0.0);
    baseline_atoms.emplace_back(std::move(baseline_center));
    ModelObject baseline_model{ std::move(baseline_atoms) };

    std::vector<std::unique_ptr<AtomObject>> hydrogen_atoms;
    auto hydrogen_center{ std::make_unique<AtomObject>() };
    hydrogen_center->SetSerialID(1);
    hydrogen_center->SetElement(Element::CARBON);
    hydrogen_center->SetPosition(0.0, 0.0, 0.0);
    hydrogen_atoms.emplace_back(std::move(hydrogen_center));
    auto hydrogen_neighbor{ std::make_unique<AtomObject>() };
    hydrogen_neighbor->SetSerialID(2);
    hydrogen_neighbor->SetElement(Element::HYDROGEN);
    hydrogen_neighbor->SetPosition(1.0, 0.0, 0.0);
    hydrogen_atoms.emplace_back(std::move(hydrogen_neighbor));
    ModelObject hydrogen_model{ std::move(hydrogen_atoms) };

    const auto baseline{ SampleAtomMapValues(map, *baseline_model.FindAtomPtr(1)) };
    const auto with_hydrogen{ SampleAtomMapValues(map, *hydrogen_model.FindAtomPtr(1)) };

    ASSERT_EQ(baseline.size(), with_hydrogen.size());
    for (std::size_t i = 0; i < baseline.size(); ++i)
        EXPECT_EQ(baseline[i].point.is_selected, with_hydrogen[i].point.is_selected);
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

TEST(MapSamplerTest, PotentialSamplingWorkflowWritesSamplesToSelectedAtoms)
{
    auto map{ MakeWorkflowMapObject() };
    auto model{ MakeLinearNeighborModel() };
    model->SelectAllAtoms();

    RunPotentialSamplingWorkflow(
        map,
        *model,
        2);

    ASSERT_EQ(2u, model->GetSelectedAtomCount());
    for (const auto * atom : model->GetSelectedAtoms())
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        EXPECT_FALSE(local_view.GetRawSamplingEntries(false).empty());
    }
}

TEST(MapSamplerTest, PotentialSamplingWorkflowUsesExplicitAtomDomain)
{
    auto map{ MakeWorkflowMapObject() };
    auto model{ MakeLinearNeighborModel() };
    model->SelectAtoms([](const AtomObject & atom) { return atom.GetSerialID() == 1; });
    const auto selected{ model->GetSelectedAtoms() };
    const std::vector<AtomObject *> atoms{ model->FindAtomPtr(2) };

    RunPotentialSamplingWorkflow(
        map,
        *model,
        atoms,
        1);

    EXPECT_EQ(model->GetSelectedAtoms(), selected);
    EXPECT_FALSE(AtomLocalPotentialView::For(*model->FindAtomPtr(1)).IsAvailable());
    const auto sampled_view{ AtomLocalPotentialView::For(*model->FindAtomPtr(2)) };
    ASSERT_TRUE(sampled_view.IsAvailable());
    EXPECT_FALSE(sampled_view.GetRawSamplingEntries(false).empty());
}

TEST(MapSamplerTest, DefaultSamplingWorkflowPreservesAtomOutputOrderAcrossWorkerCounts)
{
    auto map{ MakeWorkflowMapObject() };
    auto serial_model{ MakeLinearNeighborModel() };
    auto parallel_model{ MakeLinearNeighborModel() };
    auto four_thread_model{ MakeLinearNeighborModel() };
    serial_model->SelectAllAtoms();
    parallel_model->SelectAllAtoms();
    four_thread_model->SelectAllAtoms();
    const std::vector<AtomObject *> serial_atoms{
        serial_model->FindAtomPtr(2), serial_model->FindAtomPtr(1)
    };
    const std::vector<AtomObject *> parallel_atoms{
        parallel_model->FindAtomPtr(2), parallel_model->FindAtomPtr(1)
    };
    const std::vector<AtomObject *> four_thread_atoms{
        four_thread_model->FindAtomPtr(2), four_thread_model->FindAtomPtr(1)
    };

    RunPotentialSamplingWorkflow(map, *serial_model, serial_atoms, 1);
    RunPotentialSamplingWorkflow(map, *parallel_model, parallel_atoms, 2);
    RunPotentialSamplingWorkflow(map, *four_thread_model, four_thread_atoms, 4);

    for (const auto serial_id : { 1, 2 })
    {
        const auto serial_samples = AtomLocalPotentialView::For(
            *serial_model->FindAtomPtr(serial_id)).GetRawSamplingEntries(false);
        for (const auto * model : { parallel_model.get(), four_thread_model.get() })
        {
            const auto parallel_samples = AtomLocalPotentialView::For(
                *model->FindAtomPtr(serial_id)).GetRawSamplingEntries(false);
            ASSERT_EQ(serial_samples.size(), parallel_samples.size());
            for (std::size_t i = 0; i < serial_samples.size(); ++i)
            {
                EXPECT_DOUBLE_EQ(serial_samples[i].response, parallel_samples[i].response);
                EXPECT_EQ(serial_samples[i].point.position, parallel_samples[i].point.position);
                EXPECT_DOUBLE_EQ(serial_samples[i].point.distance, parallel_samples[i].point.distance);
                EXPECT_EQ(serial_samples[i].point.is_selected, parallel_samples[i].point.is_selected);
            }
        }
    }
}

TEST(MapSamplerTest, PotentialSamplingProgressLabelAndQuietMode)
{
    auto map{ MakeWorkflowMapObject() };
    auto model{ MakeLinearNeighborModel() };
    model->SelectAllAtoms();
    const auto atoms{ model->GetSelectedAtoms() };
    const auto previous_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    RunPotentialSamplingWorkflow(map, *model, atoms,
        1, " Joint sampling");
    const auto labelled{ testing::internal::GetCapturedStdout() };
    EXPECT_NE(labelled.find("100% Joint sampling"), std::string::npos);

    Logger::SetLogLevel(LogLevel::Debug);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    RunPotentialSamplingWorkflow(map, *model, atoms,
        1, " Joint sampling", true);
    const auto quiet{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_level);
    EXPECT_TRUE(quiet.empty());
}

TEST(MapSamplerTest, AtomSamplerRequiresAttachedAtomBeforeSampling)
{
    auto map{ MakeMapObject() };
    AtomObject detached_atom;
    detached_atom.SetPosition(0.0, 0.0, 0.0);

    EXPECT_THROW(
        (void)SampleAtomMapValues(
            map,
            detached_atom,
            SphereSamplingMethod::FibonacciDeterministic),
        std::runtime_error);
}

TEST(MapSamplerTest, SharedTricubicStencilPreservesClampAndEffectiveWeights)
{
    const auto map = MakeMapObject();
    for (const auto & position : std::vector<std::array<double, 3>>{{0, 0, 0}, {0.3, 0.6, 0.8}, {1, 1, 1}})
    {
        const auto stencil = detail::MakeTricubicStencil(map, position);
        const double value = detail::InterpolateTricubic(stencil, [&](const auto & p) {
            return map.GetMapValue(p[0], p[1], p[2]);
        });
        double reference = 0, weight_sum = 0;
        std::array<std::array<double,4>,3> weights{};
        for(std::size_t axis=0;axis<3;++axis)
        {
            const double t=stencil.local[axis], t2=t*t, t3=t2*t;
            weights[axis]={-.5*t+t2-.5*t3,1.-2.5*t2+1.5*t3,
                .5*t+2.*t2-1.5*t3,-.5*t2+.5*t3};
        }
        for(int i=-1;i<3;++i) for(int j=-1;j<3;++j) for(int k=-1;k<3;++k)
        {
            const std::array<int,3> node{std::clamp(stencil.base[0]+i,0,stencil.dimensions[0]-1),
                std::clamp(stencil.base[1]+j,0,stencil.dimensions[1]-1),
                std::clamp(stencil.base[2]+k,0,stencil.dimensions[2]-1)};
            EXPECT_EQ(stencil.Node(i,j,k),node);
            const double weight=weights[0][static_cast<std::size_t>(i+1)]*
                weights[1][static_cast<std::size_t>(j+1)]*weights[2][static_cast<std::size_t>(k+1)];
            reference+=weight*map.GetMapValue(node[0],node[1],node[2]);
            weight_sum+=weight;
        }
        EXPECT_NEAR(value, reference, 1e-13);
        EXPECT_NEAR(weight_sum, 1, 1e-14);
    }
}

TEST(MapSamplerTest, TricubicBaselineAtOffGridAndHalfVoxelBoundary)
{
    const auto map = MakeMapObject();
    const SamplingPointList points{{0., {0.3, 0.6, 0.8}}, {0., {-0.5, 0., 0.}}};
    const auto samples = detail::BuildLocalPotentialSampleList(map, points);
    ASSERT_EQ(samples.size(), 2u);
    EXPECT_NEAR(samples[0].response, 5.898, 2e-14);
    EXPECT_DOUBLE_EQ(samples[1].response, 0.9375);
}

TEST(MapSamplerTest, ExplicitMapPointsUseTheProductionInterpolationPath)
{
    MapObject map({4,4,4},{0.1,0.1,0.1},{-0.1,-0.1,-0.1});
    auto values{std::make_unique<double[]>(64)};
    for (std::size_t i=0;i<64;++i) values[i]=static_cast<double>(i);
    map.SetMapValueArray(std::move(values));
    const SamplingPointList points{{0.25,{-0.1,-0.1,-0.1},true},{1.5,{0.0,0.0,0.0},false}};
    const auto samples{detail::BuildLocalPotentialSampleList(map,points)};
    ASSERT_EQ(samples.size(), points.size());
    EXPECT_DOUBLE_EQ(samples[0].response,map.GetMapValue(0,0,0));
    EXPECT_DOUBLE_EQ(samples[1].response,map.GetMapValue(1,1,1));
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(samples[i].point.distance, points[i].distance);
        EXPECT_EQ(samples[i].point.position, points[i].position);
        EXPECT_EQ(samples[i].point.is_selected, points[i].is_selected);
    }
}
