#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include "core/detail/gaussian_fit/GaussianModelOperations.hpp"
#include "core/detail/second_stage/CouplingGraph.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/SecondStageState.hpp"
#include "core/detail/second_stage/observation/SecondStageLogging.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
using rhbm_gem::FittingStage;

using second_stage_test::BuildSamples;
using second_stage_test::Distance;
using second_stage_test::ExpectGaussianModelsNear;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeAtom;
using second_stage_test::MakeGaussianResult;
using second_stage_test::MakeSecondStageOptions;

TEST(EstimatorSecondStageDefenseTest, HydrogenNeverCountsAsSecondStageNeighbor)
{
    auto model{ second_stage_test::BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 0.8, 0.0, 0.0 },
            std::array<double, 3>{ 1.2, 0.0, 0.0 }
        },
        { Spot::CA, Spot::H, Spot::CB },
        { Element::CARBON, Element::HYDROGEN, Element::CARBON },
        {
            rg::GaussianModel3D{ 6.0, 0.48, 0.05 },
            rg::GaussianModel3D{ 2.0, 0.42, 0.01 },
            rg::GaussianModel3D{ 4.8, 0.62, -0.08 }
        },
        rg::GaussianModel3D{ 5.3, 0.57, 0.02 }) };
    model->SelectAtoms([](const rg::AtomObject & atom) { return atom.GetSerialID() == 1; });
    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);

    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    EXPECT_EQ(rg::AtomLocalPotentialView::For(*model->FindAtomPtr(1)).GetNeighborCountForPeeling(), 1);
}

void AddCouplingGraphSample(
    detail::CouplingGraphBuilder & builder,
    detail::SampleRef sample_id,
    std::vector<detail::GraphParticipant> participant_list)
{
    builder.AddSample(sample_id, participant_list);
}

bool HasCouplingNeighbor(
    const detail::GraphTopology & topology,
    std::size_t atom_index,
    std::size_t neighbor_index)
{
    const auto & neighbor_index_list{ topology.adjacency_list.at(atom_index) };
    return std::find(
        neighbor_index_list.begin(),
        neighbor_index_list.end(),
        neighbor_index) != neighbor_index_list.end();
}


} // namespace


TEST(EstimatorSecondStageDefenseTest, CouplingGraphAdaptiveHysteresisAddsAndRemovesEdges)
{
    const auto build_topology = [](
        double edge_weight,
        const detail::GraphTopology * previous_topology)
    {
        const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
        const auto self_scale{ std::sqrt(1.0 / edge_weight - 1.0) };
        detail::CouplingGraphBuilder builder{ 2 };
        AddCouplingGraphSample(
            builder,
            { 0, 0 },
            { { 0, unit }, { 1, unit } });
        AddCouplingGraphSample(
            builder,
            { 0, 1 },
            { { 0, self_scale * unit } });
        AddCouplingGraphSample(
            builder,
            { 1, 0 },
            { { 1, self_scale * unit } });
        detail::CouplingGraphOptions options;
        options.minimum_weight = 0.06;
        options.retained_edge_minimum_weight = 0.04;
        return builder.BuildTopology(
            options,
            previous_topology);
    };

    detail::GraphTopology absent_previous;
    absent_previous.adjacency_list.resize(2);
    const auto absent_midpoint{ build_topology(0.05, &absent_previous) };
    EXPECT_FALSE(HasCouplingNeighbor(absent_midpoint, 0, 1));
    const auto added{ build_topology(0.061, &absent_previous) };
    EXPECT_TRUE(HasCouplingNeighbor(added, 0, 1));
    const auto cutoff_previous{ detail::ApplyGraphAtomCutoff(added, 1) };
    EXPECT_FALSE(HasCouplingNeighbor(build_topology(0.05, &cutoff_previous), 0, 1));
    const auto retained_midpoint{ build_topology(0.05, &added) };
    EXPECT_TRUE(HasCouplingNeighbor(retained_midpoint, 0, 1));
    const auto removed{ build_topology(0.039, &retained_midpoint) };
    EXPECT_FALSE(HasCouplingNeighbor(removed, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphKeepsFormalEdgesAndAtomCutoff)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    detail::CouplingGraphBuilder builder{ 7 };
    const std::array<double, 3> edge_weight_list{ 0.06, 0.12, 0.25 };
    for (std::size_t edge_index = 0; edge_index < edge_weight_list.size(); edge_index++)
    {
        const auto left_index{ 2 * edge_index };
        const auto right_index{ left_index + 1 };
        const auto self_scale{
            std::sqrt(1.0 / edge_weight_list.at(edge_index) - 1.0)
        };
        AddCouplingGraphSample(
            builder,
            { edge_index, 0 },
            { { left_index, unit }, { right_index, unit } });
        AddCouplingGraphSample(
            builder,
            { edge_index, 1 },
            { { left_index, self_scale * unit } });
        AddCouplingGraphSample(
            builder,
            { edge_index, 2 },
            { { right_index, self_scale * unit } });
    }

    detail::CouplingGraphOptions options;
    options.maximum_atom_count = 1;
    const auto topology{
        builder.BuildTopology(options)
    };
    ASSERT_EQ(topology.retained_edge_list.size(), edge_weight_list.size());
    for (std::size_t edge_index = 0;
        edge_index < topology.retained_edge_list.size();
        edge_index++)
    {
        const auto & edge{ topology.retained_edge_list.at(edge_index) };
        EXPECT_EQ(edge.left_atom_index, 2 * edge_index);
        EXPECT_EQ(edge.right_atom_index, 2 * edge_index + 1);
    }
    EXPECT_EQ(topology.summary.retained_edge_count, 3U);
    EXPECT_EQ(topology.summary.cut_edge_count, 0U);
    const auto formal_partition{
        detail::BuildGraphPartition(
            topology,
            { 0, 1, 2, 3, 4, 5, 6 })
    };
    EXPECT_EQ(formal_partition.sample_id_list_by_key.size(), 7U);
    EXPECT_EQ(topology.summary.component_count, 7U);
    EXPECT_EQ(topology.summary.maximum_component_size, 1U);
    EXPECT_EQ(topology.atom_cutoff_summary.cut_edge_count, 3U);
    EXPECT_NEAR(topology.summary.maximum_component_ratio, 1.0 / 7.0, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, CouplingPartitionCutsWeakBridgeAndDuplicatesBoundarySample)
{
    detail::GraphTopology topology;
    topology.adjacency_list.resize(3);
    topology.adjacency_list.at(0).push_back(1);
    topology.adjacency_list.at(1).push_back(0);
    topology.sample_dependency_list = {
        { { 0, 0 }, { 0, 1 } },
        { { 1, 0 }, { 1, 2 } },
        { { 0, 1 }, { 1, 2 } }
    };

    const auto partition{
        detail::BuildGraphPartition(topology, { 0, 1, 2 })
    };
    ASSERT_EQ(partition.sample_id_list_by_key.size(), 2U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 2 }), 1U);
    ASSERT_EQ(partition.boundary_sample_dependency_list.size(), 2U);
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front().sample_id,
        (detail::SampleRef{ 0, 1 }));
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front().cluster_key_list,
        (std::vector<detail::ClusterKey>{ { 0, 1 }, { 2 } }));
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front()
            .contributor_atom_index_list,
        (std::vector<std::size_t>{ 1, 2 }));
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.back().sample_id,
        (detail::SampleRef{ 1, 0 }));
    auto contributor_changed_partition{ partition };
    contributor_changed_partition.boundary_sample_dependency_list.front()
        .contributor_atom_index_list = { 1 };
    EXPECT_NE(
        contributor_changed_partition.boundary_sample_dependency_list,
        partition.boundary_sample_dependency_list);
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 0, 1 }).size(), 3U);
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 2 }).size(), 2U);

    const auto key_list{
        detail::BuildGraphClusterKeyList(partition)
    };
    EXPECT_EQ(key_list, (std::vector<std::vector<std::size_t>>{ { 0, 1 }, { 2 } }));
    const auto affected_sample_list{
        detail::BuildGraphAffectedSampleUnion(
            partition,
            key_list)
    };
    EXPECT_EQ(affected_sample_list.size(), 3U);
    EXPECT_EQ(affected_sample_list.at(0).atom_index, 0U);
    EXPECT_EQ(affected_sample_list.at(0).sample_index, 0U);
    EXPECT_EQ(affected_sample_list.at(1).atom_index, 0U);
    EXPECT_EQ(affected_sample_list.at(1).sample_index, 1U);
    EXPECT_EQ(affected_sample_list.at(2).atom_index, 1U);
    EXPECT_EQ(affected_sample_list.at(2).sample_index, 0U);

    const auto inactive_partition{
        detail::BuildGraphPartition(topology, { 2, 0 })
    };
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 0 }), 1U);
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_TRUE(inactive_partition.boundary_sample_dependency_list.empty());
}


TEST(EstimatorSecondStageDefenseTest,
    ResidualBaselineAndOverlayPreserveFrozenBackground)
{
    detail::SecondStageContext context;
    context.atom_list.resize(1);
    context.atom_list.at(0).neighbor_atom_sample_offset_list = { 0, 0, 0 };
    const rg::GaussianModel3D previous_model{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D candidate_model{ 10.0, 0.60, 0.20 };
    constexpr double contributor_distance{ 0.25 };
    for (const auto distance : { 0.15, 0.45 })
    {
        context.atom_list.at(0).raw_sampling_entries.emplace_back(LocalPotentialSample{
            candidate_model.ResponseAtDistance(distance) + previous_model.ResponseAtDistance(contributor_distance),
            SamplingPoint{ distance } });
        context.atom_list.at(0).unselected_distance_list_by_sample.push_back({ contributor_distance });
    }
    const detail::FitState previous_state{ MakeGaussianResult(previous_model) };
    context.frozen_background = detail::BuildFrozenBackground(context, previous_state);
    ASSERT_TRUE(context.frozen_background);
    const auto baseline{ detail::BuildResidualBaseline(context, previous_state) };
    ASSERT_TRUE(baseline.sample_list.at(0).at(1).has_value());
    EXPECT_NEAR(baseline.sample_list.at(0).at(1)->adjusted_response, candidate_model.ResponseAtDistance(0.45), 1.0e-12);
    detail::FitStatePatch patch;
    patch.atom_index_list = { 0 };
    patch.mdpde_list = { MakeGaussianResult(candidate_model).mdpde };
    const detail::CandidateEvaluationOverlay overlay{
        context,
        baseline,
        previous_state,
        patch
    };
    const detail::SampleRef sample_ref{ 0, 1 };
    const auto candidate_snapshot{ detail::BuildSecondStageModelSnapshot(context, overlay.GetState()) };
    const auto direct{ detail::EvaluateResidualSample(context, sample_ref, candidate_snapshot) };
    const auto overlaid{ overlay(sample_ref) };
    ASSERT_TRUE(direct.has_value());
    ASSERT_TRUE(overlaid.has_value());
    EXPECT_NEAR(direct->residual, 0.0, 1.0e-12);
    EXPECT_DOUBLE_EQ(direct->adjusted_response, overlaid->adjusted_response);
    EXPECT_DOUBLE_EQ(direct->residual, overlaid->residual);
    const auto samples{ detail::BuildSecondStageAdjustedSamples(context, 0, candidate_snapshot) };
    EXPECT_DOUBLE_EQ(samples.at(1).response, direct->adjusted_response);
    EXPECT_EQ(overlay.GetState().size(), 1U);

    // A later refresh must not change an already captured snapshot or overlay.
    const detail::FitState next_state{ MakeGaussianResult(candidate_model) };
    context.frozen_background = detail::BuildFrozenBackground(context, next_state);
    ASSERT_TRUE(context.frozen_background);
    EXPECT_NE(context.frozen_background->response_by_atom, candidate_snapshot.frozen_background->response_by_atom);
    const auto old_response{ detail::EvaluateResidualSample(context, sample_ref, candidate_snapshot) };
    ASSERT_TRUE(old_response.has_value());
    EXPECT_DOUBLE_EQ(old_response->residual, direct->residual);
    EXPECT_DOUBLE_EQ(overlay(sample_ref)->residual, direct->residual);
    const auto refreshed{ detail::BuildResidualBaseline(context, next_state) };
    ASSERT_TRUE(refreshed.sample_list.front().back().has_value());
    EXPECT_NE(refreshed.sample_list.front().back()->residual, direct->residual);

}
