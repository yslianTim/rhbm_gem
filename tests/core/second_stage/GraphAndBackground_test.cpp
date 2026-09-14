#include "support/SecondStageNumericalProbe.hpp"
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

std::unique_ptr<rg::ModelObject> BuildUnselectedContributorDefenseModel(
    const std::array<rg::GaussianModel3D, 2> & selected_seed_model_list,
    const std::array<rg::GaussianModel3D, 7> & truth_model_list,
    bool use_alternate_group_keys = false,
    bool shared_cluster = false,
    bool shared_contributor = false,
    bool use_alternate_residue_keys = false)
{
    std::vector<std::array<double, 3>> position_list{
        { 0.0, 0.0, 0.0 },
        { 7.0, 0.0, 0.0 },
        { 0.60, 0.0, 0.0 },
        { 7.60, 0.0, 0.0 },
        { 0.90, 0.0, 0.0 },
        { 15.0, 0.0, 0.0 },
        { 22.0, 0.0, 0.0 }
    };
    if (shared_cluster)
    {
        // Selected sample responses overlap, independently of residue identity.
        position_list.at(1) = { 0.8, 0.0, 0.0 };
        position_list.at(3) = { 1.4, 0.0, 0.0 };
    }
    if (shared_contributor)
    {
        // Contributor 3 affects both targets, but selected responses do not overlap.
        position_list.at(1) = { 4.0, 0.0, 0.0 };
        position_list.at(2) = { 2.0, 0.0, 0.0 };
        position_list.at(3) = { 4.6, 0.0, 0.0 };
    }
    auto spot_list = std::vector<Spot>{
        Spot::C,
        Spot::C,
        Spot::C,
        Spot::C,
        Spot::C,
        Spot::O,
        Spot::O
    };
    auto element_list = std::vector<Element>{
        Element::CARBON,
        Element::CARBON,
        Element::CARBON,
        Element::CARBON,
        Element::HYDROGEN,
        Element::OXYGEN,
        Element::OXYGEN
    };
    if (use_alternate_group_keys)
    {
        spot_list.at(0) = Spot::N;
        spot_list.at(1) = Spot::O;
        spot_list.at(2) = Spot::N;
        spot_list.at(3) = Spot::O;
    }
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    for (std::size_t i = 0; i < position_list.size(); i++)
    {
        atom_list.emplace_back(MakeAtom(
            static_cast<int>(i + 1),
            spot_list.at(i),
            element_list.at(i),
            position_list.at(i)));
        if (use_alternate_residue_keys)
        {
            atom_list.back()->SetChainID("relabeled");
            atom_list.back()->SetSequenceID(42);
            atom_list.back()->SetResidue(Residue::GLY);
        }
    }
    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();
    for (int serial_id = 3; serial_id <= 7; serial_id++)
    {
        model->SetAtomSelected(serial_id, false);
    }

    std::vector<rg::AtomObject *> all_atoms;
    std::vector<rg::GaussianModel3D> truth_models;
    for (std::size_t atom_index = 0;
        atom_index < model->GetAtomList().size(); atom_index++)
    {
        all_atoms.emplace_back(model->GetAtomList().at(atom_index).get());
        truth_models.emplace_back(truth_model_list.at(atom_index));
    }

    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    for (std::size_t atom_index = 0;
        atom_index < selected_atoms.size();
        atom_index++)
    {
        auto * atom{ selected_atoms.at(atom_index) };
        analysis.SetAtomLocalAlphaR(FittingStage::Second, *atom, 0.0);
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
            MakeGaussianResult(selected_seed_model_list.at(atom_index)));
        analysis.SetAtomLocalRawSamplingEntries(
            *atom, BuildSamples(*atom, all_atoms, truth_models));
    }
    return model;
}

} // namespace

TEST(EstimatorSecondStageDefenseTest, AdaptiveTopologyDriftTracksReferenceState)
{
    const auto make_state = [](double log_width)
    {
        return detail::FitState{ MakeGaussianResult(
            *rg::GaussianModel3D::FromTransformedCoordinates(
                rg::GaussianModel3D::TransformedCoordinates{ 0.0, log_width, 0.0 })) };
    };
    // These widths round-trip to a transformed difference of exactly 0.10.
    const auto reference_state{ make_state(-0.097) };
    const detail::FittedGaussianSnapshot reference{
        reference_state.at(0).mdpde.GetModel() };
    const auto unchanged{
        detail::CalculateAdaptiveTopologyDrift(reference_state, reference, { 0 })
    };
    EXPECT_EQ(unchanged, 0.0);

    // More than three accepted updates stay below the original reference threshold.
    for (const auto drift : { 0.02, 0.04, 0.06, 0.08, 0.099 })
    {
        const auto maximum_transformed_drift{ detail::CalculateAdaptiveTopologyDrift(
            make_state(-0.097 + drift), reference, { 0 }) };
        EXPECT_NEAR(maximum_transformed_drift, drift, 1.0e-12);
    }
    const auto threshold_state{ make_state(-0.097 + 0.10) };
    const auto threshold{ detail::CalculateAdaptiveTopologyDrift(
        threshold_state, reference, { 0 }) };
    ASSERT_EQ(threshold, detail::kAdaptiveTopologyRebuildDriftThreshold);
    const auto above{ detail::CalculateAdaptiveTopologyDrift(
        make_state(-0.097 + 0.101), reference, { 0 }) };
    EXPECT_NEAR(above, 0.101, 1.0e-12);

    const detail::FittedGaussianSnapshot rebuilt_reference{
        threshold_state.at(0).mdpde.GetModel() };
    const auto after_rebuild{ detail::CalculateAdaptiveTopologyDrift(
        make_state(-0.097 + 0.12), rebuilt_reference, { 0 }) };
    EXPECT_NEAR(after_rebuild, 0.02, 1.0e-12);


    detail::SecondStageContext context;
    context.atom_list.resize(2);
    detail::FitState partition_state{
        MakeGaussianResult({ 4.0, 0.4, 1.0 }), MakeGaussianResult({ 8.0, 0.8, 3.0 }) };
    for (auto & atom : context.atom_list)
    {
        atom.raw_sampling_entries.resize(1);
        atom.unselected_distance_list_by_sample = { { 0.3 } };
    }
    const auto background{ detail::BuildFrozenBackground(context, partition_state) };
    ASSERT_TRUE(background);
    ExpectGaussianModelsNear(background->model_by_atom.at(0), { 6.0, 0.6, 2.0 }, 1.0e-12);
    ExpectGaussianModelsNear(background->model_by_atom.at(1), { 6.0, 0.6, 2.0 }, 1.0e-12);
    EXPECT_EQ(background->response_by_atom.at(0), background->response_by_atom.at(1));
    context.frozen_background = background;
    context.atom_list.at(1).unselected_distance_list_by_sample.front().front() = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(detail::BuildFrozenBackground(context, partition_state));
    EXPECT_EQ(context.frozen_background, background);
    context.atom_list.at(1).unselected_distance_list_by_sample.front().front() = 0.3;
    EXPECT_FALSE(detail::BuildFrozenBackground(context, {}));
    const auto empty{ detail::BuildFrozenBackground({}, {}) };
    ASSERT_TRUE(empty);
    EXPECT_TRUE(empty->model_by_atom.empty());
    EXPECT_TRUE(empty->response_by_atom.empty());
    partition_state.front().mdpde = MakeGaussianResult({ -1.0, 0.4, 1.0 }).mdpde;
    EXPECT_FALSE(detail::BuildFrozenBackground(context, partition_state));
    ExpectGaussianModelsNear(background->model_by_atom.at(0), { 6.0, 0.6, 2.0 }, 1.0e-12);

}

TEST(EstimatorSecondStageDefenseTest, AdaptiveTopologyDriftIsIntensityScaleInvariant)
{
    constexpr double intensity_scale{ 100.0 };
    const detail::FitState reference_state{
        MakeGaussianResult(rg::GaussianModel3D{ 8.0, 0.50, 0.10 })
    };
    const detail::FitState accepted_state{
        MakeGaussianResult(rg::GaussianModel3D{ 9.0, 0.60, 0.15 })
    };
    const detail::FitState scaled_reference_state{
        MakeGaussianResult(rg::GaussianModel3D{
            8.0 * intensity_scale,
            0.50,
            0.10 * intensity_scale })
    };
    const detail::FitState scaled_accepted_state{
        MakeGaussianResult(rg::GaussianModel3D{
            9.0 * intensity_scale,
            0.60,
            0.15 * intensity_scale })
    };

    const auto base{
        detail::CalculateAdaptiveTopologyDrift(
            accepted_state,
            { reference_state.at(0).mdpde.GetModel() },
            { 0 })
    };
    const auto scaled{
        detail::CalculateAdaptiveTopologyDrift(
            scaled_accepted_state,
            { scaled_reference_state.at(0).mdpde.GetModel() },
            { 0 })
    };
    EXPECT_NEAR(
        base,
        scaled,
        1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphNormalizesFullJacobianEnergy)
{
    Eigen::Vector3d jacobian;
    jacobian << 1.0, 2.0, 3.0;
    detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(builder, { 0, 0 }, { { 0, jacobian }, { 1, jacobian } });
    AddCouplingGraphSample(builder, { 0, 1 }, { { 0, 2.0 * jacobian }, { 1, 2.0 * jacobian } });
    const auto topology{ builder.BuildTopology() };
    EXPECT_TRUE(HasCouplingNeighbor(topology, 0, 1));
    EXPECT_NEAR(topology.summary.weight_median, 1.0, 1.0e-12);
    EXPECT_NEAR(topology.summary.weight_percentile_95, 1.0, 1.0e-12);
    EXPECT_NEAR(topology.summary.weight_maximum, 1.0, 1.0e-12);

    detail::CouplingGraphBuilder scaled_builder{ 2 };
    AddCouplingGraphSample(
        scaled_builder,
        { 0, 0 },
        { { 0, 5.0 * jacobian }, { 1, jacobian } });
    AddCouplingGraphSample(
        scaled_builder,
        { 0, 1 },
        { { 0, 10.0 * jacobian }, { 1, 2.0 * jacobian } });
    const auto scaled_topology{
        scaled_builder.BuildTopology()
    };
    EXPECT_TRUE(HasCouplingNeighbor(scaled_topology, 0, 1));
    EXPECT_NEAR(
        scaled_topology.summary.weight_median,
        topology.summary.weight_median,
        1.0e-12);

    detail::CouplingGraphBuilder tiny_builder{ 2 };
    AddCouplingGraphSample(
        tiny_builder,
        { 0, 0 },
        { { 0, 1.0e-100 * jacobian }, { 1, 1.0e-100 * jacobian } });
    const auto tiny_topology{
        tiny_builder.BuildTopology()
    };
    EXPECT_TRUE(HasCouplingNeighbor(tiny_topology, 0, 1));
    EXPECT_NEAR(
        tiny_topology.summary.weight_median,
        topology.summary.weight_median,
        1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphNormalizesDuplicateParticipants)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(
        builder,
        { 0, 0 },
        {
            { 1, unit },
            { 0, unit },
            { 1, unit }
        });
    AddCouplingGraphSample(builder, { 0, 1 }, { { 0, unit } });

    const auto topology{ builder.BuildTopology() };
    ASSERT_EQ(topology.sample_dependency_list.size(), 2U);
    EXPECT_EQ(
        topology.sample_dependency_list.front().contributor_atom_index_list,
        (std::vector<std::size_t>{ 0, 1 }));
    EXPECT_TRUE(HasCouplingNeighbor(topology, 0, 1));
    ASSERT_EQ(topology.retained_edge_list.size(), 1U);
    EXPECT_NEAR(
        topology.retained_edge_list.front().weight,
        1.0 / std::sqrt(2.0),
        1.0e-12);
    EXPECT_EQ(topology.summary.component_count, 1U);
    EXPECT_EQ(topology.summary.maximum_component_size, 2U);
    EXPECT_DOUBLE_EQ(topology.summary.maximum_component_ratio, 1.0);
    EXPECT_DOUBLE_EQ(topology.summary.configured_minimum_weight, 0.05);
    EXPECT_EQ(topology.atom_cutoff_summary.maximum_atom_count_limit, 100U);

    const auto previous_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    detail::LogGraphTopology(topology, false);
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_level);
    EXPECT_NE(output.find(
        "Local-fitting atom cutoff: atoms=2, limit=100, clusters=1, max-atoms=2, cutoff-edges=0."),
        std::string::npos);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphPropagatesInvalidDuplicateJacobian)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    const auto invalid{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(
        builder,
        { 0, 0 },
        {
            { 0, unit },
            { 1, unit },
            { 1, invalid }
        });

    const auto topology{ builder.BuildTopology() };
    EXPECT_FALSE(topology.summary.uses_weighted_graph);
    ASSERT_EQ(topology.sample_dependency_list.size(), 1U);
    EXPECT_EQ(
        topology.sample_dependency_list.front().contributor_atom_index_list,
        (std::vector<std::size_t>{ 0, 1 }));
    EXPECT_TRUE(HasCouplingNeighbor(topology, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphSummaryUsesOnlySelectedSampleConnectivity)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(builder, { 0, 0 }, { { 0, unit } });
    AddCouplingGraphSample(builder, { 1, 0 }, { { 1, unit } });

    const auto topology{
        builder.BuildTopology()
    };
    EXPECT_TRUE(topology.adjacency_list.at(0).empty());
    EXPECT_TRUE(topology.adjacency_list.at(1).empty());
    EXPECT_EQ(topology.summary.component_count, 2U);
    EXPECT_EQ(topology.summary.maximum_component_size, 1U);
    EXPECT_DOUBLE_EQ(topology.summary.maximum_component_ratio, 0.5);

    constexpr std::size_t selected_count{ 11 };
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    detail::SecondStageContext context;
    detail::FitState state;
    context.atom_list.resize(selected_count);
    const rg::GaussianModel3D model{ 6.0, 0.55, 0.10 };
    for (std::size_t i = 0; i < selected_count; i++)
    {
        atoms.emplace_back(MakeAtom(static_cast<int>(i + 1), Spot::C,
            Element::CARBON, { static_cast<double>(i), 0.0, 0.0 }));
        atoms.back()->SetSequenceID(1);
        context.atom_list.at(i).atom = atoms.back().get();
        state.emplace_back(MakeGaussianResult(model));
        context.atom_list.at(i).raw_sampling_entries = {
            { model.ResponseAtDistance(0.2) + model.ResponseAtDistance(0.3), SamplingPoint{ 0.2 } } };
        context.atom_list.at(i).neighbor_atom_sample_offset_list = { 0, 0 };
        context.atom_list.at(i).unselected_distance_list_by_sample = { { 0.3 } };
    }
    const auto background_topology{ detail::BuildSecondStageGraphTopology(context, state, true) };
    EXPECT_EQ(background_topology.adjacency_list.size(), selected_count);
    for (const auto & neighbors : background_topology.adjacency_list) EXPECT_TRUE(neighbors.empty());
    EXPECT_EQ(background_topology.summary.component_count, selected_count);
    EXPECT_EQ(background_topology.summary.maximum_component_size, 1U);
    const auto partition{ detail::BuildGraphPartition(
        background_topology, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }) };
    EXPECT_EQ(partition.sample_id_list_by_key.size(), selected_count);
    for (std::size_t i = 0; i < atoms.size(); i++)
    {
        atoms.at(i)->SetChainID(i % 2 == 0 ? "B" : "C");
        atoms.at(i)->SetSequenceID(static_cast<int>(100 - i));
    }
    const auto relabeled{ detail::BuildSecondStageGraphTopology(context, state, true) };
    EXPECT_EQ(relabeled.adjacency_list, background_topology.adjacency_list);
    EXPECT_EQ(relabeled.summary.component_count, selected_count);
    EXPECT_EQ(detail::BuildGraphPartition(
        relabeled, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }).sample_id_list_by_key,
        partition.sample_id_list_by_key);
    EXPECT_THROW(detail::BuildSecondStageGraphTopology(context, {}, true), std::invalid_argument);

}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphCutsWeakAndCancelledEdges)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    detail::CouplingGraphBuilder weak_builder{ 2 };
    AddCouplingGraphSample(weak_builder, { 0, 0 }, { { 0, unit }, { 1, unit } });
    AddCouplingGraphSample(weak_builder, { 0, 1 }, { { 0, 10.0 * unit } });
    AddCouplingGraphSample(weak_builder, { 1, 0 }, { { 1, 10.0 * unit } });
    const auto weak_topology{
        weak_builder.BuildTopology()
    };
    EXPECT_FALSE(HasCouplingNeighbor(weak_topology, 0, 1));
    EXPECT_EQ(weak_topology.summary.candidate_edge_count, 1U);
    EXPECT_EQ(weak_topology.summary.cut_edge_count, 1U);

    detail::CouplingGraphBuilder cancelled_builder{ 2 };
    AddCouplingGraphSample(cancelled_builder, { 0, 0 }, { { 0, unit }, { 1, unit } });
    AddCouplingGraphSample(cancelled_builder, { 0, 1 }, { { 0, unit }, { 1, -unit } });
    const auto cancelled_topology{
        cancelled_builder.BuildTopology()
    };
    EXPECT_FALSE(HasCouplingNeighbor(cancelled_topology, 0, 1));
}

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

TEST(EstimatorSecondStageDefenseTest, BoundaryReconciliationComponentsUseAcceptedSharedSamples)
{
    const detail::ClusterKey key_a{ 0 };
    const detail::ClusterKey key_b{ 1 };
    const detail::ClusterKey key_c{ 2 };
    const detail::ClusterKey key_d{ 3 };
    const detail::ClusterKey key_e{ 4 };
    const detail::ClusterKey key_f{ 5 };
    const detail::SampleRef sample_ab{ 0, 0 };
    const detail::SampleRef sample_bc{ 1, 0 };
    const detail::SampleRef sample_de{ 3, 0 };
    detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key = {
        { key_a, { sample_ab } },
        { key_b, { sample_ab, sample_bc } },
        { key_c, { sample_bc } },
        { key_d, { sample_de } },
        { key_e, { sample_de } },
        { key_f, {} }
    };
    partition.boundary_sample_dependency_list = {
        { sample_ab, { key_a, key_b }, { 0, 1 } },
        { sample_bc, { key_b, key_c }, { 1, 2 } },
        { sample_de, { key_d, key_e }, { 3, 4 } }
    };
    detail::SecondStageContext context;
    context.atom_list.resize(6);
    const auto component_list{
        detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_f, key_e, key_c, key_a, key_d, key_b })
    };
    ASSERT_EQ(component_list.size(), 2U);
    EXPECT_EQ(
        component_list.at(0).key_list,
        (std::vector<detail::ClusterKey>{ key_a, key_b, key_c }));
    EXPECT_EQ(
        component_list.at(0).affected_sample_ref_list,
        (std::vector<detail::SampleRef>{ sample_ab, sample_bc }));
    EXPECT_EQ(component_list.at(0).boundary_sample_count, 2U);
    EXPECT_EQ(
        component_list.at(0).interface_atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 2 }));
    EXPECT_TRUE(component_list.at(0).halo_atom_index_list.empty());
    EXPECT_EQ(
        component_list.at(1).key_list,
        (std::vector<detail::ClusterKey>{ key_d, key_e }));
    EXPECT_EQ(
        component_list.at(1).affected_sample_ref_list,
        (std::vector<detail::SampleRef>{ sample_de }));
    EXPECT_EQ(component_list.at(1).boundary_sample_count, 1U);
    EXPECT_EQ(
        component_list.at(1).interface_atom_index_list,
        (std::vector<std::size_t>{ 3, 4 }));
    EXPECT_TRUE(component_list.at(1).halo_atom_index_list.empty());
    for (const auto & component : component_list)
    {
        const auto expanded{ detail::ExpandBoundaryReconciliationHalo(context, component, 0) };
        EXPECT_EQ(expanded.halo_atom_index_list, component.interface_atom_index_list);
        EXPECT_EQ(expanded.interface_atom_index_list, component.interface_atom_index_list);
        EXPECT_EQ(expanded.key_list, component.key_list);
        EXPECT_EQ(expanded.affected_sample_ref_list, component.affected_sample_ref_list);
    }

    EXPECT_TRUE(
        detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_c, key_a }).empty());
    EXPECT_EQ(
        detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_c, key_b, key_a }),
        std::vector<detail::BoundaryReconciliationComponent>{
            component_list.front()
        });
}

TEST(EstimatorSecondStageDefenseTest, BoundaryPhysicalHaloRejectsOutOfRangeAtoms)
{
    const detail::ClusterKey key_a{ 0 };
    const detail::ClusterKey key_b{ 1 };
    detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key = {
        { key_a, {} },
        { key_b, {} }
    };
    partition.boundary_sample_dependency_list = {
        { { 0, 0 }, { key_a, key_b }, { 0, 1 } }
    };
    detail::SecondStageContext context;
    context.atom_list.resize(1);

    EXPECT_THROW(
        detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_a, key_b }),
        std::invalid_argument);

    partition.boundary_sample_dependency_list.front()
        .contributor_atom_index_list = { 0 };
    EXPECT_THROW(
        detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_a, key_b }),
        std::invalid_argument);

    const detail::BoundaryReconciliationComponent invalid_interface{
        .key_list = { { 0, 1 } },
        .affected_sample_ref_list = {},
        .interface_atom_index_list = { 1 },
        .halo_atom_index_list = {},
            .boundary_sample_count = 0
    };
    EXPECT_THROW(
        detail::ExpandBoundaryReconciliationHalo(
            context,
            invalid_interface,
            0),
        std::invalid_argument);

    auto invalid_component{ invalid_interface };
    invalid_component.interface_atom_index_list = { 0 };
    EXPECT_THROW(
        detail::ExpandBoundaryReconciliationHalo(
            context,
            invalid_component,
            0),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, BoundaryHaloExpandsPhysicalParticipantsByHop)
{
    detail::SecondStageContext context;
    context.atom_list.resize(5);
    for (auto & atom_context : context.atom_list)
    {
        atom_context.raw_sampling_entries.resize(1);
        atom_context.neighbor_atom_sample_offset_list = { 0, 0 };
    }
    context.atom_list.at(0).neighbor_atom_sample_list = {
        { 1, 0.5 },
        { 3, 0.5 }
    };
    context.atom_list.at(0).neighbor_atom_sample_offset_list = { 0, 2 };
    context.atom_list.at(1).neighbor_atom_sample_list = {
        { 2, 0.5 }
    };
    context.atom_list.at(1).neighbor_atom_sample_offset_list = { 0, 1 };
    context.atom_list.at(2).neighbor_atom_sample_list = {
        { 3, 0.5 },
        { 4, 0.5 }
    };
    context.atom_list.at(2).neighbor_atom_sample_offset_list = { 0, 2 };
    context.atom_list.at(2).unselected_distance_list_by_sample = { { 0.5 } };

    const detail::BoundaryReconciliationComponent component{
        .key_list = { { 0, 1 }, { 2, 3 } },
        .affected_sample_ref_list = { { 0, 0 }, { 1, 0 }, { 2, 0 } },
        .interface_atom_index_list = { 0 },
        .halo_atom_index_list = {},
        .boundary_sample_count = 1
    };
    const auto depth_zero{
        detail::ExpandBoundaryReconciliationHalo(context, component, 0)
    };
    EXPECT_EQ(depth_zero.interface_atom_index_list, (std::vector<std::size_t>{ 0 }));
    EXPECT_EQ(depth_zero.halo_atom_index_list, (std::vector<std::size_t>{ 0 }));

    const auto depth_one{
        detail::ExpandBoundaryReconciliationHalo(context, component, 1)
    };
    EXPECT_EQ(depth_one.interface_atom_index_list, (std::vector<std::size_t>{ 0 }));
    EXPECT_EQ(
        depth_one.halo_atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 3 }));

    const auto fixed_point{
        detail::ExpandBoundaryReconciliationHalo(context, component, 10)
    };
    EXPECT_EQ(
        fixed_point.halo_atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 2, 3 }));
}

TEST(EstimatorSecondStageDefenseTest, UncutDependencyPolishMergesWholeActiveClusters)
{
    detail::GraphTopology topology;
    topology.adjacency_list.resize(6);
    topology.sample_dependency_list = {
        { { 0, 0 }, { 1, 2 } },
        { { 2, 0 }, { 2, 3 } },
        { { 3, 0 }, { 3, 4 } }
    };
    detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key = {
        { { 0, 1 }, { { 0, 0 } } },
        { { 2 }, { { 2, 0 } } },
        { { 3 }, { { 3, 0 } } },
        { { 5 }, { { 5, 0 } } }
    };
    const std::vector<detail::ClusterKey> owner_key_by_atom_index{
        { 0, 1 }, { 0, 1 }, { 2 }, { 3 }, {}, { 5 }
    };
    const auto component_list{
        detail::BuildUncutDependencyPolishComponents(
            topology,
            partition,
            owner_key_by_atom_index)
    };
    ASSERT_EQ(component_list.size(), 1U);
    EXPECT_EQ(
        component_list.front().key_list,
        (std::vector<detail::ClusterKey>{ { 0, 1 }, { 2 }, { 3 } }));
    EXPECT_EQ(
        component_list.front().atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 2, 3 }));
    EXPECT_EQ(
        component_list.front().affected_sample_ref_list,
        (std::vector<detail::SampleRef>{ { 0, 0 }, { 2, 0 }, { 3, 0 } }));
}

TEST(EstimatorSecondStageDefenseTest, CouplingAtomCutoffBoundsComponentsAndPreservesDependencies)
{
    for (const std::size_t atom_count : { 100U, 101U, 102U })
    {
        detail::GraphTopology topology;
        topology.adjacency_list.resize(atom_count);
        // The 102-atom case also has an isolated atom, with no forced closure.
        const std::size_t first_connected_atom{ atom_count == 102 ? 1U : 0U };
        for (std::size_t atom_index = first_connected_atom; atom_index + 1 < atom_count; atom_index++)
        {
            topology.retained_edge_list.push_back({
                atom_index, atom_index + 1, 1.0 - 0.001 * static_cast<double>(atom_index) });
        }
        std::vector<std::size_t> active_index_list(atom_count);
        std::iota(active_index_list.begin(), active_index_list.end(), 0U);
        topology.sample_dependency_list = { { { 0, 0 }, active_index_list } };
        const auto capped_topology{ detail::ApplyGraphAtomCutoff(topology, 100) };
        EXPECT_EQ(capped_topology.adjacency_list.size(), atom_count);
        EXPECT_EQ(capped_topology.summary.maximum_component_size, 100U);
        EXPECT_EQ(capped_topology.summary.component_count, atom_count - 99);
        EXPECT_EQ(capped_topology.atom_cutoff_summary.cut_edge_count, atom_count == 100 ? 0U : 1U);
        EXPECT_EQ(capped_topology.retained_edge_list.size(), topology.retained_edge_list.size());
        ASSERT_EQ(capped_topology.sample_dependency_list.size(), 1U);
        EXPECT_EQ(capped_topology.sample_dependency_list.front().contributor_atom_index_list,
            active_index_list);
        const auto partition{ detail::BuildGraphPartition(capped_topology, active_index_list) };
        EXPECT_EQ(partition.boundary_sample_dependency_list.size(), atom_count == 100 ? 0U : 1U);
        EXPECT_EQ(partition.sample_id_list_by_key.size(), capped_topology.summary.component_count);
        std::vector<detail::ClusterKey> owner_key_by_atom_index(atom_count);
        std::size_t maximum_component_size{ 0 };
        for (const auto & [key, sample_id_list] : partition.sample_id_list_by_key)
        {
            EXPECT_EQ(sample_id_list.size(), 1U);
            EXPECT_LE(key.size(), 100U);
            maximum_component_size = std::max(maximum_component_size, key.size());
            for (const auto atom_index : key) owner_key_by_atom_index.at(atom_index) = key;
        }
        EXPECT_EQ(capped_topology.summary.maximum_component_size, maximum_component_size);
        EXPECT_DOUBLE_EQ(capped_topology.summary.maximum_component_ratio,
            static_cast<double>(maximum_component_size) / static_cast<double>(atom_count));
        if (first_connected_atom == 1)
            EXPECT_EQ(partition.sample_id_list_by_key.count({ 0 }), 1U);
        const auto polish_components{ detail::BuildUncutDependencyPolishComponents(
            capped_topology, partition, owner_key_by_atom_index) };
        ASSERT_EQ(polish_components.size(), 1U);
        EXPECT_EQ(polish_components.front().atom_index_list.size(), atom_count);

        std::reverse(active_index_list.begin(), active_index_list.end());
        const auto reversed_partition{ detail::BuildGraphPartition(capped_topology, active_index_list) };
        EXPECT_EQ(partition.sample_id_list_by_key, reversed_partition.sample_id_list_by_key);
        ASSERT_EQ(partition.boundary_sample_dependency_list.size(), reversed_partition.boundary_sample_dependency_list.size());
        if (!partition.boundary_sample_dependency_list.empty())
        {
            auto expected{ partition.boundary_sample_dependency_list.front() };
            auto actual{ reversed_partition.boundary_sample_dependency_list.front() };
            std::sort(expected.cluster_key_list.begin(), expected.cluster_key_list.end());
            std::sort(actual.cluster_key_list.begin(), actual.cluster_key_list.end());
            EXPECT_EQ(expected, actual);
        }
    }

    const auto empty{ detail::CouplingGraphBuilder{ 0 }.BuildTopology() };
    EXPECT_TRUE(empty.adjacency_list.empty());
    EXPECT_EQ(empty.summary.component_count, 0U);
    EXPECT_EQ(empty.summary.maximum_component_size, 0U);
    EXPECT_DOUBLE_EQ(empty.summary.maximum_component_ratio, 0.0);
    EXPECT_TRUE(detail::BuildGraphPartition(empty, {}).sample_id_list_by_key.empty());
    EXPECT_THROW(detail::ApplyGraphAtomCutoff({}, 0), std::invalid_argument);
    detail::CouplingGraphOptions invalid_options;
    invalid_options.maximum_atom_count = 0;
    EXPECT_THROW(detail::CouplingGraphBuilder{ 0 }.BuildTopology(invalid_options),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, CouplingAtomCutoffPrioritizesStrongEdgesAndStableTies)
{
    detail::GraphTopology topology;
    topology.adjacency_list.resize(3);
    topology.retained_edge_list = { { 1, 2, 0.80 }, { 0, 1, 0.90 }, { 0, 2, 0.70 } };
    const auto capped_topology{ detail::ApplyGraphAtomCutoff(topology, 2) };
    const auto partition{ detail::BuildGraphPartition(capped_topology, { 2, 1, 0 }) };
    EXPECT_EQ(capped_topology.summary.component_count, partition.sample_id_list_by_key.size());
    EXPECT_EQ(capped_topology.summary.maximum_component_size, 2U);
    EXPECT_DOUBLE_EQ(capped_topology.summary.maximum_component_ratio, 2.0 / 3.0);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_TRUE(HasCouplingNeighbor(capped_topology, 0, 1));
    EXPECT_FALSE(HasCouplingNeighbor(capped_topology, 1, 2));

    const auto singleton_topology{ detail::ApplyGraphAtomCutoff(topology, 1) };
    EXPECT_EQ(singleton_topology.summary.component_count, 3U);
    EXPECT_EQ(singleton_topology.summary.maximum_component_size, 1U);
    EXPECT_DOUBLE_EQ(singleton_topology.summary.maximum_component_ratio, 1.0 / 3.0);
    EXPECT_EQ(singleton_topology.atom_cutoff_summary.cut_edge_count, 3U);
    const auto whole_topology{ detail::ApplyGraphAtomCutoff(topology, 3) };
    EXPECT_EQ(whole_topology.summary.component_count, 1U);
    EXPECT_EQ(whole_topology.summary.maximum_component_size, 3U);
    EXPECT_DOUBLE_EQ(whole_topology.summary.maximum_component_ratio, 1.0);
    // All three internal edges survive, not only the union-find spanning tree.
    for (const auto & neighbors : whole_topology.adjacency_list) EXPECT_EQ(neighbors.size(), 2U);

    for (auto & edge : topology.retained_edge_list) edge.weight = 0.8;
    for (std::size_t order = 0; order < 3; order++)
    {
        std::rotate(topology.retained_edge_list.begin(),
            topology.retained_edge_list.begin() + 1, topology.retained_edge_list.end());
        for (auto & edge : topology.retained_edge_list)
            std::swap(edge.left_atom_index, edge.right_atom_index);
        const auto tied_topology{ detail::ApplyGraphAtomCutoff(topology, 2) };
        EXPECT_EQ(detail::BuildGraphPartition(tied_topology, { 0, 1, 2 }).sample_id_list_by_key,
            partition.sample_id_list_by_key);
    }

    topology.retained_edge_list.front().weight = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(detail::ApplyGraphAtomCutoff(topology, 2), std::invalid_argument);
    topology.retained_edge_list.front().weight = 1.0;
    topology.retained_edge_list.front().left_atom_index = 3;
    EXPECT_THROW(detail::ApplyGraphAtomCutoff(topology, 2), std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, CouplingPartitionKeepsStrongChainAndBinaryFallback)
{
    detail::GraphTopology strong_topology;
    strong_topology.adjacency_list = {
        { 1 },
        { 0, 2 },
        { 1 }
    };
    const auto strong_partition{
        detail::BuildGraphPartition(
            strong_topology,
            { 0, 1, 2 })
    };
    EXPECT_EQ(strong_partition.sample_id_list_by_key.count({ 0, 1, 2 }), 1U);

    detail::CouplingGraphBuilder builder{ 2 };
    const auto invalid{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    AddCouplingGraphSample(
        builder,
        { 0, 0 },
        { { 0, Eigen::Vector3d::Ones() }, { 1, invalid } });
    detail::CouplingGraphOptions fallback_options;
    const auto binary_topology{
        builder.BuildTopology(fallback_options)
    };
    EXPECT_FALSE(binary_topology.summary.uses_weighted_graph);
    const auto binary_partition{
        detail::BuildGraphPartition(
            binary_topology,
            { 0, 1 })
    };
    EXPECT_EQ(binary_partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    const auto capped_binary_topology{
        detail::ApplyGraphAtomCutoff(
            binary_topology,
            1)
    };
    const auto capped_binary_partition{
        detail::BuildGraphPartition(
            capped_binary_topology,
            { 0, 1 })
    };
    EXPECT_EQ(capped_binary_partition.sample_id_list_by_key.size(), 2U);

    detail::CouplingGraphBuilder overflow_builder{ 2 };
    const auto huge{ Eigen::Vector3d::Constant(1.0e200) };
    AddCouplingGraphSample(
        overflow_builder,
        { 0, 0 },
        { { 0, huge }, { 1, huge } });
    fallback_options.maximum_atom_count = 1;
    const auto overflow_topology{
        overflow_builder.BuildTopology(fallback_options)
    };
    EXPECT_FALSE(overflow_topology.summary.uses_weighted_graph);
    EXPECT_EQ(overflow_topology.summary.component_count, 2U);
    EXPECT_EQ(overflow_topology.atom_cutoff_summary.maximum_atom_count_limit, 1U);
    EXPECT_FALSE(HasCouplingNeighbor(overflow_topology, 0, 1));
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

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageIterationsUsesFrozenGlobalBackgroundWithoutGroupOrResidueKeys)
{
    const std::array seeds{ rg::GaussianModel3D{ 5.0, 0.50, 0.05 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.15 } };
    const std::array truth{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 }, rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 5.0, 0.48, 0.10 }, rg::GaussianModel3D{ 7.0, 0.63, 0.10 },
        rg::GaussianModel3D{ 2.5, 0.45, 0.05 }, rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 } };
    auto options{ MakeSecondStageOptions() };
    options.exclude_hydrogen = true;
    for (const auto & [shared_cluster, shared_contributor] : {
        std::pair{ false, false }, std::pair{ true, false }, std::pair{ false, true } })
    {
        std::array<rg::GaussianModel3D, 2> reference_selected;
        std::array<rg::GaussianModel3D, 2> reference_background;
        for (const double scale : { 1.0, 100.0 })
        {
            auto scaled_seeds{ seeds };
            auto scaled_truth{ truth };
            for (auto & model : scaled_seeds)
                model = { scale * model.GetAmplitude(), model.GetWidth(), scale * model.GetOffset() };
            for (auto & model : scaled_truth)
                model = { scale * model.GetAmplitude(), model.GetWidth(), scale * model.GetOffset() };
            auto serial{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, false, shared_cluster, shared_contributor) };
            auto parallel{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, true, shared_cluster, shared_contributor) };
            auto logged{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, false, shared_cluster, shared_contributor) };
            EXPECT_EQ(rg::data_internal::GetGroupKey(serial->FindAtomPtr(3)),
                rg::data_internal::GetGroupKey(serial->FindAtomPtr(4)));
            EXPECT_NE(rg::data_internal::GetGroupKey(serial->FindAtomPtr(3)),
                rg::data_internal::GetGroupKey(parallel->FindAtomPtr(3)));
            options.thread_size = 1;
            serial->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
            detail::RunSecondStageIterations(*serial, options);
            options.thread_size = 2;
            parallel->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
            detail::RunSecondStageIterations(*parallel, options);
            const auto previous_level{ Logger::GetLogLevel() };
            Logger::SetLogLevel(LogLevel::Debug);
            testing::internal::CaptureStdout();
            options.thread_size = 1;
            options.quiet_mode = false;
            logged->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
            second_stage_test::BeginNumericalCapture();
            detail::RunSecondStageIterations(*logged, options);
            const auto captured{ second_stage_test::EndNumericalCapture() };
            const auto output{ testing::internal::GetCapturedStdout() };
            auto alternate_logged{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, true, shared_cluster, shared_contributor) };
            testing::internal::CaptureStdout();
            options.thread_size = 2;
            alternate_logged->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
            detail::RunSecondStageIterations(*alternate_logged, options);
            const auto alternate_output{ testing::internal::GetCapturedStdout() };
            auto relabeled_logged{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, false, shared_cluster, shared_contributor, true) };
            testing::internal::CaptureStdout();
            relabeled_logged->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
            detail::RunSecondStageIterations(*relabeled_logged, options);
            const auto relabeled_output{ testing::internal::GetCapturedStdout() };
            Logger::SetLogLevel(previous_level);
            options.quiet_mode = true;
            const auto audit_records = [](const std::string & log)
            {
                std::vector<std::string> records;
                std::istringstream lines{ log };
                std::string line;
                while (std::getline(lines, line))
                {
                    for (const std::string marker : {
                        "Convergence safeguard audit:", "Second-stage audit terminal:",
                        "Second-stage audit terminal atom:", "Second-stage local fitting summary:",
                        "Local-fitting atom cutoff:", "Adaptive local-fitting topology rebuild:",
                        "Cluster best source:", "Cluster best publication:" })
                    {
                        const auto position{ line.find(marker) };
                        if (position != std::string::npos)
                        {
                            if (marker == "Adaptive local-fitting topology rebuild:")
                                EXPECT_NE(line.find(", trigger=drift, drift="), std::string::npos);
                            records.emplace_back(line.substr(position));
                        }
                    }
                }
                return records;
            };
            EXPECT_EQ(audit_records(output), audit_records(alternate_output));
            EXPECT_EQ(audit_records(output), audit_records(relabeled_output));

            EXPECT_NE(output.find(shared_cluster ?
                "initial components/max atoms/ratio = 1/2/1.00" :
                "initial components/max atoms/ratio = 2/1/0.50"), std::string::npos);
            if (shared_contributor)
            {
                EXPECT_NE(output.find("candidate/retained/cut edges = 0/0/0"), std::string::npos);
            }
            std::array<std::optional<rg::GaussianModel3D>, 2> first_background;
            std::array<std::optional<rg::GaussianModel3D>, 2> last_background;
            ASSERT_GT(captured.backgrounds.size(), 1U);
            for (std::size_t target = 0; target < 2; target++)
            {
                const auto & first{ captured.backgrounds.front().at(target) };
                const auto & last{ captured.backgrounds.back().at(target) };
                first_background.at(target) = rg::GaussianModel3D{first[0], first[1], first[2]};
                last_background.at(target) = rg::GaussianModel3D{last[0], last[1], last[2]};
                for (const auto & background : captured.backgrounds)
                    EXPECT_EQ(background.at(0), background.at(1));
                const auto expected_initial{
                    *detail::BuildGaussianParameterMedian({ scaled_seeds[0], scaled_seeds[1] }) };
                ExpectGaussianModelsNear(*first_background.at(target), expected_initial, 1.0e-12 * scale);
                const int serial_id{ static_cast<int>(target + 1) };
                const auto selected{ GetEstimateModel(*serial->FindAtomPtr(serial_id)) };
                ExpectGaussianModelsNear(selected, GetEstimateModel(*parallel->FindAtomPtr(serial_id)), 1.0e-10 * scale);
                ExpectGaussianModelsNear(selected, GetEstimateModel(*logged->FindAtomPtr(serial_id)), 1.0e-10 * scale);
                ExpectGaussianModelsNear(selected, GetEstimateModel(*alternate_logged->FindAtomPtr(serial_id)), 0.0);
                ExpectGaussianModelsNear(selected, GetEstimateModel(*relabeled_logged->FindAtomPtr(serial_id)), 0.0);
                const auto view{ rg::AtomLocalPotentialView::For(*serial->FindAtomPtr(serial_id)) };
                const auto raw{ view.GetRawSamplingEntries(false) };
                const auto peeled{ view.GetPeelingSamplingEntries(false) };
                const auto alternate{ rg::AtomLocalPotentialView::For(*parallel->FindAtomPtr(serial_id)).GetPeelingSamplingEntries(false) };
                const auto relabeled{ rg::AtomLocalPotentialView::For(*relabeled_logged->FindAtomPtr(serial_id)).GetPeelingSamplingEntries(false) };
                ASSERT_EQ(raw.size(), peeled.size());
                ASSERT_EQ(alternate.size(), peeled.size());
                ASSERT_EQ(relabeled.size(), peeled.size());
                for (std::size_t row = 0; row < raw.size(); row++)
                {
                    double background{ 0.0 };
                    for (const int contributor_id : { 3, 4, 6, 7 })
                    {
                        const auto distance{ Distance(raw.at(row).point.position,
                            serial->FindAtomPtr(contributor_id)->GetPosition()) };
                        if (distance <= 2.5)
                            background += last_background.at(target)->ResponseAtDistance(distance);
                    }
                    const auto * selected_neighbor{ serial->FindAtomPtr(target == 0 ? 2 : 1) };
                    const auto selected_distance{ Distance(raw.at(row).point.position, selected_neighbor->GetPosition()) };
                    if (selected_distance <= 2.5)
                        background += GetEstimateModel(*selected_neighbor).ResponseAtDistance(selected_distance);
                    EXPECT_NEAR(peeled.at(row).response, raw.at(row).response - background, 1.0e-10 * scale);
                    EXPECT_NEAR(peeled.at(row).response, alternate.at(row).response, 1.0e-10 * scale);
                    EXPECT_DOUBLE_EQ(peeled.at(row).response, relabeled.at(row).response);
                }
                EXPECT_EQ(view.GetNeighborCountForPeeling(), shared_cluster ? 3 :
                    shared_contributor && target == 1 ? 2 : 1);
                const rg::GaussianModel3D normalized_selected{
                    selected.GetAmplitude() / scale, selected.GetWidth(), selected.GetOffset() / scale };
                const rg::GaussianModel3D normalized_background{
                    last_background.at(target)->GetAmplitude() / scale, last_background.at(target)->GetWidth(),
                    last_background.at(target)->GetOffset() / scale };
                if (scale == 1.0)
                {
                    reference_selected.at(target) = normalized_selected;
                    reference_background.at(target) = normalized_background;
                }
                else
                {
                    ExpectGaussianModelsNear(normalized_selected, reference_selected.at(target), 1.0e-6);
                    ExpectGaussianModelsNear(normalized_background, reference_background.at(target), 1.0e-6);
                }
            }
            if (shared_cluster)
                ExpectGaussianModelsNear(*last_background[0], *last_background[1], 0.0);
            EXPECT_EQ(serial->GetSelectedAtomCount(), 2U);
            for (int id = 3; id <= 7; id++)
            {
                EXPECT_FALSE(rg::AtomLocalPotentialView::For(*serial->FindAtomPtr(id)).IsAvailable());
                EXPECT_FALSE(rg::AtomLocalPotentialView::For(*parallel->FindAtomPtr(id)).IsAvailable());
            }
        }
    }
    auto excluded{ BuildUnselectedContributorDefenseModel(seeds, truth) };
    auto included{ BuildUnselectedContributorDefenseModel(seeds, truth) };
    excluded->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*excluded, options);
    options.exclude_hydrogen = false;
    included->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*included, options);
    EXPECT_EQ(rg::AtomLocalPotentialView::For(*included->FindAtomPtr(1)).GetNeighborCountForPeeling(), 2);
    EXPECT_NE(rg::AtomLocalPotentialView::For(*excluded->FindAtomPtr(1)).GetPeelingSamplingEntries(false).front().response,
        rg::AtomLocalPotentialView::For(*included->FindAtomPtr(1)).GetPeelingSamplingEntries(false).front().response);

}
