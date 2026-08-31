#pragma once

#include "core/detail/SecondStageFitting.hpp"

#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::core::detail {

using ResidueKey = std::pair<std::string, int>;

struct GraphWeightedEdge
{
    std::size_t left_atom_index{ 0 };
    std::size_t right_atom_index{ 0 };
    double weight{ 0.0 };
};

struct GraphParticipant
{
    std::size_t atom_index{ 0 };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
};

struct GraphSampleDependency
{
    SampleRef sample_id{};
    std::vector<std::size_t> contributor_atom_index_list{};
};

struct CouplingGraphSummary
{
    struct ThresholdSensitivity
    {
        double minimum_weight{ 0.0 };
        std::size_t retained_edge_count{ 0 };
        std::size_t cut_edge_count{ 0 };
        std::size_t component_count{ 0 };
        std::size_t maximum_component_size{ 0 };
        double maximum_component_ratio{ 0.0 };
    };

    bool uses_weighted_graph{ false };
    std::size_t candidate_edge_count{ 0 };
    std::size_t retained_edge_count{ 0 };
    std::size_t cut_edge_count{ 0 };
    double weight_median{ 0.0 };
    double weight_percentile_95{ 0.0 };
    double weight_maximum{ 0.0 };
    std::vector<ThresholdSensitivity> threshold_sensitivity_list{};
    // The requested threshold, retained for diagnostics when fallback occurs.
    double configured_minimum_weight{ 0.0 };
    std::size_t component_count{ 0 };
    std::size_t maximum_component_size{ 0 };
    double maximum_component_ratio{ 0.0 };
};

struct GraphTopology
{
    std::vector<std::vector<std::size_t>> adjacency_list{};
    std::vector<GraphWeightedEdge> retained_edge_list{};
    std::vector<ResidueKey> residue_key_by_atom_index{};
    std::vector<GraphSampleDependency> sample_dependency_list{};
    CouplingGraphSummary summary{};
    struct ResidueCutoffSummary
    {
        std::size_t residue_count{ 0 };
        std::size_t maximum_residue_count{ 0 };
        std::size_t cluster_count{ 0 };
        std::size_t cut_edge_count{ 0 };
        std::size_t maximum_residue_count_limit{ 0 };
    } residue_cutoff_summary{};
};

struct CouplingGraphOptions
{
    double minimum_weight{ 0.05 };
    std::optional<double> retained_edge_minimum_weight{};
    std::vector<double> sensitivity_minimum_weight_list{
        0.05, 0.075, 0.10, 0.15, 0.20, 0.30
    };
    std::size_t maximum_residue_count{ 10 };
};

struct CouplingGraphPartition
{
    struct BoundarySampleDependency
    {
        SampleRef sample_id{};
        std::vector<ClusterKey> cluster_key_list{};
        std::vector<std::size_t> contributor_atom_index_list{};

        friend bool operator==(
            const BoundarySampleDependency &,
            const BoundarySampleDependency &) = default;
    };

    std::map<ClusterKey, std::vector<SampleRef>> sample_id_list_by_key{};
    std::vector<BoundarySampleDependency> boundary_sample_dependency_list{};
};

struct BoundaryReconciliationComponent
{
    std::vector<ClusterKey> key_list{};
    std::vector<SampleRef> affected_sample_ref_list{};
    std::vector<std::size_t> interface_atom_index_list{};
    std::vector<std::size_t> shape_active_atom_index_list{};
    std::vector<std::size_t> offset_closure_atom_index_list{};
    std::size_t boundary_sample_count{ 0 };

    friend bool operator==(
        const BoundaryReconciliationComponent &,
        const BoundaryReconciliationComponent &) = default;
};

struct DependencyPolishComponent
{
    std::vector<ClusterKey> key_list{};
    std::vector<SampleRef> affected_sample_ref_list{};
    std::vector<std::size_t> atom_index_list{};

    friend bool operator==(
        const DependencyPolishComponent &,
        const DependencyPolishComponent &) = default;
};

struct GraphIndexPairHash
{
    std::size_t operator()(const std::pair<std::size_t, std::size_t> & pair) const noexcept
    {
        const auto left_hash{ std::hash<std::size_t>{}(pair.first) };
        const auto right_hash{ std::hash<std::size_t>{}(pair.second) };
        return left_hash ^ (right_hash + static_cast<std::size_t>(0x9e3779b9) +
            (left_hash << 6) + (left_hash >> 2));
    }
};

class CouplingGraphBuilder
{
    using AtomPair = std::pair<std::size_t, std::size_t>;
    std::size_t m_atom_count{ 0 };
    std::vector<Eigen::Matrix3d> m_self_gram_list{};
    std::unordered_map<AtomPair, Eigen::Matrix3d, GraphIndexPairHash> m_pair_accumulator_by_pair{};
    std::vector<GraphSampleDependency> m_sample_dependency_list{};
    bool m_has_invalid_jacobian{ false };

public:
    explicit CouplingGraphBuilder(std::size_t atom_count);
    void AddSample(SampleRef sample_id, std::vector<GraphParticipant> & participant_list);
    GraphTopology BuildTopology(
        std::vector<ResidueKey> residue_key_by_atom_index,
        const CouplingGraphOptions & options = {},
        const GraphTopology * previous_topology = nullptr);

private:
    std::vector<CouplingGraphSummary::ThresholdSensitivity> BuildThresholdSensitivity(
        const std::vector<GraphWeightedEdge> & weighted_edge_list,
        const std::vector<double> & minimum_weight_list) const;

    GraphTopology BuildFromWeights(
        const std::vector<GraphWeightedEdge> & weighted_edge_list,
        const CouplingGraphOptions & options,
        const GraphTopology * previous_topology);

    GraphTopology BuildWeightedOrBinary(
        const CouplingGraphOptions & options,
        const GraphTopology * previous_topology);

    GraphTopology BuildBinary();
};

GraphTopology BuildSecondStageGraphTopology(
    const SecondStageContext & context,
    const FitState & initial_state,
    bool quiet_mode);

GraphTopology BuildAdaptiveSecondStageGraphTopology(
    const SecondStageContext & context,
    const FitState & accepted_state,
    const GraphTopology & previous_topology,
    bool quiet_mode);

void LogGraphTopology(const GraphTopology & topology, bool quiet_mode);

GraphTopology ApplyGraphResidueCutoff(GraphTopology topology, std::size_t maximum_residue_count);

CouplingGraphPartition BuildGraphPartition(
    const GraphTopology & topology,
    const std::vector<std::size_t> & active_index_list);

std::vector<ClusterKey> BuildGraphClusterKeyList(const CouplingGraphPartition & partition);

std::vector<SampleRef> BuildGraphAffectedSampleUnion(
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & key_list);

std::vector<BoundaryReconciliationComponent> BuildBoundaryReconciliationComponents(
    const SecondStageContext & context,
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & accepted_key_list);

BoundaryReconciliationComponent ExpandBoundaryReconciliationHalo(
    const SecondStageContext & context,
    BoundaryReconciliationComponent component,
    std::size_t halo_depth);

std::vector<DependencyPolishComponent> BuildUncutDependencyPolishComponents(
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & owner_key_by_atom_index);

} // namespace rhbm_gem::core::detail
