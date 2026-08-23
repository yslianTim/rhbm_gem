#pragma once

#include "core/detail/FittingModel.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::core::detail {

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

class DisjointSet
{
    std::vector<std::size_t> m_parent_list{};
    std::vector<std::size_t> m_component_size_list{};

public:
    explicit DisjointSet(std::size_t item_count);

    std::size_t Find(std::size_t index);

    void Merge(std::size_t left, std::size_t right);

    std::size_t ComponentSize(std::size_t index);
};

struct DisjointSetComponentSummary
{
    std::size_t component_count{ 0 };
    std::size_t maximum_component_size{ 0 };
};

DisjointSetComponentSummary SummarizeDisjointSetComponents(
    DisjointSet & component_set,
    std::size_t item_count);

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
    std::vector<double> sensitivity_minimum_weight_list{
        0.05, 0.075, 0.10, 0.15, 0.20, 0.30
    };
    std::size_t maximum_residue_count{ 10 };
};

struct CouplingGraphPartition
{
    std::map<ClusterKey, std::vector<SampleRef>> sample_id_list_by_key{};
    std::size_t boundary_sample_count{ 0 };
};

GraphTopology ApplyGraphResidueCutoff(GraphTopology topology, std::size_t maximum_residue_count);
void UpdateGraphComponentSummary(GraphTopology & topology);

class CouplingGraphBuilder
{
    using AtomPair = std::pair<std::size_t, std::size_t>;

    std::size_t m_atom_count{ 0 };
    std::vector<Eigen::Matrix3d> m_self_gram_list{};
    std::map<AtomPair, Eigen::Matrix3d> m_pair_accumulator_by_pair{};
    std::vector<GraphSampleDependency> m_sample_dependency_list{};
    bool m_has_invalid_jacobian{ false };

    static void NormalizeParticipantList(
        std::vector<GraphParticipant> & participant_list);

    static void ValidateBuildOptions(const CouplingGraphOptions & options);

    static double FrobeniusNorm(const Eigen::Matrix3d & matrix);

    std::vector<CouplingGraphSummary::ThresholdSensitivity> BuildThresholdSensitivity(
        const std::vector<GraphWeightedEdge> & weighted_edge_list,
        const std::vector<double> & minimum_weight_list) const;

    GraphTopology BuildFromWeights(
        const std::vector<GraphWeightedEdge> & weighted_edge_list,
        double minimum_weight);

public:
    explicit CouplingGraphBuilder(std::size_t atom_count);

    void AddSample(
        SampleRef sample_id,
        std::vector<GraphParticipant> & participant_list);

private:
    GraphTopology BuildWeightedOrBinary(
        double minimum_weight,
        const std::vector<double> & sensitivity_minimum_weight_list);

    GraphTopology BuildBinary();

public:
    GraphTopology BuildTopology(
        std::vector<ResidueKey> residue_key_by_atom_index,
        const CouplingGraphOptions & options = {});
};

Eigen::Vector3d EvaluateCouplingGraphJacobian(
    const std::optional<TransformedModelInvariants> & invariants,
    double distance);

GraphTopology BuildSecondStageGraphTopology(
    const SecondStageContext & context,
    const FitState & initial_state,
    bool quiet_mode);

void LogGraphTopology(const GraphTopology & topology, bool quiet_mode);

GraphTopology ApplyGraphResidueCutoff(
    GraphTopology topology,
    std::size_t maximum_residue_count);

void UpdateGraphComponentSummary(GraphTopology & topology);

CouplingGraphPartition BuildGraphPartition(
    const GraphTopology & topology,
    const std::vector<std::size_t> & active_index_list);

std::vector<ClusterKey> BuildGraphClusterKeyList(const CouplingGraphPartition & partition);

std::vector<SampleRef> BuildGraphAffectedSampleUnion(
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & key_list);

} // namespace rhbm_gem::core::detail
