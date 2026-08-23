#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include "core/detail/ResidualEvaluation.hpp"
#include "core/detail/SecondStageContext.hpp"
#include "core/detail/TransformedGaussianModel.hpp"

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
    explicit DisjointSet(std::size_t item_count)
        : m_parent_list(item_count),
          m_component_size_list(item_count, 1)
    {
        for (std::size_t i = 0; i < item_count; i++) m_parent_list.at(i) = i;
    }

    std::size_t Find(std::size_t index)
    {
        if (m_parent_list.at(index) == index) return index;
        m_parent_list.at(index) = Find(m_parent_list.at(index));
        return m_parent_list.at(index);
    }

    void Merge(std::size_t left, std::size_t right)
    {
        const auto left_root{ Find(left) };
        const auto right_root{ Find(right) };
        if (left_root == right_root) return;
        m_parent_list.at(right_root) = left_root;
        m_component_size_list.at(left_root) += m_component_size_list.at(right_root);
    }

    std::size_t ComponentSize(std::size_t index)
    {
        return m_component_size_list.at(Find(index));
    }
};

struct DisjointSetComponentSummary
{
    std::size_t component_count{ 0 };
    std::size_t maximum_component_size{ 0 };
};

inline DisjointSetComponentSummary SummarizeDisjointSetComponents(
    DisjointSet & component_set,
    std::size_t item_count)
{
    DisjointSetComponentSummary summary;
    for (std::size_t item_index = 0; item_index < item_count; item_index++)
    {
        if (component_set.Find(item_index) != item_index) continue;
        summary.component_count++;
        summary.maximum_component_size = std::max(
            summary.maximum_component_size,
            component_set.ComponentSize(item_index));
    }
    return summary;
}

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

inline GraphTopology ApplyGraphResidueCutoff(GraphTopology topology, std::size_t maximum_residue_count);
inline void UpdateGraphComponentSummary(GraphTopology & topology);

class CouplingGraphBuilder
{
    using AtomPair = std::pair<std::size_t, std::size_t>;

    std::size_t m_atom_count{ 0 };
    std::vector<Eigen::Matrix3d> m_self_gram_list{};
    std::map<AtomPair, Eigen::Matrix3d> m_pair_accumulator_by_pair{};
    std::vector<GraphSampleDependency> m_sample_dependency_list{};
    bool m_has_invalid_jacobian{ false };

    static void NormalizeParticipantList(std::vector<GraphParticipant> & participant_list)
    {
        std::ranges::sort(participant_list, {}, &GraphParticipant::atom_index);

        std::size_t normalized_size{ 0 };
        for (std::size_t index = 0; index < participant_list.size(); index++)
        {
            auto & participant{ participant_list.at(index) };
            if (normalized_size == 0 ||
                participant_list.at(normalized_size - 1).atom_index != participant.atom_index)
            {
                if (normalized_size != index)
                {
                    participant_list.at(normalized_size) = std::move(participant);
                }
                normalized_size++;
                continue;
            }

            auto & accumulated_jacobian{
                participant_list.at(normalized_size - 1).jacobian
            };
            if (!accumulated_jacobian.allFinite() || !participant.jacobian.allFinite())
            {
                accumulated_jacobian = Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
                continue;
            }
            accumulated_jacobian += participant.jacobian;
            if (!accumulated_jacobian.allFinite())
            {
                accumulated_jacobian = Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());
            }
        }
        participant_list.resize(normalized_size);
    }

    static void ValidateBuildOptions(const CouplingGraphOptions & options)
    {
        if (!std::isfinite(options.minimum_weight) ||
            options.minimum_weight < 0.0 || options.minimum_weight > 1.0)
        {
            throw std::invalid_argument(
                "Local fitting coupling minimum weight must be in [0, 1].");
        }
        if (options.maximum_residue_count == 0)
        {
            throw std::invalid_argument(
                "Local fitting coupling maximum residue count must be positive.");
        }
        for (const auto minimum_weight : options.sensitivity_minimum_weight_list)
        {
            if (!std::isfinite(minimum_weight) || minimum_weight < 0.0 || minimum_weight > 1.0)
            {
                throw std::invalid_argument(
                    "Local fitting coupling sensitivity minimum weight must be in [0, 1].");
            }
        }
    }

    static double FrobeniusNorm(const Eigen::Matrix3d & matrix)
    {
        double norm{ 0.0 };
        for (Eigen::Index row = 0; row < matrix.rows(); row++)
        {
            for (Eigen::Index column = 0; column < matrix.cols(); column++)
            {
                norm = std::hypot(norm, matrix(row, column));
            }
        }
        return norm;
    }

    std::vector<CouplingGraphSummary::ThresholdSensitivity> BuildThresholdSensitivity(
        const std::vector<GraphWeightedEdge> & weighted_edge_list,
        const std::vector<double> & minimum_weight_list) const
    {
        std::vector<CouplingGraphSummary::ThresholdSensitivity> sensitivity_list;
        sensitivity_list.reserve(minimum_weight_list.size());
        for (const auto minimum_weight : minimum_weight_list)
        {
            DisjointSet component_set{ m_atom_count };
            std::size_t retained_edge_count{ 0 };
            for (const auto & weighted_edge : weighted_edge_list)
            {
                const auto weight{ weighted_edge.weight };
                if (weight < minimum_weight) continue;
                retained_edge_count++;
                component_set.Merge(weighted_edge.left_atom_index, weighted_edge.right_atom_index);
            }

            const auto component_summary{
                SummarizeDisjointSetComponents(component_set, m_atom_count)
            };
            sensitivity_list.emplace_back(
                CouplingGraphSummary::ThresholdSensitivity{
                    minimum_weight,
                    retained_edge_count,
                    weighted_edge_list.size() - retained_edge_count,
                    component_summary.component_count,
                    component_summary.maximum_component_size,
                    m_atom_count == 0 ? 0.0 :
                        static_cast<double>(component_summary.maximum_component_size) / static_cast<double>(m_atom_count)
                });
        }
        return sensitivity_list;
    }

    GraphTopology BuildFromWeights(
        const std::vector<GraphWeightedEdge> & weighted_edge_list,
        double minimum_weight)
    {
        GraphTopology topology;
        topology.adjacency_list.resize(m_atom_count);
        topology.sample_dependency_list = std::move(m_sample_dependency_list);
        topology.summary.candidate_edge_count = weighted_edge_list.size();

        std::vector<double> weight_list;
        weight_list.reserve(weighted_edge_list.size());
        for (const auto & weighted_edge : weighted_edge_list)
        {
            const auto weight{ weighted_edge.weight };
            weight_list.emplace_back(weight);
            if (weight < minimum_weight) continue;

            topology.retained_edge_list.emplace_back(weighted_edge);
            topology.summary.retained_edge_count++;
        }
        topology.summary.cut_edge_count = topology.summary.candidate_edge_count - topology.summary.retained_edge_count;
        topology.summary.weight_median = array_helper::ComputePercentile(weight_list, 0.5);
        topology.summary.weight_percentile_95 = array_helper::ComputePercentile(weight_list, 0.95);
        topology.summary.weight_maximum = weight_list.empty() ? 0.0 : std::ranges::max(weight_list);
        return topology;
    }

public:
    explicit CouplingGraphBuilder(std::size_t atom_count)
        : m_atom_count{ atom_count },
          m_self_gram_list(atom_count, Eigen::Matrix3d::Zero())
    {
    }

    void AddSample(SampleRef sample_id, std::vector<GraphParticipant> & participant_list)
    {
        NormalizeParticipantList(participant_list);

        for (std::size_t i = 0; i < participant_list.size(); i++)
        {
            const auto & participant{ participant_list.at(i) };
            if (participant.atom_index >= m_atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling participant index is out of range.");
            }
            if (!participant.jacobian.allFinite())
            {
                m_has_invalid_jacobian = true;
            }
        }

        GraphSampleDependency dependency;
        dependency.sample_id = sample_id;
        dependency.contributor_atom_index_list.reserve(participant_list.size());
        for (const auto & participant : participant_list)
        {
            dependency.contributor_atom_index_list.emplace_back(participant.atom_index);
            if (!m_has_invalid_jacobian)
            {
                m_self_gram_list.at(participant.atom_index) +=
                    participant.jacobian * participant.jacobian.transpose();
            }
        }
        m_sample_dependency_list.emplace_back(std::move(dependency));

        for (std::size_t i = 0; i < participant_list.size(); i++)
        {
            for (std::size_t j = i + 1; j < participant_list.size(); j++)
            {
                const auto & left{ participant_list.at(i) };
                const auto & right{ participant_list.at(j) };
                const AtomPair pair{ left.atom_index, right.atom_index };
                auto pair_iter{ m_pair_accumulator_by_pair.find(pair) };
                if (pair_iter == m_pair_accumulator_by_pair.end())
                {
                    pair_iter = m_pair_accumulator_by_pair.emplace(pair, Eigen::Matrix3d::Zero()).first;
                }
                if (m_has_invalid_jacobian) continue;
                pair_iter->second += left.jacobian * right.jacobian.transpose();
            }
        }
    }

private:
    GraphTopology BuildWeightedOrBinary(
        double minimum_weight,
        const std::vector<double> & sensitivity_minimum_weight_list)
    {
        if (m_has_invalid_jacobian) return BuildBinary();

        std::vector<double> self_gram_norm_list(m_atom_count, 0.0);
        std::vector<char> self_gram_norm_is_built(m_atom_count, 0);
        const auto get_self_gram_norm = [&](std::size_t atom_index)
        {
            if (self_gram_norm_is_built.at(atom_index) != 0)
            {
                return self_gram_norm_list.at(atom_index);
            }
            const auto & self_gram{ m_self_gram_list.at(atom_index) };
            if (!self_gram.allFinite())
            {
                return std::numeric_limits<double>::quiet_NaN();
            }
            const auto norm{ FrobeniusNorm(self_gram) };
            self_gram_norm_is_built.at(atom_index) = 1;
            self_gram_norm_list.at(atom_index) = norm;
            return norm;
        };

        std::vector<GraphWeightedEdge> weighted_edge_list;
        weighted_edge_list.reserve(m_pair_accumulator_by_pair.size());
        for (const auto & [pair, cross_gram] : m_pair_accumulator_by_pair)
        {
            if (!cross_gram.allFinite()) return BuildBinary();
            const auto left_norm{ get_self_gram_norm(pair.first) };
            const auto right_norm{ get_self_gram_norm(pair.second) };
            if (!std::isfinite(left_norm) || !std::isfinite(right_norm)) return BuildBinary();
            double weight{ 0.0 };
            if (left_norm == 0.0 || right_norm == 0.0)
            {
                weighted_edge_list.emplace_back(
                    GraphWeightedEdge{ pair.first, pair.second, weight });
                continue;
            }
            const auto denominator{ std::sqrt(left_norm) * std::sqrt(right_norm) };
            if (!std::isfinite(denominator)) return BuildBinary();
            const auto raw_weight{ FrobeniusNorm(cross_gram) / denominator };
            if (!std::isfinite(raw_weight)) return BuildBinary();
            weight = std::clamp(raw_weight, 0.0, 1.0);
            weighted_edge_list.emplace_back(
                GraphWeightedEdge{ pair.first, pair.second, weight });
        }
        auto topology{ BuildFromWeights(weighted_edge_list, minimum_weight) };
        topology.summary.uses_weighted_graph = true;
        topology.summary.threshold_sensitivity_list =
            BuildThresholdSensitivity(weighted_edge_list, sensitivity_minimum_weight_list);
        return topology;
    }

    GraphTopology BuildBinary()
    {
        std::vector<GraphWeightedEdge> weighted_edge_list;
        weighted_edge_list.reserve(m_pair_accumulator_by_pair.size());
        for (const auto & entry : m_pair_accumulator_by_pair)
        {
            weighted_edge_list.emplace_back(
                GraphWeightedEdge{ entry.first.first, entry.first.second, 1.0 });
        }
        return BuildFromWeights(weighted_edge_list, 0.0);
    }

public:
    GraphTopology BuildTopology(
        std::vector<ResidueKey> residue_key_by_atom_index,
        const CouplingGraphOptions & options = {})
    {
        ValidateBuildOptions(options);
        if (residue_key_by_atom_index.size() != m_atom_count)
        {
            throw std::invalid_argument(
                "Local fitting coupling residue key count must match atom count.");
        }

        auto topology{ BuildWeightedOrBinary(
            options.minimum_weight,
            options.sensitivity_minimum_weight_list) };
        topology.residue_key_by_atom_index = std::move(residue_key_by_atom_index);
        topology.summary.configured_minimum_weight = options.minimum_weight;
        return ApplyGraphResidueCutoff(std::move(topology), options.maximum_residue_count);
    }
};

inline Eigen::Vector3d EvaluateCouplingGraphJacobian(
    const std::optional<TransformedModelInvariants> & invariants,
    double distance)
{
    static const Eigen::Vector3d invalid_jacobian{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    if (!invariants.has_value()) return invalid_jacobian;
    return EvaluateTransformedJacobian(*invariants, distance).value_or(invalid_jacobian);
}

inline GraphTopology BuildSecondStageGraphTopology(
    const SecondStageContext & context,
    const FitState & initial_state,
    bool quiet_mode)
{
    std::size_t total_sample_count{ 0 };
    for (const auto & atom_context : context)
    {
        total_sample_count += atom_context.raw_sampling_entries.size();
    }
    const std::size_t total_work{ total_sample_count + 1 };
    std::size_t completed_work{ 0 };
    const std::string progress_message{ " Build local-fitting coupling topology" };
    if (!quiet_mode)
    {
        Logger::ProgressPercent(completed_work, total_work, 50, progress_message);
    }

    CouplingGraphBuilder builder{ context.size() };
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, initial_state)
    };
    std::vector<std::optional<TransformedModelInvariants>> selected_model_invariants;
    selected_model_invariants.reserve(model_snapshot.selected.size());
    for (const auto & model : model_snapshot.selected)
    {
        selected_model_invariants.emplace_back(BuildTransformedModelInvariants(model));
    }
    std::vector<std::optional<TransformedModelInvariants>> unselected_model_invariants;
    unselected_model_invariants.reserve(model_snapshot.unselected.size());
    for (const auto & model : model_snapshot.unselected)
    {
        unselected_model_invariants.emplace_back(BuildTransformedModelInvariants(model));
    }
    std::vector<GraphParticipant> participant_list;
    participant_list.reserve(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto & atom_context{ context.at(i) };
        for (std::size_t j = 0; j < atom_context.raw_sampling_entries.size(); j++)
        {
            const auto & sample{ atom_context.raw_sampling_entries.at(j) };
            const auto target_jacobian{ EvaluateCouplingGraphJacobian(
                selected_model_invariants.at(i),
                static_cast<double>(sample.point.distance)) };
            participant_list.clear();
            participant_list.emplace_back(GraphParticipant{ i, target_jacobian });
            for (const auto & neighbor_atom_sample : atom_context.Neighbors(j))
            {
                const auto neighbor_jacobian{ neighbor_atom_sample.is_selected ?
                    EvaluateCouplingGraphJacobian(
                        selected_model_invariants.at(neighbor_atom_sample.atom_index),
                        neighbor_atom_sample.distance) :
                    EvaluateCouplingGraphJacobian(
                        unselected_model_invariants.at(neighbor_atom_sample.atom_index),
                        neighbor_atom_sample.distance) };
                if (neighbor_atom_sample.is_selected)
                {
                    participant_list.emplace_back(
                        GraphParticipant{
                            neighbor_atom_sample.atom_index,
                            neighbor_jacobian
                        });
                    continue;
                }
                const auto selected_group_id{
                    context.unselected_atom_list.at(neighbor_atom_sample.atom_index).selected_group_id
                };
                if (!selected_group_id.has_value()) continue;
                for (const auto selected_index :
                    context.selected_atom_index_list_by_group.at(*selected_group_id))
                {
                    participant_list.emplace_back(
                        GraphParticipant{ selected_index, neighbor_jacobian });
                }
            }
            builder.AddSample(SampleRef{ i, j }, participant_list);
            completed_work++;
            if (!quiet_mode)
            {
                Logger::ProgressPercent(completed_work, total_work, 50, progress_message);
            }
        }
    }

    std::vector<ResidueKey> residue_key_by_atom_index;
    residue_key_by_atom_index.reserve(context.size());
    for (const auto & atom_context : context)
    {
        residue_key_by_atom_index.emplace_back(
            atom_context.atom->GetChainID(),
            atom_context.atom->GetSequenceID());
    }
    const auto topology{ builder.BuildTopology(std::move(residue_key_by_atom_index)) };
    completed_work++;
    if (!quiet_mode)
    {
        Logger::ProgressPercent(completed_work, total_work, 50, progress_message);
    }
    return topology;
}

inline void LogGraphTopology(const GraphTopology & topology, bool quiet_mode)
{
    if (quiet_mode) return;

    const auto & summary{ topology.summary };
    if (!summary.uses_weighted_graph)
    {
        Logger::Log(LogLevel::Warning,
            "Weighted local-fitting coupling graph is unavailable; using binary connectivity.");
        if (Logger::GetLogLevel() >= LogLevel::Debug)
        {
            Logger::Log(LogLevel::Debug,
                "Local-fitting weighted threshold sensitivity is unavailable in binary fallback mode.");
        }
    }
    std::ostringstream message;
    message << "Local-fitting coupling graph mode = "
        << (summary.uses_weighted_graph ? "weighted" : "binary-fallback")
        << std::scientific << std::setprecision(2)
        << ", minimum weight = " << summary.configured_minimum_weight
        << ", candidate/retained/cut edges = "
        << summary.candidate_edge_count << "/"
        << summary.retained_edge_count << "/"
        << summary.cut_edge_count
        << ", weight p50/p95/max = "
        << summary.weight_median << "/"
        << summary.weight_percentile_95 << "/"
        << summary.weight_maximum
        << ", initial components/max atoms/ratio = "
        << summary.component_count << "/"
        << summary.maximum_component_size << "/"
        << std::fixed << std::setprecision(2)
        << summary.maximum_component_ratio << ".";
    Logger::Log(LogLevel::Info, message.str());

    const auto & residue_cutoff_summary{ topology.residue_cutoff_summary };
    std::ostringstream residue_cutoff_message;
    residue_cutoff_message
        << "Local-fitting residue cutoff: residues="
        << residue_cutoff_summary.residue_count
        << ", limit=" << residue_cutoff_summary.maximum_residue_count_limit
        << ", clusters=" << residue_cutoff_summary.cluster_count
        << ", max-residues=" << residue_cutoff_summary.maximum_residue_count
        << ", cutoff-edges=" << residue_cutoff_summary.cut_edge_count << ".";
    Logger::Log(LogLevel::Info, residue_cutoff_message.str());

    for (const auto & sensitivity : summary.threshold_sensitivity_list)
    {
        std::ostringstream sensitivity_message;
        sensitivity_message
            << std::scientific << std::setprecision(2)
            << "Coupling sensitivity: threshold=" << sensitivity.minimum_weight
            << ", retained/cut="
            << sensitivity.retained_edge_count << "/"
            << sensitivity.cut_edge_count
            << ", components/max-atoms/ratio="
            << sensitivity.component_count << "/"
            << sensitivity.maximum_component_size << "/"
            << std::fixed << std::setprecision(2)
            << sensitivity.maximum_component_ratio << ".";
        Logger::Log(LogLevel::Info, sensitivity_message.str());
    }
}

inline GraphTopology ApplyGraphResidueCutoff(
    GraphTopology topology,
    std::size_t maximum_residue_count)
{
    const auto atom_count{ topology.adjacency_list.size() };
    const auto & residue_key_by_atom_index{ topology.residue_key_by_atom_index };
    if (residue_key_by_atom_index.size() != atom_count)
    {
        throw std::invalid_argument(
            "Local fitting coupling residue key count must match atom count.");
    }
    if (maximum_residue_count == 0)
    {
        throw std::invalid_argument(
            "Local fitting coupling maximum residue count must be positive.");
    }

    std::map<ResidueKey, std::size_t> residue_index_by_key;
    for (const auto & residue_key : residue_key_by_atom_index)
    {
        residue_index_by_key.emplace(residue_key, 0);
    }
    std::size_t next_residue_index{ 0 };
    for (auto & entry : residue_index_by_key)
    {
        entry.second = next_residue_index++;
    }

    std::vector<std::size_t> residue_index_by_atom_index;
    residue_index_by_atom_index.reserve(atom_count);
    for (const auto & residue_key : residue_key_by_atom_index)
    {
        residue_index_by_atom_index.emplace_back(residue_index_by_key.at(residue_key));
    }

    using ResiduePair = std::pair<std::size_t, std::size_t>;
    std::map<ResiduePair, double> maximum_weight_by_residue_pair;
    for (const auto & edge : topology.retained_edge_list)
    {
        if (edge.left_atom_index >= atom_count || edge.right_atom_index >= atom_count)
        {
            throw std::invalid_argument(
                "Local fitting coupling retained edge index is out of range.");
        }
        if (!std::isfinite(edge.weight))
        {
            throw std::invalid_argument(
                "Local fitting coupling retained edge weight must be finite.");
        }
        const auto left_residue_index{
            residue_index_by_atom_index.at(edge.left_atom_index)
        };
        const auto right_residue_index{
            residue_index_by_atom_index.at(edge.right_atom_index)
        };
        if (left_residue_index == right_residue_index) continue;

        const auto residue_pair{ std::minmax(left_residue_index, right_residue_index) };
        auto weight_iter{ maximum_weight_by_residue_pair.find(residue_pair) };
        if (weight_iter == maximum_weight_by_residue_pair.end())
        {
            maximum_weight_by_residue_pair.emplace(residue_pair, edge.weight);
        }
        else
        {
            weight_iter->second = std::max(weight_iter->second, edge.weight);
        }
    }

    struct WeightedResidueEdge
    {
        ResiduePair residue_pair{};
        double weight{ 0.0 };
    };
    std::vector<WeightedResidueEdge> weighted_residue_edge_list;
    weighted_residue_edge_list.reserve(maximum_weight_by_residue_pair.size());
    for (const auto & [residue_pair, weight] : maximum_weight_by_residue_pair)
    {
        weighted_residue_edge_list.emplace_back(WeightedResidueEdge{ residue_pair, weight });
    }
    std::sort(
        weighted_residue_edge_list.begin(),
        weighted_residue_edge_list.end(),
        [](const auto & lhs, const auto & rhs)
        {
            if (lhs.weight != rhs.weight) return lhs.weight > rhs.weight;
            return lhs.residue_pair < rhs.residue_pair;
        });

    DisjointSet residue_component_set{ residue_index_by_key.size() };
    for (const auto & edge : weighted_residue_edge_list)
    {
        const auto left_root{ residue_component_set.Find(edge.residue_pair.first) };
        const auto right_root{ residue_component_set.Find(edge.residue_pair.second) };
        if (left_root == right_root) continue;
        if (residue_component_set.ComponentSize(left_root) + residue_component_set.ComponentSize(right_root) > maximum_residue_count)
        {
            continue;
        }
        residue_component_set.Merge(left_root, right_root);
    }

    topology.adjacency_list.assign(atom_count, {});
    std::size_t cut_edge_count{ 0 };
    for (const auto & edge : topology.retained_edge_list)
    {
        const auto left_residue_index{
            residue_index_by_atom_index.at(edge.left_atom_index)
        };
        const auto right_residue_index{
            residue_index_by_atom_index.at(edge.right_atom_index)
        };
        if (residue_component_set.Find(left_residue_index) !=
            residue_component_set.Find(right_residue_index))
        {
            cut_edge_count++;
            continue;
        }
        topology.adjacency_list.at(edge.left_atom_index).emplace_back(edge.right_atom_index);
        topology.adjacency_list.at(edge.right_atom_index).emplace_back(edge.left_atom_index);
    }

    const auto residue_component_summary{
        SummarizeDisjointSetComponents(residue_component_set, residue_index_by_key.size())
    };
    topology.residue_cutoff_summary =
        GraphTopology::ResidueCutoffSummary{
            residue_index_by_key.size(),
            residue_component_summary.maximum_component_size,
            residue_component_summary.component_count,
            cut_edge_count,
            maximum_residue_count
        };
    UpdateGraphComponentSummary(topology);
    return topology;
}

inline void UpdateGraphComponentSummary(GraphTopology & topology)
{
    const auto atom_count{ topology.adjacency_list.size() };
    DisjointSet component_set{ atom_count };
    if (!topology.residue_key_by_atom_index.empty())
    {
        if (topology.residue_key_by_atom_index.size() != atom_count)
        {
            throw std::invalid_argument(
                "Local fitting coupling residue key count must match atom count.");
        }
        std::map<ResidueKey, std::size_t> first_atom_index_by_residue_key;
        for (std::size_t atom_index = 0; atom_index < atom_count; atom_index++)
        {
            const auto [iter, inserted]{
                first_atom_index_by_residue_key.emplace(
                    topology.residue_key_by_atom_index.at(atom_index),
                    atom_index)
            };
            if (!inserted)
            {
                component_set.Merge(iter->second, atom_index);
            }
        }
    }
    for (std::size_t atom_index = 0; atom_index < atom_count; atom_index++)
    {
        for (const auto neighbor_index : topology.adjacency_list.at(atom_index))
        {
            if (neighbor_index >= atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling edge index is out of range.");
            }
            component_set.Merge(atom_index, neighbor_index);
        }
    }

    const auto component_summary{
        SummarizeDisjointSetComponents(component_set, atom_count)
    };
    topology.summary.component_count = component_summary.component_count;
    topology.summary.maximum_component_size = component_summary.maximum_component_size;
    topology.summary.maximum_component_ratio = atom_count == 0 ? 0.0 :
        static_cast<double>(component_summary.maximum_component_size) / static_cast<double>(atom_count);
}

inline CouplingGraphPartition BuildGraphPartition(
    const GraphTopology & topology,
    const std::vector<std::size_t> & active_index_list)
{
    const auto atom_count{ topology.adjacency_list.size() };
    const auto inactive_position{ active_index_list.size() };
    std::vector<std::size_t> active_position_by_atom_index(atom_count, inactive_position);
    for (std::size_t position = 0; position < active_index_list.size(); position++)
    {
        const auto atom_index{ active_index_list.at(position) };
        if (atom_index >= atom_count)
        {
            throw std::invalid_argument("Local fitting coupling active index is out of range.");
        }
        if (active_position_by_atom_index.at(atom_index) != inactive_position)
        {
            throw std::invalid_argument("Local fitting coupling active indexes must be unique.");
        }
        active_position_by_atom_index.at(atom_index) = position;
    }

    DisjointSet component_set{ active_index_list.size() };

    if (!topology.residue_key_by_atom_index.empty())
    {
        if (topology.residue_key_by_atom_index.size() != atom_count)
        {
            throw std::invalid_argument(
                "Local fitting coupling residue key count must match atom count.");
        }
        std::map<ResidueKey, std::size_t> first_position_by_residue_key;
        for (std::size_t position = 0; position < active_index_list.size(); position++)
        {
            const auto atom_index{ active_index_list.at(position) };
            const auto & residue_key{
                topology.residue_key_by_atom_index.at(atom_index)
            };
            const auto [iter, inserted]{
                first_position_by_residue_key.emplace(residue_key, position)
            };
            if (!inserted) component_set.Merge(iter->second, position);
        }
    }

    for (const auto atom_index : active_index_list)
    {
        const auto position{ active_position_by_atom_index.at(atom_index) };
        for (const auto neighbor_index : topology.adjacency_list.at(atom_index))
        {
            if (neighbor_index >= atom_count)
            {
                throw std::invalid_argument("Local fitting coupling edge index is out of range.");
            }
            const auto neighbor_position{
                active_position_by_atom_index.at(neighbor_index)
            };
            if (neighbor_position != inactive_position)
            {
                component_set.Merge(position, neighbor_position);
            }
        }
    }

    std::map<std::size_t, ClusterKey> key_by_root;
    for (std::size_t position = 0; position < active_index_list.size(); position++)
    {
        key_by_root[component_set.Find(position)].emplace_back(active_index_list.at(position));
    }
    for (auto & entry : key_by_root)
    {
        std::ranges::sort(entry.second);
    }

    std::map<std::size_t, std::vector<SampleRef>> sample_id_list_by_root;
    CouplingGraphPartition partition;
    for (const auto & dependency : topology.sample_dependency_list)
    {
        std::vector<std::size_t> root_list;
        for (const auto atom_index : dependency.contributor_atom_index_list)
        {
            if (atom_index >= atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling sample contributor is out of range.");
            }
            const auto position{ active_position_by_atom_index.at(atom_index) };
            if (position == inactive_position) continue;
            root_list.emplace_back(component_set.Find(position));
        }
        std::ranges::sort(root_list);
        root_list.erase(std::ranges::unique(root_list).begin(), root_list.end());
        if (root_list.size() > 1) partition.boundary_sample_count++;
        for (const auto root : root_list)
        {
            sample_id_list_by_root[root].emplace_back(dependency.sample_id);
        }
    }

    for (auto & [root, key] : key_by_root)
    {
        partition.sample_id_list_by_key.emplace(std::move(key), std::move(sample_id_list_by_root[root]));
    }
    return partition;
}

inline std::vector<ClusterKey> BuildGraphClusterKeyList(const CouplingGraphPartition & partition)
{
    const auto key_view{ partition.sample_id_list_by_key | std::views::keys };
    return { key_view.begin(), key_view.end() };
}

inline std::vector<SampleRef> BuildGraphAffectedSampleUnion(
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & key_list)
{
    std::vector<SampleRef> sample_id_list;
    for (const auto & key : key_list)
    {
        const auto iter{ partition.sample_id_list_by_key.find(key) };
        if (iter == partition.sample_id_list_by_key.end()) continue;
        sample_id_list.insert(
            sample_id_list.end(),
            iter->second.begin(),
            iter->second.end());
    }
    std::ranges::sort(sample_id_list);
    sample_id_list.erase(
        std::ranges::unique(sample_id_list).begin(),
        sample_id_list.end());
    return sample_id_list;
}

} // namespace rhbm_gem::core::detail
