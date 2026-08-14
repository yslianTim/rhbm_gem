#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include "core/detail/SecondStageIdentifiers.hpp"

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
    // The threshold applied to retain edges; binary fallback uses zero.
    double minimum_weight{ 0.0 };
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

inline GraphTopology ApplyGraphResidueCutoff(
    GraphTopology topology,
    std::vector<ResidueKey> residue_key_by_atom_index,
    std::size_t maximum_residue_count);

inline void UpdateGraphComponentSummary(GraphTopology & topology);

class CouplingGraphBuilder
{
    using AtomPair = std::pair<std::size_t, std::size_t>;

    struct WeightedPair
    {
        AtomPair pair{};
        double weight{ 0.0 };
    };

    std::size_t m_atom_count{ 0 };
    std::vector<Eigen::Matrix3d> m_self_gram_list{};
    std::map<AtomPair, Eigen::Matrix3d> m_pair_accumulator_by_pair{};
    std::vector<GraphSampleDependency> m_sample_dependency_list{};
    bool m_has_invalid_jacobian{ false };

    static void NormalizeParticipantList(
        std::vector<GraphParticipant> & participant_list)
    {
        std::sort(
            participant_list.begin(),
            participant_list.end(),
            [](const auto & lhs, const auto & rhs)
            {
                return lhs.atom_index < rhs.atom_index;
            });

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

    std::vector<CouplingGraphSummary::ThresholdSensitivity>
    BuildThresholdSensitivity(
        const std::vector<WeightedPair> & weighted_pair_list,
        const std::vector<double> & minimum_weight_list) const
    {
        std::vector<CouplingGraphSummary::ThresholdSensitivity> sensitivity_list;
        sensitivity_list.reserve(minimum_weight_list.size());
        for (const auto minimum_weight : minimum_weight_list)
        {
            if (!std::isfinite(minimum_weight) || minimum_weight < 0.0 || minimum_weight > 1.0)
            {
                throw std::invalid_argument(
                    "Local fitting coupling sensitivity minimum weight must be in [0, 1].");
            }

            DisjointSet component_set{ m_atom_count };

            std::size_t retained_edge_count{ 0 };
            for (const auto & weighted_pair : weighted_pair_list)
            {
                const auto weight{ weighted_pair.weight };
                if (weight < minimum_weight) continue;
                retained_edge_count++;
                component_set.Merge(weighted_pair.pair.first, weighted_pair.pair.second);
            }

            std::size_t component_count{ 0 };
            std::size_t maximum_component_size{ 0 };
            for (std::size_t atom_index = 0; atom_index < m_atom_count; atom_index++)
            {
                if (component_set.Find(atom_index) != atom_index) continue;
                component_count++;
                maximum_component_size = std::max(
                    maximum_component_size,
                    component_set.ComponentSize(atom_index));
            }
            sensitivity_list.emplace_back(
                CouplingGraphSummary::ThresholdSensitivity{
                    minimum_weight,
                    retained_edge_count,
                    weighted_pair_list.size() - retained_edge_count,
                    component_count,
                    maximum_component_size,
                    m_atom_count == 0 ? 0.0 :
                        static_cast<double>(maximum_component_size) / static_cast<double>(m_atom_count)
                });
        }
        return sensitivity_list;
    }

    GraphTopology BuildFromWeights(
        const std::vector<WeightedPair> & weighted_pair_list,
        double minimum_weight,
        bool uses_weighted_graph)
    {
        GraphTopology topology;
        topology.adjacency_list.resize(m_atom_count);
        topology.sample_dependency_list = std::move(m_sample_dependency_list);
        topology.summary.uses_weighted_graph = uses_weighted_graph;
        topology.summary.minimum_weight = minimum_weight;
        topology.summary.candidate_edge_count = weighted_pair_list.size();

        std::vector<double> weight_list;
        weight_list.reserve(weighted_pair_list.size());
        for (const auto & weighted_pair : weighted_pair_list)
        {
            const auto weight{ weighted_pair.weight };
            weight_list.emplace_back(weight);
            if (weight < minimum_weight) continue;

            topology.adjacency_list.at(weighted_pair.pair.first).emplace_back(weighted_pair.pair.second);
            topology.adjacency_list.at(weighted_pair.pair.second).emplace_back(weighted_pair.pair.first);
            topology.retained_edge_list.emplace_back(
                GraphWeightedEdge{
                    weighted_pair.pair.first,
                    weighted_pair.pair.second,
                    weight
                });
            topology.summary.retained_edge_count++;
        }
        topology.summary.cut_edge_count = topology.summary.candidate_edge_count - topology.summary.retained_edge_count;
        topology.summary.weight_median = array_helper::ComputePercentile(weight_list, 0.5);
        topology.summary.weight_percentile_95 = array_helper::ComputePercentile(weight_list, 0.95);
        topology.summary.weight_maximum = weight_list.empty() ? 0.0 :
            *std::max_element(weight_list.begin(), weight_list.end());
        return topology;
    }

    void AddSampleData(
        SampleRef sample_id,
        const std::vector<GraphParticipant> & participant_list)
    {
        bool sample_has_invalid_jacobian{ false };
        for (std::size_t i = 0; i < participant_list.size(); i++)
        {
            const auto & participant{ participant_list.at(i) };
            if (participant.atom_index >= m_atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling participant index is out of range.");
            }
            if (i > 0 && participant_list.at(i - 1).atom_index == participant.atom_index)
            {
                throw std::invalid_argument(
                    "Local fitting coupling sample participants must be unique.");
            }
            if (!participant.jacobian.allFinite())
            {
                sample_has_invalid_jacobian = true;
            }
        }
        m_has_invalid_jacobian = m_has_invalid_jacobian || sample_has_invalid_jacobian;

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
                    pair_iter = m_pair_accumulator_by_pair.emplace(
                        pair,
                        Eigen::Matrix3d::Zero()).first;
                }
                if (m_has_invalid_jacobian) continue;
                pair_iter->second += left.jacobian * right.jacobian.transpose();
            }
        }
    }

public:
    explicit CouplingGraphBuilder(std::size_t atom_count)
        : m_atom_count{ atom_count },
          m_self_gram_list(atom_count, Eigen::Matrix3d::Zero())
    {
    }

    void AddSample(
        SampleRef sample_id,
        std::vector<GraphParticipant> & participant_list)
    {
        NormalizeParticipantList(participant_list);
        AddSampleData(sample_id, participant_list);
    }

    void AddSample(
        SampleRef sample_id,
        std::vector<GraphParticipant> && participant_list)
    {
        NormalizeParticipantList(participant_list);
        AddSampleData(sample_id, participant_list);
    }

    void AddSortedSample(
        SampleRef sample_id,
        const std::vector<GraphParticipant> & participant_list)
    {
        auto normalized_list{ participant_list };
        NormalizeParticipantList(normalized_list);
        AddSampleData(sample_id, normalized_list);
    }

private:
    std::optional<GraphTopology> BuildWeighted(
        double minimum_weight,
        const std::vector<double> & sensitivity_minimum_weight_list = {})
    {
        if (!std::isfinite(minimum_weight) || minimum_weight < 0.0 || minimum_weight > 1.0)
        {
            throw std::invalid_argument(
                "Local fitting coupling minimum weight must be in [0, 1].");
        }
        if (m_has_invalid_jacobian) return std::nullopt;

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

        std::vector<WeightedPair> weighted_pair_list;
        weighted_pair_list.reserve(m_pair_accumulator_by_pair.size());
        for (const auto & [pair, cross_gram] : m_pair_accumulator_by_pair)
        {
            if (!cross_gram.allFinite())
            {
                return std::nullopt;
            }
            const auto left_norm{ get_self_gram_norm(pair.first) };
            const auto right_norm{ get_self_gram_norm(pair.second) };
            if (!std::isfinite(left_norm) || !std::isfinite(right_norm))
            {
                return std::nullopt;
            }
            double weight{ 0.0 };
            if (left_norm == 0.0 || right_norm == 0.0)
            {
                weighted_pair_list.emplace_back(WeightedPair{ pair, weight });
                continue;
            }
            const auto denominator{
                std::sqrt(left_norm) * std::sqrt(right_norm)
            };
            if (!std::isfinite(denominator)) return std::nullopt;
            const auto raw_weight{
                FrobeniusNorm(cross_gram) / denominator
            };
            if (!std::isfinite(raw_weight)) return std::nullopt;
            weight = std::clamp(raw_weight, 0.0, 1.0);
            weighted_pair_list.emplace_back(WeightedPair{ pair, weight });
        }
        auto sensitivity_list{ BuildThresholdSensitivity(weighted_pair_list, sensitivity_minimum_weight_list) };
        auto topology{ BuildFromWeights(weighted_pair_list, minimum_weight, true) };
        topology.summary.threshold_sensitivity_list = std::move(sensitivity_list);
        return topology;
    }

    GraphTopology BuildBinary()
    {
        std::vector<WeightedPair> weighted_pair_list;
        weighted_pair_list.reserve(m_pair_accumulator_by_pair.size());
        for (const auto & [pair, cross_gram] : m_pair_accumulator_by_pair)
        {
            static_cast<void>(cross_gram);
            weighted_pair_list.emplace_back(WeightedPair{ pair, 1.0 });
        }
        return BuildFromWeights(weighted_pair_list, 0.0, false);
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

        auto weighted_topology{ BuildWeighted(
            options.minimum_weight,
            options.sensitivity_minimum_weight_list) };
        auto topology{
            weighted_topology.has_value() ? std::move(*weighted_topology) : BuildBinary()
        };
        topology.summary.configured_minimum_weight = options.minimum_weight;
        topology = ApplyGraphResidueCutoff(
            std::move(topology),
            std::move(residue_key_by_atom_index),
            options.maximum_residue_count);
        return topology;
    }
};

inline GraphTopology ApplyGraphResidueCutoff(
    GraphTopology topology,
    std::vector<ResidueKey> residue_key_by_atom_index,
    std::size_t maximum_residue_count)
{
    const auto atom_count{ topology.adjacency_list.size() };
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
    for (auto & [residue_key, residue_index] : residue_index_by_key)
    {
        static_cast<void>(residue_key);
        residue_index = next_residue_index++;
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

        const auto residue_pair{
            std::minmax(left_residue_index, right_residue_index)
        };
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

    std::size_t cluster_count{ 0 };
    std::size_t observed_maximum_residue_count{ 0 };
    for (std::size_t residue_index = 0; residue_index < residue_index_by_key.size(); residue_index++)
    {
        if (residue_component_set.Find(residue_index) != residue_index) continue;
        cluster_count++;
        observed_maximum_residue_count = std::max(
            observed_maximum_residue_count,
            residue_component_set.ComponentSize(residue_index));
    }
    topology.residue_key_by_atom_index = std::move(residue_key_by_atom_index);
    topology.residue_cutoff_summary =
        GraphTopology::ResidueCutoffSummary{
            residue_index_by_key.size(),
            observed_maximum_residue_count,
            cluster_count,
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

    std::size_t component_count{ 0 };
    std::size_t maximum_component_size{ 0 };
    for (std::size_t atom_index = 0; atom_index < atom_count; atom_index++)
    {
        if (component_set.Find(atom_index) != atom_index) continue;
        component_count++;
        maximum_component_size = std::max(
            maximum_component_size,
            component_set.ComponentSize(atom_index));
    }
    topology.summary.component_count = component_count;
    topology.summary.maximum_component_size = maximum_component_size;
    topology.summary.maximum_component_ratio = atom_count == 0 ? 0.0 :
        static_cast<double>(maximum_component_size) / static_cast<double>(atom_count);
}

inline CouplingGraphPartition BuildGraphPartition(
    const GraphTopology & topology,
    const std::vector<std::size_t> & active_index_list)
{
    const auto atom_count{ topology.adjacency_list.size() };
    std::vector<std::optional<std::size_t>> active_position_by_atom_index(atom_count);
    for (std::size_t position = 0; position < active_index_list.size(); position++)
    {
        const auto atom_index{ active_index_list.at(position) };
        if (atom_index >= atom_count)
        {
            throw std::invalid_argument(
                "Local fitting coupling active index is out of range.");
        }
        if (active_position_by_atom_index.at(atom_index).has_value())
        {
            throw std::invalid_argument(
                "Local fitting coupling active indexes must be unique.");
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
        const auto position{ *active_position_by_atom_index.at(atom_index) };
        for (const auto neighbor_index : topology.adjacency_list.at(atom_index))
        {
            if (neighbor_index >= atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling edge index is out of range.");
            }
            const auto neighbor_position{
                active_position_by_atom_index.at(neighbor_index)
            };
            if (neighbor_position.has_value())
            {
                component_set.Merge(position, *neighbor_position);
            }
        }
    }

    std::map<std::size_t, ClusterKey> key_by_root;
    for (std::size_t position = 0; position < active_index_list.size(); position++)
    {
        key_by_root[component_set.Find(position)].emplace_back(active_index_list.at(position));
    }
    for (auto & [root, key] : key_by_root)
    {
        static_cast<void>(root);
        std::sort(key.begin(), key.end());
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
            if (!position.has_value()) continue;
            root_list.emplace_back(component_set.Find(*position));
        }
        std::sort(root_list.begin(), root_list.end());
        root_list.erase(std::unique(root_list.begin(), root_list.end()), root_list.end());
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

inline std::vector<ClusterKey>
BuildGraphClusterKeyList(
    const CouplingGraphPartition & partition)
{
    std::vector<ClusterKey> cluster_key_list;
    cluster_key_list.reserve(partition.sample_id_list_by_key.size());
    for (const auto & [key, sample_id_list] : partition.sample_id_list_by_key)
    {
        static_cast<void>(sample_id_list);
        cluster_key_list.emplace_back(key);
    }
    return cluster_key_list;
}

inline std::vector<SampleRef>
BuildGraphAffectedSampleUnion(
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
    std::sort(
        sample_id_list.begin(),
        sample_id_list.end(),
        [](const auto & lhs, const auto & rhs)
        {
            if (lhs.atom_index != rhs.atom_index)
            {
                return lhs.atom_index < rhs.atom_index;
            }
            return lhs.sample_index < rhs.sample_index;
        });
    sample_id_list.erase(
        std::unique(
            sample_id_list.begin(),
            sample_id_list.end(),
            [](const auto & lhs, const auto & rhs)
            {
                return lhs.atom_index == rhs.atom_index &&
                    lhs.sample_index == rhs.sample_index;
            }),
        sample_id_list.end());
    return sample_id_list;
}

} // namespace rhbm_gem::core::detail
