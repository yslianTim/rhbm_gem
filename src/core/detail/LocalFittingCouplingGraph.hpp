#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/algorithm/ClusteredAndersonAcceleration.hpp>

namespace rhbm_gem::core::detail {

struct LocalFittingCouplingSampleId
{
    std::size_t atom_index{ 0 };
    std::size_t sample_index{ 0 };
};

struct LocalFittingCouplingParticipant
{
    std::size_t atom_index{ 0 };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
};

struct LocalFittingCouplingSampleDependency
{
    LocalFittingCouplingSampleId sample_id{};
    std::vector<std::size_t> contributor_atom_index_list{};
};

struct LocalFittingCouplingEdge
{
    std::size_t neighbor_index{ 0 };
    double weight{ 0.0 };
};

struct LocalFittingCouplingGraphSummary
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
};

struct LocalFittingCouplingTopology
{
    std::vector<std::vector<LocalFittingCouplingEdge>> adjacency_list{};
    std::vector<LocalFittingCouplingSampleDependency> sample_dependency_list{};
    LocalFittingCouplingGraphSummary summary{};
};

struct LocalFittingCouplingPartition
{
    std::map<algorithm::ClusterKey, std::vector<LocalFittingCouplingSampleId>>
        sample_id_list_by_key{};
    std::size_t boundary_sample_count{ 0 };
};

class LocalFittingCouplingGraphBuilder
{
    using AtomPair = std::pair<std::size_t, std::size_t>;

    std::size_t m_atom_count{ 0 };
    std::vector<Eigen::Matrix3d> m_self_gram_list{};
    std::map<AtomPair, Eigen::Matrix3d> m_cross_gram_by_pair{};
    std::set<AtomPair> m_candidate_pair_set{};
    std::vector<LocalFittingCouplingSampleDependency> m_sample_dependency_list{};
    bool m_has_invalid_jacobian{ false };

    static double Percentile(std::vector<double> values, double quantile)
    {
        if (values.empty()) return 0.0;
        std::sort(values.begin(), values.end());
        const auto position{
            quantile * static_cast<double>(values.size() - 1)
        };
        const auto lower_index{ static_cast<std::size_t>(std::floor(position)) };
        const auto upper_index{ static_cast<std::size_t>(std::ceil(position)) };
        const auto fraction{ position - static_cast<double>(lower_index) };
        return values.at(lower_index) * (1.0 - fraction) +
            values.at(upper_index) * fraction;
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

    std::vector<LocalFittingCouplingGraphSummary::ThresholdSensitivity>
    BuildThresholdSensitivity(
        const std::map<AtomPair, double> & weight_by_pair,
        const std::vector<double> & minimum_weight_list) const
    {
        std::vector<LocalFittingCouplingGraphSummary::ThresholdSensitivity> sensitivity_list;
        sensitivity_list.reserve(minimum_weight_list.size());
        for (const auto minimum_weight : minimum_weight_list)
        {
            if (!std::isfinite(minimum_weight) ||
                minimum_weight < 0.0 || minimum_weight > 1.0)
            {
                throw std::invalid_argument(
                    "Local fitting coupling sensitivity minimum weight must be in [0, 1].");
            }

            std::vector<std::size_t> parent_list(m_atom_count);
            std::vector<std::size_t> component_size_list(m_atom_count, 1);
            for (std::size_t i = 0; i < parent_list.size(); i++) parent_list.at(i) = i;
            const auto find_root = [&](std::size_t index, auto && self) -> std::size_t
            {
                if (parent_list.at(index) == index) return index;
                parent_list.at(index) = self(parent_list.at(index), self);
                return parent_list.at(index);
            };
            const auto merge = [&](std::size_t left, std::size_t right)
            {
                const auto left_root{ find_root(left, find_root) };
                const auto right_root{ find_root(right, find_root) };
                if (left_root == right_root) return;
                parent_list.at(right_root) = left_root;
                component_size_list.at(left_root) += component_size_list.at(right_root);
            };

            std::size_t retained_edge_count{ 0 };
            for (const auto & pair : m_candidate_pair_set)
            {
                const auto weight_iter{ weight_by_pair.find(pair) };
                const auto weight{
                    weight_iter == weight_by_pair.end() ? 0.0 : weight_iter->second
                };
                if (weight < minimum_weight) continue;
                retained_edge_count++;
                merge(pair.first, pair.second);
            }

            std::size_t component_count{ 0 };
            std::size_t maximum_component_size{ 0 };
            for (std::size_t atom_index = 0; atom_index < m_atom_count; atom_index++)
            {
                if (find_root(atom_index, find_root) != atom_index) continue;
                component_count++;
                maximum_component_size = std::max(
                    maximum_component_size,
                    component_size_list.at(atom_index));
            }
            sensitivity_list.emplace_back(
                LocalFittingCouplingGraphSummary::ThresholdSensitivity{
                    minimum_weight,
                    retained_edge_count,
                    m_candidate_pair_set.size() - retained_edge_count,
                    component_count,
                    maximum_component_size,
                    m_atom_count == 0 ? 0.0 :
                        static_cast<double>(maximum_component_size) /
                            static_cast<double>(m_atom_count)
                });
        }
        return sensitivity_list;
    }

    LocalFittingCouplingTopology BuildFromWeights(
        const std::map<AtomPair, double> & weight_by_pair,
        double minimum_weight,
        bool uses_weighted_graph) const
    {
        LocalFittingCouplingTopology topology;
        topology.adjacency_list.resize(m_atom_count);
        topology.sample_dependency_list = m_sample_dependency_list;
        topology.summary.uses_weighted_graph = uses_weighted_graph;
        topology.summary.candidate_edge_count = m_candidate_pair_set.size();

        std::vector<double> weight_list;
        weight_list.reserve(m_candidate_pair_set.size());
        for (const auto & pair : m_candidate_pair_set)
        {
            const auto weight_iter{ weight_by_pair.find(pair) };
            const auto weight{
                weight_iter == weight_by_pair.end() ? 0.0 : weight_iter->second
            };
            weight_list.emplace_back(weight);
            if (weight < minimum_weight) continue;

            topology.adjacency_list.at(pair.first).emplace_back(
                LocalFittingCouplingEdge{ pair.second, weight });
            topology.adjacency_list.at(pair.second).emplace_back(
                LocalFittingCouplingEdge{ pair.first, weight });
            topology.summary.retained_edge_count++;
        }
        topology.summary.cut_edge_count =
            topology.summary.candidate_edge_count -
            topology.summary.retained_edge_count;
        topology.summary.weight_median = Percentile(weight_list, 0.5);
        topology.summary.weight_percentile_95 = Percentile(weight_list, 0.95);
        topology.summary.weight_maximum = weight_list.empty() ? 0.0 :
            *std::max_element(weight_list.begin(), weight_list.end());
        return topology;
    }

public:
    explicit LocalFittingCouplingGraphBuilder(std::size_t atom_count)
        : m_atom_count{ atom_count },
          m_self_gram_list(atom_count, Eigen::Matrix3d::Zero())
    {
    }

    void AddSample(
        LocalFittingCouplingSampleId sample_id,
        std::vector<LocalFittingCouplingParticipant> participant_list)
    {
        std::sort(
            participant_list.begin(),
            participant_list.end(),
            [](const auto & lhs, const auto & rhs)
            {
                return lhs.atom_index < rhs.atom_index;
            });
        for (std::size_t i = 0; i < participant_list.size(); i++)
        {
            const auto & participant{ participant_list.at(i) };
            if (participant.atom_index >= m_atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling participant index is out of range.");
            }
            if (i > 0 &&
                participant_list.at(i - 1).atom_index == participant.atom_index)
            {
                throw std::invalid_argument(
                    "Local fitting coupling sample participants must be unique.");
            }
            if (!participant.jacobian.allFinite())
            {
                m_has_invalid_jacobian = true;
            }
        }

        LocalFittingCouplingSampleDependency dependency;
        dependency.sample_id = sample_id;
        dependency.contributor_atom_index_list.reserve(participant_list.size());
        for (const auto & participant : participant_list)
        {
            dependency.contributor_atom_index_list.emplace_back(participant.atom_index);
            if (participant.jacobian.allFinite())
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
                m_candidate_pair_set.emplace(pair);
                if (!left.jacobian.allFinite() || !right.jacobian.allFinite()) continue;
                auto iter{ m_cross_gram_by_pair.find(pair) };
                if (iter == m_cross_gram_by_pair.end())
                {
                    iter = m_cross_gram_by_pair.emplace(
                        pair,
                        Eigen::Matrix3d::Zero()).first;
                }
                iter->second += left.jacobian * right.jacobian.transpose();
            }
        }
    }

    std::optional<LocalFittingCouplingTopology> BuildWeighted(
        double minimum_weight,
        const std::vector<double> & sensitivity_minimum_weight_list = {}) const
    {
        if (!std::isfinite(minimum_weight) ||
            minimum_weight < 0.0 || minimum_weight > 1.0)
        {
            throw std::invalid_argument(
                "Local fitting coupling minimum weight must be in [0, 1].");
        }
        if (m_has_invalid_jacobian) return std::nullopt;

        std::map<AtomPair, double> weight_by_pair;
        for (const auto & pair : m_candidate_pair_set)
        {
            const auto cross_iter{ m_cross_gram_by_pair.find(pair) };
            if (cross_iter == m_cross_gram_by_pair.end()) continue;
            const auto & left_self_gram{ m_self_gram_list.at(pair.first) };
            const auto & right_self_gram{ m_self_gram_list.at(pair.second) };
            if (!left_self_gram.allFinite() ||
                !right_self_gram.allFinite() ||
                !cross_iter->second.allFinite())
            {
                return std::nullopt;
            }
            const auto left_norm{ FrobeniusNorm(left_self_gram) };
            const auto right_norm{ FrobeniusNorm(right_self_gram) };
            if (left_norm == 0.0 || right_norm == 0.0) continue;
            const auto denominator{
                std::sqrt(left_norm) * std::sqrt(right_norm)
            };
            if (!std::isfinite(denominator)) return std::nullopt;
            const auto raw_weight{
                FrobeniusNorm(cross_iter->second) / denominator
            };
            if (!std::isfinite(raw_weight)) return std::nullopt;
            weight_by_pair.emplace(pair, std::clamp(raw_weight, 0.0, 1.0));
        }
        auto topology{ BuildFromWeights(weight_by_pair, minimum_weight, true) };
        topology.summary.threshold_sensitivity_list = BuildThresholdSensitivity(
            weight_by_pair,
            sensitivity_minimum_weight_list);
        return topology;
    }

    LocalFittingCouplingTopology BuildBinary() const
    {
        std::map<AtomPair, double> weight_by_pair;
        for (const auto & pair : m_candidate_pair_set)
        {
            weight_by_pair.emplace(pair, 1.0);
        }
        return BuildFromWeights(weight_by_pair, 0.0, false);
    }
};

inline LocalFittingCouplingPartition BuildLocalFittingCouplingPartition(
    const LocalFittingCouplingTopology & topology,
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

    std::vector<std::size_t> parent_list(active_index_list.size());
    for (std::size_t i = 0; i < parent_list.size(); i++) parent_list.at(i) = i;
    const auto find_root = [&](std::size_t index, auto && self) -> std::size_t
    {
        if (parent_list.at(index) == index) return index;
        parent_list.at(index) = self(parent_list.at(index), self);
        return parent_list.at(index);
    };
    const auto merge = [&](std::size_t left, std::size_t right)
    {
        const auto left_root{ find_root(left, find_root) };
        const auto right_root{ find_root(right, find_root) };
        if (left_root != right_root) parent_list.at(right_root) = left_root;
    };

    for (const auto atom_index : active_index_list)
    {
        const auto position{ *active_position_by_atom_index.at(atom_index) };
        for (const auto & edge : topology.adjacency_list.at(atom_index))
        {
            if (edge.neighbor_index >= atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling edge index is out of range.");
            }
            const auto neighbor_position{
                active_position_by_atom_index.at(edge.neighbor_index)
            };
            if (neighbor_position.has_value()) merge(position, *neighbor_position);
        }
    }

    std::map<std::size_t, algorithm::ClusterKey> key_by_root;
    for (std::size_t position = 0; position < active_index_list.size(); position++)
    {
        key_by_root[find_root(position, find_root)].emplace_back(active_index_list.at(position));
    }
    for (auto & [root, key] : key_by_root)
    {
        static_cast<void>(root);
        std::sort(key.begin(), key.end());
    }

    std::map<std::size_t, std::vector<LocalFittingCouplingSampleId>>
        sample_id_list_by_root;
    LocalFittingCouplingPartition partition;
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
            root_list.emplace_back(find_root(*position, find_root));
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
        partition.sample_id_list_by_key.emplace(
            std::move(key),
            std::move(sample_id_list_by_root[root]));
    }
    return partition;
}

} // namespace rhbm_gem::core::detail
