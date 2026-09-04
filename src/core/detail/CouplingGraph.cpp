#include "core/detail/CouplingGraph.hpp"

#include "core/detail/GaussianModelOperations.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

using ResiduePair = std::pair<std::size_t, std::size_t>;

class DisjointSet
{
    std::vector<std::size_t> m_parent_list{};
    std::vector<std::size_t> m_component_size_list{};

public:
    explicit DisjointSet(std::size_t item_count)
        : m_parent_list(item_count), m_component_size_list(item_count, 1)
    {
        for (std::size_t i = 0; i < item_count; i++)
        {
            m_parent_list.at(i) = i;
        }
    }

    std::size_t Find(std::size_t index)
    {
        if (m_parent_list.at(index) == index) return index;
        m_parent_list.at(index) = Find(m_parent_list.at(index));
        return m_parent_list.at(index);
    }

    void Merge(std::size_t left, std::size_t right)
    {
        auto left_root{ Find(left) };
        auto right_root{ Find(right) };
        if (left_root == right_root) return;
        if (m_component_size_list.at(left_root) <
            m_component_size_list.at(right_root))
        {
            std::swap(left_root, right_root);
        }
        m_parent_list.at(right_root) = left_root;
        m_component_size_list.at(left_root) += m_component_size_list.at(right_root);
    }

    std::size_t ComponentSize(std::size_t index) { return m_component_size_list.at(Find(index)); }
};

struct DisjointSetComponentSummary
{
    std::size_t component_count{ 0 };
    std::size_t maximum_component_size{ 0 };
};

DisjointSetComponentSummary SummarizeDisjointSetComponents(
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

Eigen::Vector3d EvaluateCouplingGraphJacobian(
    const std::optional<TransformedModelInvariants> & invariants,
    double distance)
{
    static const Eigen::Vector3d invalid_jacobian{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    if (!invariants.has_value()) return invalid_jacobian;
    return EvaluateTransformedJacobian(*invariants, distance).value_or(invalid_jacobian);
}

void UpdateGraphComponentSummary(GraphTopology & topology)
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
                throw std::invalid_argument("Local fitting coupling edge index is out of range.");
            }
            if (atom_index >= neighbor_index) continue;
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

void NormalizeParticipantList(std::vector<GraphParticipant> & participant_list)
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
            accumulated_jacobian = Eigen::Vector3d::Constant(
                std::numeric_limits<double>::quiet_NaN());
            continue;
        }
        accumulated_jacobian += participant.jacobian;
        if (!accumulated_jacobian.allFinite())
        {
            accumulated_jacobian = Eigen::Vector3d::Constant(
                std::numeric_limits<double>::quiet_NaN());
        }
    }
    participant_list.resize(normalized_size);
}

void ValidateBuildOptions(const CouplingGraphOptions & options)
{
    if (!std::isfinite(options.minimum_weight) ||
        options.minimum_weight < 0.0 || options.minimum_weight > 1.0)
    {
        throw std::invalid_argument("Local fitting coupling minimum weight must be in [0, 1].");
    }
    if (options.retained_edge_minimum_weight.has_value() &&
        (!std::isfinite(*options.retained_edge_minimum_weight) ||
         *options.retained_edge_minimum_weight < 0.0 ||
         *options.retained_edge_minimum_weight > options.minimum_weight))
    {
        throw std::invalid_argument(
            "Local fitting coupling retained-edge minimum weight must be in [0, minimum weight].");
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

double FrobeniusNorm(const Eigen::Matrix3d & matrix)
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

void SortGraphWeightedEdges(std::vector<GraphWeightedEdge> & weighted_edge_list)
{
    std::ranges::sort(
        weighted_edge_list,
        {},
        [](const GraphWeightedEdge & edge)
        {
            return std::pair{ edge.left_atom_index, edge.right_atom_index };
        });
}

} // namespace

std::vector<CouplingGraphSummary::ThresholdSensitivity>
CouplingGraphBuilder::BuildThresholdSensitivity(
    const std::vector<GraphWeightedEdge> & weighted_edge_list,
    const std::vector<double> & minimum_weight_list) const
{
    std::vector<CouplingGraphSummary::ThresholdSensitivity> sensitivity_list;
    sensitivity_list.reserve(minimum_weight_list.size());
    if (minimum_weight_list.empty()) return sensitivity_list;

    struct ThresholdEntry
    {
        double minimum_weight{ 0.0 };
        std::size_t original_index{ 0 };
    };

    std::vector<ThresholdEntry> threshold_entry_list;
    threshold_entry_list.reserve(minimum_weight_list.size());
    for (std::size_t index = 0; index < minimum_weight_list.size(); index++)
    {
        threshold_entry_list.emplace_back(
            ThresholdEntry{ minimum_weight_list.at(index), index });
    }
    std::sort(
        threshold_entry_list.begin(),
        threshold_entry_list.end(),
        [](const auto & lhs, const auto & rhs)
        {
            if (lhs.minimum_weight != rhs.minimum_weight)
            {
                return lhs.minimum_weight > rhs.minimum_weight;
            }
            return lhs.original_index < rhs.original_index;
        });

    std::vector<double> threshold_list;
    std::vector<std::size_t> bucket_index_by_original_index(minimum_weight_list.size());
    for (const auto & entry : threshold_entry_list)
    {
        if (threshold_list.empty() || threshold_list.back() != entry.minimum_weight)
        {
            threshold_list.emplace_back(entry.minimum_weight);
        }
        bucket_index_by_original_index.at(entry.original_index) = threshold_list.size() - 1;
    }

    std::vector<std::vector<std::size_t>> edge_index_list_by_bucket(threshold_list.size());
    for (std::size_t edge_index = 0; edge_index < weighted_edge_list.size(); edge_index++)
    {
        const auto bucket_iter{ std::lower_bound(
            threshold_list.begin(),
            threshold_list.end(),
            weighted_edge_list.at(edge_index).weight,
            [](const double threshold, const double weight)
            {
                return threshold > weight;
            }) };
        if (bucket_iter == threshold_list.end()) continue;
        edge_index_list_by_bucket.at(
            static_cast<std::size_t>(bucket_iter - threshold_list.begin()))
            .emplace_back(edge_index);
    }

    struct ThresholdSnapshot
    {
        std::size_t retained_edge_count{ 0 };
        DisjointSetComponentSummary component_summary{};
    };
    std::vector<ThresholdSnapshot> snapshot_list(threshold_list.size());
    DisjointSet component_set{ m_atom_count };
    std::size_t retained_edge_count{ 0 };
    for (std::size_t bucket_index = 0; bucket_index < threshold_list.size(); bucket_index++)
    {
        for (const auto edge_index : edge_index_list_by_bucket.at(bucket_index))
        {
            const auto & weighted_edge{ weighted_edge_list.at(edge_index) };
            component_set.Merge(weighted_edge.left_atom_index, weighted_edge.right_atom_index);
        }
        retained_edge_count += edge_index_list_by_bucket.at(bucket_index).size();
        snapshot_list.at(bucket_index) = ThresholdSnapshot{
            retained_edge_count,
            SummarizeDisjointSetComponents(component_set, m_atom_count)
        };
    }

    for (std::size_t original_index = 0; original_index < minimum_weight_list.size(); original_index++)
    {
        const auto & snapshot{
            snapshot_list.at(bucket_index_by_original_index.at(original_index))
        };
        sensitivity_list.emplace_back(
            CouplingGraphSummary::ThresholdSensitivity{
                minimum_weight_list.at(original_index),
                snapshot.retained_edge_count,
                weighted_edge_list.size() - snapshot.retained_edge_count,
                snapshot.component_summary.component_count,
                snapshot.component_summary.maximum_component_size,
                m_atom_count == 0 ? 0.0 :
                    static_cast<double>(snapshot.component_summary.maximum_component_size) /
                    static_cast<double>(m_atom_count)
            });
    }
    return sensitivity_list;
}

GraphTopology CouplingGraphBuilder::BuildFromWeights(
    const std::vector<GraphWeightedEdge> & weighted_edge_list,
    const CouplingGraphOptions & options,
    const GraphTopology * previous_topology)
{
    if (previous_topology != nullptr && previous_topology->adjacency_list.size() != m_atom_count)
    {
        throw std::invalid_argument(
            "Previous local fitting coupling topology has an inconsistent atom count.");
    }

    std::unordered_set<AtomPair, GraphIndexPairHash> previous_edge_set;
    const bool use_previous_edge_set{
        previous_topology != nullptr &&
        options.retained_edge_minimum_weight.has_value()
    };
    if (use_previous_edge_set)
    {
        std::size_t previous_edge_count{ 0 };
        for (std::size_t atom_index = 0;
            atom_index < previous_topology->adjacency_list.size();
            atom_index++)
        {
            for (const auto neighbor_index :
                previous_topology->adjacency_list.at(atom_index))
            {
                if (atom_index < neighbor_index) previous_edge_count++;
            }
        }
        previous_edge_set.reserve(previous_edge_count);
        for (std::size_t atom_index = 0;
            atom_index < previous_topology->adjacency_list.size();
            atom_index++)
        {
            for (const auto neighbor_index :
                previous_topology->adjacency_list.at(atom_index))
            {
                if (atom_index >= neighbor_index) continue;
                previous_edge_set.emplace(atom_index, neighbor_index);
            }
        }
    }

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
        auto minimum_weight{ options.minimum_weight };
        if (use_previous_edge_set &&
            previous_edge_set.contains(AtomPair{
                weighted_edge.left_atom_index,
                weighted_edge.right_atom_index
            }))
        {
            minimum_weight = *options.retained_edge_minimum_weight;
        }
        if (weight < minimum_weight) continue;

        topology.retained_edge_list.emplace_back(weighted_edge);
    }
    topology.summary.retained_edge_count = topology.retained_edge_list.size();
    topology.summary.cut_edge_count = topology.summary.candidate_edge_count - topology.summary.retained_edge_count;
    topology.summary.weight_median = array_helper::ComputePercentile(weight_list, 0.5);
    topology.summary.weight_percentile_95 = array_helper::ComputePercentile(weight_list, 0.95);
    topology.summary.weight_maximum = weight_list.empty() ? 0.0 : std::ranges::max(weight_list);
    return topology;
}

CouplingGraphBuilder::CouplingGraphBuilder(std::size_t atom_count)
    : m_atom_count{ atom_count },
      m_self_gram_list(atom_count, Eigen::Matrix3d::Zero())
{
    m_pair_accumulator_by_pair.reserve(atom_count);
}

void CouplingGraphBuilder::AddSample(
    SampleRef sample_id,
    std::vector<GraphParticipant> & participant_list)
{
    NormalizeParticipantList(participant_list);

    for (const auto & participant : participant_list)
    {
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
        dependency.contributor_atom_index_list.emplace_back(
            participant.atom_index);
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
            auto [pair_iter, inserted]{ m_pair_accumulator_by_pair.try_emplace(pair) };
            if (inserted) pair_iter->second.setZero();
            if (m_has_invalid_jacobian) continue;
            pair_iter->second +=
                left.jacobian * right.jacobian.transpose();
        }
    }
}

GraphTopology CouplingGraphBuilder::BuildWeightedOrBinary(
    const CouplingGraphOptions & options,
    const GraphTopology * previous_topology)
{
    if (m_has_invalid_jacobian) return BuildBinary();

    std::vector<double> self_gram_scale_list(m_atom_count, 0.0);
    std::vector<char> self_gram_scale_is_built(m_atom_count, 0);
    const auto get_self_gram_scale = [&](std::size_t atom_index)
    {
        if (self_gram_scale_is_built.at(atom_index) != 0)
        {
            return self_gram_scale_list.at(atom_index);
        }
        const auto & self_gram{ m_self_gram_list.at(atom_index) };
        if (!self_gram.allFinite())
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const auto norm{ FrobeniusNorm(self_gram) };
        const auto scale{ std::sqrt(norm) };
        self_gram_scale_is_built.at(atom_index) = 1;
        self_gram_scale_list.at(atom_index) = scale;
        return scale;
    };

    std::vector<GraphWeightedEdge> weighted_edge_list;
    weighted_edge_list.reserve(m_pair_accumulator_by_pair.size());
    for (const auto & [pair, cross_gram] : m_pair_accumulator_by_pair)
    {
        if (!cross_gram.allFinite()) return BuildBinary();
        const auto left_scale{ get_self_gram_scale(pair.first) };
        const auto right_scale{ get_self_gram_scale(pair.second) };
        if (!std::isfinite(left_scale) || !std::isfinite(right_scale))
        {
            return BuildBinary();
        }
        if (left_scale == 0.0 || right_scale == 0.0)
        {
            weighted_edge_list.emplace_back(
                GraphWeightedEdge{ pair.first, pair.second, 0.0 });
            continue;
        }
        const auto denominator{ left_scale * right_scale };
        if (!std::isfinite(denominator)) return BuildBinary();
        const auto raw_weight{ FrobeniusNorm(cross_gram) / denominator };
        if (!std::isfinite(raw_weight)) return BuildBinary();
        const auto weight{ std::clamp(raw_weight, 0.0, 1.0) };
        weighted_edge_list.emplace_back(
            GraphWeightedEdge{ pair.first, pair.second, weight });
    }
    SortGraphWeightedEdges(weighted_edge_list);
    auto topology{ BuildFromWeights(weighted_edge_list, options, previous_topology) };
    topology.summary.uses_weighted_graph = true;
    topology.summary.threshold_sensitivity_list = BuildThresholdSensitivity(
        weighted_edge_list,
        options.sensitivity_minimum_weight_list);
    return topology;
}

GraphTopology CouplingGraphBuilder::BuildBinary()
{
    std::vector<GraphWeightedEdge> weighted_edge_list;
    weighted_edge_list.reserve(m_pair_accumulator_by_pair.size());
    for (const auto & entry : m_pair_accumulator_by_pair)
    {
        weighted_edge_list.emplace_back(GraphWeightedEdge{
            entry.first.first,
            entry.first.second,
            1.0
        });
    }
    SortGraphWeightedEdges(weighted_edge_list);
    CouplingGraphOptions options;
    options.minimum_weight = 0.0;
    return BuildFromWeights(weighted_edge_list, options, nullptr);
}

GraphTopology CouplingGraphBuilder::BuildTopology(
    std::vector<ResidueKey> residue_key_by_atom_index,
    const CouplingGraphOptions & options,
    const GraphTopology * previous_topology)
{
    ValidateBuildOptions(options);
    if (residue_key_by_atom_index.size() != m_atom_count)
    {
        throw std::invalid_argument(
            "Local fitting coupling residue key count must match atom count.");
    }

    auto topology{ BuildWeightedOrBinary(options, previous_topology) };
    topology.residue_key_by_atom_index = std::move(residue_key_by_atom_index);
    topology.summary.configured_minimum_weight = options.minimum_weight;
    return ApplyGraphResidueCutoff(
        std::move(topology),
        options.maximum_residue_count);
}

static GraphTopology BuildSecondStageGraphTopologyImpl(
    const SecondStageContext & context,
    const FitState & state,
    bool quiet_mode,
    const CouplingGraphOptions & options,
    const GraphTopology * previous_topology)
{
    std::size_t total_sample_count{ 0 };
    for (const auto & atom_context : context)
    {
        total_sample_count += atom_context.raw_sampling_entries.size();
    }
    const std::size_t total_work{ total_sample_count + 1 };
    std::size_t completed_work{ 0 };
    const std::string progress_message{
        previous_topology == nullptr ?
            " Build local-fitting coupling topology" :
            " Rebuild local-fitting coupling topology"
    };
    if (!quiet_mode)
    {
        Logger::ProgressPercent(completed_work, total_work, 50, progress_message);
    }

    CouplingGraphBuilder builder{ context.size() };
    if (state.size() != context.size())
    {
        throw std::invalid_argument("Second-stage node snapshot size is inconsistent.");
    }
    std::vector<std::optional<TransformedModelInvariants>> model_invariants;
    model_invariants.reserve(state.size());
    for (const auto & result : state)
    {
        model_invariants.emplace_back(BuildTransformedModelInvariants(result.mdpde.GetModel()));
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
                model_invariants.at(i),
                sample.point.distance) };
            participant_list.clear();
            participant_list.emplace_back(GraphParticipant{ i, target_jacobian });
            for (const auto & neighbor_atom_sample : atom_context.Neighbors(j))
            {
                const auto neighbor_jacobian{ EvaluateCouplingGraphJacobian(
                    model_invariants.at(neighbor_atom_sample.atom_index),
                    neighbor_atom_sample.distance) };
                participant_list.emplace_back(
                    GraphParticipant{
                        neighbor_atom_sample.atom_index,
                        neighbor_jacobian
                    });
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
    const auto topology{
        builder.BuildTopology(
            std::move(residue_key_by_atom_index),
            options,
            previous_topology)
    };
    completed_work++;
    if (!quiet_mode)
    {
        Logger::ProgressPercent(completed_work, total_work, 50, progress_message);
    }
    return topology;
}

GraphTopology BuildSecondStageGraphTopology(
    const SecondStageContext & context,
    const FitState & initial_state,
    bool quiet_mode)
{
    return BuildSecondStageGraphTopologyImpl(
        context,
        initial_state,
        quiet_mode,
        CouplingGraphOptions{},
        nullptr);
}

GraphTopology BuildAdaptiveSecondStageGraphTopology(
    const SecondStageContext & context,
    const FitState & accepted_state,
    const GraphTopology & previous_topology,
    bool quiet_mode)
{
    CouplingGraphOptions options;
    options.minimum_weight = 0.06;
    options.retained_edge_minimum_weight = 0.04;
    return BuildSecondStageGraphTopologyImpl(
        context,
        accepted_state,
        quiet_mode,
        options,
        &previous_topology);
}

void LogGraphTopology(const GraphTopology & topology, bool quiet_mode)
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

GraphTopology ApplyGraphResidueCutoff(GraphTopology topology, std::size_t maximum_residue_count)
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

    std::unordered_map<ResiduePair, double, GraphIndexPairHash> maximum_weight_by_residue_pair;
    std::size_t maximum_residue_pair_count{ 0 };
    if (residue_index_by_key.size() > 1)
    {
        const auto residue_count{ residue_index_by_key.size() };
        const auto residue_factor{ residue_count - 1 };
        maximum_residue_pair_count =
            residue_count > std::numeric_limits<std::size_t>::max() / residue_factor ?
                std::numeric_limits<std::size_t>::max() :
                residue_count * residue_factor / 2;
    }
    maximum_weight_by_residue_pair.reserve(std::min(
        topology.retained_edge_list.size(),
        maximum_residue_pair_count));
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
    std::vector<std::size_t> adjacency_degree_list(atom_count, 0);
    std::vector<const GraphWeightedEdge *> connected_edge_list;
    connected_edge_list.reserve(topology.retained_edge_list.size());
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
        connected_edge_list.emplace_back(&edge);
        adjacency_degree_list.at(edge.left_atom_index)++;
        adjacency_degree_list.at(edge.right_atom_index)++;
    }
    for (std::size_t atom_index = 0; atom_index < atom_count; atom_index++)
    {
        topology.adjacency_list.at(atom_index).reserve(adjacency_degree_list.at(atom_index));
    }
    for (const auto * edge : connected_edge_list)
    {
        topology.adjacency_list.at(edge->left_atom_index).emplace_back(edge->right_atom_index);
        topology.adjacency_list.at(edge->right_atom_index).emplace_back(edge->left_atom_index);
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

CouplingGraphPartition BuildGraphPartition(
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
        std::vector<std::size_t> active_contributor_atom_index_list;
        for (const auto atom_index : dependency.contributor_atom_index_list)
        {
            if (atom_index >= atom_count)
            {
                throw std::invalid_argument(
                    "Local fitting coupling sample contributor is out of range.");
            }
            const auto position{ active_position_by_atom_index.at(atom_index) };
            if (position == inactive_position) continue;
            active_contributor_atom_index_list.emplace_back(atom_index);
            root_list.emplace_back(component_set.Find(position));
        }
        std::ranges::sort(active_contributor_atom_index_list);
        active_contributor_atom_index_list.erase(
            std::ranges::unique(active_contributor_atom_index_list).begin(),
            active_contributor_atom_index_list.end());
        std::ranges::sort(root_list);
        root_list.erase(std::ranges::unique(root_list).begin(), root_list.end());
        if (root_list.size() > 1)
        {
            CouplingGraphPartition::BoundarySampleDependency boundary_dependency{
                dependency.sample_id,
                {},
                std::move(active_contributor_atom_index_list)
            };
            boundary_dependency.cluster_key_list.reserve(root_list.size());
            for (const auto root : root_list)
            {
                boundary_dependency.cluster_key_list.emplace_back(key_by_root.at(root));
            }
            partition.boundary_sample_dependency_list.emplace_back(std::move(boundary_dependency));
        }
        for (const auto root : root_list)
        {
            sample_id_list_by_root[root].emplace_back(dependency.sample_id);
        }
    }

    std::ranges::sort(
        partition.boundary_sample_dependency_list,
        {},
        &CouplingGraphPartition::BoundarySampleDependency::sample_id);

    for (auto & [root, key] : key_by_root)
    {
        partition.sample_id_list_by_key.emplace(std::move(key), std::move(sample_id_list_by_root[root]));
    }
    return partition;
}

std::vector<ClusterKey> BuildGraphClusterKeyList(const CouplingGraphPartition & partition)
{
    const auto key_view{ partition.sample_id_list_by_key | std::views::keys };
    return { key_view.begin(), key_view.end() };
}

std::vector<SampleRef> BuildGraphAffectedSampleUnion(
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

std::vector<BoundaryReconciliationComponent> BuildBoundaryReconciliationComponents(
    const SecondStageContext & context,
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & accepted_key_list)
{
    auto sorted_accepted_key_list{ accepted_key_list };
    std::ranges::sort(sorted_accepted_key_list);
    sorted_accepted_key_list.erase(
        std::ranges::unique(sorted_accepted_key_list).begin(),
        sorted_accepted_key_list.end());
    if (sorted_accepted_key_list.size() < 2 || partition.boundary_sample_dependency_list.empty())
    {
        return {};
    }

    std::map<ClusterKey, std::size_t> accepted_position_by_key;
    for (std::size_t position = 0; position < sorted_accepted_key_list.size(); position++)
    {
        const auto & key{ sorted_accepted_key_list.at(position) };
        if (!partition.sample_id_list_by_key.contains(key))
        {
            throw std::invalid_argument(
                "Boundary reconciliation accepted key is absent from the graph partition.");
        }
        accepted_position_by_key.emplace(key, position);
    }

    DisjointSet component_set{ sorted_accepted_key_list.size() };
    for (const auto & dependency : partition.boundary_sample_dependency_list)
    {
        std::vector<std::size_t> accepted_position_list;
        for (const auto & key : dependency.cluster_key_list)
        {
            const auto iter{ accepted_position_by_key.find(key) };
            if (iter != accepted_position_by_key.end())
            {
                accepted_position_list.emplace_back(iter->second);
            }
        }
        if (accepted_position_list.size() < 2) continue;
        const auto first_position{ accepted_position_list.front() };
        for (std::size_t i = 1; i < accepted_position_list.size(); i++)
        {
            component_set.Merge(first_position, accepted_position_list.at(i));
        }
    }

    std::map<std::size_t, std::vector<ClusterKey>> key_list_by_root;
    for (std::size_t position = 0; position < sorted_accepted_key_list.size(); position++)
    {
        if (component_set.ComponentSize(position) < 2) continue;
        key_list_by_root[component_set.Find(position)].emplace_back(sorted_accepted_key_list.at(position));
    }

    std::vector<BoundaryReconciliationComponent> component_list;
    component_list.reserve(key_list_by_root.size());
    for (auto & key_list : key_list_by_root | std::views::values)
    {
        std::size_t boundary_sample_count{ 0 };
        std::vector<std::size_t> interface_atom_index_list;
        std::vector<std::size_t> component_atom_index_list;
        for (const auto & key : key_list)
        {
            component_atom_index_list.insert(
                component_atom_index_list.end(),
                key.begin(),
                key.end());
        }
        std::ranges::sort(component_atom_index_list);
        for (const auto & dependency : partition.boundary_sample_dependency_list)
        {
            const auto accepted_contributor_count{
                std::ranges::count_if(
                    dependency.cluster_key_list,
                    [&](const auto & key)
                    {
                        return std::ranges::find(key_list, key) != key_list.end();
                    })
            };
            if (accepted_contributor_count <= 1) continue;
            boundary_sample_count++;
            for (const auto atom_index : dependency.contributor_atom_index_list)
            {
                if (std::ranges::binary_search(component_atom_index_list, atom_index))
                {
                    interface_atom_index_list.emplace_back(atom_index);
                }
            }
        }
        std::ranges::sort(interface_atom_index_list);
        interface_atom_index_list.erase(
            std::ranges::unique(interface_atom_index_list).begin(),
            interface_atom_index_list.end());
        for (const auto atom_index : component_atom_index_list)
        {
            if (atom_index >= context.size())
            {
                throw std::invalid_argument("Boundary reconciliation component atom is out of range.");
            }
        }
        component_list.emplace_back(BoundaryReconciliationComponent{
            key_list,
            BuildGraphAffectedSampleUnion(partition, key_list),
            interface_atom_index_list,
            std::move(interface_atom_index_list),
            boundary_sample_count
        });
    }
    std::ranges::sort(component_list, {}, &BoundaryReconciliationComponent::key_list);
    return component_list;
}

BoundaryReconciliationComponent ExpandBoundaryReconciliationHalo(
    const SecondStageContext & context,
    BoundaryReconciliationComponent component,
    std::size_t halo_depth)
{
    std::vector<std::size_t> component_atom_index_list;
    for (const auto & key : component.key_list)
    {
        component_atom_index_list.insert(
            component_atom_index_list.end(),
            key.begin(),
            key.end());
    }
    std::ranges::sort(component_atom_index_list);
    component_atom_index_list.erase(
        std::ranges::unique(component_atom_index_list).begin(),
        component_atom_index_list.end());

    for (const auto atom_index : component_atom_index_list)
    {
        if (atom_index >= context.size())
        {
            throw std::invalid_argument("Boundary halo component atom is out of range.");
        }
    }
    auto halo_atom_index_list{ component.interface_atom_index_list };
    std::ranges::sort(halo_atom_index_list);
    halo_atom_index_list.erase(
        std::ranges::unique(halo_atom_index_list).begin(),
        halo_atom_index_list.end());
    for (const auto atom_index : halo_atom_index_list)
    {
        if (!std::ranges::binary_search(component_atom_index_list, atom_index))
        {
            throw std::invalid_argument(
                "Boundary halo interface atom is absent from its component.");
        }
    }

    for (std::size_t depth = 0; depth < halo_depth; depth++)
    {
        const auto previous_size{ halo_atom_index_list.size() };
        auto expanded_atom_index_list{ halo_atom_index_list };
        for (const auto & sample_ref : component.affected_sample_ref_list)
        {
            if (sample_ref.atom_index >= context.size())
            {
                throw std::invalid_argument("Boundary halo sample owner is out of range.");
            }
            const auto & atom_context{ context.at(sample_ref.atom_index) };
            if (sample_ref.sample_index >= atom_context.raw_sampling_entries.size())
            {
                throw std::invalid_argument("Boundary halo sample index is out of range.");
            }
            std::vector<std::size_t> direct_participant_list;
            if (std::ranges::binary_search(component_atom_index_list, sample_ref.atom_index))
            {
                direct_participant_list.emplace_back(sample_ref.atom_index);
            }
            for (const auto & neighbor : atom_context.Neighbors(sample_ref.sample_index))
            {
                if (!std::ranges::binary_search(
                        component_atom_index_list,
                        neighbor.atom_index))
                {
                    continue;
                }
                direct_participant_list.emplace_back(neighbor.atom_index);
            }
            std::ranges::sort(direct_participant_list);
            direct_participant_list.erase(
                std::ranges::unique(direct_participant_list).begin(),
                direct_participant_list.end());
            const auto touches_halo{
                std::ranges::any_of(
                    direct_participant_list,
                    [&](const auto atom_index)
                    {
                        return std::ranges::binary_search(halo_atom_index_list, atom_index);
                    })
            };
            if (!touches_halo) continue;
            expanded_atom_index_list.insert(
                expanded_atom_index_list.end(),
                direct_participant_list.begin(),
                direct_participant_list.end());
        }
        std::ranges::sort(expanded_atom_index_list);
        expanded_atom_index_list.erase(
            std::ranges::unique(expanded_atom_index_list).begin(),
            expanded_atom_index_list.end());
        halo_atom_index_list = std::move(expanded_atom_index_list);
        if (halo_atom_index_list.size() == previous_size) break;
    }

    component.halo_atom_index_list = std::move(halo_atom_index_list);
    return component;
}

std::vector<DependencyPolishComponent> BuildUncutDependencyPolishComponents(
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const std::vector<ClusterKey> & owner_key_by_atom_index)
{
    const auto key_list{ BuildGraphClusterKeyList(partition) };
    if (key_list.empty()) return {};
    if (owner_key_by_atom_index.size() != topology.adjacency_list.size())
    {
        throw std::invalid_argument("Dependency polish owner mapping has an inconsistent atom count.");
    }

    const auto inactive_position{ key_list.size() };
    std::vector<std::size_t> key_position_by_atom_index(
        topology.adjacency_list.size(),
        inactive_position);
    for (std::size_t key_position = 0; key_position < key_list.size(); key_position++)
    {
        for (const auto atom_index : key_list.at(key_position))
        {
            if (atom_index >= key_position_by_atom_index.size())
            {
                throw std::invalid_argument("Dependency polish cluster atom is out of range.");
            }
            key_position_by_atom_index.at(atom_index) = key_position;
        }
    }

    DisjointSet component_set{ key_list.size() };
    for (const auto & dependency : topology.sample_dependency_list)
    {
        std::vector<std::size_t> participant_key_position_list;
        for (const auto atom_index : dependency.contributor_atom_index_list)
        {
            if (atom_index >= key_position_by_atom_index.size())
            {
                throw std::invalid_argument("Dependency polish sample contributor is out of range.");
            }
            const auto key_position{ key_position_by_atom_index.at(atom_index) };
            if (key_position == inactive_position || owner_key_by_atom_index.at(atom_index).empty())
            {
                continue;
            }
            participant_key_position_list.emplace_back(key_position);
        }
        std::ranges::sort(participant_key_position_list);
        participant_key_position_list.erase(
            std::ranges::unique(participant_key_position_list).begin(),
            participant_key_position_list.end());
        if (participant_key_position_list.empty()) continue;
        const auto first_position{ participant_key_position_list.front() };
        for (std::size_t i = 1; i < participant_key_position_list.size(); i++)
        {
            component_set.Merge(first_position, participant_key_position_list.at(i));
        }
    }

    std::map<std::size_t, std::vector<ClusterKey>> key_list_by_root;
    for (std::size_t key_position = 0; key_position < key_list.size(); key_position++)
    {
        key_list_by_root[component_set.Find(key_position)].emplace_back(key_list.at(key_position));
    }

    std::vector<DependencyPolishComponent> component_list;
    for (auto & component_key_list : key_list_by_root | std::views::values)
    {
        std::vector<std::size_t> atom_index_list;
        for (const auto & key : component_key_list)
        {
            atom_index_list.insert(
                atom_index_list.end(),
                key.begin(),
                key.end());
        }
        std::ranges::sort(atom_index_list);
        if (atom_index_list.size() < 2) continue;

        auto affected_sample_ref_list{
            BuildGraphAffectedSampleUnion(partition, component_key_list)
        };
        affected_sample_ref_list.erase(
            std::remove_if(
                affected_sample_ref_list.begin(),
                affected_sample_ref_list.end(),
                [&](const auto & sample_ref)
                {
                    return sample_ref.atom_index >= owner_key_by_atom_index.size() ||
                        owner_key_by_atom_index.at(sample_ref.atom_index).empty();
                }),
            affected_sample_ref_list.end());
        if (affected_sample_ref_list.empty()) continue;
        component_list.emplace_back(DependencyPolishComponent{
            std::move(component_key_list),
            std::move(affected_sample_ref_list),
            std::move(atom_index_list)
        });
    }
    std::ranges::sort(component_list, {}, &DependencyPolishComponent::key_list);
    return component_list;
}

} // namespace rhbm_gem::core::detail
