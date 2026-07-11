#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/algorithm/ClusteredAndersonAcceleration.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>
#include <rhbm_gem/utils/algorithm/ScaleReferenceTracker.hpp>

namespace rhbm_gem::algorithm {

struct ClusteredFittingQualityOptions
{
    std::size_t scale_warmup_count{ 1 };
    double objective_relative_tolerance{ 1.0e-3 };
    double objective_tie_relative_tolerance{ 1.0e-8 };
    double ridge_multiplier_min{ 1.0 };
    double ridge_multiplier_max{ 10.0 };
    double ridge_growth{ 2.0 };
    double ridge_shrink{ 0.8 };
};

template <typename ObjectiveSamples>
struct ClusteredFittingQualityTrackedCandidate
{
    FittingQualityCandidateStats candidate_stats{};
    std::optional<ObjectiveSamples> objective_samples{};
};

template <typename ObjectiveSamples>
struct ClusteredFittingQualityInitialState
{
    std::optional<double> initial_scale_sample{};
    FittingQualityCandidateStats candidate_stats{};
    std::optional<ObjectiveSamples> objective_samples{};
};

template <typename ObjectiveSamples>
struct ClusteredFittingQualityCandidateScore
{
    bool has_objective_reference{ false };
    FittingQualityCandidateStats candidate_stats{};
    std::optional<double> committed_quality_objective{};
    std::optional<FittingQualityCandidateStats> best_candidate_stats{};
    std::optional<ObjectiveSamples> objective_samples{};
    std::optional<double> objective_scale_sample{};
};

template <typename ObjectiveSamples>
class ClusteredFittingQualityStateSet
{
public:
    using InitialState = ClusteredFittingQualityInitialState<ObjectiveSamples>;
    using CandidateScore = ClusteredFittingQualityCandidateScore<ObjectiveSamples>;
    using TrackedCandidate = ClusteredFittingQualityTrackedCandidate<ObjectiveSamples>;

private:
    struct ClusterState
    {
        ScaleReferenceTracker objective_scale_tracker;
        FittingQualityCandidateStats previous_candidate_stats{};
        std::optional<ObjectiveSamples> previous_objective_samples{};
        std::optional<TrackedCandidate> best_candidate{};
        double objective_ridge_multiplier{ 1.0 };

        ClusterState(
            const ClusteredFittingQualityOptions & options,
            InitialState initial_state)
            : objective_scale_tracker{
                  options.scale_warmup_count,
                  initial_state.initial_scale_sample
              },
              previous_candidate_stats{ std::move(initial_state.candidate_stats) },
              previous_objective_samples{ std::move(initial_state.objective_samples) },
              objective_ridge_multiplier{ options.ridge_multiplier_min }
        {
            if (previous_candidate_stats.quality_objective.has_value())
            {
                best_candidate = TrackedCandidate{
                    previous_candidate_stats,
                    previous_objective_samples
                };
            }
        }
    };

    ClusteredFittingQualityOptions m_options{};
    std::map<ClusterKey, ClusterState> m_state_by_key{};

    static void ValidateOptions(const ClusteredFittingQualityOptions & options)
    {
        if (options.scale_warmup_count == 0)
        {
            throw std::invalid_argument("Clustered fitting quality scale warmup count must be positive.");
        }
        if (!std::isfinite(options.objective_relative_tolerance) ||
            options.objective_relative_tolerance < 0.0)
        {
            throw std::invalid_argument("Clustered fitting quality objective tolerance must be finite and non-negative.");
        }
        if (!std::isfinite(options.objective_tie_relative_tolerance) ||
            options.objective_tie_relative_tolerance < 0.0)
        {
            throw std::invalid_argument("Clustered fitting quality tie tolerance must be finite and non-negative.");
        }
        if (!std::isfinite(options.ridge_multiplier_min) ||
            options.ridge_multiplier_min <= 0.0)
        {
            throw std::invalid_argument("Clustered fitting quality minimum ridge multiplier must be finite and positive.");
        }
        if (!std::isfinite(options.ridge_multiplier_max) ||
            options.ridge_multiplier_max < options.ridge_multiplier_min)
        {
            throw std::invalid_argument("Clustered fitting quality maximum ridge multiplier is invalid.");
        }
        if (!std::isfinite(options.ridge_growth) || options.ridge_growth < 1.0)
        {
            throw std::invalid_argument("Clustered fitting quality ridge growth must be finite and at least one.");
        }
        if (!std::isfinite(options.ridge_shrink) ||
            options.ridge_shrink <= 0.0 ||
            options.ridge_shrink > 1.0)
        {
            throw std::invalid_argument("Clustered fitting quality ridge shrink must be finite and in (0, 1].");
        }
    }

    double IncreaseObjectiveRidgeMultiplier(double multiplier) const
    {
        return std::min(
            m_options.ridge_multiplier_max,
            std::max(m_options.ridge_multiplier_min, multiplier) * m_options.ridge_growth);
    }

    double DecreaseObjectiveRidgeMultiplier(double multiplier) const
    {
        return std::max(
            m_options.ridge_multiplier_min,
            multiplier * m_options.ridge_shrink);
    }

public:
    explicit ClusteredFittingQualityStateSet(ClusteredFittingQualityOptions options = {})
        : m_options{ options }
    {
        ValidateOptions(m_options);
    }

    template <typename Initializer>
    void Reconcile(
        const std::vector<ClusterKey> & key_list,
        Initializer initializer)
    {
        std::map<ClusterKey, ClusterState> next_state_by_key;
        for (const auto & key : key_list)
        {
            auto iter{ m_state_by_key.find(key) };
            if (iter != m_state_by_key.end())
            {
                next_state_by_key.emplace(key, std::move(iter->second));
                continue;
            }

            next_state_by_key.emplace(
                key,
                ClusterState{ m_options, initializer(key) });
        }
        m_state_by_key = std::move(next_state_by_key);
    }

    template <typename Scorer>
    bool TryCommitCandidate(
        const ClusterKey & key,
        Scorer scorer)
    {
        auto iter{ m_state_by_key.find(key) };
        if (iter == m_state_by_key.end())
        {
            throw std::invalid_argument("Clustered fitting quality state is missing.");
        }
        auto & state{ iter->second };

        auto previous_candidate_stats_for_backtracking{ state.previous_candidate_stats };
        auto score{
            scorer(
                state.objective_scale_tracker,
                previous_candidate_stats_for_backtracking,
                state.previous_objective_samples,
                state.best_candidate)
        };

        if (score.has_objective_reference)
        {
            const auto * best_candidate_stats_for_backtracking{
                score.best_candidate_stats.has_value() ?
                    &*score.best_candidate_stats :
                    nullptr
            };
            if (!IsFittingQualityAcceptableForProgress(
                    score.candidate_stats,
                    previous_candidate_stats_for_backtracking,
                    best_candidate_stats_for_backtracking,
                    m_options.objective_relative_tolerance))
            {
                return false;
            }
        }

        if (state.best_candidate.has_value() &&
            score.best_candidate_stats.has_value())
        {
            state.best_candidate->candidate_stats = *score.best_candidate_stats;
        }
        auto committed_candidate_stats{ score.candidate_stats };
        committed_candidate_stats.quality_objective = score.committed_quality_objective;
        if (committed_candidate_stats.quality_objective.has_value() &&
            score.objective_scale_sample.has_value())
        {
            state.objective_scale_tracker.CommitScaleSample(
                *score.objective_scale_sample);
        }
        state.previous_candidate_stats = committed_candidate_stats;
        state.previous_objective_samples = std::move(score.objective_samples);
        if (committed_candidate_stats.quality_objective.has_value() &&
            (!state.best_candidate.has_value() ||
                IsBetterFittingQualityCandidate(
                    committed_candidate_stats,
                    state.best_candidate->candidate_stats,
                    m_options.objective_tie_relative_tolerance)))
        {
            state.best_candidate = TrackedCandidate{
                committed_candidate_stats,
                state.previous_objective_samples
            };
        }
        return true;
    }

    bool IncreaseObjectiveRidge(const std::vector<ClusterKey> & key_list)
    {
        bool increased{ false };
        for (const auto & key : key_list)
        {
            auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end()) continue;

            const auto previous_multiplier{ iter->second.objective_ridge_multiplier };
            iter->second.objective_ridge_multiplier =
                IncreaseObjectiveRidgeMultiplier(previous_multiplier);
            if (iter->second.objective_ridge_multiplier > previous_multiplier)
            {
                increased = true;
            }
        }
        return increased;
    }

    void DecreaseObjectiveRidge(const std::vector<ClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end()) continue;

            iter->second.objective_ridge_multiplier =
                DecreaseObjectiveRidgeMultiplier(iter->second.objective_ridge_multiplier);
        }
    }

    std::vector<double> BuildObjectiveRidgeMultiplierList(std::size_t atom_size) const
    {
        std::vector<double> multiplier_list(atom_size, m_options.ridge_multiplier_min);
        for (const auto & [key, state] : m_state_by_key)
        {
            for (const auto active_index : key)
            {
                if (active_index >= multiplier_list.size())
                {
                    throw std::invalid_argument("Clustered fitting quality atom index is out of range.");
                }
                multiplier_list.at(active_index) =
                    std::max(multiplier_list.at(active_index), state.objective_ridge_multiplier);
            }
        }
        return multiplier_list;
    }

    bool AllActiveReferencesLocked(const std::vector<ClusterKey> & key_list) const
    {
        for (const auto & key : key_list)
        {
            const auto iter{ m_state_by_key.find(key) };
            if (iter == m_state_by_key.end())
            {
                throw std::invalid_argument("Clustered fitting quality state is missing.");
            }
            const auto & tracker{ iter->second.objective_scale_tracker };
            if (tracker.GetCommittedReference().has_value() && !tracker.IsLocked()) return false;
        }
        return true;
    }
};

} // namespace rhbm_gem::algorithm
