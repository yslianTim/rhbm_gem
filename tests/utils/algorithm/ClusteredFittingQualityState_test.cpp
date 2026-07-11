#include <gtest/gtest.h>

#include <rhbm_gem/utils/algorithm/ClusteredFittingQualityState.hpp>

#include <optional>
#include <vector>

namespace {
namespace alg = rhbm_gem::algorithm;

using Samples = std::vector<double>;
using StateSet = alg::ClusteredFittingQualityStateSet<Samples>;
using CandidateScore = alg::ClusteredFittingQualityCandidateScore<Samples>;
using InitialState = alg::ClusteredFittingQualityInitialState<Samples>;
using TrackedCandidate = alg::ClusteredFittingQualityTrackedCandidate<Samples>;

alg::ClusteredFittingQualityOptions MakeOptions()
{
    return alg::ClusteredFittingQualityOptions{
        2,
        1.0e-3,
        1.0e-8,
        1.0,
        4.0,
        2.0,
        0.5
    };
}

alg::FittingQualityCandidateStats MakeStats(double objective, double change)
{
    return alg::FittingQualityCandidateStats{
        objective,
        alg::ParameterChangeStats{ std::vector<double>{ change } }
    };
}

InitialState MakeInitialState(double objective)
{
    return InitialState{
        1.0,
        MakeStats(objective, 0.0),
        Samples{ objective }
    };
}

CandidateScore MakeScore(
    double objective,
    double scale_sample,
    const std::optional<TrackedCandidate> & best_candidate)
{
    CandidateScore score;
    score.has_objective_reference = true;
    score.candidate_stats = MakeStats(objective, 0.1);
    score.committed_quality_objective = objective;
    if (best_candidate.has_value())
    {
        score.best_candidate_stats = best_candidate->candidate_stats;
    }
    score.objective_samples = Samples{ objective };
    score.objective_scale_sample = scale_sample;
    return score;
}

} // namespace

TEST(ClusteredFittingQualityStateTest, AcceptsImprovedCandidateAndLocksReferenceAfterWarmup)
{
    StateSet state_set{ MakeOptions() };
    const std::vector<alg::ClusterKey> key_list{ { 0 } };
    state_set.Reconcile(
        key_list,
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(10.0);
        });

    EXPECT_FALSE(state_set.AllActiveReferencesLocked(key_list));

    EXPECT_TRUE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> & best_candidate)
            {
                return MakeScore(9.0, 1.1, best_candidate);
            }));

    EXPECT_TRUE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> & best_candidate)
            {
                return MakeScore(8.0, 1.2, best_candidate);
            }));

    EXPECT_TRUE(state_set.AllActiveReferencesLocked(key_list));
}

TEST(ClusteredFittingQualityStateTest, AcceptedCandidateCanCommitDifferentReferenceStats)
{
    StateSet state_set{ MakeOptions() };
    state_set.Reconcile(
        { { 0 } },
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(10.0);
        });

    EXPECT_TRUE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> & best_candidate)
            {
                auto score{ MakeScore(9.0, 1.0, best_candidate) };
                score.committed_quality_objective = 5.0;
                return score;
            }));

    EXPECT_FALSE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> & best_candidate)
            {
                return MakeScore(7.0, 1.0, best_candidate);
            }));
}

TEST(ClusteredFittingQualityStateTest, RejectsDeterioratedCandidateWithoutChangingState)
{
    StateSet state_set{ MakeOptions() };
    state_set.Reconcile(
        { { 0 } },
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(10.0);
        });

    EXPECT_FALSE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> & best_candidate)
            {
                return MakeScore(20.0, 1.0, best_candidate);
            }));

    EXPECT_TRUE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> & best_candidate)
            {
                return MakeScore(9.0, 1.0, best_candidate);
            }));
}

TEST(ClusteredFittingQualityStateTest, AcceptsCandidateWhenNoScaleReferenceExists)
{
    StateSet state_set{ MakeOptions() };
    state_set.Reconcile(
        { { 0 } },
        [](const alg::ClusterKey &)
        {
            return InitialState{
                std::nullopt,
                alg::FittingQualityCandidateStats{
                    std::nullopt,
                    alg::ParameterChangeStats{ std::vector<double>{ 0.0 } }
                },
                std::nullopt
            };
        });

    EXPECT_TRUE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> &)
            {
                CandidateScore score;
                score.candidate_stats.parameter_change_stats =
                    alg::ParameterChangeStats{ std::vector<double>{ 0.1 } };
                return score;
            }));
}

TEST(ClusteredFittingQualityStateTest, RejectsMissingCandidateObjectiveWhenReferenceExists)
{
    StateSet state_set{ MakeOptions() };
    state_set.Reconcile(
        { { 0 } },
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(10.0);
        });

    EXPECT_FALSE(
        state_set.TryCommitCandidate(
            { 0 },
            [](const alg::ScaleReferenceTracker &,
                alg::FittingQualityCandidateStats &,
                const std::optional<Samples> &,
                const std::optional<TrackedCandidate> &)
            {
                CandidateScore score;
                score.has_objective_reference = true;
                score.candidate_stats.parameter_change_stats =
                    alg::ParameterChangeStats{ std::vector<double>{ 0.1 } };
                return score;
            }));
}

TEST(ClusteredFittingQualityStateTest, ObjectiveRidgeIncreasesSaturatesAndDecreasesPerCluster)
{
    StateSet state_set{ MakeOptions() };
    state_set.Reconcile(
        { { 0 }, { 2 } },
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(10.0);
        });

    EXPECT_TRUE(state_set.IncreaseObjectiveRidge({ { 0 } }));
    {
        const auto multiplier_list{ state_set.BuildObjectiveRidgeMultiplierList(3) };
        EXPECT_DOUBLE_EQ(2.0, multiplier_list.at(0));
        EXPECT_DOUBLE_EQ(1.0, multiplier_list.at(1));
        EXPECT_DOUBLE_EQ(1.0, multiplier_list.at(2));
    }

    EXPECT_TRUE(state_set.IncreaseObjectiveRidge({ { 0 } }));
    EXPECT_FALSE(state_set.IncreaseObjectiveRidge({ { 0 } }));
    {
        const auto multiplier_list{ state_set.BuildObjectiveRidgeMultiplierList(3) };
        EXPECT_DOUBLE_EQ(4.0, multiplier_list.at(0));
    }

    state_set.DecreaseObjectiveRidge({ { 0 } });
    {
        const auto multiplier_list{ state_set.BuildObjectiveRidgeMultiplierList(3) };
        EXPECT_DOUBLE_EQ(2.0, multiplier_list.at(0));
    }
}

TEST(ClusteredFittingQualityStateTest, UpdatesAcceptedAndRejectedRidgesIndependently)
{
    StateSet state_set{ MakeOptions() };
    state_set.Reconcile(
        { { 0 }, { 1 } },
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(10.0);
        });
    ASSERT_TRUE(state_set.IncreaseObjectiveRidge({ { 0 }, { 1 } }));

    state_set.DecreaseObjectiveRidge({ { 0 } });
    EXPECT_TRUE(state_set.IncreaseObjectiveRidge({ { 1 } }));

    const auto multiplier_list{ state_set.BuildObjectiveRidgeMultiplierList(2) };
    EXPECT_DOUBLE_EQ(1.0, multiplier_list.at(0));
    EXPECT_DOUBLE_EQ(4.0, multiplier_list.at(1));
}

TEST(ClusteredFittingQualityStateTest, ReconcilePreservesUnchangedClusterAndDropsRemovedCluster)
{
    StateSet state_set{ MakeOptions() };
    state_set.Reconcile(
        { { 0 }, { 1 } },
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(10.0);
        });
    state_set.IncreaseObjectiveRidge({ { 0 }, { 1 } });
    state_set.Reconcile(
        { { 1 } },
        [](const alg::ClusterKey &)
        {
            return MakeInitialState(5.0);
        });

    const auto multiplier_list{ state_set.BuildObjectiveRidgeMultiplierList(2) };
    EXPECT_DOUBLE_EQ(1.0, multiplier_list.at(0));
    EXPECT_DOUBLE_EQ(2.0, multiplier_list.at(1));
}
