#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include <rhbm_gem/utils/algorithm/AdaptiveRelaxationController.hpp>
#include <rhbm_gem/utils/algorithm/ConvergenceFreezeTracker.hpp>
#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

alg::ParameterChange MakeChange(double value)
{
    return alg::ParameterChange{ std::vector<double>{ value } };
}

alg::FittingQualityCandidateStats MakeQualityCandidate(
    bool has_quality_objective,
    double quality_objective,
    double change = 0.1)
{
    return alg::FittingQualityCandidateStats{
        has_quality_objective,
        quality_objective,
        alg::ParameterChangeStats{ std::vector<double>{ change } }
    };
}

} // namespace

TEST(ConvergenceAlgorithmTest, SummarizesParameterChangePercentiles)
{
    const std::vector<alg::ParameterChange> change_list{
        alg::ParameterChange{ std::vector<double>{ 1.0, 10.0 } },
        alg::ParameterChange{ std::vector<double>{ 3.0, 30.0 } },
        alg::ParameterChange{ std::vector<double>{ 5.0, 50.0 } }
    };

    const auto stats{ alg::SummarizeParameterChangeStats(change_list, { 0, 2 }, 0.5) };

    ASSERT_EQ(2, stats.percentile_list.size());
    EXPECT_DOUBLE_EQ(3.0, stats.percentile_list.at(0));
    EXPECT_DOUBLE_EQ(30.0, stats.percentile_list.at(1));
}

TEST(ConvergenceAlgorithmTest, FreezeTrackerFreezesStableEntries)
{
    alg::ConvergenceFreezeTracker tracker{
        2,
        1.0e-4,
        0.5,
        2
    };
    const std::vector<alg::ParameterChange> change_list{
        MakeChange(0.001),
        MakeChange(0.1)
    };

    tracker.Update(change_list, tracker.BuildActiveIndexList());
    EXPECT_EQ(0, tracker.GetFrozenCount());

    tracker.Update(change_list, tracker.BuildActiveIndexList());

    EXPECT_EQ(1, tracker.GetFrozenCount());
    EXPECT_EQ(1, tracker.GetActiveCount());
    const auto active_index_list{ tracker.BuildActiveIndexList() };
    ASSERT_EQ(1, active_index_list.size());
    EXPECT_EQ(1, active_index_list.at(0));
}

TEST(ConvergenceAlgorithmTest, FreezeTrackerThawsFrozenEntryAndResetsStableCount)
{
    alg::ConvergenceFreezeTracker tracker{
        1,
        1.0e-4,
        0.5,
        2
    };
    const std::vector<alg::ParameterChange> stable_change_list{
        MakeChange(0.001)
    };

    tracker.Update(stable_change_list, tracker.BuildActiveIndexList());
    tracker.Update(stable_change_list, tracker.BuildActiveIndexList());
    ASSERT_TRUE(tracker.IsFrozen(0));
    ASSERT_EQ(0u, tracker.GetActiveCount());

    EXPECT_TRUE(tracker.Thaw(0));
    EXPECT_FALSE(tracker.IsFrozen(0));
    EXPECT_EQ(1u, tracker.GetActiveCount());
    const auto active_index_list{ tracker.BuildActiveIndexList() };
    ASSERT_EQ(1u, active_index_list.size());
    EXPECT_EQ(0u, active_index_list.front());

    tracker.Update(stable_change_list, tracker.BuildActiveIndexList());

    EXPECT_FALSE(tracker.IsFrozen(0));
}

TEST(ConvergenceAlgorithmTest, AdaptiveRelaxationGrowsShrinksAndClamps)
{
    alg::AdaptiveRelaxationController controller{
        2.0,
        0.1,
        1.0,
        2.0,
        0.5,
        0.01,
        2
    };

    EXPECT_DOUBLE_EQ(1.0, controller.GetBeta());
    controller.Update(10.0);
    controller.Update(9.0);
    controller.Update(8.0);
    EXPECT_DOUBLE_EQ(1.0, controller.GetBeta());

    controller.Update(20.0);

    EXPECT_DOUBLE_EQ(0.5, controller.GetBeta());
}

TEST(ConvergenceAlgorithmTest, AdaptiveRelaxationManualShrinkClampsAtMinimum)
{
    alg::AdaptiveRelaxationController controller{
        0.8,
        0.1,
        1.0,
        2.0,
        0.5,
        0.01,
        2
    };

    EXPECT_FALSE(controller.IsAtMinimum());
    EXPECT_DOUBLE_EQ(0.4, controller.Shrink());
    EXPECT_FALSE(controller.IsAtMinimum());
    EXPECT_DOUBLE_EQ(0.2, controller.Shrink());
    EXPECT_DOUBLE_EQ(0.1, controller.Shrink());
    EXPECT_TRUE(controller.IsAtMinimum());
    EXPECT_DOUBLE_EQ(0.1, controller.Shrink());
}

TEST(ConvergenceAlgorithmTest, NormalizedVectorChangeHandlesLargeScaleParameters)
{
    Eigen::VectorXd previous{ Eigen::VectorXd::Constant(1, 1000.0) };
    Eigen::VectorXd current{ Eigen::VectorXd::Constant(1, 1001.0) };

    EXPECT_NEAR(
        alg::CalculateMaximumNormalizedVectorChange(current, previous, 1.0e-2),
        1.0 / 1001.0,
        1.0e-12);
}

TEST(ConvergenceAlgorithmTest, NormalizedChangeKeepsLargeAbsoluteMovementScaleRelative)
{
    EXPECT_LT(
        alg::CalculateNormalizedChange(1001.0, 1000.0, 1.0),
        1.0e-3);
    EXPECT_GE(
        alg::CalculateNormalizedChange(2.0, 1.0, 1.0),
        1.0e-3);
}

TEST(ConvergenceAlgorithmTest, SummarizesPerParameterNormalizedChangePercentiles)
{
    const std::vector<alg::ParameterChange> change_list{
        alg::ParameterChange{ std::vector<double>{ 0.0005, 0.0004, 0.0003 } },
        alg::ParameterChange{ std::vector<double>{ 0.0007, 0.0015, 0.0004 } },
        alg::ParameterChange{ std::vector<double>{ 0.0006, 0.0005, 0.0020 } }
    };

    const auto stats{ alg::SummarizeParameterChangeStats(change_list, { 0, 1, 2 }, 0.5) };

    ASSERT_EQ(3, stats.percentile_list.size());
    EXPECT_LT(stats.percentile_list.at(0), 1.0e-3);
    EXPECT_LT(stats.percentile_list.at(1), 1.0e-3);
    EXPECT_LT(stats.percentile_list.at(2), 1.0e-3);
}

TEST(ConvergenceAlgorithmTest, NormalizedVectorChangeKeepsSmallScaleMovementVisible)
{
    Eigen::VectorXd previous{ Eigen::VectorXd::Constant(1, 0.01) };
    Eigen::VectorXd current{ Eigen::VectorXd::Constant(1, 0.011) };

    EXPECT_NEAR(
        alg::CalculateMaximumNormalizedVectorChange(current, previous, 1.0e-2),
        0.001 / 0.011,
        1.0e-12);
}

TEST(ConvergenceAlgorithmTest, NormalizedVectorChangeUsesFloorNearZero)
{
    Eigen::VectorXd previous{ Eigen::VectorXd::Zero(1) };
    Eigen::VectorXd current{ Eigen::VectorXd::Constant(1, 1.0e-5) };

    EXPECT_NEAR(
        alg::CalculateMaximumNormalizedVectorChange(current, previous, 1.0e-2),
        1.0e-3,
        1.0e-12);
}

TEST(ConvergenceAlgorithmTest, FittingQualityCandidateRankingUsesChangeAsTieBreaker)
{
    alg::FittingQualityCandidateStats lower_change{
        true,
        10.0,
        alg::ParameterChangeStats{ std::vector<double>{ 0.1 } }
    };
    alg::FittingQualityCandidateStats better_quality{
        true,
        9.0,
        alg::ParameterChangeStats{ std::vector<double>{ 1.0 } }
    };
    alg::FittingQualityCandidateStats tied_quality_higher_change{
        true,
        10.0 + 1.0e-10,
        alg::ParameterChangeStats{ std::vector<double>{ 0.2 } }
    };

    EXPECT_TRUE(alg::IsBetterFittingQualityCandidate(
        better_quality,
        lower_change,
        1.0e-8));
    EXPECT_FALSE(alg::IsBetterFittingQualityCandidate(
        tied_quality_higher_change,
        lower_change,
        1.0e-8));
}

TEST(ConvergenceAlgorithmTest, FittingQualityBacktrackingAcceptsNonDeterioratedCandidate)
{
    const auto candidate{ MakeQualityCandidate(true, 9.5) };
    const auto previous{ MakeQualityCandidate(true, 10.0) };
    const auto best{ MakeQualityCandidate(true, 9.6) };

    const auto decision{ alg::EvaluateFittingQualityBacktracking(
        candidate,
        previous,
        true,
        best,
        1.0e-3,
        0,
        3) };

    EXPECT_TRUE(decision.accepted);
    EXPECT_FALSE(decision.should_shrink_beta);
    EXPECT_FALSE(decision.reached_retry_limit);
}

TEST(ConvergenceAlgorithmTest, FittingQualityBacktrackingRejectsPreviousDeterioration)
{
    const auto candidate{ MakeQualityCandidate(true, 10.2) };
    const auto previous{ MakeQualityCandidate(true, 10.0) };
    const auto best{ MakeQualityCandidate(true, 9.5) };

    const auto decision{ alg::EvaluateFittingQualityBacktracking(
        candidate,
        previous,
        true,
        best,
        1.0e-3,
        0,
        3) };

    EXPECT_FALSE(decision.accepted);
    EXPECT_TRUE(decision.should_shrink_beta);
    EXPECT_FALSE(decision.reached_retry_limit);
}

TEST(ConvergenceAlgorithmTest, FittingQualityBacktrackingRejectsBestDeterioration)
{
    const auto candidate{ MakeQualityCandidate(true, 9.5) };
    const auto previous{ MakeQualityCandidate(true, 10.0) };
    const auto best{ MakeQualityCandidate(true, 9.0) };

    const auto decision{ alg::EvaluateFittingQualityBacktracking(
        candidate,
        previous,
        true,
        best,
        1.0e-3,
        0,
        3) };

    EXPECT_FALSE(decision.accepted);
    EXPECT_TRUE(decision.should_shrink_beta);
    EXPECT_FALSE(decision.reached_retry_limit);
}

TEST(ConvergenceAlgorithmTest, FittingQualityBacktrackingReportsRetryLimit)
{
    const auto candidate{ MakeQualityCandidate(true, 10.2) };
    const auto previous{ MakeQualityCandidate(true, 10.0) };
    const auto best{ MakeQualityCandidate(true, 9.5) };

    const auto decision{ alg::EvaluateFittingQualityBacktracking(
        candidate,
        previous,
        true,
        best,
        1.0e-3,
        2,
        3) };

    EXPECT_FALSE(decision.accepted);
    EXPECT_FALSE(decision.should_shrink_beta);
    EXPECT_TRUE(decision.reached_retry_limit);
}

TEST(ConvergenceAlgorithmTest, FittingQualityBacktrackingAcceptsWhenNoObjectiveReferenceExists)
{
    const auto candidate{ MakeQualityCandidate(false, 0.0) };
    const auto previous{ MakeQualityCandidate(false, 0.0) };
    const auto best{ MakeQualityCandidate(false, 0.0) };

    const auto decision{ alg::EvaluateFittingQualityBacktracking(
        candidate,
        previous,
        false,
        best,
        1.0e-3,
        0,
        3) };

    EXPECT_TRUE(decision.accepted);
    EXPECT_FALSE(decision.should_shrink_beta);
    EXPECT_FALSE(decision.reached_retry_limit);
}

TEST(ConvergenceAlgorithmTest, FittingQualityBacktrackingRejectsMissingCandidateObjective)
{
    const auto candidate{ MakeQualityCandidate(false, 0.0) };
    const auto previous{ MakeQualityCandidate(true, 10.0) };
    const auto best{ MakeQualityCandidate(true, 9.5) };

    const auto decision{ alg::EvaluateFittingQualityBacktracking(
        candidate,
        previous,
        true,
        best,
        1.0e-3,
        0,
        3) };

    EXPECT_FALSE(decision.accepted);
    EXPECT_TRUE(decision.should_shrink_beta);
    EXPECT_FALSE(decision.reached_retry_limit);
}
