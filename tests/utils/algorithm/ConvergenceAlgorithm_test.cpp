#include <gtest/gtest.h>

#include <cstddef>
#include <optional>
#include <vector>

#include <rhbm_gem/utils/algorithm/AdaptiveRelaxationController.hpp>
#include <rhbm_gem/utils/algorithm/ConvergenceFreezeTracker.hpp>
#include <rhbm_gem/utils/algorithm/DependencyThawHysteresisTracker.hpp>
#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

alg::ParameterChange MakeChange(double value)
{
    return alg::ParameterChange{ std::vector<double>{ value } };
}

alg::FittingQualityCandidateStats MakeQualityCandidate(
    std::optional<double> quality_objective,
    double change = 0.1)
{
    return alg::FittingQualityCandidateStats{
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

TEST(ConvergenceAlgorithmTest, FreezeTrackerUpdatesOnlyActiveIndexes)
{
    alg::ConvergenceFreezeTracker tracker{
        2,
        1.0e-4,
        0.5,
        1
    };
    const std::vector<alg::ParameterChange> change_list{
        MakeChange(0.0),
        MakeChange(0.0)
    };

    tracker.Update(change_list, { 0 });

    EXPECT_TRUE(tracker.IsFrozen(0));
    EXPECT_FALSE(tracker.IsFrozen(1));
    const auto active_index_list{ tracker.BuildActiveIndexList() };
    ASSERT_EQ(1u, active_index_list.size());
    EXPECT_EQ(1u, active_index_list.front());
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

TEST(ConvergenceAlgorithmTest, DependencyThawHysteresisTracksThresholdGrowthAndDecay)
{
    alg::DependencyThawHysteresisTracker tracker{
        1,
        2.0,
        4.0,
        0.5,
        3
    };

    EXPECT_DOUBLE_EQ(10.0, tracker.GetThreshold(0, 10.0));
    EXPECT_FALSE(tracker.ShouldThaw(0, 9.9, 10.0));
    EXPECT_TRUE(tracker.ShouldThaw(0, 10.0, 10.0));

    tracker.RecordDependencyThaw(0);
    EXPECT_DOUBLE_EQ(20.0, tracker.GetThreshold(0, 10.0));
    tracker.RecordDependencyThaw(0);
    EXPECT_DOUBLE_EQ(40.0, tracker.GetThreshold(0, 10.0));
    tracker.RecordDependencyThaw(0);
    EXPECT_DOUBLE_EQ(40.0, tracker.GetThreshold(0, 10.0));

    tracker.DecayFrozen(0);
    EXPECT_DOUBLE_EQ(20.0, tracker.GetThreshold(0, 10.0));
    tracker.DecayFrozen(0);
    EXPECT_DOUBLE_EQ(10.0, tracker.GetThreshold(0, 10.0));
    tracker.DecayFrozen(0);
    EXPECT_DOUBLE_EQ(10.0, tracker.GetThreshold(0, 10.0));
}

TEST(ConvergenceAlgorithmTest, DependencyThawHysteresisCapsAndLocksThawCount)
{
    alg::DependencyThawHysteresisTracker tracker{
        1,
        2.0,
        4.0,
        0.5,
        2
    };

    EXPECT_TRUE(tracker.CanDependencyThaw(0));
    tracker.RecordDependencyThaw(0);
    EXPECT_TRUE(tracker.CanDependencyThaw(0));
    tracker.RecordDependencyThaw(0);

    EXPECT_FALSE(tracker.CanDependencyThaw(0));
    EXPECT_FALSE(tracker.CanDependencyThaw(0));
}

TEST(ConvergenceAlgorithmTest, DependencyThawHysteresisRejectsInvalidSettings)
{
    EXPECT_THROW(
        alg::DependencyThawHysteresisTracker(1, 0.9, 2.0, 0.5, 1),
        std::invalid_argument);
    EXPECT_THROW(
        alg::DependencyThawHysteresisTracker(1, 2.0, 0.9, 0.5, 1),
        std::invalid_argument);
    EXPECT_THROW(
        alg::DependencyThawHysteresisTracker(1, 2.0, 4.0, -0.1, 1),
        std::invalid_argument);
    EXPECT_THROW(
        alg::DependencyThawHysteresisTracker(1, 2.0, 4.0, 1.1, 1),
        std::invalid_argument);
    EXPECT_THROW(
        alg::DependencyThawHysteresisTracker(1, 2.0, 4.0, 0.5, -1),
        std::invalid_argument);
}

TEST(ConvergenceAlgorithmTest, DependencyThawHysteresisRejectsOutOfRangeIndex)
{
    alg::DependencyThawHysteresisTracker tracker{
        1,
        2.0,
        4.0,
        0.5,
        1
    };

    EXPECT_THROW(tracker.GetThreshold(1, 10.0), std::invalid_argument);
    EXPECT_THROW(tracker.ShouldThaw(1, 10.0, 10.0), std::invalid_argument);
    EXPECT_THROW(tracker.CanDependencyThaw(1), std::invalid_argument);
    EXPECT_THROW(tracker.RecordDependencyThaw(1), std::invalid_argument);
    EXPECT_THROW(tracker.DecayFrozen(1), std::invalid_argument);
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
        10.0,
        alg::ParameterChangeStats{ std::vector<double>{ 0.1 } }
    };
    alg::FittingQualityCandidateStats better_quality{
        9.0,
        alg::ParameterChangeStats{ std::vector<double>{ 1.0 } }
    };
    alg::FittingQualityCandidateStats tied_quality_higher_change{
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

TEST(ConvergenceAlgorithmTest, FittingQualityCandidateRankingPrefersAvailableObjective)
{
    const auto with_objective{ MakeQualityCandidate(10.0, 1.0) };
    const auto without_objective{ MakeQualityCandidate(std::nullopt, 0.1) };

    EXPECT_TRUE(alg::IsBetterFittingQualityCandidate(
        with_objective,
        without_objective,
        1.0e-8));
    EXPECT_FALSE(alg::IsBetterFittingQualityCandidate(
        without_objective,
        with_objective,
        1.0e-8));
}

TEST(ConvergenceAlgorithmTest, FittingQualityAcceptanceAcceptsNonDeterioratedCandidate)
{
    const auto candidate{ MakeQualityCandidate(9.5) };
    const auto previous{ MakeQualityCandidate(10.0) };
    const auto best{ MakeQualityCandidate(9.6) };

    EXPECT_TRUE(alg::IsFittingQualityAcceptableForProgress(
        candidate,
        previous,
        &best,
        1.0e-3));
}

TEST(ConvergenceAlgorithmTest, FittingQualityAcceptanceRejectsPreviousDeterioration)
{
    const auto candidate{ MakeQualityCandidate(10.2) };
    const auto previous{ MakeQualityCandidate(10.0) };
    const auto best{ MakeQualityCandidate(9.5) };

    EXPECT_FALSE(alg::IsFittingQualityAcceptableForProgress(
        candidate,
        previous,
        &best,
        1.0e-3));
}

TEST(ConvergenceAlgorithmTest, FittingQualityAcceptanceRejectsBestDeterioration)
{
    const auto candidate{ MakeQualityCandidate(9.5) };
    const auto previous{ MakeQualityCandidate(10.0) };
    const auto best{ MakeQualityCandidate(9.0) };

    EXPECT_FALSE(alg::IsFittingQualityAcceptableForProgress(
        candidate,
        previous,
        &best,
        1.0e-3));
}

TEST(ConvergenceAlgorithmTest, FittingQualityAcceptanceAcceptsWhenNoObjectiveReferenceExists)
{
    const auto candidate{ MakeQualityCandidate(std::nullopt) };
    const auto previous{ MakeQualityCandidate(std::nullopt) };

    EXPECT_TRUE(alg::IsFittingQualityAcceptableForProgress(
        candidate,
        previous,
        nullptr,
        1.0e-3));
}

TEST(ConvergenceAlgorithmTest, FittingQualityAcceptanceRejectsMissingCandidateObjective)
{
    const auto candidate{ MakeQualityCandidate(std::nullopt) };
    const auto previous{ MakeQualityCandidate(10.0) };
    const auto best{ MakeQualityCandidate(9.5) };

    EXPECT_FALSE(alg::IsFittingQualityAcceptableForProgress(
        candidate,
        previous,
        &best,
        1.0e-3));
}
