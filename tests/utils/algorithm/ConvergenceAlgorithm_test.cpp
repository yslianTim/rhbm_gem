#include <gtest/gtest.h>

#include <optional>
#include <vector>

#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

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
