#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include <rhbm_gem/utils/algorithm/AdaptiveRelaxationController.hpp>
#include <rhbm_gem/utils/algorithm/ConvergenceFreezeTracker.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

alg::ParameterChange MakeChange(double value)
{
    return alg::ParameterChange{ std::vector<double>{ value } };
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
