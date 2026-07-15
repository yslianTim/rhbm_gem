#include <gtest/gtest.h>

#include <vector>

#include <rhbm_gem/utils/algorithm/Convergence.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

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
