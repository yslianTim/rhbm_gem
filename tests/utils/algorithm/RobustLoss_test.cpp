#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;
}

TEST(RobustLossTest, HuberLossAndWeightMatchPiecewiseFormula)
{
    constexpr double cutoff{ 2.0 };

    EXPECT_DOUBLE_EQ(
        0.5,
        alg::CalculateRobustLoss(alg::RobustLossKind::Huber, 1.0, cutoff));
    EXPECT_DOUBLE_EQ(
        8.0,
        alg::CalculateRobustLoss(alg::RobustLossKind::Huber, 5.0, cutoff));

    EXPECT_DOUBLE_EQ(
        1.0,
        alg::CalculateRobustWeight(alg::RobustLossKind::Huber, 1.0, 1.0, cutoff));
    EXPECT_DOUBLE_EQ(
        0.4,
        alg::CalculateRobustWeight(alg::RobustLossKind::Huber, 5.0, 1.0, cutoff));
}

TEST(RobustLossTest, CauchyLossAndWeightMatchRedescendingFormula)
{
    constexpr double cutoff{ 2.0 };
    constexpr double residual{ 3.0 };
    const auto normalized_residual{ residual / cutoff };

    EXPECT_DOUBLE_EQ(
        0.5 * cutoff * cutoff * std::log1p(normalized_residual * normalized_residual),
        alg::CalculateRobustLoss(alg::RobustLossKind::Cauchy, residual, cutoff));
    EXPECT_DOUBLE_EQ(
        1.0 / (1.0 + normalized_residual * normalized_residual),
        alg::CalculateRobustWeight(alg::RobustLossKind::Cauchy, residual, 1.0, cutoff));

    EXPECT_GT(
        alg::CalculateRobustWeight(alg::RobustLossKind::Cauchy, 2.0, 1.0, cutoff),
        alg::CalculateRobustWeight(alg::RobustLossKind::Cauchy, 20.0, 1.0, cutoff));
}

TEST(RobustLossTest, CauchyHandlesVeryLargeResidualWithoutOverflow)
{
    const auto residual{ std::numeric_limits<double>::max() };

    EXPECT_TRUE(std::isfinite(
        alg::CalculateRobustLoss(alg::RobustLossKind::Cauchy, residual, 1.0)));
    EXPECT_DOUBLE_EQ(
        0.0,
        alg::CalculateRobustWeight(alg::RobustLossKind::Cauchy, residual, 1.0, 1.0));
}
