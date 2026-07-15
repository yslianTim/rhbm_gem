#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;
}

TEST(RobustLossTest, CauchyLossAndWeightMatchRedescendingFormula)
{
    constexpr double cutoff{ 2.0 };
    constexpr double residual{ 3.0 };
    const auto normalized_residual{ residual / cutoff };

    EXPECT_DOUBLE_EQ(
        0.5 * cutoff * cutoff * std::log1p(normalized_residual * normalized_residual),
        alg::CalculateCauchyLoss(residual, cutoff));
    EXPECT_DOUBLE_EQ(
        1.0 / (1.0 + normalized_residual * normalized_residual),
        alg::CalculateCauchyWeight(residual, 1.0, cutoff));

    EXPECT_GT(
        alg::CalculateCauchyWeight(2.0, 1.0, cutoff),
        alg::CalculateCauchyWeight(20.0, 1.0, cutoff));
}

TEST(RobustLossTest, CauchyHandlesVeryLargeResidualWithoutOverflow)
{
    const auto residual{ std::numeric_limits<double>::max() };

    EXPECT_TRUE(std::isfinite(
        alg::CalculateCauchyLoss(residual, 1.0)));
    EXPECT_DOUBLE_EQ(
        0.0,
        alg::CalculateCauchyWeight(residual, 1.0, 1.0));
}

TEST(RobustLossTest, CauchyRejectsInvalidScalesAndHandlesNonFiniteResidual)
{
    const auto infinity{ std::numeric_limits<double>::infinity() };

    EXPECT_THROW(alg::CalculateCauchyLoss(1.0, 0.0), std::invalid_argument);
    EXPECT_THROW(alg::CalculateCauchyWeight(1.0, 0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(alg::CalculateCauchyWeight(1.0, 1.0, infinity), std::invalid_argument);
    EXPECT_TRUE(std::isinf(alg::CalculateCauchyLoss(infinity, 1.0)));
    EXPECT_DOUBLE_EQ(0.0, alg::CalculateCauchyWeight(infinity, 1.0, 1.0));
}
