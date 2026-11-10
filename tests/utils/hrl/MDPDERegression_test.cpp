#include <gtest/gtest.h>

#include "utils/hrl/MDPDEEquations.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {
struct ShapeFixture
{
    rhbm_gem::RHBMMemberDataset dataset;
    double alpha{ 0.1 };
    rhbm_gem::RHBMExecutionOptions options;
    rhbm_gem::RHBMBetaEstimateResult expected;
};

ShapeFixture MakeFixture()
{
    ShapeFixture f;
    f.dataset.X.resize(40, 2);
    f.dataset.y.resize(40);
    for (int i = 0; i < 40; ++i)
    {
        const double r{ static_cast<double>(i) / 40.0 };
        f.dataset.X.row(i) << 1.0, -0.5 * r * r;
        f.dataset.y(i) = 1.2 - 2.0 * r * r + 0.01 * std::sin(3.0 * i + 0.2);
    }
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha, f.dataset, f.options);
    return f;
}
} // namespace

TEST(MDPDENumericsTest, FreshRobustEquationsMatchIndependentScalarCalculation)
{
    const auto f{ MakeFixture() };
    Eigen::Vector2d beta;
    beta << 1.19, 3.95;
    constexpr double variance{ 0.0001 }, floor{ 0.01 }, alpha{ 0.5 };
    const auto equations{ rhbm_gem::mdpde_detail::EvaluateMDPDEEquations(f.dataset, alpha, beta, variance, floor) };
    Eigen::Vector3d expected{ Eigen::Vector3d::Zero() };
    double denominator{};
    for (Eigen::Index i = 0; i < f.dataset.y.size(); ++i)
    {
        const double residual{ f.dataset.y(i) - f.dataset.X.row(i).dot(beta) };
        const double weight{ std::max(floor, std::exp(-alpha * residual * residual / (2 * variance))) };
        expected.head<2>() += f.dataset.X.row(i).transpose() * weight * residual;
        expected(2) += weight * (residual * residual / variance - 1) /
            static_cast<double>(f.dataset.y.size());
        denominator += weight;
        EXPECT_DOUBLE_EQ(equations.weights(i), weight);
    }
    const double correction{ alpha / std::pow(1 + alpha, 1.5) };
    expected(2) += correction;
    denominator -= static_cast<double>(f.dataset.y.size()) * correction;
    EXPECT_TRUE(equations.raw.isApprox(expected, 1e-12));
    EXPECT_NEAR(equations.denominator, denominator, 1e-12);
}

TEST(MDPDENumericsTest, BoundaryRankFloorAndDenominatorAreExplicit)
{
    auto f{ MakeFixture() };
    const auto beta{ f.expected.beta_ols };
    EXPECT_EQ(rhbm_gem::mdpde_detail::EvaluateMDPDEEquations(
        f.dataset, 0.1, beta, 0.0, 1e-8).reason, "variance-boundary");
    f.dataset.X.col(1) = f.dataset.X.col(0);
    EXPECT_EQ(rhbm_gem::mdpde_detail::EvaluateMDPDEEquations(
        f.dataset, 0.1, beta, 1.0, 1e-8).reason, "rank-deficient");
    f = MakeFixture();
    f.dataset.y.array() += 100.0;
    const auto equations{ rhbm_gem::mdpde_detail::EvaluateMDPDEEquations(
        f.dataset, 0.1, beta, 1e-6, 1e-8) };
    EXPECT_EQ(equations.floor_count, 40);
    EXPECT_EQ(equations.reason, "invalid-denominator");
}
