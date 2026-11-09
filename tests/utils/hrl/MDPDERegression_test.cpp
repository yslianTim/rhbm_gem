#include <gtest/gtest.h>

#include "utils/hrl/MDPDEEndpointRefinement.hpp"
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

TEST(EndpointRefinementTest, FailedOnlyLeavesNativeSuccessUntouchedDespiteFreshResidual)
{
    auto f{ MakeFixture() };
    f.options.tolerance = 0.1;
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha, f.dataset, f.options);
    ASSERT_EQ(f.expected.status, rhbm_gem::RHBMEstimationStatus::SUCCESS);
    const auto equations{ rhbm_gem::mdpde_detail::EvaluateMDPDEEquations(
        f.dataset, f.alpha, f.expected.beta_mdpde, f.expected.sigma_square,
        f.options.data_weight_min) };
    ASSERT_GT(equations.scaled.lpNorm<Eigen::Infinity>(), 1e-8);
    const auto result{ rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(
        f.dataset, f.alpha, f.options, f.expected) };
    EXPECT_FALSE(result.refinement);
    EXPECT_EQ(result.status, f.expected.status);
    EXPECT_EQ(result.Qualification(), rhbm_gem::RHBMSolveQualification::NativeSuccess);
    EXPECT_TRUE((result.beta_mdpde.array() == f.expected.beta_mdpde.array()).all());
    EXPECT_DOUBLE_EQ(result.sigma_square, f.expected.sigma_square);
    EXPECT_EQ(result.diagnostics.iterations, f.expected.diagnostics.iterations);
    EXPECT_TRUE((result.data_weight.diagonal().array() == f.expected.data_weight.diagonal().array()).all());
    EXPECT_TRUE((result.data_covariance.diagonal().array() == f.expected.data_covariance.diagonal().array()).all());
}

TEST(EndpointRefinementTest, FailedOnlyDoesNotPromoteAnIneligibleRankDeficientEndpoint)
{
    auto f{ MakeFixture() };
    f.dataset.X.col(1) = f.dataset.X.col(0);
    auto endpoint{ f.expected };
    endpoint.status = rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    const auto result{ rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(
        f.dataset, f.alpha, f.options, endpoint) };
    ASSERT_TRUE(result.refinement);
    EXPECT_FALSE(result.refinement->accepted);
    EXPECT_EQ(result.refinement->reason, "rank-deficient");
    EXPECT_EQ(result.refinement->reference_updates, 0);
    EXPECT_EQ(result.status, endpoint.status);
    EXPECT_EQ(result.Qualification(), rhbm_gem::RHBMSolveQualification::Unqualified);
    EXPECT_TRUE((result.beta_mdpde.array() == endpoint.beta_mdpde.array()).all());
    EXPECT_DOUBLE_EQ(result.sigma_square, endpoint.sigma_square);
}

TEST(EndpointRefinementTest, FailedOnlyRefinementQualifiesGeneratedFailure)
{
    auto f{ MakeFixture() };
    f.options.max_iterations = 1;
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha, f.dataset, f.options);
    ASSERT_EQ(f.expected.status, rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED);
    const auto result{ rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(
        f.dataset, f.alpha, f.options, f.expected) };
    ASSERT_TRUE(result.refinement);
    ASSERT_TRUE(result.refinement->accepted) << result.refinement->reason;
    EXPECT_EQ(result.status, f.expected.status);
    EXPECT_EQ(result.Qualification(), rhbm_gem::RHBMSolveQualification::RefinedSuccess);
    EXPECT_LE(result.refinement->candidate_equation_evaluations,
        rhbm_gem::mdpde_detail::kEndpointEquationBudget);
    EXPECT_LE(*result.refinement->candidate_residual, 1e-8);
    EXPECT_LE(*result.refinement->reference_residual, 1e-10);
    EXPECT_TRUE(*result.refinement->floor_masks_equal);
}
