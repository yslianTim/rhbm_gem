#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

alg::WeightedRidgeSystem MakeSingleParameterSystem(
    const std::vector<double> & basis_list,
    const std::vector<double> & response_list,
    double ridge,
    double previous_parameter = 0.0)
{
    alg::WeightedRidgeSystem system;
    system.design_matrix.resize(
        static_cast<Eigen::Index>(basis_list.size()),
        1);
    std::vector<Eigen::Triplet<double>> triplet_list;
    triplet_list.reserve(basis_list.size());
    for (std::size_t i = 0; i < basis_list.size(); i++)
    {
        triplet_list.emplace_back(static_cast<Eigen::Index>(i), 0, basis_list.at(i));
    }
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    system.response.resize(static_cast<Eigen::Index>(response_list.size()));
    for (std::size_t i = 0; i < response_list.size(); i++)
    {
        system.response(static_cast<Eigen::Index>(i)) = response_list.at(i);
    }
    system.previous_parameter = Eigen::VectorXd::Constant(1, previous_parameter);
    system.ridge_diagonal = Eigen::VectorXd::Constant(1, ridge);
    return system;
}

alg::WeightedRidgeSystem MakeTwoParameterSystem(
    double second_column_delta,
    double ridge)
{
    alg::WeightedRidgeSystem system;
    system.design_matrix.resize(2, 2);
    std::vector<Eigen::Triplet<double>> triplet_list{
        { 0, 0, 1.0 },
        { 0, 1, 1.0 },
        { 1, 0, 1.0 },
        { 1, 1, 1.0 + second_column_delta }
    };
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    system.response.resize(2);
    system.response << 0.0, 1.0;
    system.previous_parameter = Eigen::VectorXd::Zero(2);
    system.ridge_diagonal = Eigen::VectorXd::Constant(2, ridge);
    return system;
}

} // namespace

TEST(WeightedRidgeSolverTest, SolvesSingleParameterLeastSquares)
{
    const auto system{ MakeSingleParameterSystem({ 1.0, 2.0 }, { 2.0, 4.0 }, 0.0) };
    const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    alg::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd parameter;

    ASSERT_TRUE(solver.Solve(system, weight, parameter));

    ASSERT_EQ(1, parameter.size());
    EXPECT_NEAR(2.0, parameter(0), 1.0e-12);
}

TEST(WeightedRidgeSolverTest, WeightChangesSolution)
{
    const auto system{ MakeSingleParameterSystem({ 1.0, 1.0 }, { 0.0, 10.0 }, 0.0) };
    alg::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd first_weight(2);
    first_weight << 10.0, 1.0;
    Eigen::VectorXd second_weight(2);
    second_weight << 1.0, 10.0;
    Eigen::VectorXd first_parameter;
    Eigen::VectorXd second_parameter;

    ASSERT_TRUE(solver.Solve(system, first_weight, first_parameter));
    ASSERT_TRUE(solver.Solve(system, second_weight, second_parameter));

    EXPECT_NEAR(10.0 / 11.0, first_parameter(0), 1.0e-12);
    EXPECT_NEAR(100.0 / 11.0, second_parameter(0), 1.0e-12);
}

TEST(WeightedRidgeSolverTest, LargerRidgeKeepsSolutionCloserToPreviousParameter)
{
    const auto weak_ridge_system{ MakeSingleParameterSystem({ 1.0 }, { 10.0 }, 0.1, 2.0) };
    const auto strong_ridge_system{ MakeSingleParameterSystem({ 1.0 }, { 10.0 }, 10.0, 2.0) };
    const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(1) };
    alg::WeightedRidgeSolver solver{ weak_ridge_system };
    Eigen::VectorXd weak_ridge_parameter;
    Eigen::VectorXd strong_ridge_parameter;

    ASSERT_TRUE(solver.Solve(weak_ridge_system, weight, weak_ridge_parameter));
    ASSERT_TRUE(solver.Solve(strong_ridge_system, weight, strong_ridge_parameter));

    EXPECT_LT(
        std::abs(strong_ridge_parameter(0) - 2.0),
        std::abs(weak_ridge_parameter(0) - 2.0));
    EXPECT_NEAR((10.0 + 0.1 * 2.0) / 1.1, weak_ridge_parameter(0), 1.0e-12);
    EXPECT_NEAR((10.0 + 10.0 * 2.0) / 11.0, strong_ridge_parameter(0), 1.0e-12);
}

TEST(WeightedRidgeSolverTest, LargerRidgeSuppressesNearCollinearParameterMovement)
{
    const auto weak_ridge_system{ MakeTwoParameterSystem(1.0e-6, 1.0e-12) };
    const auto strong_ridge_system{ MakeTwoParameterSystem(1.0e-6, 10.0) };
    const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(2) };
    alg::WeightedRidgeSolver solver{ weak_ridge_system };
    Eigen::VectorXd weak_ridge_parameter;
    Eigen::VectorXd strong_ridge_parameter;

    ASSERT_TRUE(solver.Solve(weak_ridge_system, weight, weak_ridge_parameter));
    ASSERT_TRUE(solver.Solve(strong_ridge_system, weight, strong_ridge_parameter));

    EXPECT_GT(weak_ridge_parameter.norm(), 1.0e5);
    EXPECT_LT(strong_ridge_parameter.norm(), 0.2);
}

TEST(WeightedRidgeSolverTest, ReturnsFalseForSingularSystem)
{
    const auto system{ MakeSingleParameterSystem({ 0.0 }, { 1.0 }, 0.0) };
    const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    alg::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd parameter;

    EXPECT_FALSE(solver.Solve(system, weight, parameter));
}
