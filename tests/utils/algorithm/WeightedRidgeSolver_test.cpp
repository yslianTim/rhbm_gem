#include <gtest/gtest.h>

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSystem.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

alg::WeightedRidgeSystem MakeSingleParameterSystem(
    const std::vector<double> & basis_list,
    const std::vector<double> & response_list,
    double ridge)
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
    system.previous_parameter = Eigen::VectorXd::Zero(1);
    system.ridge_diagonal = Eigen::VectorXd::Constant(1, ridge);
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

TEST(WeightedRidgeSolverTest, ReturnsFalseForSingularSystem)
{
    const auto system{ MakeSingleParameterSystem({ 0.0 }, { 1.0 }, 0.0) };
    const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    alg::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd parameter;

    EXPECT_FALSE(solver.Solve(system, weight, parameter));
}
