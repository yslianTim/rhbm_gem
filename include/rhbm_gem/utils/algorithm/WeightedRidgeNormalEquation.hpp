#pragma once

#include <cmath>
#include <utility>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/WeightedRidgeSystem.hpp>

namespace rhbm_gem::algorithm {

struct WeightedRidgeNormalEquation
{
    Eigen::SparseMatrix<double> normal_matrix;
    Eigen::VectorXd right_hand_side;
};

inline WeightedRidgeNormalEquation BuildWeightedRidgeNormalEquation(
    const WeightedRidgeSystem & system,
    const Eigen::VectorXd & weight)
{
    auto weighted_design{ system.design_matrix };
    for (Eigen::Index column = 0; column < weighted_design.outerSize(); column++)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator iter(weighted_design, column);
             iter;
             ++iter)
        {
            iter.valueRef() *= std::sqrt(weight(iter.row()));
        }
    }
    const auto weighted_response{ weight.array().sqrt().matrix().cwiseProduct(system.response) };
    Eigen::SparseMatrix<double> normal_matrix{
        weighted_design.transpose() * weighted_design
    };
    const auto column_count{ normal_matrix.cols() };
    for (Eigen::Index i = 0; i < column_count; i++)
    {
        normal_matrix.coeffRef(i, i) += system.ridge_diagonal(i);
    }
    normal_matrix.makeCompressed();
    Eigen::VectorXd right_hand_side{
        weighted_design.transpose() * weighted_response +
        system.ridge_diagonal.cwiseProduct(system.previous_parameter)
    };

    return { std::move(normal_matrix), std::move(right_hand_side) };
}

} // namespace rhbm_gem::algorithm
