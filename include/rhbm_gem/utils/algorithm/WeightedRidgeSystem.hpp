#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace rhbm_gem::algorithm {

struct WeightedRidgeSystem
{
    Eigen::SparseMatrix<double> design_matrix;
    Eigen::VectorXd response;
    Eigen::VectorXd previous_parameter;
    Eigen::VectorXd ridge_diagonal;
};

} // namespace rhbm_gem::algorithm
