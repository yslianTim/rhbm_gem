#pragma once

#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::GaussianResponseMath
{

double GetGaussianResponseAtDistance(double distance, double width, int dimension = 3);

double GetGaussianResponseAtPoint(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    double width);

double GetGaussianResponseAtPointWithNeighborhood(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    const std::vector<Eigen::VectorXd> & neighbor_center_list,
    double width);

} // namespace rhbm_gem::GaussianResponseMath
