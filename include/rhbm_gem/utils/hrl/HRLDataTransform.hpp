#pragma once

#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

namespace HRLDataTransform
{
HRLMemberDataset BuildMemberDataset(
    const std::vector<Eigen::VectorXd> & data_vector,
    bool quiet_mode = false
);

Eigen::MatrixXd BuildBetaMatrix(
    const std::vector<Eigen::VectorXd> & beta_list,
    bool quiet_mode = false
);

HRLGroupEstimationInput BuildGroupInput(
    int basis_size,
    const std::vector<std::vector<Eigen::VectorXd>> & data_entries,
    const std::vector<Eigen::VectorXd> & beta_list,
    const std::vector<double> & sigma_square_list,
    const std::vector<HRLDiagonalMatrix> & weight_list,
    const std::vector<HRLDiagonalMatrix> & covariance_list
);
} // namespace HRLDataTransform
