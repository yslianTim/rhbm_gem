#pragma once

#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

namespace HRLModelAlgorithms
{
HRLBetaEstimateResult EstimateBetaMDPDE(
    double alpha_r,
    const HRLMemberDataset & dataset,
    const HRLExecutionOptions & options = {}
);

HRLMuEstimateResult EstimateMuMDPDE(
    double alpha_g,
    const Eigen::MatrixXd & beta_array,
    const HRLExecutionOptions & options = {}
);

HRLWebEstimateResult EstimateWEB(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLDiagonalMatrix> & capital_sigma_list,
    const Eigen::VectorXd & mu_mdpde,
    const std::vector<Eigen::MatrixXd> & member_capital_lambda_list,
    const HRLExecutionOptions & options = {}
);

Eigen::ArrayXd CalculateMemberStatisticalDistance(
    const Eigen::VectorXd & mu_prior,
    const Eigen::MatrixXd & capital_lambda,
    const Eigen::MatrixXd & beta_posterior_array
);

Eigen::Array<bool, Eigen::Dynamic, 1> CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array
);
} // namespace HRLModelAlgorithms
