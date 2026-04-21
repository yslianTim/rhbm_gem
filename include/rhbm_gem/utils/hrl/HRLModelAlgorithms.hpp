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
    const HRLBetaMatrix & beta_array,
    const HRLExecutionOptions & options = {}
);

HRLWebEstimateResult EstimateWEB(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLDiagonalMatrix> & capital_sigma_list,
    const HRLMuVector & mu_mdpde,
    const std::vector<HRLMemberCovarianceMatrix> & member_capital_lambda_list,
    const HRLExecutionOptions & options = {}
);

Eigen::ArrayXd CalculateMemberStatisticalDistance(
    const HRLMuVector & mu_prior,
    const HRLGroupCovarianceMatrix & capital_lambda,
    const HRLBetaMatrix & beta_posterior_array
);

Eigen::Array<bool, Eigen::Dynamic, 1> CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array
);
} // namespace HRLModelAlgorithms
