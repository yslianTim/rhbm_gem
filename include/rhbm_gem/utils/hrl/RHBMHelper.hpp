#pragma once

#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::rhbm_helper
{
HRLMemberDataset BuildMemberDataset(const SeriesPointList & series_point_list);

HRLBetaMatrix BuildBetaMatrix(const std::vector<HRLBetaVector> & beta_list);

HRLGroupEstimationInput BuildGroupInput(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLBetaEstimateResult> & member_fit_results
);

HRLBetaEstimateResult EstimateBetaMDPDE(
    double alpha_r,
    const HRLMemberDataset & dataset,
    const HRLExecutionOptions & options = {}
);

HRLMuEstimateResult EstimateMuMDPDE(
    double alpha_g,
    const HRLBetaMatrix & beta_matrix,
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
    const HRLBetaPosteriorMatrix & beta_posterior_matrix
);

Eigen::Array<bool, Eigen::Dynamic, 1> CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array
);
} // namespace rhbm_gem::rhbm_helper
