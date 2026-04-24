#pragma once

#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::rhbm_helper
{
RHBMMemberDataset BuildMemberDataset(const SeriesPointList & series_point_list);

RHBMBetaMatrix BuildBetaMatrix(const std::vector<RHBMBetaVector> & beta_list);

RHBMGroupEstimationInput BuildGroupInput(
    const std::vector<RHBMMemberDataset> & member_datasets,
    const std::vector<RHBMBetaEstimateResult> & member_fit_results
);

RHBMBetaEstimateResult EstimateBetaMDPDE(
    double alpha_r,
    const RHBMMemberDataset & dataset,
    const RHBMExecutionOptions & options = {}
);

RHBMMuEstimateResult EstimateMuMDPDE(
    double alpha_g,
    const RHBMBetaMatrix & beta_matrix,
    const RHBMExecutionOptions & options = {}
);

RHBMWebEstimateResult EstimateWEB(
    const std::vector<RHBMMemberDataset> & member_datasets,
    const std::vector<RHBMDiagonalMatrix> & capital_sigma_list,
    const RHBMMuVector & mu_mdpde,
    const std::vector<RHBMMemberCovarianceMatrix> & member_capital_lambda_list,
    const RHBMExecutionOptions & options = {}
);

RHBMGroupEstimationResult EstimateGroup(
    double alpha_g,
    const RHBMGroupEstimationInput & input,
    const RHBMExecutionOptions & options = {}
);

Eigen::ArrayXd CalculateMemberStatisticalDistance(
    const RHBMMuVector & mu_prior,
    const RHBMGroupCovarianceMatrix & capital_lambda,
    const RHBMBetaPosteriorMatrix & beta_posterior_matrix
);

Eigen::Array<bool, Eigen::Dynamic, 1> CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array
);
} // namespace rhbm_gem::rhbm_helper
