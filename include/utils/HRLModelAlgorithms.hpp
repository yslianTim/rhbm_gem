#pragma once

#include <vector>

#include <Eigen/Dense>

#include "HRLModelTypes.hpp"

namespace HRLModelAlgorithms
{
HRLBetaEstimateResult EstimateBetaMDPDE(
    double alpha_r,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
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

Eigen::VectorXd CalculateBetaByOLS(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y
);

Eigen::VectorXd CalculateBetaByMDPDE(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const HRLDiagonalMatrix & W
);

Eigen::VectorXd CalculateMuByMedian(
    const Eigen::MatrixXd & beta_array
);

Eigen::VectorXd CalculateMuByMDPDE(
    const Eigen::MatrixXd & beta_array,
    const Eigen::ArrayXd & omega_array,
    double omega_sum
);

HRLDiagonalMatrix CalculateDataWeight(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const Eigen::VectorXd & beta,
    double sigma_square,
    double weight_min = HRLExecutionOptions{}.data_weight_min
);

double CalculateDataVarianceSquare(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const HRLDiagonalMatrix & W,
    const Eigen::VectorXd & beta
);

HRLDiagonalMatrix CalculateDataCovariance(
    double sigma_square,
    const HRLDiagonalMatrix & W
);

Eigen::ArrayXd CalculateMemberWeight(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::MatrixXd & capital_lambda,
    double weight_min = HRLExecutionOptions{}.member_weight_min
);

Eigen::MatrixXd CalculateMemberCovariance(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum
);

std::vector<Eigen::MatrixXd> CalculateWeightedMemberCovariance(
    const Eigen::MatrixXd & capital_lambda,
    const Eigen::ArrayXd & omega_array
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
