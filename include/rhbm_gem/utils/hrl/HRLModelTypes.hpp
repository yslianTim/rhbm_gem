#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <Eigen/Dense>

using HRLDesignMatrix = Eigen::MatrixXd;
using HRLResponseVector = Eigen::VectorXd;
using HRLScoreVector = Eigen::VectorXd;
using HRLBetaVector = Eigen::VectorXd;
using HRLMuVector = Eigen::VectorXd;
using HRLBetaMatrix = Eigen::MatrixXd;
using HRLGroupCovarianceMatrix = Eigen::MatrixXd;
using HRLMemberCovarianceMatrix = Eigen::MatrixXd;
using HRLBetaPosteriorMatrix = Eigen::MatrixXd;
using HRLPosteriorCovarianceMatrix = Eigen::MatrixXd;
using HRLDiagonalMatrix = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;

enum class HRLEstimationStatus
{
    SUCCESS,
    SINGLE_MEMBER,
    INSUFFICIENT_DATA,
    MAX_ITERATIONS_REACHED,
    NUMERICAL_FALLBACK
};

struct HRLExecutionOptions
{
    bool quiet_mode{ false };
    int thread_size{ 1 };
    int max_iterations{ 100 };
    double tolerance{ 1.0e-5 };
    double data_weight_min{ 1.0e-8 };
    double member_weight_min{ 1.0e-2 };
    std::optional<std::uint32_t> random_seed{};
};

struct HRLMemberDataset
{
    HRLDesignMatrix X;     // member basis/design matrix
    HRLResponseVector y;   // regression response vector
    HRLScoreVector score;  // member score/weight vector
};

struct HRLBetaEstimateResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    HRLBetaVector beta_ols;
    HRLBetaVector beta_mdpde;
    double sigma_square{ 0.0 };
    HRLDiagonalMatrix data_weight;
    HRLDiagonalMatrix data_covariance;
};

struct HRLGroupEstimationInput
{
    int basis_size{ 0 };
    std::vector<HRLMemberDataset> member_datasets;
    std::vector<HRLBetaEstimateResult> member_fit_results;
};

struct HRLMuEstimateResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    HRLMuVector mu_mean;
    HRLMuVector mu_mdpde;
    Eigen::ArrayXd omega_array;
    double omega_sum{ 0.0 };
    HRLGroupCovarianceMatrix capital_lambda;
    std::vector<HRLMemberCovarianceMatrix> member_capital_lambda_list;
};

struct HRLWebEstimateResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    HRLMuVector mu_prior;
    HRLBetaPosteriorMatrix beta_posterior_matrix;
    std::vector<HRLPosteriorCovarianceMatrix> capital_sigma_posterior_list;
};

struct HRLGroupEstimationResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    HRLMuVector mu_mean;
    HRLMuVector mu_mdpde;
    HRLMuVector mu_prior;
    HRLGroupCovarianceMatrix capital_lambda;
    HRLBetaPosteriorMatrix beta_posterior_matrix;
    std::vector<HRLPosteriorCovarianceMatrix> capital_sigma_posterior_list;
    Eigen::ArrayXd omega_array;
    Eigen::ArrayXd statistical_distance_array;
    Eigen::Array<bool, Eigen::Dynamic, 1> outlier_flag_array;
};
