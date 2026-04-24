#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <Eigen/Dense>

using RHBMDesignMatrix = Eigen::MatrixXd;
using RHBMResponseVector = Eigen::VectorXd;
using RHBMScoreVector = Eigen::VectorXd;
using RHBMBetaVector = Eigen::VectorXd;
using RHBMMuVector = Eigen::VectorXd;
using RHBMBetaMatrix = Eigen::MatrixXd;
using RHBMGroupCovarianceMatrix = Eigen::MatrixXd;
using RHBMMemberCovarianceMatrix = Eigen::MatrixXd;
using RHBMBetaPosteriorMatrix = Eigen::MatrixXd;
using RHBMPosteriorCovarianceMatrix = Eigen::MatrixXd;
using RHBMDiagonalMatrix = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;

enum class RHBMEstimationStatus
{
    SUCCESS,
    SINGLE_MEMBER,
    INSUFFICIENT_DATA,
    MAX_ITERATIONS_REACHED,
    NUMERICAL_FALLBACK
};

struct RHBMExecutionOptions
{
    bool quiet_mode{ false };
    int thread_size{ 1 };
    int max_iterations{ 100 };
    double tolerance{ 1.0e-5 };
    double data_weight_min{ 1.0e-8 };
    double member_weight_min{ 1.0e-2 };
    std::optional<std::uint32_t> random_seed{};
};

struct RHBMMemberDataset
{
    RHBMDesignMatrix X;     // member basis/design matrix
    RHBMResponseVector y;   // regression response vector
    RHBMScoreVector score;  // member score/weight vector
};

struct RHBMBetaEstimateResult
{
    RHBMEstimationStatus status{ RHBMEstimationStatus::SUCCESS };
    RHBMBetaVector beta_ols;
    RHBMBetaVector beta_mdpde;
    double sigma_square{ 0.0 };
    RHBMDiagonalMatrix data_weight;
    RHBMDiagonalMatrix data_covariance;
};

struct RHBMGroupEstimationInput
{
    int basis_size{ 0 };
    std::vector<RHBMMemberDataset> member_datasets;
    std::vector<RHBMBetaEstimateResult> member_fit_results;
};

struct RHBMMuEstimateResult
{
    RHBMEstimationStatus status{ RHBMEstimationStatus::SUCCESS };
    RHBMMuVector mu_mean;
    RHBMMuVector mu_mdpde;
    Eigen::ArrayXd omega_array;
    double omega_sum{ 0.0 };
    RHBMGroupCovarianceMatrix capital_lambda;
    std::vector<RHBMMemberCovarianceMatrix> member_capital_lambda_list;
};

struct RHBMWebEstimateResult
{
    RHBMEstimationStatus status{ RHBMEstimationStatus::SUCCESS };
    RHBMMuVector mu_prior;
    RHBMBetaPosteriorMatrix beta_posterior_matrix;
    std::vector<RHBMPosteriorCovarianceMatrix> capital_sigma_posterior_list;
};

struct RHBMGroupEstimationResult
{
    RHBMEstimationStatus status{ RHBMEstimationStatus::SUCCESS };
    RHBMMuVector mu_mean;
    RHBMMuVector mu_mdpde;
    RHBMMuVector mu_prior;
    RHBMGroupCovarianceMatrix capital_lambda;
    RHBMBetaPosteriorMatrix beta_posterior_matrix;
    std::vector<RHBMPosteriorCovarianceMatrix> capital_sigma_posterior_list;
    Eigen::ArrayXd omega_array;
    Eigen::ArrayXd statistical_distance_array;
    Eigen::Array<bool, Eigen::Dynamic, 1> outlier_flag_array;
};
