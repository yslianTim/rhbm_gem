#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <Eigen/Dense>

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
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
};

struct HRLMemberLocalEstimate
{
    Eigen::VectorXd beta_mdpde;
    double sigma_square{ 0.0 };
    HRLDiagonalMatrix data_weight;
    HRLDiagonalMatrix data_covariance;
};

struct HRLGroupEstimationInput
{
    int basis_size{ 0 };
    std::vector<HRLMemberDataset> member_datasets;
    std::vector<HRLMemberLocalEstimate> member_estimates;
};

struct HRLBetaEstimateResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    Eigen::VectorXd beta_ols;
    Eigen::VectorXd beta_mdpde;
    double sigma_square{ 0.0 };
    HRLDiagonalMatrix data_weight;
    HRLDiagonalMatrix data_covariance;
};

struct HRLMuEstimateResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    Eigen::VectorXd mu_mean;
    Eigen::VectorXd mu_mdpde;
    Eigen::ArrayXd omega_array;
    double omega_sum{ 0.0 };
    Eigen::MatrixXd capital_lambda;
    std::vector<Eigen::MatrixXd> member_capital_lambda_list;
};

struct HRLWebEstimateResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    Eigen::VectorXd mu_prior;
    Eigen::MatrixXd beta_posterior_array;
    std::vector<Eigen::MatrixXd> capital_sigma_posterior_list;
};

struct HRLGroupEstimationResult
{
    HRLEstimationStatus status{ HRLEstimationStatus::SUCCESS };
    Eigen::VectorXd mu_mean;
    Eigen::VectorXd mu_mdpde;
    Eigen::VectorXd mu_prior;
    Eigen::MatrixXd capital_lambda;
    Eigen::MatrixXd beta_posterior_array;
    std::vector<Eigen::MatrixXd> capital_sigma_posterior_list;
    Eigen::ArrayXd omega_array;
    Eigen::ArrayXd statistical_distance_array;
    Eigen::Array<bool, Eigen::Dynamic, 1> outlier_flag_array;
};
