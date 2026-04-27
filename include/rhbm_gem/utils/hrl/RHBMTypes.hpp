#pragma once

#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>

namespace rhbm_gem {

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

// GaussianParameterVector stays in Gaussian-model space; its dimension depends on the model.
// RHBMBetaVector stays in HRL linearized coefficient space.
using GaussianParameterVector = Eigen::VectorXd;

// Shared Gaussian value types used by HRL linearization and downstream consumers.
struct GaussianEstimate
{
    double amplitude{ 0.0 };
    double width{ 0.0 };

    double GetParameter(int par_id) const
    {
        switch (par_id)
        {
        case 0:
            return amplitude;
        case 1:
            return width;
        case 2:
            return Intensity();
        default:
            throw std::out_of_range("GaussianEstimate parameter index is out of range.");
        }
    }

    double Intensity() const
    {
        if (width == 0.0)
        {
            return 0.0;
        }
        return amplitude * std::pow(Constants::two_pi * width * width, -1.5);
    }

};

struct GaussianParameterUncertainty
{
    double amplitude{ 0.0 };
    double width{ 0.0 };
};

struct GaussianEstimateWithUncertainty
{
    GaussianEstimate estimate{};
    GaussianParameterUncertainty standard_deviation{};

    double GetEstimate(int par_id) const
    {
        return estimate.GetParameter(par_id);
    }

    double GetStandardDeviation(int par_id) const
    {
        switch (par_id)
        {
        case 0:
            return standard_deviation.amplitude;
        case 1:
            return standard_deviation.width;
        case 2:
            return IntensityStandardDeviation();
        default:
            throw std::out_of_range(
                "GaussianEstimateWithUncertainty standard deviation index is out of range.");
        }
    }

    double IntensityStandardDeviation() const
    {
        const auto sigma_amplitude{ standard_deviation.amplitude };
        const auto sigma_width{ standard_deviation.width };
        const auto amplitude{ estimate.amplitude };
        const auto width{ estimate.width };
        return std::sqrt(
            std::pow(std::pow(Constants::two_pi * width * width, -1.5) * sigma_amplitude, 2) +
            std::pow(
                -3.0 * amplitude * std::pow(Constants::two_pi, -1.5)
                    * std::pow(width, -4) * sigma_width,
                2));
    }
};

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

} // namespace rhbm_gem
