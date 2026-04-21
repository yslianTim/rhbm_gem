#include <rhbm_gem/utils/hrl/HRLGroupEstimator.hpp>

#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <stdexcept>
#include <utility>
#include <vector>

HRLGroupEstimator::HRLGroupEstimator(HRLExecutionOptions options) :
    m_options{ std::move(options) }
{
}

HRLGroupEstimationResult HRLGroupEstimator::Estimate(
    const HRLGroupEstimationInput & input,
    double alpha_g) const
{
    rhbm_gem::NumericValidation::RequirePositive(input.basis_size, "basis_size");
    const auto validated_input{
        HRLDataTransform::BuildGroupInput(input.member_datasets, input.member_fit_results)
    };
    if (validated_input.basis_size != input.basis_size)
    {
        throw std::invalid_argument("basis_size is inconsistent with member datasets.");
    }

    std::vector<HRLBetaVector> beta_list;
    beta_list.reserve(validated_input.member_fit_results.size());
    for (const auto & fit_result : validated_input.member_fit_results)
    {
        beta_list.emplace_back(fit_result.beta_mdpde);
    }

    const auto beta_matrix{ HRLDataTransform::BuildBetaMatrix(beta_list) };
    auto mu_result{ HRLModelAlgorithms::EstimateMuMDPDE(alpha_g, beta_matrix, m_options) };
    if (mu_result.status == HRLEstimationStatus::SINGLE_MEMBER)
    {
        if (!m_options.quiet_mode)
        {
            Logger::Log(LogLevel::Debug,
                "HRLGroupEstimator::Estimate : single-member group, using fallback result.");
        }
        return BuildFallbackResult(validated_input, mu_result);
    }

    std::vector<HRLDiagonalMatrix> capital_sigma_list;
    capital_sigma_list.reserve(validated_input.member_fit_results.size());
    for (const auto & fit_result : validated_input.member_fit_results)
    {
        capital_sigma_list.emplace_back(fit_result.data_covariance);
    }
    auto web_result{
        HRLModelAlgorithms::EstimateWEB(
            validated_input.member_datasets,
            capital_sigma_list,
            mu_result.mu_mdpde,
            mu_result.member_capital_lambda_list,
            m_options
        )
    };
    if (web_result.status != HRLEstimationStatus::SUCCESS)
    {
        if (!m_options.quiet_mode)
        {
            Logger::Log(LogLevel::Debug,
                "HRLGroupEstimator::Estimate : WEB estimation skipped, using fallback result.");
        }
        return BuildFallbackResult(validated_input, mu_result);
    }

    HRLGroupEstimationResult result;
    result.status = mu_result.status;
    result.mu_mean = mu_result.mu_mean;
    result.mu_mdpde = mu_result.mu_mdpde;
    result.mu_prior = web_result.mu_prior;
    result.capital_lambda = mu_result.capital_lambda;
    result.beta_posterior_matrix = web_result.beta_posterior_matrix;
    result.capital_sigma_posterior_list = std::move(web_result.capital_sigma_posterior_list);
    result.omega_array = mu_result.omega_array;
    result.statistical_distance_array = HRLModelAlgorithms::CalculateMemberStatisticalDistance(
        result.mu_prior,
        result.capital_lambda,
        result.beta_posterior_matrix
    );
    result.outlier_flag_array = HRLModelAlgorithms::CalculateOutlierMemberFlag(
        validated_input.basis_size,
        result.statistical_distance_array
    );
    if (mu_result.status == HRLEstimationStatus::MAX_ITERATIONS_REACHED ||
        web_result.status == HRLEstimationStatus::MAX_ITERATIONS_REACHED)
    {
        result.status = HRLEstimationStatus::MAX_ITERATIONS_REACHED;
    }
    return result;
}

HRLGroupEstimationResult HRLGroupEstimator::BuildFallbackResult(
    const HRLGroupEstimationInput & input,
    const HRLMuEstimateResult & mu_result)
{
    std::vector<HRLBetaVector> beta_list;
    beta_list.reserve(input.member_fit_results.size());
    for (const auto & fit_result : input.member_fit_results)
    {
        beta_list.emplace_back(fit_result.beta_mdpde);
    }
    const auto beta_matrix{ HRLDataTransform::BuildBetaMatrix(beta_list) };

    HRLGroupEstimationResult result;
    result.status = (mu_result.status == HRLEstimationStatus::SUCCESS)
        ? HRLEstimationStatus::NUMERICAL_FALLBACK
        : mu_result.status;
    result.mu_mean = mu_result.mu_mean;
    result.mu_mdpde = mu_result.mu_mdpde;
    result.mu_prior = mu_result.mu_mdpde;
    result.capital_lambda = mu_result.capital_lambda;
    result.beta_posterior_matrix = beta_matrix;
    result.omega_array = mu_result.omega_array;
    result.capital_sigma_posterior_list.assign(
        input.member_fit_results.size(),
        HRLPosteriorCovarianceMatrix::Zero(input.basis_size, input.basis_size)
    );
    result.statistical_distance_array = HRLModelAlgorithms::CalculateMemberStatisticalDistance(
        result.mu_prior,
        result.capital_lambda,
        result.beta_posterior_matrix
    );
    result.outlier_flag_array = HRLModelAlgorithms::CalculateOutlierMemberFlag(
        input.basis_size,
        result.statistical_distance_array
    );
    return result;
}
