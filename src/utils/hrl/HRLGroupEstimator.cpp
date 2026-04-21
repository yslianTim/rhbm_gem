#include <rhbm_gem/utils/hrl/HRLGroupEstimator.hpp>

#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
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
    ValidateInput(input);

    std::vector<Eigen::VectorXd> beta_list;
    beta_list.reserve(input.member_estimates.size());
    for (const auto & estimate : input.member_estimates)
    {
        beta_list.emplace_back(estimate.beta_mdpde);
    }
    const auto beta_matrix{ HRLDataTransform::BuildBetaMatrix(beta_list, true) };

    auto mu_result{ HRLModelAlgorithms::EstimateMuMDPDE(alpha_g, beta_matrix, m_options) };
    if (mu_result.status == HRLEstimationStatus::SINGLE_MEMBER)
    {
        if (!m_options.quiet_mode)
        {
            Logger::Log(LogLevel::Debug,
                "HRLGroupEstimator::Estimate : single-member group, using fallback result.");
        }
        return BuildFallbackResult(input, mu_result);
    }

    std::vector<HRLDiagonalMatrix> capital_sigma_list;
    capital_sigma_list.reserve(input.member_estimates.size());
    for (const auto & estimate : input.member_estimates)
    {
        capital_sigma_list.emplace_back(estimate.data_covariance);
    }
    auto web_result{
        HRLModelAlgorithms::EstimateWEB(
            input.member_datasets,
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
        return BuildFallbackResult(input, mu_result);
    }

    HRLGroupEstimationResult result;
    result.status = mu_result.status;
    result.mu_mean = mu_result.mu_mean;
    result.mu_mdpde = mu_result.mu_mdpde;
    result.mu_prior = web_result.mu_prior;
    result.capital_lambda = mu_result.capital_lambda;
    result.beta_posterior_array = web_result.beta_posterior_array;
    result.capital_sigma_posterior_list = std::move(web_result.capital_sigma_posterior_list);
    result.omega_array = mu_result.omega_array;
    result.statistical_distance_array = HRLModelAlgorithms::CalculateMemberStatisticalDistance(
        result.mu_prior,
        result.capital_lambda,
        result.beta_posterior_array
    );
    result.outlier_flag_array = HRLModelAlgorithms::CalculateOutlierMemberFlag(
        input.basis_size,
        result.statistical_distance_array
    );
    if (mu_result.status == HRLEstimationStatus::MAX_ITERATIONS_REACHED ||
        web_result.status == HRLEstimationStatus::MAX_ITERATIONS_REACHED)
    {
        result.status = HRLEstimationStatus::MAX_ITERATIONS_REACHED;
    }
    return result;
}

void HRLGroupEstimator::ValidateInput(const HRLGroupEstimationInput & input)
{
    rhbm_gem::NumericValidation::RequirePositive(input.basis_size, "basis_size");
    if (input.member_datasets.empty())
    {
        throw std::invalid_argument("member_datasets must not be empty.");
    }
    if (input.member_datasets.size() != input.member_estimates.size())
    {
        throw std::invalid_argument("member_datasets and member_estimates sizes are inconsistent.");
    }

    for (std::size_t i = 0; i < input.member_datasets.size(); i++)
    {
        const auto & dataset{ input.member_datasets.at(i) };
        const auto & estimate{ input.member_estimates.at(i) };
        if (dataset.X.cols() != input.basis_size || dataset.X.rows() != dataset.y.rows())
        {
            throw std::invalid_argument("member dataset shape is inconsistent.");
        }
        if (estimate.beta_mdpde.rows() != input.basis_size)
        {
            throw std::invalid_argument("member beta basis size is inconsistent.");
        }
        try
        {
            rhbm_gem::EigenValidation::RequireVectorSize(
                estimate.data_weight.diagonal(),
                dataset.y.size(),
                "member weight");
            rhbm_gem::EigenValidation::RequireVectorSize(
                estimate.data_covariance.diagonal(),
                dataset.y.size(),
                "member covariance");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("member weight or covariance size is inconsistent.");
        }
    }
}

HRLGroupEstimationResult HRLGroupEstimator::BuildFallbackResult(
    const HRLGroupEstimationInput & input,
    const HRLMuEstimateResult & mu_result)
{
    std::vector<Eigen::VectorXd> beta_list;
    beta_list.reserve(input.member_estimates.size());
    for (const auto & estimate : input.member_estimates)
    {
        beta_list.emplace_back(estimate.beta_mdpde);
    }
    const auto beta_matrix{ HRLDataTransform::BuildBetaMatrix(beta_list, true) };

    HRLGroupEstimationResult result;
    result.status = (mu_result.status == HRLEstimationStatus::SUCCESS)
        ? HRLEstimationStatus::NUMERICAL_FALLBACK
        : mu_result.status;
    result.mu_mean = mu_result.mu_mean;
    result.mu_mdpde = mu_result.mu_mdpde;
    result.mu_prior = mu_result.mu_mdpde;
    result.capital_lambda = mu_result.capital_lambda;
    result.beta_posterior_array = beta_matrix;
    result.omega_array = mu_result.omega_array;
    result.capital_sigma_posterior_list.assign(
        input.member_estimates.size(),
        Eigen::MatrixXd::Zero(input.basis_size, input.basis_size)
    );
    result.statistical_distance_array = HRLModelAlgorithms::CalculateMemberStatisticalDistance(
        result.mu_prior,
        result.capital_lambda,
        result.beta_posterior_array
    );
    result.outlier_flag_array = HRLModelAlgorithms::CalculateOutlierMemberFlag(
        input.basis_size,
        result.statistical_distance_array
    );
    return result;
}
