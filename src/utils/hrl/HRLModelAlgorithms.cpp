#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>

#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/EigenMatrixUtility.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "utils/hrl/detail/ScopedEigenThreadCount.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <boost/math/distributions/chi_squared.hpp>

namespace
{
void ValidateMatrixVectorShape(
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector)
{
    if (design_matrix.rows() != response_vector.rows())
    {
        throw std::invalid_argument("X and y sizes are inconsistent.");
    }
    if (design_matrix.cols() <= 0)
    {
        throw std::invalid_argument("X must contain at least one basis column.");
    }
    try
    {
        rhbm_gem::EigenValidation::RequireFinite(design_matrix, "X");
        rhbm_gem::EigenValidation::RequireFinite(response_vector, "y");
    }
    catch (const std::invalid_argument &)
    {
        throw std::invalid_argument("X and y must contain only finite values.");
    }
}

void ValidateDiagonalSize(
    const HRLDiagonalMatrix & matrix,
    Eigen::Index expected_size,
    const char * name)
{
    try
    {
        rhbm_gem::EigenValidation::RequireVectorSize(
            matrix.diagonal(),
            expected_size,
            name);
    }
    catch (const std::invalid_argument &)
    {
        throw std::invalid_argument(std::string(name) + " size is inconsistent.");
    }
}

double CalculateChiSquareQuantile(int df)
{
    const boost::math::chi_squared_distribution<double> distribution(static_cast<double>(df));
    return boost::math::quantile(distribution, 0.99);
}

HRLBetaVector CalculateBetaByOLS(
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector);

HRLBetaVector CalculateBetaWithSelectedDataByOLS(
    const HRLMemberDataset & dataset);

HRLBetaVector CalculateBetaByMDPDE(
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector,
    const HRLDiagonalMatrix & W);

HRLMuVector CalculateMuByMedian(
    const HRLBetaMatrix & beta_matrix);

HRLMuVector CalculateMuByMDPDE(
    const HRLBetaMatrix & beta_matrix,
    const Eigen::ArrayXd & omega_array,
    double omega_sum);

HRLDiagonalMatrix CalculateDataWeight(
    double alpha,
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector,
    const HRLBetaVector & beta,
    double sigma_square,
    double weight_min);

double CalculateDataVarianceSquare(
    double alpha,
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector,
    const HRLDiagonalMatrix & W,
    const HRLBetaVector & beta);

HRLDiagonalMatrix CalculateDataCovariance(
    double sigma_square,
    const HRLDiagonalMatrix & W);

Eigen::ArrayXd CalculateMemberWeight(
    double alpha,
    const HRLBetaMatrix & beta_matrix,
    const HRLMuVector & mu,
    const HRLGroupCovarianceMatrix & capital_lambda,
    double weight_min);

HRLGroupCovarianceMatrix CalculateMemberCovariance(
    double alpha,
    const HRLBetaMatrix & beta_matrix,
    const HRLMuVector & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum);

std::vector<HRLMemberCovarianceMatrix> CalculateWeightedMemberCovariance(
    const HRLGroupCovarianceMatrix & capital_lambda,
    const Eigen::ArrayXd & omega_array);
} // namespace

HRLBetaEstimateResult HRLModelAlgorithms::EstimateBetaMDPDE(
    double alpha_r,
    const HRLMemberDataset & dataset,
    const HRLExecutionOptions & options)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(alpha_r, "alpha");
    rhbm_gem::NumericValidation::RequirePositive(options.max_iterations, "max_iterations");
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(options.tolerance, "tolerance");
    rhbm_gem::NumericValidation::RequireFinitePositive(options.data_weight_min, "data_weight_min");
    const auto & design_matrix{ dataset.X };
    const auto & response_vector{ dataset.y };
    ValidateMatrixVectorShape(design_matrix, response_vector);

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    HRLBetaEstimateResult result;
    const auto data_size{ static_cast<int>(response_vector.size()) };
    if (data_size <= 1)
    {
        result.status = HRLEstimationStatus::INSUFFICIENT_DATA;
        result.beta_ols = HRLBetaVector::Zero(design_matrix.cols());
        result.beta_mdpde = HRLBetaVector::Zero(design_matrix.cols());
        result.sigma_square = std::numeric_limits<double>::max();
        result.data_weight.resize(1);
        result.data_weight.diagonal().setOnes();
        result.data_covariance.resize(1);
        result.data_covariance.diagonal().setOnes();
        return result;
    }

    result.beta_ols = CalculateBetaByOLS(design_matrix, response_vector);
    result.beta_mdpde = result.beta_ols;
    //result.beta_mdpde = CalculateBetaWithSelectedDataByOLS(dataset); // TEST
    result.sigma_square =
        (response_vector - (design_matrix * result.beta_mdpde)).squaredNorm() /
        static_cast<double>(data_size - 1);
    result.data_weight = HRLResponseVector::Ones(data_size).asDiagonal();

    auto beta_in_previous_iter{ result.beta_mdpde };
    bool converged{ false };
    for (int t = 0; t < options.max_iterations; t++)
    {
        result.data_weight = CalculateDataWeight(
            alpha_r,
            design_matrix,
            response_vector,
            result.beta_mdpde,
            result.sigma_square,
            options.data_weight_min
        );
        result.beta_mdpde = CalculateBetaByMDPDE(
            design_matrix,
            response_vector,
            result.data_weight);
        result.sigma_square = CalculateDataVarianceSquare(
            alpha_r,
            design_matrix,
            response_vector,
            result.data_weight,
            result.beta_mdpde
        );
        if ((result.beta_mdpde - beta_in_previous_iter).squaredNorm() < options.tolerance)
        {
            converged = true;
            break;
        }
        beta_in_previous_iter = result.beta_mdpde;
    }

    result.data_covariance = CalculateDataCovariance(result.sigma_square, result.data_weight);
    if (!converged)
    {
        result.status = HRLEstimationStatus::MAX_ITERATIONS_REACHED;
    }
    if (result.sigma_square == std::numeric_limits<double>::max())
    {
        result.status = HRLEstimationStatus::NUMERICAL_FALLBACK;
    }
    return result;
}

HRLMuEstimateResult HRLModelAlgorithms::EstimateMuMDPDE(
    double alpha_g,
    const HRLBetaMatrix & beta_matrix,
    const HRLExecutionOptions & options)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(alpha_g, "alpha");
    rhbm_gem::NumericValidation::RequirePositive(options.max_iterations, "max_iterations");
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(options.tolerance, "tolerance");
    rhbm_gem::NumericValidation::RequireFinitePositive(options.member_weight_min, "member_weight_min");
    rhbm_gem::EigenValidation::RequireNonEmpty(beta_matrix, "beta_matrix");

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    HRLMuEstimateResult result;
    const auto basis_size{ static_cast<int>(beta_matrix.rows()) };
    const auto member_size{ static_cast<int>(beta_matrix.cols()) };
    if (member_size == 1)
    {
        result.status = HRLEstimationStatus::SINGLE_MEMBER;
        result.mu_mean = beta_matrix.rowwise().mean();
        result.mu_mdpde = result.mu_mean;
        result.omega_array = Eigen::ArrayXd::Ones(member_size);
        result.omega_sum = result.omega_array.sum();
        result.capital_lambda = HRLGroupCovarianceMatrix::Identity(basis_size, basis_size);
        result.member_capital_lambda_list.assign(
            static_cast<std::size_t>(member_size),
            HRLMemberCovarianceMatrix::Identity(basis_size, basis_size)
        );
        return result;
    }

    result.mu_mean = beta_matrix.rowwise().mean();
    result.mu_mdpde = CalculateMuByMedian(beta_matrix);
    result.capital_lambda = HRLGroupCovarianceMatrix::Identity(basis_size, basis_size);
    auto mu_in_previous_iter{ result.mu_mdpde };
    bool converged{ false };
    for (int t = 0; t < options.max_iterations; t++)
    {
        result.omega_array = CalculateMemberWeight(
            alpha_g,
            beta_matrix,
            result.mu_mdpde,
            result.capital_lambda,
            options.member_weight_min
        );
        result.omega_sum = result.omega_array.sum();
        result.mu_mdpde = CalculateMuByMDPDE(beta_matrix, result.omega_array, result.omega_sum);
        result.capital_lambda = CalculateMemberCovariance(
            alpha_g,
            beta_matrix,
            result.mu_mdpde,
            result.omega_array,
            result.omega_sum
        );
        if ((result.mu_mdpde - mu_in_previous_iter).squaredNorm() < options.tolerance)
        {
            converged = true;
            break;
        }
        mu_in_previous_iter = result.mu_mdpde;
    }

    result.member_capital_lambda_list = CalculateWeightedMemberCovariance(
        result.capital_lambda,
        result.omega_array
    );
    if (!converged)
    {
        result.status = HRLEstimationStatus::MAX_ITERATIONS_REACHED;
    }
    try
    {
        rhbm_gem::EigenValidation::RequireFinite(result.capital_lambda, "capital_lambda");
    }
    catch (const std::invalid_argument &)
    {
        result.status = HRLEstimationStatus::NUMERICAL_FALLBACK;
    }
    return result;
}

HRLWebEstimateResult HRLModelAlgorithms::EstimateWEB(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLDiagonalMatrix> & capital_sigma_list,
    const HRLMuVector & mu_mdpde,
    const std::vector<HRLMemberCovarianceMatrix> & member_capital_lambda_list,
    const HRLExecutionOptions & options)
{
    if (member_datasets.empty())
    {
        throw std::invalid_argument("member_datasets must not be empty.");
    }
    if (capital_sigma_list.size() != member_datasets.size() ||
        member_capital_lambda_list.size() != member_datasets.size())
    {
        throw std::invalid_argument("WEB inputs must have consistent member counts.");
    }
    rhbm_gem::EigenValidation::RequireNonEmpty(mu_mdpde, "mu_mdpde");

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    HRLWebEstimateResult result;
    const auto basis_size{ static_cast<int>(mu_mdpde.rows()) };
    const auto member_size{ static_cast<int>(member_datasets.size()) };
    result.mu_prior = HRLMuVector::Zero(basis_size);
    result.beta_posterior_matrix = HRLBetaPosteriorMatrix::Zero(basis_size, member_size);
    result.capital_sigma_posterior_list.clear();
    result.capital_sigma_posterior_list.reserve(static_cast<std::size_t>(member_size));
    if (member_size == 1)
    {
        result.status = HRLEstimationStatus::SINGLE_MEMBER;
        return result;
    }

    HRLMuVector numerator{ HRLMuVector::Zero(basis_size) };
    HRLGroupCovarianceMatrix denominator{
        HRLGroupCovarianceMatrix::Zero(basis_size, basis_size)
    };
    for (std::size_t i = 0; i < static_cast<std::size_t>(member_size); i++)
    {
        const auto & dataset{ member_datasets.at(i) };
        const auto & capital_sigma{ capital_sigma_list.at(i) };
        const auto & member_capital_lambda{ member_capital_lambda_list.at(i) };
        ValidateMatrixVectorShape(dataset.X, dataset.y);
        ValidateDiagonalSize(capital_sigma, dataset.y.size(), "capital_sigma");
        if (member_capital_lambda.rows() != basis_size || member_capital_lambda.cols() != basis_size)
        {
            throw std::invalid_argument("member_capital_lambda size is inconsistent.");
        }

        const auto inv_capital_sigma{
            EigenMatrixUtility::GetInverseDiagonalMatrix(capital_sigma)
        };
        const auto inv_member_capital_lambda{
            EigenMatrixUtility::GetInverseMatrix(member_capital_lambda)
        };
        const HRLGroupCovarianceMatrix gram_matrix{
            dataset.X.transpose() * inv_capital_sigma * dataset.X
        };
        const HRLMuVector moment_matrix{
            dataset.X.transpose() * inv_capital_sigma * dataset.y
        };
        const HRLGroupCovarianceMatrix inv_capital_sigma_posterior{
            gram_matrix + inv_member_capital_lambda
        };
        const HRLPosteriorCovarianceMatrix capital_sigma_posterior{
            EigenMatrixUtility::GetInverseMatrix(inv_capital_sigma_posterior)
        };

        result.capital_sigma_posterior_list.emplace_back(capital_sigma_posterior);
        result.beta_posterior_matrix.col(static_cast<Eigen::Index>(i)) =
            capital_sigma_posterior *
            (moment_matrix + inv_member_capital_lambda * mu_mdpde);
        numerator += inv_member_capital_lambda * capital_sigma_posterior * moment_matrix;
        denominator += inv_member_capital_lambda * capital_sigma_posterior * gram_matrix;
    }

    result.mu_prior = EigenMatrixUtility::GetInverseMatrix(denominator) * numerator;
    if (member_size == 2)
    {
        result.mu_prior = mu_mdpde;
    }
    return result;
}

namespace
{
HRLBetaVector CalculateBetaByOLS(
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector)
{
    const HRLGroupCovarianceMatrix gram_matrix{ design_matrix.transpose() * design_matrix };
    const HRLGroupCovarianceMatrix inverse_gram_matrix{
        EigenMatrixUtility::GetInverseMatrix(gram_matrix)
    };
    return inverse_gram_matrix * (design_matrix.transpose() * response_vector);
}

HRLBetaVector CalculateBetaWithSelectedDataByOLS(
    const HRLMemberDataset & dataset)
{
    ValidateMatrixVectorShape(dataset.X, dataset.y);
    if (dataset.score.size() != dataset.y.size())
    {
        throw std::invalid_argument("dataset score size is inconsistent.");
    }

    std::vector<Eigen::Index> selected_row_indices;
    selected_row_indices.reserve(static_cast<std::size_t>(dataset.y.size()));
    for (Eigen::Index row = 0; row < dataset.score.size(); row++)
    {
        if (dataset.score(row) != 0.0)
        {
            selected_row_indices.emplace_back(row);
        }
    }

    if (selected_row_indices.size() <= 1U)
    {
        return CalculateBetaByOLS(dataset.X, dataset.y);
    }

    const auto selected_row_count{ static_cast<Eigen::Index>(selected_row_indices.size()) };
    HRLDesignMatrix selected_design_matrix{
        HRLDesignMatrix::Zero(selected_row_count, dataset.X.cols())
    };
    HRLResponseVector selected_response_vector{ HRLResponseVector::Zero(selected_row_count) };
    for (Eigen::Index i = 0; i < selected_row_count; i++)
    {
        const auto row{ selected_row_indices.at(static_cast<std::size_t>(i)) };
        selected_design_matrix.row(i) = dataset.X.row(row);
        selected_response_vector(i) = dataset.y(row);
    }

    return CalculateBetaByOLS(selected_design_matrix, selected_response_vector);
}

HRLBetaVector CalculateBetaByMDPDE(
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector,
    const HRLDiagonalMatrix & W)
{
    ValidateMatrixVectorShape(design_matrix, response_vector);
    ValidateDiagonalSize(W, response_vector.size(), "W");
    const HRLGroupCovarianceMatrix gram_matrix{
        design_matrix.transpose() * W * design_matrix
    };
    const auto inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
    return inverse_gram_matrix * (design_matrix.transpose() * W * response_vector);
}

HRLMuVector CalculateMuByMedian(const HRLBetaMatrix & beta_matrix)
{
    rhbm_gem::EigenValidation::RequireNonEmpty(beta_matrix, "beta_matrix");
    const auto basis_size{ static_cast<int>(beta_matrix.rows()) };
    HRLMuVector mu{ HRLMuVector::Zero(basis_size) };
    for (int b = 0; b < basis_size; b++)
    {
        mu(b) = EigenMatrixUtility::GetMedian(beta_matrix.row(b));
    }
    return mu;
}

HRLMuVector CalculateMuByMDPDE(
    const HRLBetaMatrix & beta_matrix,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    if (beta_matrix.cols() != omega_array.rows())
    {
        throw std::invalid_argument("beta_matrix and omega_array sizes are inconsistent.");
    }
    rhbm_gem::NumericValidation::RequireFinitePositive(omega_sum, "omega_sum");
    HRLBetaMatrix numerator{ beta_matrix.array() / omega_sum };
    for (int i = 0; i < numerator.cols(); i++)
    {
        numerator.col(i) *= omega_array(i);
    }
    return numerator.rowwise().sum();
}

HRLDiagonalMatrix CalculateDataWeight(
    double alpha,
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector,
    const HRLBetaVector & beta,
    double sigma_square,
    double weight_min)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(alpha, "alpha");
    rhbm_gem::NumericValidation::RequireFinitePositive(weight_min, "weight_min");
    ValidateMatrixVectorShape(design_matrix, response_vector);
    if (beta.rows() != design_matrix.cols())
    {
        throw std::invalid_argument("beta size is inconsistent with X.");
    }

    if (alpha == 0.0)
    {
        return HRLResponseVector::Ones(response_vector.size()).asDiagonal();
    }
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0)
    {
        sigma_square = weight_min;
    }

    const auto max_log{ std::log(std::numeric_limits<double>::max()) };
    Eigen::ArrayXd exponent{
        -0.5 * alpha * (response_vector - (design_matrix * beta)).array().square() / sigma_square
    };
    exponent = exponent.cwiseMin(max_log);
    Eigen::ArrayXd W{ exponent.exp() };
    W = W.unaryExpr([weight_min](double w) {
        return std::isfinite(w) ? w : weight_min;
    });
    W = W.cwiseMax(weight_min);
    return W.matrix().asDiagonal();
}

double CalculateDataVarianceSquare(
    double alpha,
    const HRLDesignMatrix & design_matrix,
    const HRLResponseVector & response_vector,
    const HRLDiagonalMatrix & W,
    const HRLBetaVector & beta)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(alpha, "alpha");
    ValidateMatrixVectorShape(design_matrix, response_vector);
    ValidateDiagonalSize(W, response_vector.size(), "W");
    if (beta.rows() != design_matrix.cols())
    {
        throw std::invalid_argument("beta size is inconsistent with X.");
    }

    const auto n{ static_cast<double>(response_vector.size()) };
    const HRLResponseVector residual{ response_vector - (design_matrix * beta) };
    const auto numerator{ static_cast<double>(residual.transpose() * W * residual) };
    auto denominator{ W.diagonal().sum() - n * alpha * std::pow(1.0 + alpha, -1.5) };
    if (denominator <= 0.0)
    {
        denominator = std::numeric_limits<double>::min();
    }
    auto sigma_square{ numerator / denominator };
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0)
    {
        sigma_square = std::numeric_limits<double>::max();
    }
    return sigma_square;
}

HRLDiagonalMatrix CalculateDataCovariance(
    double sigma_square,
    const HRLDiagonalMatrix & W)
{
    const HRLResponseVector data_weight_array{ W.diagonal() };
    const auto n{ static_cast<int>(data_weight_array.size()) };
    const auto W_inverse_trace{
        EigenMatrixUtility::GetInverseDiagonalMatrix(W).diagonal().sum()
    };
    if (!std::isfinite(W_inverse_trace) || W_inverse_trace <= 0.0)
    {
        return HRLResponseVector::Ones(n).asDiagonal();
    }

    HRLResponseVector capital_sigma{ HRLResponseVector::Zero(n) };
    for (int j = 0; j < n; j++)
    {
        if (data_weight_array(j) == 0.0 || W_inverse_trace == 0.0)
        {
            continue;
        }
        capital_sigma(j) = n * sigma_square / data_weight_array(j) / W_inverse_trace;
        if (!std::isfinite(capital_sigma(j)) || capital_sigma(j) <= 0.0)
        {
            capital_sigma(j) = 1.0;
        }
    }
    return capital_sigma.asDiagonal();
}

Eigen::ArrayXd CalculateMemberWeight(
    double alpha,
    const HRLBetaMatrix & beta_matrix,
    const HRLMuVector & mu,
    const HRLGroupCovarianceMatrix & capital_lambda,
    double weight_min)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(alpha, "alpha");
    rhbm_gem::NumericValidation::RequireFinitePositive(weight_min, "weight_min");
    if (beta_matrix.rows() != mu.rows())
    {
        throw std::invalid_argument("beta_matrix and mu sizes are inconsistent.");
    }
    if (capital_lambda.rows() != capital_lambda.cols() ||
        capital_lambda.rows() != beta_matrix.rows())
    {
        throw std::invalid_argument("capital_lambda size is inconsistent.");
    }

    const auto member_size{ static_cast<int>(beta_matrix.cols()) };
    const auto weight_member_min{ weight_min / static_cast<double>(member_size) };
    const auto inverse_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(capital_lambda) };
    const HRLBetaMatrix residual_matrix{ beta_matrix.colwise() - mu };
    Eigen::ArrayXd omega_array{ Eigen::ArrayXd::Zero(member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const HRLBetaVector residual{ residual_matrix.col(i) };
        const auto exp_index{
            static_cast<double>(residual.transpose() * inverse_capital_lambda * residual)
        };
        if (!std::isfinite(exp_index))
        {
            omega_array(i) = weight_member_min;
        }
        else
        {
            omega_array(i) = std::exp(-0.5 * alpha * exp_index);
            if (omega_array(i) < weight_member_min)
            {
                omega_array(i) = weight_member_min;
            }
        }
    }
    return omega_array;
}

HRLGroupCovarianceMatrix CalculateMemberCovariance(
    double alpha,
    const HRLBetaMatrix & beta_matrix,
    const HRLMuVector & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    rhbm_gem::NumericValidation::RequireFiniteNonNegative(alpha, "alpha");
    if (beta_matrix.rows() != mu.rows() || beta_matrix.cols() != omega_array.rows())
    {
        throw std::invalid_argument("Member covariance inputs are inconsistent.");
    }
    const auto basis_size{ static_cast<int>(beta_matrix.rows()) };
    const auto member_size{ static_cast<int>(beta_matrix.cols()) };
    HRLGroupCovarianceMatrix numerator{
        HRLGroupCovarianceMatrix::Zero(basis_size, basis_size)
    };
    const double denominator{
        omega_sum - member_size * alpha * std::pow(1.0 + alpha, -1.0 - 0.5 * basis_size)
    };
    if (denominator <= 0.0)
    {
        return HRLGroupCovarianceMatrix::Identity(basis_size, basis_size);
    }

    const HRLBetaMatrix residual_matrix{ beta_matrix.colwise() - mu };
    for (int i = 0; i < member_size; i++)
    {
        const HRLBetaVector residual{ residual_matrix.col(i) };
        numerator += omega_array(i) * (residual * residual.transpose());
    }

    HRLGroupCovarianceMatrix capital_lambda{ numerator / denominator };
    try
    {
        rhbm_gem::EigenValidation::RequireFinite(capital_lambda, "capital_lambda");
    }
    catch (const std::invalid_argument &)
    {
        capital_lambda = HRLGroupCovarianceMatrix::Identity(basis_size, basis_size);
    }
    return capital_lambda;
}

std::vector<HRLMemberCovarianceMatrix> CalculateWeightedMemberCovariance(
    const HRLGroupCovarianceMatrix & capital_lambda,
    const Eigen::ArrayXd & omega_array)
{
    try
    {
        rhbm_gem::EigenValidation::RequireSquare(capital_lambda, "capital_lambda");
    }
    catch (const std::invalid_argument &)
    {
        throw std::invalid_argument("capital_lambda must be square.");
    }
    const auto member_size{ static_cast<int>(omega_array.size()) };
    if (member_size <= 0)
    {
        throw std::invalid_argument("omega_array must not be empty.");
    }

    auto omega_inverse_sum{ 0.0 };
    for (int i = 0; i < member_size; i++)
    {
        omega_inverse_sum += (omega_array(i) == 0.0) ? 0.0 : 1.0 / omega_array(i);
    }
    if (omega_inverse_sum <= 0.0)
    {
        return std::vector<HRLMemberCovarianceMatrix>(
            static_cast<std::size_t>(member_size),
            capital_lambda
        );
    }

    const auto omega_h{ 1.0 / omega_inverse_sum };
    std::vector<HRLMemberCovarianceMatrix> member_capital_lambda_list(
        static_cast<std::size_t>(member_size));
    for (int i = 0; i < member_size; i++)
    {
        member_capital_lambda_list.at(static_cast<std::size_t>(i)) =
            (member_size * omega_h / omega_array(i)) * capital_lambda;
    }
    return member_capital_lambda_list;
}
} // namespace

Eigen::ArrayXd HRLModelAlgorithms::CalculateMemberStatisticalDistance(
    const HRLMuVector & mu_prior,
    const HRLGroupCovarianceMatrix & capital_lambda,
    const HRLBetaPosteriorMatrix & beta_posterior_matrix)
{
    if (beta_posterior_matrix.rows() != mu_prior.rows())
    {
        throw std::invalid_argument("beta_posterior_matrix and mu_prior sizes are inconsistent.");
    }
    if (capital_lambda.rows() != capital_lambda.cols() ||
        capital_lambda.rows() != mu_prior.rows())
    {
        throw std::invalid_argument("capital_lambda size is inconsistent.");
    }

    const auto member_size{ static_cast<int>(beta_posterior_matrix.cols()) };
    Eigen::ArrayXd statistical_distance_array{ Eigen::ArrayXd::Zero(member_size) };
    const HRLBetaPosteriorMatrix error_matrix{ beta_posterior_matrix.colwise() - mu_prior };
    const auto inv_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(capital_lambda) };
    for (int i = 0; i < member_size; i++)
    {
        statistical_distance_array(i) =
            error_matrix.col(i).transpose() * inv_capital_lambda * error_matrix.col(i);
    }
    return statistical_distance_array;
}

Eigen::Array<bool, Eigen::Dynamic, 1> HRLModelAlgorithms::CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array)
{
    rhbm_gem::NumericValidation::RequirePositive(basis_size, "basis_size");
    return statistical_distance_array > CalculateChiSquareQuantile(basis_size);
}
