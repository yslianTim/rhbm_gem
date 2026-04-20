#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>

#include <rhbm_gem/utils/math/EigenMatrixUtility.hpp>

#include "utils/hrl/detail/ScopedEigenThreadCount.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <boost/math/distributions/chi_squared.hpp>

namespace
{
void ValidateAlpha(double alpha)
{
    if (!std::isfinite(alpha) || alpha < 0.0)
    {
        throw std::invalid_argument("Alpha parameters must be finite and non-negative.");
    }
}

void ValidateMaximumIteration(int size)
{
    if (size <= 0)
    {
        throw std::invalid_argument("size must be greater than 0");
    }
}

void ValidateTolerance(double value)
{
    if (!std::isfinite(value))
    {
        throw std::invalid_argument("tolerance must be finite");
    }
    if (value < 0.0)
    {
        throw std::invalid_argument("tolerance must be non-negative");
    }
}

void ValidateWeightMinimum(double value, const char * name)
{
    if (!std::isfinite(value) || value <= 0.0)
    {
        throw std::invalid_argument(std::string(name) + " must be finite and positive.");
    }
}

void ValidateMatrixVectorShape(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y)
{
    if (X.rows() != y.rows())
    {
        throw std::invalid_argument("X and y sizes are inconsistent.");
    }
    if (X.cols() <= 0)
    {
        throw std::invalid_argument("X must contain at least one basis column.");
    }
    if (!X.array().allFinite() || !y.array().allFinite())
    {
        throw std::invalid_argument("X and y must contain only finite values.");
    }
}

void ValidateDiagonalSize(
    const HRLDiagonalMatrix & matrix,
    Eigen::Index expected_size,
    const char * name)
{
    if (matrix.diagonal().size() != expected_size)
    {
        throw std::invalid_argument(std::string(name) + " size is inconsistent.");
    }
}

double CalculateChiSquareQuantile(int df)
{
    const boost::math::chi_squared_distribution<double> distribution(static_cast<double>(df));
    return boost::math::quantile(distribution, 0.99);
}

Eigen::VectorXd CalculateBetaByOLS(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y);

Eigen::VectorXd CalculateBetaWithSelectedDataByOLS(
    const HRLMemberDataset & dataset);

Eigen::VectorXd CalculateBetaByMDPDE(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const HRLDiagonalMatrix & W);

Eigen::VectorXd CalculateMuByMedian(
    const Eigen::MatrixXd & beta_array);

Eigen::VectorXd CalculateMuByMDPDE(
    const Eigen::MatrixXd & beta_array,
    const Eigen::ArrayXd & omega_array,
    double omega_sum);

HRLDiagonalMatrix CalculateDataWeight(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const Eigen::VectorXd & beta,
    double sigma_square,
    double weight_min);

double CalculateDataVarianceSquare(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const HRLDiagonalMatrix & W,
    const Eigen::VectorXd & beta);

HRLDiagonalMatrix CalculateDataCovariance(
    double sigma_square,
    const HRLDiagonalMatrix & W);

Eigen::ArrayXd CalculateMemberWeight(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::MatrixXd & capital_lambda,
    double weight_min);

Eigen::MatrixXd CalculateMemberCovariance(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum);

std::vector<Eigen::MatrixXd> CalculateWeightedMemberCovariance(
    const Eigen::MatrixXd & capital_lambda,
    const Eigen::ArrayXd & omega_array);
} // namespace

HRLBetaEstimateResult HRLModelAlgorithms::EstimateBetaMDPDE(
    double alpha_r,
    const HRLMemberDataset & dataset,
    const HRLExecutionOptions & options)
{
    ValidateAlpha(alpha_r);
    ValidateMaximumIteration(options.max_iterations);
    ValidateTolerance(options.tolerance);
    ValidateWeightMinimum(options.data_weight_min, "data_weight_min");
    const auto & X{ dataset.X };
    const auto & y{ dataset.y };
    ValidateMatrixVectorShape(X, y);

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    HRLBetaEstimateResult result;
    const auto data_size{ static_cast<int>(y.size()) };
    if (data_size <= 1)
    {
        result.status = HRLEstimationStatus::INSUFFICIENT_DATA;
        result.beta_ols = Eigen::VectorXd::Zero(X.cols());
        result.beta_mdpde = Eigen::VectorXd::Zero(X.cols());
        result.sigma_square = std::numeric_limits<double>::max();
        result.data_weight.resize(1);
        result.data_weight.diagonal().setOnes();
        result.data_covariance.resize(1);
        result.data_covariance.diagonal().setOnes();
        return result;
    }

    result.beta_ols = CalculateBetaByOLS(X, y);
    result.beta_mdpde = result.beta_ols;
    //result.beta_mdpde = CalculateBetaWithSelectedDataByOLS(dataset); // TEST
    result.sigma_square =
        (y - (X * result.beta_mdpde)).squaredNorm() / static_cast<double>(data_size - 1);
    result.data_weight = Eigen::VectorXd::Ones(data_size).asDiagonal();

    auto beta_in_previous_iter{ result.beta_mdpde };
    bool converged{ false };
    for (int t = 0; t < options.max_iterations; t++)
    {
        result.data_weight = CalculateDataWeight(
            alpha_r,
            X,
            y,
            result.beta_mdpde,
            result.sigma_square,
            options.data_weight_min
        );
        result.beta_mdpde = CalculateBetaByMDPDE(X, y, result.data_weight);
        result.sigma_square = CalculateDataVarianceSquare(
            alpha_r,
            X,
            y,
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
    const Eigen::MatrixXd & beta_array,
    const HRLExecutionOptions & options)
{
    ValidateAlpha(alpha_g);
    ValidateMaximumIteration(options.max_iterations);
    ValidateTolerance(options.tolerance);
    ValidateWeightMinimum(options.member_weight_min, "member_weight_min");
    if (beta_array.rows() <= 0 || beta_array.cols() <= 0)
    {
        throw std::invalid_argument("beta_array must not be empty.");
    }

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    HRLMuEstimateResult result;
    const auto basis_size{ static_cast<int>(beta_array.rows()) };
    const auto member_size{ static_cast<int>(beta_array.cols()) };
    if (member_size == 1)
    {
        result.status = HRLEstimationStatus::SINGLE_MEMBER;
        result.mu_mean = beta_array.rowwise().mean();
        result.mu_mdpde = result.mu_mean;
        result.omega_array = Eigen::ArrayXd::Ones(member_size);
        result.omega_sum = result.omega_array.sum();
        result.capital_lambda = Eigen::MatrixXd::Identity(basis_size, basis_size);
        result.member_capital_lambda_list.assign(
            static_cast<std::size_t>(member_size),
            Eigen::MatrixXd::Identity(basis_size, basis_size)
        );
        return result;
    }

    result.mu_mean = beta_array.rowwise().mean();
    result.mu_mdpde = CalculateMuByMedian(beta_array);
    result.capital_lambda = Eigen::MatrixXd::Identity(basis_size, basis_size);
    auto mu_in_previous_iter{ result.mu_mdpde };
    bool converged{ false };
    for (int t = 0; t < options.max_iterations; t++)
    {
        result.omega_array = CalculateMemberWeight(
            alpha_g,
            beta_array,
            result.mu_mdpde,
            result.capital_lambda,
            options.member_weight_min
        );
        result.omega_sum = result.omega_array.sum();
        result.mu_mdpde = CalculateMuByMDPDE(beta_array, result.omega_array, result.omega_sum);
        result.capital_lambda = CalculateMemberCovariance(
            alpha_g,
            beta_array,
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
    if (!result.capital_lambda.array().allFinite())
    {
        result.status = HRLEstimationStatus::NUMERICAL_FALLBACK;
    }
    return result;
}

HRLWebEstimateResult HRLModelAlgorithms::EstimateWEB(
    const std::vector<HRLMemberDataset> & member_datasets,
    const std::vector<HRLDiagonalMatrix> & capital_sigma_list,
    const Eigen::VectorXd & mu_mdpde,
    const std::vector<Eigen::MatrixXd> & member_capital_lambda_list,
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
    if (mu_mdpde.rows() <= 0)
    {
        throw std::invalid_argument("mu_mdpde must not be empty.");
    }

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    HRLWebEstimateResult result;
    const auto basis_size{ static_cast<int>(mu_mdpde.rows()) };
    const auto member_size{ static_cast<int>(member_datasets.size()) };
    result.mu_prior = Eigen::VectorXd::Zero(basis_size);
    result.beta_posterior_array = Eigen::MatrixXd::Zero(basis_size, member_size);
    result.capital_sigma_posterior_list.clear();
    result.capital_sigma_posterior_list.reserve(static_cast<std::size_t>(member_size));
    if (member_size == 1)
    {
        result.status = HRLEstimationStatus::SINGLE_MEMBER;
        return result;
    }

    Eigen::VectorXd numerator{ Eigen::VectorXd::Zero(basis_size) };
    Eigen::MatrixXd denominator{ Eigen::MatrixXd::Zero(basis_size, basis_size) };
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
        const Eigen::MatrixXd gram_matrix{ dataset.X.transpose() * inv_capital_sigma * dataset.X };
        const Eigen::VectorXd moment_matrix{ dataset.X.transpose() * inv_capital_sigma * dataset.y };
        const Eigen::MatrixXd inv_capital_sigma_posterior{
            gram_matrix + inv_member_capital_lambda
        };
        const Eigen::MatrixXd capital_sigma_posterior{
            EigenMatrixUtility::GetInverseMatrix(inv_capital_sigma_posterior)
        };

        result.capital_sigma_posterior_list.emplace_back(capital_sigma_posterior);
        result.beta_posterior_array.col(static_cast<Eigen::Index>(i)) =
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
Eigen::VectorXd CalculateBetaByOLS(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y)
{
    const Eigen::MatrixXd gram_matrix{ X.transpose() * X };
    const Eigen::MatrixXd inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
    return inverse_gram_matrix * (X.transpose() * y);
}

Eigen::VectorXd CalculateBetaWithSelectedDataByOLS(
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
    Eigen::MatrixXd X_selected{ Eigen::MatrixXd::Zero(selected_row_count, dataset.X.cols()) };
    Eigen::VectorXd y_selected{ Eigen::VectorXd::Zero(selected_row_count) };
    for (Eigen::Index i = 0; i < selected_row_count; i++)
    {
        const auto row{ selected_row_indices.at(static_cast<std::size_t>(i)) };
        X_selected.row(i) = dataset.X.row(row);
        y_selected(i) = dataset.y(row);
    }

    return CalculateBetaByOLS(X_selected, y_selected);
}

Eigen::VectorXd CalculateBetaByMDPDE(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const HRLDiagonalMatrix & W)
{
    ValidateMatrixVectorShape(X, y);
    ValidateDiagonalSize(W, y.size(), "W");
    const Eigen::MatrixXd gram_matrix{ X.transpose() * W * X };
    const auto inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
    return inverse_gram_matrix * (X.transpose() * W * y);
}

Eigen::VectorXd CalculateMuByMedian(
    const Eigen::MatrixXd & beta_array)
{
    if (beta_array.rows() <= 0 || beta_array.cols() <= 0)
    {
        throw std::invalid_argument("beta_array must not be empty.");
    }
    const auto basis_size{ static_cast<int>(beta_array.rows()) };
    Eigen::VectorXd mu{ Eigen::VectorXd::Zero(basis_size) };
    for (int b = 0; b < basis_size; b++)
    {
        mu(b) = EigenMatrixUtility::GetMedian(beta_array.row(b));
    }
    return mu;
}

Eigen::VectorXd CalculateMuByMDPDE(
    const Eigen::MatrixXd & beta_array,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    if (beta_array.cols() != omega_array.rows())
    {
        throw std::invalid_argument("beta_array and omega_array sizes are inconsistent.");
    }
    if (!std::isfinite(omega_sum) || omega_sum <= 0.0)
    {
        throw std::invalid_argument("omega_sum must be finite and positive.");
    }
    Eigen::MatrixXd numerator{ beta_array.array() / omega_sum };
    for (int i = 0; i < numerator.cols(); i++)
    {
        numerator.col(i) *= omega_array(i);
    }
    return numerator.rowwise().sum();
}

HRLDiagonalMatrix CalculateDataWeight(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const Eigen::VectorXd & beta,
    double sigma_square,
    double weight_min)
{
    ValidateAlpha(alpha);
    ValidateWeightMinimum(weight_min, "weight_min");
    ValidateMatrixVectorShape(X, y);
    if (beta.rows() != X.cols())
    {
        throw std::invalid_argument("beta size is inconsistent with X.");
    }

    if (alpha == 0.0)
    {
        return Eigen::VectorXd::Ones(y.size()).asDiagonal();
    }
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0)
    {
        sigma_square = weight_min;
    }

    const auto max_log{ std::log(std::numeric_limits<double>::max()) };
    Eigen::ArrayXd exponent{
        -0.5 * alpha * (y - (X * beta)).array().square() / sigma_square
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
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const HRLDiagonalMatrix & W,
    const Eigen::VectorXd & beta)
{
    ValidateAlpha(alpha);
    ValidateMatrixVectorShape(X, y);
    ValidateDiagonalSize(W, y.size(), "W");
    if (beta.rows() != X.cols())
    {
        throw std::invalid_argument("beta size is inconsistent with X.");
    }

    const auto n{ static_cast<double>(y.size()) };
    const Eigen::VectorXd residual{ y - (X * beta) };
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
    const Eigen::VectorXd data_weight_array{ W.diagonal() };
    const auto n{ static_cast<int>(data_weight_array.size()) };
    const auto W_inverse_trace{
        EigenMatrixUtility::GetInverseDiagonalMatrix(W).diagonal().sum()
    };
    if (!std::isfinite(W_inverse_trace) || W_inverse_trace <= 0.0)
    {
        return Eigen::VectorXd::Ones(n).asDiagonal();
    }

    Eigen::VectorXd capital_sigma{ Eigen::VectorXd::Zero(n) };
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
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::MatrixXd & capital_lambda,
    double weight_min)
{
    ValidateAlpha(alpha);
    ValidateWeightMinimum(weight_min, "weight_min");
    if (beta_array.rows() != mu.rows())
    {
        throw std::invalid_argument("beta_array and mu sizes are inconsistent.");
    }
    if (capital_lambda.rows() != capital_lambda.cols() ||
        capital_lambda.rows() != beta_array.rows())
    {
        throw std::invalid_argument("capital_lambda size is inconsistent.");
    }

    const auto member_size{ static_cast<int>(beta_array.cols()) };
    const auto weight_member_min{ weight_min / static_cast<double>(member_size) };
    const auto inverse_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(capital_lambda) };
    const Eigen::MatrixXd residual_array{ beta_array.colwise() - mu };
    Eigen::ArrayXd omega_array{ Eigen::ArrayXd::Zero(member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const Eigen::VectorXd residual{ residual_array.col(i) };
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

Eigen::MatrixXd CalculateMemberCovariance(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    ValidateAlpha(alpha);
    if (beta_array.rows() != mu.rows() || beta_array.cols() != omega_array.rows())
    {
        throw std::invalid_argument("Member covariance inputs are inconsistent.");
    }
    const auto basis_size{ static_cast<int>(beta_array.rows()) };
    const auto member_size{ static_cast<int>(beta_array.cols()) };
    Eigen::MatrixXd numerator{ Eigen::MatrixXd::Zero(basis_size, basis_size) };
    const double denominator{
        omega_sum - member_size * alpha * std::pow(1.0 + alpha, -1.0 - 0.5 * basis_size)
    };
    if (denominator <= 0.0)
    {
        return Eigen::MatrixXd::Identity(basis_size, basis_size);
    }

    const Eigen::MatrixXd residual_array{ beta_array.colwise() - mu };
    for (int i = 0; i < member_size; i++)
    {
        const Eigen::VectorXd residual{ residual_array.col(i) };
        numerator += omega_array(i) * (residual * residual.transpose());
    }

    Eigen::MatrixXd capital_lambda{ numerator / denominator };
    if (!capital_lambda.array().allFinite())
    {
        capital_lambda = Eigen::MatrixXd::Identity(basis_size, basis_size);
    }
    return capital_lambda;
}

std::vector<Eigen::MatrixXd> CalculateWeightedMemberCovariance(
    const Eigen::MatrixXd & capital_lambda,
    const Eigen::ArrayXd & omega_array)
{
    if (capital_lambda.rows() != capital_lambda.cols())
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
        return std::vector<Eigen::MatrixXd>(
            static_cast<std::size_t>(member_size),
            capital_lambda
        );
    }

    const auto omega_h{ 1.0 / omega_inverse_sum };
    std::vector<Eigen::MatrixXd> member_capital_lambda_list(static_cast<std::size_t>(member_size));
    for (int i = 0; i < member_size; i++)
    {
        member_capital_lambda_list.at(static_cast<std::size_t>(i)) =
            (member_size * omega_h / omega_array(i)) * capital_lambda;
    }
    return member_capital_lambda_list;
}
} // namespace

Eigen::ArrayXd HRLModelAlgorithms::CalculateMemberStatisticalDistance(
    const Eigen::VectorXd & mu_prior,
    const Eigen::MatrixXd & capital_lambda,
    const Eigen::MatrixXd & beta_posterior_array)
{
    if (beta_posterior_array.rows() != mu_prior.rows())
    {
        throw std::invalid_argument("beta_posterior_array and mu_prior sizes are inconsistent.");
    }
    if (capital_lambda.rows() != capital_lambda.cols() ||
        capital_lambda.rows() != mu_prior.rows())
    {
        throw std::invalid_argument("capital_lambda size is inconsistent.");
    }

    const auto member_size{ static_cast<int>(beta_posterior_array.cols()) };
    Eigen::ArrayXd statistical_distance_array{ Eigen::ArrayXd::Zero(member_size) };
    const Eigen::MatrixXd error_array{ beta_posterior_array.colwise() - mu_prior };
    const auto inv_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(capital_lambda) };
    for (int i = 0; i < member_size; i++)
    {
        statistical_distance_array(i) =
            error_array.col(i).transpose() * inv_capital_lambda * error_array.col(i);
    }
    return statistical_distance_array;
}

Eigen::Array<bool, Eigen::Dynamic, 1> HRLModelAlgorithms::CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array)
{
    if (basis_size <= 0)
    {
        throw std::invalid_argument("basis_size must be positive.");
    }
    return statistical_distance_array > CalculateChiSquareQuantile(basis_size);
}
