#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>

#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "utils/hrl/detail/ScopedEigenThreadCount.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

#include <boost/math/distributions/chi_squared.hpp>

RHBMMemberDataset rhbm_gem::rhbm_helper::BuildMemberDataset(
    const SeriesPointList & series_point_list)
{
    if (series_point_list.empty())
    {
        throw std::invalid_argument("series_point_list must not be empty.");
    }

    const auto basis_size{ static_cast<int>(series_point_list.front().GetBasisSize()) };
    rhbm_gem::numeric_validation::RequirePositive(basis_size, "basis_size");
    rhbm_gem::numeric_validation::RequireAtMost(
        series_point_list.size(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "series_point_list size");
    const auto data_size{ static_cast<int>(series_point_list.size()) };

    RHBMMemberDataset dataset;
    dataset.X = RHBMDesignMatrix::Zero(data_size, basis_size);
    dataset.y = RHBMResponseVector::Zero(data_size);
    dataset.score = RHBMScoreVector::Zero(data_size);
    for (int i = 0; i < data_size; i++)
    {
        const auto & point{ series_point_list.at(static_cast<std::size_t>(i)) };
        if (static_cast<int>(point.GetBasisSize()) != basis_size)
        {
            throw std::invalid_argument("All data entries must share the same basis size.");
        }
        rhbm_gem::numeric_validation::RequireFinite(
            point.response,
            "response",
            "Member dataset contains non-finite value.");
        rhbm_gem::numeric_validation::RequireFinite(
            point.score,
            "score",
            "Member dataset contains non-finite value.");
        rhbm_gem::numeric_validation::RequireAllFinite(
            point.basis_list,
            "basis_list",
            "Member dataset contains non-finite value.");

        dataset.X.row(i) = Eigen::Map<const Eigen::RowVectorXd>(
            point.basis_list.data(), basis_size);
        dataset.y(i) = point.response;
        dataset.score(i) = point.score;
    }

    return dataset;
}

RHBMBetaMatrix rhbm_gem::rhbm_helper::BuildBetaMatrix(
    const std::vector<RHBMBetaVector> & beta_list)
{
    if (beta_list.empty())
    {
        throw std::invalid_argument("beta_list must not be empty.");
    }

    const auto basis_size{ static_cast<int>(beta_list.front().rows()) };
    rhbm_gem::numeric_validation::RequirePositive(basis_size, "basis_size");
    rhbm_gem::numeric_validation::RequireAtMost(
        beta_list.size(),
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "beta_list size");
    const auto member_size{ static_cast<int>(beta_list.size()) };
    RHBMBetaMatrix beta_matrix{ RHBMBetaMatrix::Zero(basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const auto & beta_vector{ beta_list.at(static_cast<std::size_t>(i)) };
        rhbm_gem::eigen_validation::RequireVectorSize(
            beta_vector,
            basis_size,
            "beta",
            "All beta vectors must share the same basis size.");
        rhbm_gem::eigen_validation::RequireFinite(
            beta_vector,
            "beta",
            "beta_list contains non-finite value.");
        beta_matrix.col(i) = beta_vector;
    }
    return beta_matrix;
}

RHBMGroupEstimationInput rhbm_gem::rhbm_helper::BuildGroupInput(
    const std::vector<RHBMMemberDataset> & member_datasets,
    const std::vector<RHBMBetaEstimateResult> & member_fit_results)
{
    const auto member_size{ member_datasets.size() };
    if (member_size == 0)
    {
        throw std::invalid_argument("member_datasets must not be empty.");
    }
    if (member_fit_results.size() != member_size)
    {
        throw std::invalid_argument("Group estimation inputs must have consistent member counts.");
    }

    const auto basis_size{ static_cast<int>(member_datasets.front().X.cols()) };
    rhbm_gem::numeric_validation::RequirePositive(basis_size, "basis_size");

    RHBMGroupEstimationInput input;
    input.basis_size = basis_size;
    input.member_datasets.reserve(member_size);
    input.member_fit_results.reserve(member_size);

    for (std::size_t i = 0; i < member_size; i++)
    {
        const auto & dataset{ member_datasets.at(i) };
        const auto & fit_result{ member_fit_results.at(i) };

        if (dataset.X.cols() != basis_size)
        {
            throw std::invalid_argument("Member dataset basis size is inconsistent.");
        }
        rhbm_gem::eigen_validation::RequireVectorSize(
            dataset.y,
            dataset.X.rows(),
            "Member dataset response",
            "Member dataset shape is inconsistent.");
        rhbm_gem::eigen_validation::RequireSameSize(
            dataset.score,
            dataset.y,
            "Member dataset shape",
            "Member dataset shape is inconsistent.");
        rhbm_gem::eigen_validation::RequireVectorSize(
            fit_result.beta_mdpde,
            basis_size,
            "Member fit beta",
            "Member beta basis size is inconsistent.");
        rhbm_gem::eigen_validation::RequireSameSize(
            fit_result.data_weight.diagonal(),
            dataset.y,
            "Member fit data weight",
            "Member covariance or weight size is inconsistent.");
        rhbm_gem::eigen_validation::RequireSameSize(
            fit_result.data_covariance.diagonal(),
            dataset.y,
            "Member fit data covariance",
            "Member covariance or weight size is inconsistent.");

        input.member_datasets.emplace_back(dataset);
        input.member_fit_results.emplace_back(fit_result);
    }

    return input;
}

namespace
{
void ValidateDatasetShape(
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector,
    std::string_view context = {})
{
    rhbm_gem::eigen_validation::RequireRows(
        design_matrix,
        response_vector.rows(),
        "X",
        context);
    rhbm_gem::numeric_validation::RequirePositive(
        design_matrix.cols(),
        "X column count",
        context);
    rhbm_gem::eigen_validation::RequireFinite(design_matrix, "X", context);
    rhbm_gem::eigen_validation::RequireFinite(response_vector, "y", context);
}

void ValidateDiagonalAgainstSize(
    const RHBMDiagonalMatrix & matrix,
    Eigen::Index expected_size,
    const char * name,
    std::string_view context = {})
{
    rhbm_gem::eigen_validation::RequireVectorSize(
        matrix.diagonal(),
        expected_size,
        name,
        context);
}

template <typename Derived>
void ValidateSquareBasisMatrix(
    const Eigen::DenseBase<Derived> & matrix,
    Eigen::Index basis_size,
    const char * name,
    std::string_view context = {})
{
    rhbm_gem::eigen_validation::RequireShape(
        matrix,
        basis_size,
        basis_size,
        name,
        context);
}

double CalculateChiSquareQuantile(int df)
{
    const boost::math::chi_squared_distribution<double> distribution(static_cast<double>(df));
    return boost::math::quantile(distribution, 0.99);
}

RHBMBetaVector CalculateBetaByOLS(
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector);

[[maybe_unused]] RHBMBetaVector CalculateBetaWithSelectedDataByOLS(
    const RHBMMemberDataset & dataset);

RHBMBetaVector CalculateBetaByMDPDE(
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector,
    const RHBMDiagonalMatrix & W);

RHBMMuVector CalculateMuByMedian(
    const RHBMBetaMatrix & beta_matrix);

RHBMMuVector CalculateMuByMDPDE(
    const RHBMBetaMatrix & beta_matrix,
    const Eigen::ArrayXd & omega_array,
    double omega_sum);

RHBMDiagonalMatrix CalculateDataWeight(
    double alpha,
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector,
    const RHBMBetaVector & beta,
    double sigma_square,
    double weight_min);

double CalculateDataVarianceSquare(
    double alpha,
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector,
    const RHBMDiagonalMatrix & W,
    const RHBMBetaVector & beta);

RHBMDiagonalMatrix CalculateDataCovariance(
    double sigma_square,
    const RHBMDiagonalMatrix & W);

Eigen::ArrayXd CalculateMemberWeight(
    double alpha,
    const RHBMBetaMatrix & beta_matrix,
    const RHBMMuVector & mu,
    const RHBMGroupCovarianceMatrix & capital_lambda,
    double weight_min);

RHBMGroupCovarianceMatrix CalculateMemberCovariance(
    double alpha,
    const RHBMBetaMatrix & beta_matrix,
    const RHBMMuVector & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum);

std::vector<RHBMMemberCovarianceMatrix> CalculateWeightedMemberCovariance(
    const RHBMGroupCovarianceMatrix & capital_lambda,
    const Eigen::ArrayXd & omega_array);
} // namespace

RHBMBetaEstimateResult rhbm_gem::rhbm_helper::EstimateBetaMDPDE(
    double alpha_r,
    const RHBMMemberDataset & dataset,
    const RHBMExecutionOptions & options)
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha");
    rhbm_gem::numeric_validation::RequirePositive(options.max_iterations, "max_iterations");
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(options.tolerance, "tolerance");
    rhbm_gem::numeric_validation::RequireFinitePositive(options.data_weight_min, "data_weight_min");
    const auto & design_matrix{ dataset.X };
    const auto & response_vector{ dataset.y };
    ValidateDatasetShape(
        design_matrix,
        response_vector,
        "EstimateBetaMDPDE dataset is invalid.");

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    RHBMBetaEstimateResult result;
    const auto data_size{ static_cast<int>(response_vector.size()) };
    if (data_size <= 1)
    {
        result.status = RHBMEstimationStatus::INSUFFICIENT_DATA;
        result.beta_ols = RHBMBetaVector::Zero(design_matrix.cols());
        result.beta_mdpde = RHBMBetaVector::Zero(design_matrix.cols());
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
    result.data_weight = RHBMResponseVector::Ones(data_size).asDiagonal();

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
        result.status = RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    }
    if (result.sigma_square == std::numeric_limits<double>::max())
    {
        result.status = RHBMEstimationStatus::NUMERICAL_FALLBACK;
    }
    return result;
}

RHBMMuEstimateResult rhbm_gem::rhbm_helper::EstimateMuMDPDE(
    double alpha_g,
    const RHBMBetaMatrix & beta_matrix,
    const RHBMExecutionOptions & options)
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha");
    rhbm_gem::numeric_validation::RequirePositive(options.max_iterations, "max_iterations");
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(options.tolerance, "tolerance");
    rhbm_gem::numeric_validation::RequireFinitePositive(options.member_weight_min, "member_weight_min");
    rhbm_gem::eigen_validation::RequireNonEmpty(beta_matrix, "beta_matrix");

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    RHBMMuEstimateResult result;
    const auto basis_size{ static_cast<int>(beta_matrix.rows()) };
    const auto member_size{ static_cast<int>(beta_matrix.cols()) };
    if (member_size == 1)
    {
        result.status = RHBMEstimationStatus::SINGLE_MEMBER;
        result.mu_mean = beta_matrix.rowwise().mean();
        result.mu_mdpde = result.mu_mean;
        result.omega_array = Eigen::ArrayXd::Ones(member_size);
        result.omega_sum = result.omega_array.sum();
        result.capital_lambda = RHBMGroupCovarianceMatrix::Identity(basis_size, basis_size);
        result.member_capital_lambda_list.assign(
            static_cast<std::size_t>(member_size),
            RHBMMemberCovarianceMatrix::Identity(basis_size, basis_size)
        );
        return result;
    }

    result.mu_mean = beta_matrix.rowwise().mean();
    result.mu_mdpde = CalculateMuByMedian(beta_matrix);
    result.capital_lambda = RHBMGroupCovarianceMatrix::Identity(basis_size, basis_size);
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
        result.status = RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    }
    try
    {
        rhbm_gem::eigen_validation::RequireFinite(result.capital_lambda, "capital_lambda");
    }
    catch (const std::invalid_argument &)
    {
        result.status = RHBMEstimationStatus::NUMERICAL_FALLBACK;
    }
    return result;
}

RHBMWebEstimateResult rhbm_gem::rhbm_helper::EstimateWEB(
    const std::vector<RHBMMemberDataset> & member_datasets,
    const std::vector<RHBMDiagonalMatrix> & capital_sigma_list,
    const RHBMMuVector & mu_mdpde,
    const std::vector<RHBMMemberCovarianceMatrix> & member_capital_lambda_list,
    const RHBMExecutionOptions & options)
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
    rhbm_gem::eigen_validation::RequireNonEmpty(mu_mdpde, "mu_mdpde");

    detail::ScopedEigenThreadCount thread_guard(options.thread_size);
    (void)thread_guard;

    RHBMWebEstimateResult result;
    const auto basis_size{ static_cast<int>(mu_mdpde.rows()) };
    const auto member_size{ static_cast<int>(member_datasets.size()) };
    result.mu_prior = RHBMMuVector::Zero(basis_size);
    result.beta_posterior_matrix = RHBMBetaPosteriorMatrix::Zero(basis_size, member_size);
    result.capital_sigma_posterior_list.clear();
    result.capital_sigma_posterior_list.reserve(static_cast<std::size_t>(member_size));
    if (member_size == 1)
    {
        result.status = RHBMEstimationStatus::SINGLE_MEMBER;
        return result;
    }

    RHBMMuVector numerator{ RHBMMuVector::Zero(basis_size) };
    RHBMGroupCovarianceMatrix denominator{
        RHBMGroupCovarianceMatrix::Zero(basis_size, basis_size)
    };
    for (std::size_t i = 0; i < static_cast<std::size_t>(member_size); i++)
    {
        const auto & dataset{ member_datasets.at(i) };
        const auto & capital_sigma{ capital_sigma_list.at(i) };
        const auto & member_capital_lambda{ member_capital_lambda_list.at(i) };
        ValidateDatasetShape(
            dataset.X,
            dataset.y,
            "EstimateWEB dataset is invalid.");
        ValidateDiagonalAgainstSize(
            capital_sigma,
            dataset.y.size(),
            "capital_sigma",
            "EstimateWEB covariance input is invalid.");
        ValidateSquareBasisMatrix(
            member_capital_lambda,
            basis_size,
            "member_capital_lambda",
            "EstimateWEB member covariance input is invalid.");

        const auto inv_capital_sigma{
            rhbm_gem::eigen_helper::GetInverseDiagonalMatrix(capital_sigma)
        };
        const auto inv_member_capital_lambda{
            rhbm_gem::eigen_helper::GetInverseMatrix(member_capital_lambda)
        };
        const RHBMGroupCovarianceMatrix gram_matrix{
            dataset.X.transpose() * inv_capital_sigma * dataset.X
        };
        const RHBMMuVector moment_matrix{
            dataset.X.transpose() * inv_capital_sigma * dataset.y
        };
        const RHBMGroupCovarianceMatrix inv_capital_sigma_posterior{
            gram_matrix + inv_member_capital_lambda
        };
        const RHBMPosteriorCovarianceMatrix capital_sigma_posterior{
            rhbm_gem::eigen_helper::GetInverseMatrix(inv_capital_sigma_posterior)
        };

        result.capital_sigma_posterior_list.emplace_back(capital_sigma_posterior);
        result.beta_posterior_matrix.col(static_cast<Eigen::Index>(i)) =
            capital_sigma_posterior * (moment_matrix + inv_member_capital_lambda * mu_mdpde);
        numerator += inv_member_capital_lambda * capital_sigma_posterior * moment_matrix;
        denominator += inv_member_capital_lambda * capital_sigma_posterior * gram_matrix;
    }

    result.mu_prior = rhbm_gem::eigen_helper::GetInverseMatrix(denominator) * numerator;
    if (member_size == 2)
    {
        result.mu_prior = mu_mdpde;
    }
    return result;
}

namespace
{
RHBMBetaVector CalculateBetaByOLS(
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector)
{
    const RHBMGroupCovarianceMatrix gram_matrix{ design_matrix.transpose() * design_matrix };
    const RHBMGroupCovarianceMatrix inverse_gram_matrix{
        rhbm_gem::eigen_helper::GetInverseMatrix(gram_matrix)
    };
    return inverse_gram_matrix * (design_matrix.transpose() * response_vector);
}

[[maybe_unused]] RHBMBetaVector CalculateBetaWithSelectedDataByOLS(
    const RHBMMemberDataset & dataset)
{
    ValidateDatasetShape(
        dataset.X,
        dataset.y,
        "Selected-data OLS input is invalid.");
    rhbm_gem::eigen_validation::RequireSameSize(
        dataset.score,
        dataset.y,
        "dataset score",
        "Selected-data OLS input is invalid.");

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
    RHBMDesignMatrix selected_design_matrix{
        RHBMDesignMatrix::Zero(selected_row_count, dataset.X.cols())
    };
    RHBMResponseVector selected_response_vector{ RHBMResponseVector::Zero(selected_row_count) };
    for (Eigen::Index i = 0; i < selected_row_count; i++)
    {
        const auto row{ selected_row_indices.at(static_cast<std::size_t>(i)) };
        selected_design_matrix.row(i) = dataset.X.row(row);
        selected_response_vector(i) = dataset.y(row);
    }

    return CalculateBetaByOLS(selected_design_matrix, selected_response_vector);
}

RHBMBetaVector CalculateBetaByMDPDE(
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector,
    const RHBMDiagonalMatrix & W)
{
    ValidateDatasetShape(
        design_matrix,
        response_vector,
        "MDPDE beta update input is invalid.");
    ValidateDiagonalAgainstSize(
        W,
        response_vector.size(),
        "W",
        "MDPDE beta update input is invalid.");
    const RHBMGroupCovarianceMatrix gram_matrix{
        design_matrix.transpose() * W * design_matrix
    };
    const auto inverse_gram_matrix{ rhbm_gem::eigen_helper::GetInverseMatrix(gram_matrix) };
    return inverse_gram_matrix * (design_matrix.transpose() * W * response_vector);
}

RHBMMuVector CalculateMuByMedian(const RHBMBetaMatrix & beta_matrix)
{
    rhbm_gem::eigen_validation::RequireNonEmpty(beta_matrix, "beta_matrix");
    const auto basis_size{ static_cast<int>(beta_matrix.rows()) };
    RHBMMuVector mu{ RHBMMuVector::Zero(basis_size) };
    for (int b = 0; b < basis_size; b++)
    {
        mu(b) = rhbm_gem::eigen_helper::GetMedian(beta_matrix.row(b));
    }
    return mu;
}

RHBMMuVector CalculateMuByMDPDE(
    const RHBMBetaMatrix & beta_matrix,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    rhbm_gem::eigen_validation::RequireCols(
        beta_matrix,
        omega_array.rows(),
        "beta_matrix",
        "Mu estimation input is invalid.");
    rhbm_gem::numeric_validation::RequireFinitePositive(omega_sum, "omega_sum");
    RHBMBetaMatrix numerator{ beta_matrix.array() / omega_sum };
    for (int i = 0; i < numerator.cols(); i++)
    {
        numerator.col(i) *= omega_array(i);
    }
    return numerator.rowwise().sum();
}

RHBMDiagonalMatrix CalculateDataWeight(
    double alpha,
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector,
    const RHBMBetaVector & beta,
    double sigma_square,
    double weight_min)
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(alpha, "alpha");
    rhbm_gem::numeric_validation::RequireFinitePositive(weight_min, "weight_min");
    ValidateDatasetShape(
        design_matrix,
        response_vector,
        "Data weight input is invalid.");
    rhbm_gem::eigen_validation::RequireVectorSize(
        beta,
        design_matrix.cols(),
        "beta",
        "Data weight input is invalid.");

    if (alpha == 0.0)
    {
        return RHBMResponseVector::Ones(response_vector.size()).asDiagonal();
    }
    if (!rhbm_gem::numeric_validation::IsFinitePositive(sigma_square))
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
        return rhbm_gem::numeric_validation::IsFinite(w) ? w : weight_min;
    });
    W = W.cwiseMax(weight_min);
    return W.matrix().asDiagonal();
}

double CalculateDataVarianceSquare(
    double alpha,
    const RHBMDesignMatrix & design_matrix,
    const RHBMResponseVector & response_vector,
    const RHBMDiagonalMatrix & W,
    const RHBMBetaVector & beta)
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(alpha, "alpha");
    ValidateDatasetShape(
        design_matrix,
        response_vector,
        "Data variance input is invalid.");
    ValidateDiagonalAgainstSize(
        W,
        response_vector.size(),
        "W",
        "Data variance input is invalid.");
    rhbm_gem::eigen_validation::RequireVectorSize(
        beta,
        design_matrix.cols(),
        "beta",
        "Data variance input is invalid.");

    const auto n{ static_cast<double>(response_vector.size()) };
    const RHBMResponseVector residual{ response_vector - (design_matrix * beta) };
    const auto numerator{ static_cast<double>(residual.transpose() * W * residual) };
    auto denominator{ W.diagonal().sum() - n * alpha * std::pow(1.0 + alpha, -1.5) };
    if (denominator <= 0.0)
    {
        denominator = std::numeric_limits<double>::min();
    }
    auto sigma_square{ numerator / denominator };
    if (!rhbm_gem::numeric_validation::IsFinitePositive(sigma_square))
    {
        sigma_square = std::numeric_limits<double>::max();
    }
    return sigma_square;
}

RHBMDiagonalMatrix CalculateDataCovariance(
    double sigma_square,
    const RHBMDiagonalMatrix & W)
{
    const RHBMResponseVector data_weight_array{ W.diagonal() };
    const auto n{ static_cast<int>(data_weight_array.size()) };
    const auto W_inverse_trace{
        rhbm_gem::eigen_helper::GetInverseDiagonalMatrix(W).diagonal().sum()
    };
    if (!rhbm_gem::numeric_validation::IsFinitePositive(W_inverse_trace))
    {
        return RHBMResponseVector::Ones(n).asDiagonal();
    }

    RHBMResponseVector capital_sigma{ RHBMResponseVector::Zero(n) };
    for (int j = 0; j < n; j++)
    {
        if (data_weight_array(j) == 0.0 || W_inverse_trace == 0.0)
        {
            continue;
        }
        capital_sigma(j) = n * sigma_square / data_weight_array(j) / W_inverse_trace;
        if (!rhbm_gem::numeric_validation::IsFinitePositive(capital_sigma(j)))
        {
            capital_sigma(j) = 1.0;
        }
    }
    return capital_sigma.asDiagonal();
}

Eigen::ArrayXd CalculateMemberWeight(
    double alpha,
    const RHBMBetaMatrix & beta_matrix,
    const RHBMMuVector & mu,
    const RHBMGroupCovarianceMatrix & capital_lambda,
    double weight_min)
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(alpha, "alpha");
    rhbm_gem::numeric_validation::RequireFinitePositive(weight_min, "weight_min");
    rhbm_gem::eigen_validation::RequireRows(
        beta_matrix,
        mu.rows(),
        "beta_matrix",
        "Member weight input is invalid.");
    ValidateSquareBasisMatrix(
        capital_lambda,
        beta_matrix.rows(),
        "capital_lambda",
        "Member weight input is invalid.");

    const auto member_size{ static_cast<int>(beta_matrix.cols()) };
    const auto weight_member_min{ weight_min / static_cast<double>(member_size) };
    const auto inverse_capital_lambda{ rhbm_gem::eigen_helper::GetInverseMatrix(capital_lambda) };
    const RHBMBetaMatrix residual_matrix{ beta_matrix.colwise() - mu };
    Eigen::ArrayXd omega_array{ Eigen::ArrayXd::Zero(member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const RHBMBetaVector residual{ residual_matrix.col(i) };
        const auto exp_index{
            static_cast<double>(residual.transpose() * inverse_capital_lambda * residual)
        };
        if (!rhbm_gem::numeric_validation::IsFinite(exp_index))
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

RHBMGroupCovarianceMatrix CalculateMemberCovariance(
    double alpha,
    const RHBMBetaMatrix & beta_matrix,
    const RHBMMuVector & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    rhbm_gem::numeric_validation::RequireFiniteNonNegative(alpha, "alpha");
    rhbm_gem::eigen_validation::RequireRows(
        beta_matrix,
        mu.rows(),
        "beta_matrix",
        "Member covariance input is invalid.");
    rhbm_gem::eigen_validation::RequireCols(
        beta_matrix,
        omega_array.rows(),
        "beta_matrix",
        "Member covariance input is invalid.");
    const auto basis_size{ static_cast<int>(beta_matrix.rows()) };
    const auto member_size{ static_cast<int>(beta_matrix.cols()) };
    RHBMGroupCovarianceMatrix numerator{
        RHBMGroupCovarianceMatrix::Zero(basis_size, basis_size)
    };
    const double denominator{
        omega_sum - member_size * alpha * std::pow(1.0 + alpha, -1.0 - 0.5 * basis_size)
    };
    if (denominator <= 0.0)
    {
        return RHBMGroupCovarianceMatrix::Identity(basis_size, basis_size);
    }

    const RHBMBetaMatrix residual_matrix{ beta_matrix.colwise() - mu };
    for (int i = 0; i < member_size; i++)
    {
        const RHBMBetaVector residual{ residual_matrix.col(i) };
        numerator += omega_array(i) * (residual * residual.transpose());
    }

    RHBMGroupCovarianceMatrix capital_lambda{ numerator / denominator };
    try
    {
        rhbm_gem::eigen_validation::RequireFinite(capital_lambda, "capital_lambda");
    }
    catch (const std::invalid_argument &)
    {
        capital_lambda = RHBMGroupCovarianceMatrix::Identity(basis_size, basis_size);
    }
    return capital_lambda;
}

std::vector<RHBMMemberCovarianceMatrix> CalculateWeightedMemberCovariance(
    const RHBMGroupCovarianceMatrix & capital_lambda,
    const Eigen::ArrayXd & omega_array)
{
    rhbm_gem::eigen_validation::RequireSquare(
        capital_lambda,
        "capital_lambda",
        "Weighted member covariance input is invalid.");
    rhbm_gem::eigen_validation::RequireNonEmpty(
        omega_array,
        "omega_array",
        "Weighted member covariance input is invalid.");
    const auto member_size{ static_cast<int>(omega_array.size()) };

    auto omega_inverse_sum{ 0.0 };
    for (int i = 0; i < member_size; i++)
    {
        omega_inverse_sum += (omega_array(i) == 0.0) ? 0.0 : 1.0 / omega_array(i);
    }
    if (omega_inverse_sum <= 0.0)
    {
        return std::vector<RHBMMemberCovarianceMatrix>(
            static_cast<std::size_t>(member_size),
            capital_lambda
        );
    }

    const auto omega_h{ 1.0 / omega_inverse_sum };
    std::vector<RHBMMemberCovarianceMatrix> member_capital_lambda_list(
        static_cast<std::size_t>(member_size));
    for (int i = 0; i < member_size; i++)
    {
        member_capital_lambda_list.at(static_cast<std::size_t>(i)) =
            (member_size * omega_h / omega_array(i)) * capital_lambda;
    }
    return member_capital_lambda_list;
}
} // namespace

Eigen::ArrayXd rhbm_gem::rhbm_helper::CalculateMemberStatisticalDistance(
    const RHBMMuVector & mu_prior,
    const RHBMGroupCovarianceMatrix & capital_lambda,
    const RHBMBetaPosteriorMatrix & beta_posterior_matrix)
{
    rhbm_gem::eigen_validation::RequireRows(
        beta_posterior_matrix,
        mu_prior.rows(),
        "beta_posterior_matrix",
        "Statistical distance input is invalid.");
    ValidateSquareBasisMatrix(
        capital_lambda,
        mu_prior.rows(),
        "capital_lambda",
        "Statistical distance input is invalid.");

    const auto member_size{ static_cast<int>(beta_posterior_matrix.cols()) };
    Eigen::ArrayXd statistical_distance_array{ Eigen::ArrayXd::Zero(member_size) };
    const RHBMBetaPosteriorMatrix error_matrix{ beta_posterior_matrix.colwise() - mu_prior };
    const auto inv_capital_lambda{ rhbm_gem::eigen_helper::GetInverseMatrix(capital_lambda) };
    for (int i = 0; i < member_size; i++)
    {
        statistical_distance_array(i) =
            error_matrix.col(i).transpose() * inv_capital_lambda * error_matrix.col(i);
    }
    return statistical_distance_array;
}

Eigen::Array<bool, Eigen::Dynamic, 1> rhbm_gem::rhbm_helper::CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array)
{
    rhbm_gem::numeric_validation::RequirePositive(basis_size, "basis_size");
    return statistical_distance_array > CalculateChiSquareQuantile(basis_size);
}
