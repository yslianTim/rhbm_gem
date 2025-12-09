#include "HRLModelHelper.hpp"
#include "EigenMatrixUtility.hpp"
#include "Logger.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <stdexcept>
#include <sstream>
#include <random>

using std::string;
using std::vector;
using Eigen::ArrayXf;
using Eigen::ArrayXd;
using Eigen::ArrayXXf;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using DMatrixXd = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;
using ArrayXb = Eigen::Array<bool, Eigen::Dynamic, 1>;

namespace
{
int ValidatePositive(int value, const char * name)
{
    if (value <= 0)
    {
        throw std::invalid_argument(std::string(name) + " must be positive value");
    }
    return value;
}

template <typename SizeType>
int CheckedCastToInt(SizeType value)
{
    if (value > static_cast<SizeType>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("data_array size exceeds maximum int");
    }
    return static_cast<int>(value);
}
} // namespace

HRLModelHelper::HRLModelHelper(int basis_size, int member_size) :
    m_basis_size{ ValidatePositive(basis_size, "basis_size") },
    m_member_size{ ValidatePositive(member_size, "member_size") },
    m_sigma_square_array{ ArrayXd::Zero(m_member_size) },
    m_omega_array{ ArrayXd::Ones(m_member_size) },
    m_statistical_distance_array{ ArrayXd::Zero(m_member_size) },
    m_outlier_flag_array{ ArrayXb::Constant(m_member_size, false) },
    m_capital_lambda{ MatrixXd::Identity(m_basis_size, m_basis_size) },
    m_mu_MDPDE{ VectorXd::Zero(m_basis_size) },
    m_mu_prior{ VectorXd::Zero(m_basis_size) },
    m_mu_mean{ VectorXd::Zero(m_basis_size) },
    m_beta_MDPDE_array{ MatrixXd::Zero(m_basis_size, m_member_size) },
    m_beta_posterior_array{ MatrixXd::Zero(m_basis_size, m_member_size) }
{
    auto size{ static_cast<size_t>(m_member_size) };
    m_X_list.reserve(size);
    m_y_list.reserve(size);
    m_W_list.reserve(size);
    m_capital_sigma_list.reserve(size);
    m_capital_sigma_posterior_list.assign(size, MatrixXd::Identity(m_basis_size, m_basis_size));
    m_capital_lambda_list.assign(size, MatrixXd::Identity(m_basis_size, m_basis_size));
}

HRLModelHelper::~HRLModelHelper()
{
}

void HRLModelHelper::SetThreadSize(int thread_size)
{
    if (Eigen::nbThreads() == thread_size) return;
    if (thread_size <= 0)
    {
        Logger::Log(LogLevel::Warning, "thread_size must be positive value, resetting to 1");
        thread_size = 1;
    }
    Eigen::setNbThreads(thread_size);
    Logger::Log(LogLevel::Debug, "Thread size = " + std::to_string(Eigen::nbThreads()));
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd> HRLModelHelper::BuildBasisVectorAndResponseArray(
    const std::vector<Eigen::VectorXd> & data_vector,
    bool quiet_mode)
{
    auto data_size{ CheckedCastToInt(data_vector.size()) };
    if (data_size == 0 && quiet_mode == false)
    {
        Logger::Log(LogLevel::Warning,
            "HRLModelHelper::BuildBasisVectorAndResponseArray : dataset is empty.");
    }
    
    auto basis_size{ static_cast<int>(data_vector[0].rows() - 1) };
    MatrixXd x_data_matrix{ MatrixXd::Zero(data_size, basis_size) };
    VectorXd y_data_vector{ VectorXd::Zero(data_size) };
    for (int i = 0; i < data_size; i++)
    {
        const auto & data{ data_vector.at(static_cast<size_t>(i)) };
        if (!data.array().allFinite())
        {
            throw std::invalid_argument("Member dataset contains non-finite value.");
        }
        x_data_matrix.row(i) = data.head(basis_size);
        y_data_vector(i) = data(basis_size);
    }

    return std::make_tuple(std::move(x_data_matrix), std::move(y_data_vector));
}

Eigen::MatrixXd HRLModelHelper::ConvertBetaListToMatrix(
    const std::vector<Eigen::VectorXd> & beta_list,
    bool quiet_mode)
{
    auto basis_size{ static_cast<int>(beta_list.front().rows()) };
    auto member_size{ CheckedCastToInt(beta_list.size()) };
    if (member_size == 0 && quiet_mode == false)
    {
        Logger::Log(LogLevel::Warning,
            "HRLModelHelper::BuildBetaArray : dataset is empty.");
    }

    MatrixXd beta_array{ MatrixXd::Zero(basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const auto & beta{ beta_list.at(static_cast<size_t>(i)) };
        beta_array.col(i) = beta;
    }
    return beta_array;
}

void HRLModelHelper::SetMemberDataEntriesList(
    const std::vector<std::vector<Eigen::VectorXd>> & data_list)
{
    auto member_size{ CheckedCastToInt(data_list.size()) };
    ValidateMemberSize(member_size);
    m_X_list.clear();
    m_y_list.clear();
    for (auto & data_entry : data_list)
    {
        ValidateBasisSize(static_cast<int>(data_entry[0].rows() - 1));
        auto data{ BuildBasisVectorAndResponseArray(data_entry) };
        m_X_list.emplace_back(std::get<0>(data));
        m_y_list.emplace_back(std::get<1>(data));
    }
}

void HRLModelHelper::SetMemberBetaMDPDEList(
    const std::vector<Eigen::VectorXd> & beta_list,
    const std::vector<double> & sigma_square_list,
    const std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> & W_list,
    const std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> & capital_sigma_list)
{
    auto member_size{ CheckedCastToInt(beta_list.size()) };
    ValidateMemberSize(member_size);
    m_beta_MDPDE_array.setZero();
    m_sigma_square_array.setZero();
    m_W_list.clear();
    m_capital_sigma_list.clear();
    for (int i = 0; i < member_size; i++)
    {
        ValidateBasisSize(static_cast<int>(beta_list.at(static_cast<size_t>(i)).rows()));
        m_beta_MDPDE_array.col(i) = beta_list.at(static_cast<size_t>(i));
        m_sigma_square_array(i) = sigma_square_list.at(static_cast<size_t>(i));
        m_W_list.emplace_back(W_list.at(static_cast<size_t>(i)));
        m_capital_sigma_list.emplace_back(capital_sigma_list.at(static_cast<size_t>(i)));
    }
}

void HRLModelHelper::RunGroupEstimation(
    double alpha_g, bool quiet_mode, int thread_size, int iteration, double tolerance)
{
    VectorXd mu_median;
    VectorXd mu_MDPDE;
    ArrayXd omega_array;
    double omega_sum;
    MatrixXd capital_lambda;
    std::vector<MatrixXd> member_capital_lambda_list;
    AlgorithmMuMDPDE(
        alpha_g, m_beta_MDPDE_array, mu_median, mu_MDPDE,
        omega_array, omega_sum, capital_lambda, member_capital_lambda_list,
        quiet_mode, thread_size, iteration, tolerance
    );

    m_mu_mean = mu_median;
    m_mu_MDPDE = mu_MDPDE;
    m_omega_array = omega_array;
    m_capital_lambda = capital_lambda;
    m_capital_lambda_list = std::move(member_capital_lambda_list);

    Eigen::VectorXd mu_prior;
    Eigen::MatrixXd beta_posterior_array;
    std::vector<Eigen::MatrixXd> capital_sigma_posterior_list;
    auto web_estimator{ AlgorithmWEB(
        m_X_list, m_y_list, m_capital_sigma_list, m_mu_MDPDE, m_capital_lambda_list,
        mu_prior, beta_posterior_array, capital_sigma_posterior_list,
        quiet_mode, thread_size)
    };

    if (web_estimator == false)
    {
        mu_prior = mu_MDPDE;
        beta_posterior_array = m_beta_MDPDE_array;
        for (int i = 0; i < m_member_size; i++)
        {
            capital_sigma_posterior_list.emplace_back(
                Eigen::MatrixXd::Zero(m_basis_size, m_basis_size)
            );
        }
    }

    m_mu_prior = mu_prior;
    m_beta_posterior_array = beta_posterior_array;
    m_capital_sigma_posterior_list = std::move(capital_sigma_posterior_list);

    m_statistical_distance_array = CalculateMemberStatisticalDistance(
        m_mu_prior, m_capital_lambda, m_beta_posterior_array
    );
    m_outlier_flag_array = CalculateOutlierMemberFlag(m_basis_size, m_statistical_distance_array);
}

Eigen::VectorXd HRLModelHelper::RunAlphaRTraining(
    const std::vector<Eigen::VectorXd> & data_list,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    std::vector<std::vector<Eigen::VectorXd>> data_subset_list;
    std::vector<std::vector<Eigen::VectorXd>> data_set_test;
    std::vector<std::vector<Eigen::VectorXd>> data_set_training;
    data_subset_list.resize(subset_size);
    size_t total_entry_size{ data_list.size() };
    size_t entries_in_subset_size{ total_entry_size / subset_size + 1 };
    for (size_t i = 0; i < subset_size; i++)
    {
        data_subset_list[i].reserve(entries_in_subset_size);
    }

    size_t count{ 0 };
    for (const auto & entry : data_list)
    {
        data_subset_list[count % subset_size].emplace_back(entry);
        count++;
    }

    data_set_test.resize(subset_size);
    data_set_training.resize(subset_size);
    for (size_t i = 0; i < subset_size; i++)
    {
        size_t test_set_size{ data_subset_list[i].size() };
        size_t training_set_size{ data_list.size() - test_set_size };
        data_set_test[i].reserve(test_set_size);
        data_set_training[i].reserve(training_set_size);
        data_set_test[i].insert(
            data_set_test[i].end(),
            data_subset_list[i].begin(),
            data_subset_list[i].end());
        for (size_t j = 0; j < subset_size; j++)
        {
            if (i == j) continue;
            data_set_training[i].insert(
                data_set_training[i].end(),
                data_subset_list[j].begin(),
                data_subset_list[j].end());
        }
    }

    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    Eigen::VectorXd error_sum_list{ Eigen::VectorXd::Zero(alpha_size) };
    for (int p = 0; p < alpha_size; p++)
    {
        auto alpha{ alpha_list.at(static_cast<size_t>(p)) };
        auto beta_error_sum{ 0.0 };
        for (size_t i = 0; i < subset_size; i++)
        {
            auto data_array_test{ BuildBasisVectorAndResponseArray(data_set_test.at(i)) };
            const auto & X_test{ std::get<0>(data_array_test) };
            const auto & y_test{ std::get<1>(data_array_test) };
            Eigen::VectorXd beta_ols_test;
            Eigen::VectorXd beta_mdpde_test;
            double sigma_square_test;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> W_test;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> capital_sigma_test;
            AlgorithmBetaMDPDE(
                alpha, X_test, y_test, beta_ols_test, beta_mdpde_test,
                sigma_square_test, W_test, capital_sigma_test, true
            );

            auto data_array_training{ BuildBasisVectorAndResponseArray(data_set_training.at(i)) };
            const auto & X_training{ std::get<0>(data_array_training) };
            const auto & y_training{ std::get<1>(data_array_training) };
            Eigen::VectorXd beta_ols_training;
            Eigen::VectorXd beta_mdpde_training;
            double sigma_square_training;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> W_training;
            Eigen::DiagonalMatrix<double, Eigen::Dynamic> capital_sigma_training;
            AlgorithmBetaMDPDE(
                alpha, X_training, y_training, beta_ols_training, beta_mdpde_training,
                sigma_square_training, W_training, capital_sigma_training
            );

            beta_error_sum += (beta_mdpde_test - beta_mdpde_training).norm();
        }
        error_sum_list(p) = beta_error_sum;
    }
    return error_sum_list;
}

Eigen::VectorXd HRLModelHelper::RunAlphaGTraining(
    const std::vector<Eigen::VectorXd> & data_list,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    auto data_size_in_half{ data_list.size() / 2 };
    std::vector<std::vector<Eigen::VectorXd>> data_set_test(subset_size);
    std::vector<std::vector<Eigen::VectorXd>> data_set_training(subset_size);
    for (size_t i = 0; i < subset_size; i++)
    {
        // Randomly pick the half of data into test set and training set for each group
        data_set_test[i].reserve(data_size_in_half);
        data_set_training[i].reserve(data_size_in_half);
        std::vector<Eigen::VectorXd> data_1, data_2;
        std::vector<Eigen::VectorXd> shuffled_data{ data_list };
        std::shuffle(shuffled_data.begin(), shuffled_data.end(), std::mt19937{std::random_device{}()});
        auto diff{ static_cast<std::vector<Eigen::VectorXd>::difference_type>(data_size_in_half) };
        data_1.assign(shuffled_data.begin(), shuffled_data.begin() + diff);
        data_2.assign(shuffled_data.begin() + diff, shuffled_data.end());
        data_set_test[i] = std::move(data_1);
        data_set_training[i] = std::move(data_2);
    }

    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    Eigen::VectorXd error_sum_list{ Eigen::VectorXd::Zero(alpha_size) };
    for (int p = 0; p < alpha_size; p++)
    {
        auto alpha{ alpha_list.at(static_cast<size_t>(p)) };
        auto mu_error_sum{ 0.0 };
        for (size_t i = 0; i < subset_size; i++)
        {
            auto beta_matrix_test{ ConvertBetaListToMatrix(data_set_test.at(i), true) };
            Eigen::VectorXd mu_median_test;
            Eigen::VectorXd mu_mdpde_test;
            Eigen::ArrayXd omega_array_test;
            double omega_sum_test;
            Eigen::MatrixXd capital_lambda_test;
            std::vector<Eigen::MatrixXd> member_capital_lambda_list_test;
            AlgorithmMuMDPDE(
                alpha, beta_matrix_test, mu_median_test, mu_mdpde_test,
                omega_array_test, omega_sum_test, capital_lambda_test,
                member_capital_lambda_list_test, true
            );

            auto beta_matrix_training{ ConvertBetaListToMatrix(data_set_training.at(i), true) };
            Eigen::VectorXd mu_median_training;
            Eigen::VectorXd mu_mdpde_training;
            Eigen::ArrayXd omega_array_training;
            double omega_sum_training;
            Eigen::MatrixXd capital_lambda_training;
            std::vector<Eigen::MatrixXd> member_capital_lambda_list_training;
            AlgorithmMuMDPDE(
                alpha, beta_matrix_training, mu_median_training, mu_mdpde_training,
                omega_array_training, omega_sum_training, capital_lambda_training,
                member_capital_lambda_list_training, true
            );

            mu_error_sum += (mu_mdpde_test - mu_mdpde_training).norm();
        }
        error_sum_list(p) = mu_error_sum;
    }

    return error_sum_list;
}

bool HRLModelHelper::AlgorithmBetaMDPDE(
    double alpha_r,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    Eigen::VectorXd & beta_OLS,
    Eigen::VectorXd & beta,
    double & sigma_square,
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W,
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> & capital_sigma,
    bool quiet_mode,
    int thread_size,
    int iteration,
    double tolerance,
    double weight_min)
{
    SetThreadSize(thread_size);
    ValidateMaximumIteration(iteration);
    ValidateTolerance(tolerance);
    ValidateAlpha(alpha_r);
    auto data_size{ static_cast<int>(y.size()) };
    if (data_size <= 1)
    {
        if (quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::AlgorithmBetaMDPDE : "
                "Dataset is empty, skipping Beta MDPDE estimation.");
        }
        beta_OLS = VectorXd::Zero(X.cols());
        beta = VectorXd::Zero(X.cols());
        sigma_square = std::numeric_limits<double>::max();
        W = Eigen::DiagonalMatrix<double, Eigen::Dynamic>::Identity(1);
        capital_sigma = Eigen::DiagonalMatrix<double, Eigen::Dynamic>::Identity(1);
        return false;
    }

    beta_OLS = CalculateBetaByOLS(X, y);
    beta = beta_OLS;
    sigma_square = (y - (X * beta)).squaredNorm() / (data_size - 1);
    W = VectorXd::Ones(data_size).asDiagonal();
    VectorXd beta_in_previous_iter{ beta };
    for (int t = 0; t < iteration; t++)
    { //=== Begin of iteration loop
        W = CalculateDataWeight(alpha_r, X, y, beta, sigma_square, weight_min);
        beta = CalculateBetaByMDPDE(X, y, W);
        sigma_square = CalculateDataVarianceSquare(alpha_r, X, y, W, beta);
        if ((beta - beta_in_previous_iter).squaredNorm() < tolerance) break;
        beta_in_previous_iter = beta;
        if (t == iteration - 1 && quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::AlgorithmBetaMDPDE : "
                "Reach maximum iterations before achieving tolerance.");
        }
    } //=== End of iteration loop

    capital_sigma = CalculateDataCovariance(sigma_square, W);
    return true;
}

bool HRLModelHelper::AlgorithmMuMDPDE(
    double alpha_g,
    const Eigen::MatrixXd & beta_array,
    Eigen::VectorXd & mu_median,
    Eigen::VectorXd & mu,
    Eigen::ArrayXd & omega_array,
    double & omega_sum,
    Eigen::MatrixXd & capital_lambda,
    std::vector<Eigen::MatrixXd> & member_capital_lambda_list,
    bool quiet_mode,
    int thread_size,
    int iteration,
    double tolerance,
    double weight_min)
{
    SetThreadSize(thread_size);
    ValidateMaximumIteration(iteration);
    ValidateTolerance(tolerance);
    ValidateAlpha(alpha_g);
    auto basis_size{ static_cast<int>(beta_array.rows()) };
    auto member_size{ static_cast<int>(beta_array.cols()) };
    
    if (member_size == 1)
    {
        if (quiet_mode == false)
        {
            Logger::Log(LogLevel::Debug,
                "HRLModelHelper::AlgorithmMuMDPDE : "
                "Only one member is present, skipping Mu MDPDE estimation.");
        }
        mu_median = beta_array.rowwise().mean();
        mu = mu_median;
        omega_array = ArrayXd::Ones(member_size);
        omega_sum = omega_array.sum(); 
        capital_lambda = MatrixXd::Identity(basis_size, basis_size);
        member_capital_lambda_list.assign(
            static_cast<size_t>(member_size), MatrixXd::Identity(member_size, member_size)
        );
        return false;
    }

    mu_median = CalculateMuByMedian(beta_array);
    mu = mu_median;
    capital_lambda = MatrixXd::Identity(basis_size, basis_size);
    VectorXd mu_in_previous_iter{ mu };
    for (int t = 0; t < iteration; t++)
    { //=== Begin of iteration loop
        omega_array = CalculateMemberWeight(alpha_g, beta_array, mu, capital_lambda, weight_min);
        omega_sum = omega_array.sum();
        mu = CalculateMuByMDPDE(beta_array, omega_array, omega_sum);
        capital_lambda = CalculateMemberCovariance(alpha_g, beta_array, mu, omega_array, omega_sum);
        if ((mu - mu_in_previous_iter).squaredNorm() < tolerance) break;
        mu_in_previous_iter = mu;
        if (t == iteration - 1 && quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::AlgorithmMuMDPDE : "
                "Reach maximum iterations before achieving tolerance.");
        }
    } //=== End of iteration loop

    member_capital_lambda_list = CalculateWeightedMemberCovariance(capital_lambda, omega_array);
    return true;
}

bool HRLModelHelper::AlgorithmWEB(
    const std::vector<Eigen::MatrixXd> & X_list,
    const std::vector<Eigen::VectorXd> & y_list,
    const std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> & capital_sigma_list,
    const Eigen::VectorXd & mu_MDPDE,
    const std::vector<Eigen::MatrixXd> & member_capital_lambda_list,
    Eigen::VectorXd & mu_prior,
    Eigen::MatrixXd & beta_posterior_array,
    std::vector<Eigen::MatrixXd> & capital_sigma_posterior_list,
    bool quiet_mode,
    int thread_size)
{
    SetThreadSize(thread_size);
    auto basis_size{ static_cast<int>(mu_MDPDE.rows()) };
    auto member_size{ static_cast<int>(member_capital_lambda_list.size()) };
    mu_prior = VectorXd::Zero(basis_size);
    beta_posterior_array = MatrixXd::Zero(basis_size, member_size);
    capital_sigma_posterior_list.clear();
    capital_sigma_posterior_list.reserve(static_cast<size_t>(member_size));

    if (member_size == 1)
    {
        if (quiet_mode == false)
        {
            Logger::Log(LogLevel::Debug,
                "HRLModelHelper::AlgorithmWEB : "
                "Only one member is present, WEB estimation is skipped.");
        }
        return false;
    }
    
    VectorXd numerator{ VectorXd::Zero(basis_size) };
    MatrixXd denominator{ MatrixXd::Zero(basis_size, basis_size) };
    for (size_t i = 0; i < static_cast<size_t>(member_size); i++)
    {
        const auto & X{ X_list.at(i) };
        const auto & y{ y_list.at(i) };
        const auto & capital_sigma{ capital_sigma_list.at(i) };
        const auto & member_capital_lambda{ member_capital_lambda_list.at(i) };
        auto inv_capital_sigma{ EigenMatrixUtility::GetInverseDiagonalMatrix(capital_sigma) };
        auto inv_member_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(member_capital_lambda) };
        MatrixXd gram_matrix{ X.transpose() * inv_capital_sigma * X };
        VectorXd moment_matrix{ X.transpose() * inv_capital_sigma * y };
        MatrixXd inv_capital_sigma_posterior{ gram_matrix + inv_member_capital_lambda };
        MatrixXd capital_sigma_posterior{
            EigenMatrixUtility::GetInverseMatrix(inv_capital_sigma_posterior)
        };
        capital_sigma_posterior_list.emplace_back(capital_sigma_posterior);
        beta_posterior_array.col(static_cast<int>(i)) =
            capital_sigma_posterior * (moment_matrix + inv_member_capital_lambda * mu_MDPDE);
        numerator += inv_member_capital_lambda * capital_sigma_posterior * moment_matrix;
        denominator += inv_member_capital_lambda * capital_sigma_posterior * gram_matrix;
    }
    
    MatrixXd inv_denominator{ EigenMatrixUtility::GetInverseMatrix(denominator) };
    mu_prior = inv_denominator * numerator;
    if (member_size == 2) mu_prior = mu_MDPDE;
    return true;
}

double HRLModelHelper::CalculateDataVarianceSquare(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W,
    const Eigen::VectorXd & beta)
{
    auto n{ static_cast<double>(y.size()) };
    VectorXd residual{ y - (X * beta) };
    auto numerator{ static_cast<double>(residual.transpose() * W * residual) };
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

Eigen::DiagonalMatrix<double, Eigen::Dynamic> HRLModelHelper::CalculateDataCovariance(
    double sigma_square,
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W)
{
    VectorXd data_weight_array{ W.diagonal() };
    const auto n{ static_cast<int>(data_weight_array.size()) };
    const auto W_inverse_trace{ EigenMatrixUtility::GetInverseDiagonalMatrix(W).diagonal().sum() };
    if (!std::isfinite(W_inverse_trace) || W_inverse_trace <= 0.0)
    {
        return VectorXd::Ones(n).asDiagonal();
    }
    
    VectorXd capital_sigma{ VectorXd::Zero(n) };
    for (int j = 0; j < n; j++)
    {
        if (data_weight_array(j) == 0.0 || W_inverse_trace == 0.0) continue;
        capital_sigma(j) = n * sigma_square / data_weight_array(j) / W_inverse_trace;
        if (!std::isfinite(capital_sigma(j)) || capital_sigma(j) <= 0.0)
        {
            capital_sigma(j) = 1.0;
        }
    }
    return capital_sigma.asDiagonal();
}

Eigen::MatrixXd HRLModelHelper::CalculateMemberCovariance(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    auto basis_size{ static_cast<int>(beta_array.rows()) };
    auto member_size{ static_cast<int>(beta_array.cols()) };
    MatrixXd numerator{ MatrixXd::Zero(basis_size, basis_size) };
    double denominator{ omega_sum - member_size * alpha * std::pow(1.0+alpha, -1.0 - 0.5*basis_size) };
    if (denominator <= 0.0)
    {
        return MatrixXd::Identity(basis_size, basis_size);
    }
    MatrixXd residual_array{ beta_array.colwise() - mu };
    for (int i = 0; i < member_size; i++)
    {
        VectorXd residual{ residual_array.col(i) };
        numerator += omega_array(i) * (residual * residual.transpose());
    }

    MatrixXd capital_lambda{ numerator / denominator };
    if (!capital_lambda.array().allFinite())
    {
        capital_lambda = MatrixXd::Identity(basis_size, basis_size);
    }
    return capital_lambda;
}

std::vector<Eigen::MatrixXd> HRLModelHelper::CalculateWeightedMemberCovariance(
    const Eigen::MatrixXd & capital_lambda,
    const Eigen::ArrayXd & omega_array)
{
    auto member_size{ static_cast<int>(omega_array.size()) };
    auto omega_inverse_sum{ 0.0 };
    for (int i = 0; i < member_size; i++)
    {
        omega_inverse_sum += (omega_array(i) == 0.0) ? 0.0 : 1.0 / omega_array(i);
    }
    auto omega_h{ 1.0 / omega_inverse_sum };

    std::vector<MatrixXd> member_capital_lambda_list(static_cast<size_t>(member_size));
    for (int i = 0; i < member_size; i++)
    {
        member_capital_lambda_list.at(static_cast<size_t>(i)) =
            (member_size * omega_h / omega_array(i)) * capital_lambda;
    }
    return member_capital_lambda_list;
}

Eigen::VectorXd HRLModelHelper::CalculateBetaByOLS(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y)
{
    MatrixXd gram_matrix{ X.transpose() * X };
    MatrixXd inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
    return inverse_gram_matrix * (X.transpose() * y);
}

Eigen::VectorXd HRLModelHelper::CalculateMuByMedian(const Eigen::MatrixXd & beta_array)
{
    auto basis_size{ static_cast<int>(beta_array.rows()) };
    VectorXd mu{ VectorXd::Zero(basis_size) };
    for (int b = 0; b < basis_size; b++)
    {
        mu(b) = EigenMatrixUtility::GetMedian(beta_array.row(b));
    }
    return mu;
}

Eigen::VectorXd HRLModelHelper::CalculateBetaByMDPDE(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W)
{
    MatrixXd gram_matrix{ (X.transpose() * W * X) };
    auto inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
    return inverse_gram_matrix * (X.transpose() * W * y);
}

Eigen::VectorXd HRLModelHelper::CalculateMuByMDPDE(
    const Eigen::MatrixXd & beta_array,
    const Eigen::ArrayXd & omega_array,
    double omega_sum)
{
    MatrixXd numerator{ beta_array.array() / omega_sum };
    for (int i = 0; i < numerator.cols(); i++) numerator.col(i) *= omega_array(i);
    return numerator.rowwise().sum();
}

Eigen::ArrayXd HRLModelHelper::CalculateMemberWeight(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::MatrixXd & capital_lambda,
    double weight_min)
{
    auto member_size{ static_cast<int>(beta_array.cols()) };
    auto weight_member_min{ weight_min / static_cast<double>(member_size) };
    auto inverse_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(capital_lambda) };
    MatrixXd residual_array{ beta_array.colwise() - mu };
    Eigen::ArrayXd omega_array{ Eigen::ArrayXd::Zero(member_size) };
    for (int i = 0; i < member_size; i++)
    {
        VectorXd residual{ residual_array.col(i) };
        auto exp_index{ static_cast<double>(residual.transpose() * inverse_capital_lambda * residual) };
        if (!std::isfinite(exp_index))
        {
            omega_array(i) = weight_member_min;
        }
        else
        {
            omega_array(i) = std::exp(-0.5 * alpha * exp_index);
            if (omega_array(i) < weight_member_min) omega_array(i) = weight_member_min;
        }
    }
    return omega_array;
}

DMatrixXd HRLModelHelper::CalculateDataWeight(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const VectorXd & beta,
    double sigma_square,
    double weight_min)
{
    if (alpha == 0.0) return VectorXd::Ones(y.size()).asDiagonal();
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0) sigma_square = weight_min;

    const auto max_log{ std::log(std::numeric_limits<double>::max()) };
    ArrayXd exponent{ -0.5 * alpha * (y - (X * beta)).array().square() / sigma_square };
    exponent = exponent.cwiseMin(max_log);
    ArrayXd W{ exponent.exp() };
    W = W.unaryExpr([&](double w) { return (std::isfinite(w)) ? w : weight_min; });
    W = W.cwiseMax(weight_min);
    return W.matrix().asDiagonal();
}

Eigen::ArrayXd HRLModelHelper::CalculateMemberStatisticalDistance(
    const Eigen::VectorXd & mu_prior,
    const Eigen::MatrixXd & capital_lambda,
    const Eigen::MatrixXd & beta_posterior_array)
{
    auto member_size{ static_cast<int>(beta_posterior_array.cols()) };
    Eigen::ArrayXd statistical_distance_array{ Eigen::ArrayXd::Zero(member_size) };
    MatrixXd error_array{ beta_posterior_array.colwise() - mu_prior };
    auto inv_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(capital_lambda) };
    for (int i = 0; i < member_size; i++)
    {
        statistical_distance_array(i) =
            error_array.col(i).transpose() * inv_capital_lambda * error_array.col(i);
    }
    return statistical_distance_array;
}

Eigen::Array<bool, Eigen::Dynamic, 1> HRLModelHelper::CalculateOutlierMemberFlag(
    int basis_size,
    const Eigen::ArrayXd & statistical_distance_array)
{
    auto chi_square_quantile = [](int df)
    {
        switch (df)
        {
            case 1: return 6.63;  // 99th percentile
            case 2: return 9.21;
            case 3: return 11.34;
            case 4: return 13.28;
            case 5: return 15.09;
            case 6: return 16.81;
            case 7: return 18.48;
            case 8: return 20.09;
            case 9: return 21.67;
            case 10: return 23.21;
            default: // Wilson-Hilferty approximation
                const double z{ 2.3263478740408408 }; // standard normal 0.99 quantile
                const double d{ static_cast<double>(df) };
                const double term{ 1.0 - 2.0 / (9.0 * d) + z * std::sqrt(2.0 / (9.0 * d)) };
                return d * term * term * term;
        }
    };

    return (statistical_distance_array > chi_square_quantile(basis_size));
    //return ((m_omega_array/m_omega_sum) < 0.05/static_cast<double>(member_size));
}

void HRLModelHelper::ValidateBasisSize(int size) const
{
    if (size != m_basis_size)
    {
        throw std::invalid_argument(
            "The input size of data list isn't consistent with basis size : "
            + std::to_string(size) + " -> Expected: " + std::to_string(m_basis_size)
        );
    }
}

void HRLModelHelper::ValidateMemberSize(int size) const
{
    if (size != m_member_size)
    {
        throw std::invalid_argument(
            "The input size of data list isn't consistent with member size : "
            + std::to_string(size) + " -> Expected: " + std::to_string(m_member_size)
        );
    }
}

void HRLModelHelper::ValidateMemberId(int id) const
{
    if (id < 0 || id >= m_member_size)
    {
        throw std::out_of_range("HRLModelHelper::ValidateMemberId : member ID is out of range");
    }
}

void HRLModelHelper::ValidateAlpha(double alpha)
{
    if (!std::isfinite(alpha) || alpha < 0.0)
    {
        throw std::invalid_argument("Alpha parameters must be finite and non-negative.");
    }
}

void HRLModelHelper::ValidateMaximumIteration(int size)
{
    if (size <= 0)
    {
        throw std::invalid_argument("size must be greater than 0");
    }
}

void HRLModelHelper::ValidateTolerance(double value)
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

bool HRLModelHelper::GetOutlierFlag(int id) const
{
    ValidateMemberId(id);
    return m_outlier_flag_array(id); 
}

double HRLModelHelper::GetStatisticalDistance(int id) const
{
    ValidateMemberId(id);
    return m_statistical_distance_array(id);
}

Eigen::Ref<const Eigen::VectorXd> HRLModelHelper::GetBetaPosterior(int id) const
{
    ValidateMemberId(id);
    return Eigen::Ref<const Eigen::VectorXd>(m_beta_posterior_array.col(static_cast<Eigen::Index>(id)));
}

const Eigen::MatrixXd & HRLModelHelper::GetCapitalSigmaMatrixPosterior(int id) const
{
    ValidateMemberId(id);
    return m_capital_sigma_posterior_list[static_cast<size_t>(id)];
}
