#include "HRLModelHelper.hpp"
#include "EigenMatrixUtility.hpp"
#include "Logger.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <stdexcept>
#include <sstream>

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
    m_quiet_mode{ false },
    m_basis_size{ ValidatePositive(basis_size, "basis_size") },
    m_member_size{ ValidatePositive(member_size, "member_size") },
    m_maximum_iteration{ DEFAULT_MAXIMUM_ITERATION }, m_tolerance{ DEFAULT_TOLERANCE },
    m_omega_sum{ 0.0 },
    m_weight_data_min{ DEFAULT_WEIGHT_DATA_MIN },
    m_weight_member_min{ DEFAULT_WEIGHT_MEMBER_MIN / m_member_size },
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
    Finalization();
}

void HRLModelHelper::SetThreadSize(int thread_size)
{
    if (thread_size <= 0)
    {
        Logger::Log(LogLevel::Warning, "thread_size must be positive value, resetting to 1");
        thread_size = 1;
    }
    Eigen::setNbThreads(thread_size);
    Logger::Log(LogLevel::Debug, "Thread size = " + std::to_string(Eigen::nbThreads()));
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd> HRLModelHelper::BuildBasisVectorAndResponseArray(
    const std::vector<Eigen::VectorXd> & data_vector)
{
    auto data_size{ CheckedCastToInt(data_vector.size()) };
    if (data_size == 0 && m_quiet_mode == false)
    {
        Logger::Log(LogLevel::Warning,
            "HRLModelHelper::BuildBasisVectorAndResponseArray : dataset is empty.");
    }
    
    auto basis_size{ (data_vector[0].rows() - 1) };
    if (basis_size != m_basis_size)
    {
        Logger::Log(LogLevel::Warning,"basis size wrong : "+ std::to_string(basis_size));
    }
    MatrixXd x_data_matrix{ MatrixXd::Zero(data_size, basis_size) };
    VectorXd y_data_vector{ VectorXd::Zero(data_size) };
    for (int i = 0; i < data_size; i++)
    {
        const auto & sample{ data_vector.at(static_cast<size_t>(i)) };
        if (sample.size() != basis_size + 1)
        {
            throw std::invalid_argument("The input data size isn't consistent with basis size.");
        }
        if (!sample.array().allFinite())
        {
            throw std::invalid_argument("Member dataset contains non-finite value.");
        }
        x_data_matrix.row(i) = sample.head(basis_size);
        y_data_vector(i) = sample(basis_size);
    }

    return std::make_tuple(std::move(x_data_matrix), std::move(y_data_vector));
}

Eigen::MatrixXd HRLModelHelper::BuildBetaArray(const std::vector<Eigen::VectorXd> & beta_vector)
{
    auto member_size{ CheckedCastToInt(beta_vector.size()) };
    if (member_size == 0 && m_quiet_mode == false)
    {
        Logger::Log(LogLevel::Warning,
            "HRLModelHelper::BuildBetaArray : dataset is empty.");
    }

    MatrixXd beta_array{ MatrixXd::Zero(m_basis_size, member_size) };
    for (int i = 0; i < member_size; i++)
    {
        const auto & beta{ beta_vector.at(static_cast<size_t>(i)) };
        if (beta.size() != m_basis_size)
        {
            throw std::invalid_argument("The input data size isn't consistent with basis size.");
        }
        beta_array.col(i) = beta;
    }
    return beta_array;
}

void HRLModelHelper::SetMemberDataEntriesList(const std::vector<std::vector<Eigen::VectorXd>> & data_list)
{
    auto member_size{ CheckedCastToInt(data_list.size()) };
    if (member_size != m_member_size)
    {
        throw std::invalid_argument("The input size of data list isn't consistent with member size.");
    }
    m_X_list.clear();
    m_y_list.clear();
    for (auto & data_entry : data_list)
    {
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
    if (member_size != m_member_size)
    {
        throw std::invalid_argument("The input size of data list isn't consistent with member size.");
    }
    m_beta_MDPDE_array.setZero();
    m_sigma_square_array.setZero();
    m_W_list.clear();
    m_capital_sigma_list.clear();
    for (int i = 0; i < member_size; i++)
    {
        m_beta_MDPDE_array.col(i) = beta_list.at(static_cast<size_t>(i));
        m_sigma_square_array(i) = sigma_square_list.at(static_cast<size_t>(i));
        m_W_list.emplace_back(W_list.at(static_cast<size_t>(i)));
        m_capital_sigma_list.emplace_back(capital_sigma_list.at(static_cast<size_t>(i)));
    }
}

void HRLModelHelper::SetMaximumIteration(int size)
{
    if (size <= 0)
    {
        throw std::invalid_argument("size must be greater than 0");
    }
    m_maximum_iteration = size;
}

void HRLModelHelper::SetTolerance(double value)
{
    if (!std::isfinite(value))
    {
        throw std::invalid_argument("tolerance must be finite");
    }
    if (value < 0.0)
    {
        throw std::invalid_argument("tolerance must be non-negative");
    }
    m_tolerance = value;
}

void HRLModelHelper::RunGroupEstimation(double alpha_g)
{
    ValidateAlpha(alpha_g);
    RunAlgorithmMuMDPDE(alpha_g);
    RunAlgorithmWEB();
    m_statistical_distance_array =
        CalculateMemberStatisticalDistance(m_mu_prior, m_capital_lambda, m_beta_posterior_array);
    m_outlier_flag_array = CalculateOutlierMemberFlag(m_basis_size, m_statistical_distance_array);
}

std::tuple<
    Eigen::VectorXd,
    double,
    Eigen::DiagonalMatrix<double, Eigen::Dynamic>,
    Eigen::DiagonalMatrix<double, Eigen::Dynamic>
> HRLModelHelper::RunBetaMDPDE(
    const std::vector<Eigen::VectorXd> & data_vector,
    double alpha_r,
    Eigen::VectorXd & beta_OLS)
{
    ValidateAlpha(alpha_r);
    auto data{ BuildBasisVectorAndResponseArray(data_vector) };
    const auto & X{ std::get<0>(data) };
    const auto & y{ std::get<1>(data) };
    beta_OLS = CalculateBetaByOLS(X, y);
    Eigen::VectorXd beta{ beta_OLS };
    double sigma_square;
    DMatrixXd W;
    DMatrixXd capital_sigma;
    AlgorithmBetaMDPDE(alpha_r, X, y, beta, sigma_square, W, capital_sigma);
    return std::make_tuple(beta, sigma_square, W, capital_sigma);
}

Eigen::VectorXd HRLModelHelper::RunMuMDPDE(
    const std::vector<Eigen::VectorXd> & beta_vector, double alpha_g)
{
    ValidateAlpha(alpha_g);
    auto data{ BuildBetaArray(beta_vector) };
    VectorXd mu{ VectorXd::Zero(m_basis_size) };
    ArrayXd omega_array;
    double omega_sum;
    MatrixXd capital_lambda;
    std::vector<MatrixXd> member_capital_lambda_list;
    AlgorithmMuMDPDE(
        alpha_g, data, mu,
        omega_array, omega_sum, capital_lambda, member_capital_lambda_list);
    return mu;
}

void HRLModelHelper::RunAlgorithmWEB(void)
{
    if (m_member_size == 1)
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Debug,
                "HRLModelHelper::RunAlgorithmWEB : "
                "Only one member is present, using MDPDE estimate for WEB.");
        }
        m_beta_posterior_array = m_beta_MDPDE_array;
        m_capital_sigma_posterior_list[0] = Eigen::MatrixXd::Zero(m_basis_size, m_basis_size);
        m_mu_prior = m_mu_MDPDE;
        return;
    }

    Eigen::MatrixXd beta_posterior_array;
    std::vector<Eigen::MatrixXd> capital_sigma_posterior_list;
    Eigen::VectorXd mu_prior{ AlgorithmWEB(
        m_X_list, m_y_list, m_capital_sigma_list, m_mu_MDPDE, m_capital_lambda_list,
        beta_posterior_array, capital_sigma_posterior_list)
    };
    m_beta_posterior_array = beta_posterior_array;
    m_capital_sigma_posterior_list = std::move(capital_sigma_posterior_list);
    m_mu_prior = mu_prior;
}

void HRLModelHelper::AlgorithmBetaMDPDE(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    Eigen::VectorXd & beta,
    double & sigma_square,
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W,
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> & capital_sigma)
{
    int n{ static_cast<int>(y.size()) };
    sigma_square = (n < 2) ?
        std::numeric_limits<double>::max() : (y - (X * beta)).squaredNorm() / (n-1);
    W = VectorXd::Ones(n).asDiagonal();
    VectorXd beta_in_previous_iter{ beta };
    for (int t = 0; t < m_maximum_iteration; t++)
    { //=== Begin of iteration loop
        W = CalculateDataWeight(alpha, X, y, beta, sigma_square);
        beta = CalculateBetaByMDPDE(X, y, W);
        sigma_square = CalculateDataVarianceSquare(alpha, X, y, W, beta);
        if ((beta - beta_in_previous_iter).squaredNorm() < m_tolerance) break;
        beta_in_previous_iter = beta;
        if (t == m_maximum_iteration - 1 && m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::AlgorithmBetaMDPDE : "
                "Reach maximum iterations before achieving tolerance.");
        }
    } //=== End of iteration loop

    capital_sigma = CalculateDataCovariance(sigma_square, W);
}

void HRLModelHelper::RunAlgorithmBetaMDPDE(double alpha_r)
{
    for (int i = 0; i < m_member_size; i++)
    { //=== Begin of member ID loop (0 ... I-1)
        const auto & X{ m_X_list.at(static_cast<size_t>(i)) };
        const auto & y{ m_y_list.at(static_cast<size_t>(i)) };
        VectorXd beta_OLS{ CalculateBetaByOLS(X, y) };
        VectorXd beta{ beta_OLS };
        double sigma_square;
        DMatrixXd W;
        DMatrixXd capital_sigma;
        AlgorithmBetaMDPDE(alpha_r, X, y, beta, sigma_square, W, capital_sigma);
        m_mu_mean += beta_OLS;
        m_beta_MDPDE_array.col(i) = beta;
        m_W_list.at(static_cast<size_t>(i)) = W;
        m_capital_sigma_list.at(static_cast<size_t>(i)) = capital_sigma;
    } //=== End of member ID loop

    m_mu_mean /= m_member_size;
}

void HRLModelHelper::AlgorithmMuMDPDE(
    double alpha_g,
    const Eigen::MatrixXd & beta_array,
    Eigen::VectorXd & mu,
    Eigen::ArrayXd & omega_array,
    double & omega_sum,
    Eigen::MatrixXd & capital_lambda,
    std::vector<Eigen::MatrixXd> & member_capital_lambda_list)
{
    int basis_size{ static_cast<int>(beta_array.rows()) };
    int member_size{ static_cast<int>(beta_array.cols()) };
    omega_array = ArrayXd::Ones(member_size);
    capital_lambda = MatrixXd::Identity(basis_size, basis_size);
    member_capital_lambda_list.reserve(static_cast<size_t>(member_size));
    member_capital_lambda_list.assign(
        static_cast<size_t>(member_size), MatrixXd::Identity(basis_size, basis_size));

    if (member_size == 1)
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Debug,
                "HRLModelHelper::AlgorithmMuMDPDE : "
                "Only one member is present, using beta MDPDE for Mu.");
        }
        mu = beta_array.col(0);
        return;
    }

    for (int b = 0; b < basis_size; b++)
    {
        mu(b) = EigenMatrixUtility::GetMedian(beta_array.row(b));
    }

    VectorXd mu_in_previous_iter;
    for (int t = 0; t < m_maximum_iteration; t++)
    { //=== Begin of iteration loop
        mu_in_previous_iter = mu;
        omega_array = CalculateMemberWeight(alpha_g, beta_array, mu, capital_lambda);
        omega_sum = omega_array.sum();
        mu = CalculateMuByMDPDE(beta_array, omega_array, omega_sum);
        capital_lambda = CalculateMemberCovariance(alpha_g, beta_array, mu, omega_array, omega_sum);
        if ((mu - mu_in_previous_iter).squaredNorm() < m_tolerance) break;
        if (t == m_maximum_iteration - 1 && m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "Reach maximum iterations (Mu MDPDE) before achieving tolerance.");
        }
    } //=== End of iteration loop

    auto omega_h{ CalculateInverseMemberWeightSum(omega_array) };
    for (int i = 0; i < member_size; i++)
    {
        member_capital_lambda_list.at(static_cast<size_t>(i)) =
            (member_size * omega_h / omega_array(i)) * capital_lambda;
    }
}

void HRLModelHelper::RunAlgorithmMuMDPDE(double alpha_g)
{
    VectorXd mu{ VectorXd::Zero(m_basis_size) };
    ArrayXd omega_array;
    double omega_sum;
    MatrixXd capital_lambda;
    std::vector<MatrixXd> member_capital_lambda_list;
    AlgorithmMuMDPDE(
        alpha_g, m_beta_MDPDE_array, mu,
        omega_array, omega_sum, capital_lambda, member_capital_lambda_list);
    m_mu_MDPDE = mu;
    m_omega_array = omega_array;
    m_omega_sum = omega_sum;
    m_capital_lambda = capital_lambda;
    m_capital_lambda_list = std::move(member_capital_lambda_list);
}

Eigen::VectorXd HRLModelHelper::AlgorithmWEB(
    const std::vector<Eigen::MatrixXd> & X_list,
    const std::vector<Eigen::VectorXd> & y_list,
    const std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> & capital_sigma_list,
    const Eigen::VectorXd & mu_MDPDE,
    const std::vector<Eigen::MatrixXd> & member_capital_lambda_list,
    Eigen::MatrixXd & beta_posterior_array,
    std::vector<Eigen::MatrixXd> & capital_sigma_posterior_list)
{
    auto basis_size{ static_cast<int>(mu_MDPDE.rows()) };
    auto member_size{ static_cast<int>(member_capital_lambda_list.size()) };
    beta_posterior_array = MatrixXd::Zero(basis_size, member_size);

    capital_sigma_posterior_list.clear();
    capital_sigma_posterior_list.reserve(static_cast<size_t>(member_size));

    VectorXd mu_prior{ VectorXd::Zero(basis_size) };
    VectorXd numerator{ VectorXd::Zero(basis_size) };
    MatrixXd denominator{ MatrixXd::Zero(basis_size, basis_size) };
    for (size_t i = 0; i < static_cast<size_t>(member_size); i++)
    {
        const auto & X{ X_list.at(i) };
        const auto & y{ y_list.at(i) };
        const auto & capital_sigma{ capital_sigma_list.at(i) };
        const auto & member_capital_lambda{ member_capital_lambda_list.at(i) };
        auto inv_capital_sigma{ EigenMatrixUtility::GetInverseDiagonalMatrix(capital_sigma) };
        MatrixXd inv_member_capital_lambda;
        if (member_size == 1)
        {
            inv_member_capital_lambda = MatrixXd::Zero(basis_size, basis_size);
        }
        else
        {
            inv_member_capital_lambda = EigenMatrixUtility::GetInverseMatrix(member_capital_lambda);
        }
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
    if (m_member_size == 2) mu_prior = mu_MDPDE;

    return mu_prior;
}

double HRLModelHelper::CalculateDataVarianceSquare(
    double alpha,
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y,
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W,
    const Eigen::VectorXd & beta)
{
    double alpha_pow{ std::pow(1.0 + alpha, -1.5) };
    if (!std::isfinite(alpha_pow))
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning, "non-finite alpha power term; using 0.0");
        }
        alpha_pow = 0.0;
    }
    auto n{ static_cast<double>(y.size()) };
    VectorXd residual{ y - (X * beta) };
    auto numerator{ static_cast<double>(residual.transpose() * W * residual) };
    auto denominator{ W.diagonal().sum() - n * alpha * alpha_pow };
    if (denominator <= 0.0)
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(
                LogLevel::Warning,
                "Non-positive denominator in CalculateDataVarianceSquare,"
                " using small positive value");
        }
        denominator = 1.0e-10; // avoid division by zero
    }
    auto sigma_square{ numerator / denominator };
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0)
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(
                LogLevel::Warning,
                "Non-finite variance in CalculateDataVarianceSquare,"
                " using fallback value");
        }
        sigma_square = (std::isfinite(m_weight_data_min) && m_weight_data_min > 0.0)
            ? m_weight_data_min
            : std::numeric_limits<double>::max();
    }
    return sigma_square;
}

Eigen::DiagonalMatrix<double, Eigen::Dynamic> HRLModelHelper::CalculateDataCovariance(
    double sigma_square,
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & W)
{
    VectorXd data_weight_array{ W.diagonal() };
    const auto data_size{ static_cast<int>(data_weight_array.size()) };
    const auto W_inverse_trace{ EigenMatrixUtility::GetInverseDiagonalMatrix(W).diagonal().sum() };
    if (!std::isfinite(W_inverse_trace) || W_inverse_trace <= 0.0)
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::CalculateDataCovariance : "
                "degenerate weights; using fallback covariance");
        }
        return VectorXd::Ones(data_size).asDiagonal();
    }
    
    VectorXd capital_sigma{ VectorXd::Zero(data_size) };
    for (int j = 0; j < data_size; j++)
    {
        if (data_weight_array(j) == 0.0 || W_inverse_trace == 0.0) continue;
        capital_sigma(j) = data_size * sigma_square / data_weight_array(j) / W_inverse_trace;
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
    double alpha_pow{ std::pow(1.0 + alpha, -1.0 - 0.5 * basis_size) };
    if (!std::isfinite(alpha_pow))
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::CalculateMemberCovariance : "
                "non-finite alpha_g power term; using 0.0");
        }
        alpha_pow = 0.0;
    }

    MatrixXd numerator{ MatrixXd::Zero(basis_size, basis_size) };
    double denominator{ omega_sum - member_size * alpha * alpha_pow };
    if (denominator <= 0.0)
    {
        if (m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::CalculateMemberCovariance : "
                "Member covariance denominator is non-positive; using previous value.");
        }
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
        if (m_quiet_mode == false)
        {
            Logger::Log(LogLevel::Warning,
                "HRLModelHelper::CalculateMemberCovariance : "
                "Resulting covariance has non-finite entries; resetting to identity.");
        }
        capital_lambda = MatrixXd::Identity(basis_size, basis_size);
    }
    return capital_lambda;
}

double HRLModelHelper::CalculateInverseMemberWeightSum(const Eigen::ArrayXd & omega_array)
{
    auto member_size{ static_cast<int>(omega_array.size()) };
    auto omega_inverse_sum{ 0.0 };
    for (int i = 0; i < member_size; i++)
    {
        omega_inverse_sum += (omega_array(i) == 0.0) ? 0.0 : 1.0 / omega_array(i);
    }
    return 1.0 / omega_inverse_sum;
}

Eigen::VectorXd HRLModelHelper::CalculateBetaByOLS(
    const Eigen::MatrixXd & X,
    const Eigen::VectorXd & y)
{
    MatrixXd gram_matrix{ X.transpose() * X };
    MatrixXd inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
    return inverse_gram_matrix * (X.transpose() * y);
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
    // Method 1
    MatrixXd numerator{ beta_array.array() / omega_sum };
    for (int i = 0; i < numerator.cols(); i++) numerator.col(i) *= omega_array(i);
    return numerator.rowwise().sum();

    // Method 2 (To be checked the consistency with Method 1)
    //return (beta_array * omega_array.matrix()) / omega_sum;
}

Eigen::ArrayXd HRLModelHelper::CalculateMemberWeight(
    double alpha,
    const Eigen::MatrixXd & beta_array,
    const Eigen::VectorXd & mu,
    const Eigen::MatrixXd & capital_lambda)
{
    auto member_size{ static_cast<int>(beta_array.cols()) };
    auto inverse_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(capital_lambda) };
    MatrixXd residual_array{ beta_array.colwise() - mu };
    Eigen::ArrayXd omega_array{ Eigen::ArrayXd::Zero(member_size) };
    for (int i = 0; i < member_size; i++)
    {
        VectorXd residual{ residual_array.col(i) };
        auto exp_index{ static_cast<double>(residual.transpose() * inverse_capital_lambda * residual) };
        if (!std::isfinite(exp_index))
        {
            if (m_quiet_mode == false)
            {
                Logger::Log(LogLevel::Warning,
                    "HRLModelHelper::CalculateMemberWeight : non-finite exponent index");
            }
            omega_array(i) = m_weight_member_min;
        }
        else
        {
            omega_array(i) = std::exp(-0.5 * alpha * exp_index);
            if (omega_array(i) < m_weight_member_min) omega_array(i) = m_weight_member_min;
        }
    }
    return omega_array;
}

DMatrixXd HRLModelHelper::CalculateDataWeight(
    double alpha, const Eigen::MatrixXd & X, const Eigen::VectorXd & y,
    const VectorXd & beta, double sigma_square)
{
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0)
    {
        if (!std::isfinite(m_weight_data_min) || m_weight_data_min <= 0.0)
        {
            throw std::runtime_error(
                "HRLModelHelper::CalculateDataWeight : invalid sigma square and weight min");
        }
        sigma_square = m_weight_data_min;
    }

    if (alpha == 0.0)
    {
        const auto data_size{ static_cast<int>(y.size()) };
        return VectorXd::Ones(data_size).asDiagonal();
    }

    const double max_log{ std::log(std::numeric_limits<double>::max()) };
    const double fallback_weight{ m_weight_data_min };
    ArrayXd exponent{ -0.5 * alpha * (y - (X * beta)).array().square() / sigma_square };
    exponent = exponent.cwiseMin(max_log);
    ArrayXd W{ exponent.exp() };
    W = W.unaryExpr([&](double w) { return (std::isfinite(w)) ? w : fallback_weight; });
    W = W.cwiseMax(m_weight_data_min);
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

void HRLModelHelper::Finalization(void)
{
    m_X_list.clear();
    m_y_list.clear();
}

void HRLModelHelper::ValidateMemberId(int id) const
{
    if (id < 0 || id >= m_member_size)
    {
        throw std::out_of_range("HRLModelHelper::ValidateMemberId : member ID is out of range");
    }
}

void HRLModelHelper::ValidateAlpha(double alpha) const
{
    if (!std::isfinite(alpha) || alpha < 0.0)
    {
        throw std::invalid_argument("Alpha parameters must be finite and non-negative.");
    }

    if (alpha > ALPHA_LIMIT)
    {
        throw std::overflow_error("Alpha parameters are too large and may overflow calculations.");
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

double HRLModelHelper::GetMemberWeight(int id) const
{
    ValidateMemberId(id);
    return m_omega_array(id);
}

const Eigen::DiagonalMatrix<double, Eigen::Dynamic> &
HRLModelHelper::GetDataWeightMatrix(int id) const
{
    ValidateMemberId(id);
    return m_W_list[static_cast<size_t>(id)];
}

const Eigen::DiagonalMatrix<double, Eigen::Dynamic> &
HRLModelHelper::GetDataCovarianceMatrix(int id) const
{
    ValidateMemberId(id);
    return m_capital_sigma_list[static_cast<size_t>(id)];
}

const Eigen::MatrixXd & HRLModelHelper::GetCapitalSigmaMatrixPosterior(int id) const
{
    ValidateMemberId(id);
    return m_capital_sigma_posterior_list[static_cast<size_t>(id)];
}

Eigen::Ref<const Eigen::VectorXd> HRLModelHelper::GetBetaMatrixPosterior(int id) const
{
    ValidateMemberId(id);
    return Eigen::Ref<const Eigen::VectorXd>(m_beta_posterior_array.col(static_cast<Eigen::Index>(id)));
}

Eigen::Ref<const Eigen::VectorXd> HRLModelHelper::GetBetaMatrixMDPDE(int id) const
{
    ValidateMemberId(id);
    return Eigen::Ref<const Eigen::VectorXd>(m_beta_MDPDE_array.col(static_cast<Eigen::Index>(id)));
}
