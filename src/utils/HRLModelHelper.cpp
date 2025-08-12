#include "HRLModelHelper.hpp"
#include "EigenMatrixUtility.hpp"
#include "Logger.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

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
} // namespace

HRLModelHelper::HRLModelHelper(int basis_size, int member_size) :
    m_basis_size{ ValidatePositive(basis_size, "basis_size") },
    m_member_size{ ValidatePositive(member_size, "member_size") },
    m_maximum_iteration{ DEFAULT_MAXIMUM_ITERATION }, m_tolerance{ DEFAULT_TOLERANCE },
    m_omega_sum{ 0.0 }, m_omega_h{ 0.0 },
    m_weight_data_min{ DEFAULT_WEIGHT_DATA_MIN },
    m_weight_member_min{ DEFAULT_WEIGHT_MEMBER_MIN / m_member_size },
    m_sigma_square_array{ ArrayXd::Ones(m_member_size) },
    m_omega_array{ ArrayXd::Ones(m_member_size) },
    m_statistical_distance_array{ ArrayXd::Zero(m_member_size) },
    m_outlier_flag_array{ ArrayXb::Constant(m_member_size, false) },
    m_capital_lambda{ MatrixXd::Identity(m_basis_size, m_basis_size) },
    m_mu_iter{ VectorXd::Zero(m_basis_size) },
    m_mu_MDPDE{ VectorXd::Zero(m_basis_size) },
    m_mu_prior{ VectorXd::Zero(m_basis_size) },
    m_mu_mean{ VectorXd::Zero(m_basis_size) },
    m_beta_iter_array{ MatrixXd::Zero(m_basis_size, m_member_size) },
    m_beta_OLS_array{ MatrixXd::Zero(m_basis_size, m_member_size) },
    m_beta_MDPDE_array{ MatrixXd::Zero(m_basis_size, m_member_size) },
    m_beta_posterior_array{ MatrixXd::Zero(m_basis_size, m_member_size) }
{
    m_data_size_list.reserve(static_cast<size_t>(m_member_size));
    m_member_info_list.reserve(static_cast<size_t>(m_member_size));
    m_X_list.reserve(static_cast<size_t>(m_member_size));
    m_y_list.reserve(static_cast<size_t>(m_member_size));
    m_W_list.reserve(static_cast<size_t>(m_member_size));
    m_capital_sigma_list.reserve(static_cast<size_t>(m_member_size));
    m_capital_sigma_posterior_list.reserve(static_cast<size_t>(m_member_size));
    m_capital_lambda_list.reserve(static_cast<size_t>(m_member_size));

    for (int i = 0; i < m_member_size; i++)
    {
        m_capital_sigma_posterior_list.emplace_back(MatrixXd::Identity(m_basis_size, m_basis_size));
        m_capital_lambda_list.emplace_back(MatrixXd::Identity(m_basis_size, m_basis_size));
    }
}

void HRLModelHelper::SetDataArray(
    const std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> & data_array)
{
    if (data_array.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("data_array size exceeds maximum int");
    }
    int data_array_size{ static_cast<int>(data_array.size()) };
    if (data_array_size != m_member_size)
    {
        throw std::invalid_argument("The input size of data list isn't consistent with member size.");
    }

    // Build new data structures locally so existing member data remains
    // untouched if validation fails midway through the process.
    std::vector<int> data_size_list;
    std::vector<std::string> member_info_list;
    std::vector<MatrixXd> X_list;
    std::vector<VectorXd> y_list;
    data_size_list.reserve(static_cast<size_t>(m_member_size));
    member_info_list.reserve(static_cast<size_t>(m_member_size));
    X_list.reserve(static_cast<size_t>(m_member_size));
    y_list.reserve(static_cast<size_t>(m_member_size));

    for (const auto & [member_data, member_info] : data_array)
    {
        if (member_data.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("member_data size exceeds maximum int");
        }
        auto data_size{ static_cast<int>(member_data.size()) };
        if (data_size == 0)
        {
            throw std::invalid_argument("Member dataset must not be empty.");
        }
        MatrixXd x_data_matrix{ MatrixXd::Zero(data_size, m_basis_size) };
        VectorXd y_data_vector{ VectorXd::Zero(data_size) };
        for (int i = 0; i < data_size; i++)
        {
            const auto & sample{ member_data.at(static_cast<size_t>(i)) };
            if (sample.size() != m_basis_size + 1)
            {
                throw std::invalid_argument("The input data size isn't consistent with basis size.");
            }
            if (!sample.array().allFinite())
            {
                throw std::invalid_argument("Member dataset contains non-finite value.");
            }
            x_data_matrix.row(i) = sample.head(m_basis_size);
            y_data_vector(i) = sample(m_basis_size);
        }
        data_size_list.emplace_back(data_size);
        X_list.emplace_back(std::move(x_data_matrix));
        y_list.emplace_back(std::move(y_data_vector));
        member_info_list.emplace_back(member_info);
    }

    // All validation and transformations succeeded, swap into member variables.
    std::swap(m_data_size_list, data_size_list);
    std::swap(m_X_list, X_list);
    std::swap(m_y_list, y_list);
    std::swap(m_member_info_list, member_info_list);
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

void HRLModelHelper::RunEstimation(double alpha_r, double alpha_g)
{
    if (!std::isfinite(alpha_r) || !std::isfinite(alpha_g) || alpha_r < 0.0 || alpha_g < 0.0)
    {
        throw std::invalid_argument("Alpha parameters must be finite and non-negative.");
    }

    const double alpha_limit{ std::sqrt(std::numeric_limits<double>::max()) };
    if (alpha_r > alpha_limit || alpha_g > alpha_limit)
    {
        throw std::overflow_error("Alpha parameters are too large and may overflow calculations.");
    }

    if (m_data_size_list.size() != static_cast<size_t>(m_member_size))
    {
        throw std::runtime_error("Data array not set");
    }

    Initialization();
    AlgorithmBetaMDPDE(alpha_r);
    AlgorithmMuMDPDE(alpha_g);
    AlgorithmWEB();
    CalculateStatisticalDistance();
    LabelOutlierMember();
}

void HRLModelHelper::Initialization(void)
{
    // Reset scalars and arrays to their default values
    m_capital_lambda.setIdentity();
    m_omega_sum = 0.0;
    m_omega_h = 0.0;
    m_sigma_square_array.setOnes();
    m_omega_array.setOnes();
    m_statistical_distance_array.setZero();
    m_outlier_flag_array.setConstant(false);

    // Reset member-wise parameter estimates
    m_beta_iter_array.setZero();
    m_beta_OLS_array.setZero();
    m_beta_MDPDE_array.setZero();
    m_beta_posterior_array.setZero();
    m_mu_iter.setZero();
    m_mu_MDPDE.setZero();
    m_mu_prior.setZero();
    m_mu_mean.setZero();

    // Reset the containers of matrix objects
    m_W_list.clear();
    m_capital_sigma_list.clear();
    m_capital_sigma_posterior_list.assign(static_cast<size_t>(m_member_size),
                                          MatrixXd::Identity(m_basis_size, m_basis_size));
    m_capital_lambda_list.assign(static_cast<size_t>(m_member_size),
                                 MatrixXd::Identity(m_basis_size, m_basis_size));

    for (int i = 0; i < m_member_size; i++)
    { //=== Begin of member ID loop (0 ... I-1)
        const auto & data_size{ m_data_size_list.at(static_cast<size_t>(i)) };
        DMatrixXd w{ data_size };
        w.setIdentity();
        m_W_list.emplace_back(std::move(w));
        DMatrixXd sigma{ data_size };
        sigma.setIdentity();
        m_capital_sigma_list.emplace_back(std::move(sigma));
    } //=== End of member ID loop
}

void HRLModelHelper::AlgorithmBetaMDPDE(double alpha_r)
{
    for (int i = 0; i < m_member_size; i++)
    { //=== Begin of member ID loop (0 ... I-1)
        const auto & X{ m_X_list.at(static_cast<size_t>(i)) };
        const auto & y{ m_y_list.at(static_cast<size_t>(i)) };
        const auto & n{ m_data_size_list.at(static_cast<size_t>(i)) };
        MatrixXd gram_matrix{ X.transpose() * X };
        auto inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
        m_beta_OLS_array.col(i) = inverse_gram_matrix * (X.transpose() * y);
        m_beta_iter_array.col(i) = m_beta_OLS_array.col(i);
        m_sigma_square_array(i) = (n < 2) ?
            std::numeric_limits<double>::max() :
            (y - (X * m_beta_iter_array.col(i))).squaredNorm() / (n-1);
        VectorXd beta_in_previous_iter;
        for (int iter = 0; iter < m_maximum_iteration; iter++)
        { //=== Begin of iteration loop
            beta_in_previous_iter = m_beta_iter_array.col(i);
            CalculateDataWeight(i, alpha_r);
            CalculateBetaByMDPDE(i);
            CalculateDataVarianceSquare(i, alpha_r);
            if ((m_beta_iter_array.col(i) - beta_in_previous_iter).squaredNorm() < m_tolerance) break;
            if (iter == m_maximum_iteration - 1)
            {
                Logger::Log(LogLevel::Warning,
                            "HRLModelHelper::AlgorithmBetaMDPDE : "
                            "Reach maximum iterations before achieving tolerance for member -> "
                            + m_member_info_list.at(static_cast<size_t>(i)));
            }
        } //=== End of iteration loop
        m_beta_MDPDE_array.col(i) = m_beta_iter_array.col(i);
        CalculateDataCovariance(i);
    } //=== End of member ID loop

    m_mu_mean = m_beta_OLS_array.rowwise().mean();
}

void HRLModelHelper::AlgorithmMuMDPDE(double alpha_g)
{
    if (m_member_size == 1)
    {
        Logger::Log(LogLevel::Debug,
                    "HRLModelHelper::AlgorithmMuMDPDE : "
                    "Only one member is present, using OLS estimate for Mu MDPDE.");
        m_mu_MDPDE = m_beta_MDPDE_array.col(0);
        m_capital_lambda = MatrixXd::Zero(m_basis_size, m_basis_size);
        m_capital_lambda_list.at(0) = m_capital_lambda;
        return;
    }

    for (int b = 0; b < m_basis_size; b++)
    {
        m_mu_iter(b) = EigenMatrixUtility::GetMedian(m_beta_MDPDE_array.row(b));
    }

    VectorXd mu_in_previous_iter;
    for (int iter = 0; iter < m_maximum_iteration; iter++)
    { //=== Begin of iteration loop
        mu_in_previous_iter = m_mu_iter;
        CalculateMemberWeight(alpha_g);
        CalculateMuByMDPDE();
        CalculateMemberCovariance(alpha_g);
        if ((m_mu_iter - mu_in_previous_iter).squaredNorm() < m_tolerance) break;
        if (iter == m_maximum_iteration - 1)
        {
            Logger::Log(LogLevel::Warning,
                        "Reach maximum iterations (Mu MDPDE) before achieving tolerance for member: "
                        + m_member_info_list.at(0));
        }
    } //=== End of iteration loop

    m_mu_MDPDE = m_mu_iter;
    for (int i = 0; i < m_member_size; i++)
    {
        m_capital_lambda_list.at(static_cast<size_t>(i)) =
            (m_member_size * m_omega_h / m_omega_array(i)) * m_capital_lambda;
    }

}

void HRLModelHelper::AlgorithmWEB(void)
{
    VectorXd numerator{ VectorXd::Zero(m_basis_size, 1) };
    MatrixXd denominator{ MatrixXd::Zero(m_basis_size, m_basis_size) };
    for (size_t i = 0; i < static_cast<size_t>(m_member_size); i++)
    {
        const auto & X{ m_X_list.at(i) };
        const auto & y{ m_y_list.at(i) };
        auto inv_capital_sigma{
            EigenMatrixUtility::GetInverseDiagonalMatrix(m_capital_sigma_list.at(i))
        };
        MatrixXd inv_capital_lambda;
        if (m_member_size == 1)
        {
            inv_capital_lambda = MatrixXd::Zero(m_basis_size, m_basis_size);
        }
        else
        {
            inv_capital_lambda = EigenMatrixUtility::GetInverseMatrix(m_capital_lambda_list.at(i));
        }
        MatrixXd gram_matrix{ X.transpose() * inv_capital_sigma * X };
        VectorXd moment_matrix{ X.transpose() * inv_capital_sigma * y };
        MatrixXd inv_capital_sigma_posterior{ gram_matrix + inv_capital_lambda };
        m_capital_sigma_posterior_list.at(i) =
            EigenMatrixUtility::GetInverseMatrix(inv_capital_sigma_posterior);
        const auto & capital_sigma_posterior{ m_capital_sigma_posterior_list.at(i) };
        m_beta_posterior_array.col(static_cast<int>(i)) =
            capital_sigma_posterior * (moment_matrix + inv_capital_lambda * m_mu_MDPDE);
        numerator += inv_capital_lambda * capital_sigma_posterior * moment_matrix;
        denominator += inv_capital_lambda * capital_sigma_posterior * gram_matrix;
    }
    
    if (m_member_size == 1)
    {
        Logger::Log(LogLevel::Debug,
                    "HRLModelHelper::AlgorithmWEB : "
                    "Only one member is present, using MDPDE estimate for Mu prior.");
        m_mu_prior = m_beta_MDPDE_array.col(0);
    }
    else
    {
        MatrixXd inv_denominator{ EigenMatrixUtility::GetInverseMatrix(denominator) };
        m_mu_prior = inv_denominator * numerator;
        if (m_member_size == 2) m_mu_prior = m_mu_MDPDE;
    }
}

void HRLModelHelper::CalculateDataVarianceSquare(int member_id, double alpha_r)
{
    const auto & X{ m_X_list.at(static_cast<size_t>(member_id)) };
    const auto & y{ m_y_list.at(static_cast<size_t>(member_id)) };
    const auto & W{ m_W_list.at(static_cast<size_t>(member_id)) };
    const auto & n{ m_data_size_list.at(static_cast<size_t>(member_id)) };
    VectorXd beta{ m_beta_iter_array.col(member_id) };
    VectorXd residual{ y - (X * beta) };
    auto numerator{ static_cast<double>(residual.transpose() * W * residual) };
    auto denominator{ W.diagonal().sum() - n * alpha_r * pow(1.0 + alpha_r, -1.5) };
    if (denominator <= 0.0)
    {
        Logger::Log(
            LogLevel::Warning,
            "Non-positive denominator in CalculateDataVarianceSquare for member -> "
            + m_member_info_list.at(static_cast<size_t>(member_id))
            + ", using small positive value");
        denominator = 1.0e-10; // avoid division by zero
    }
    auto sigma_square{ numerator / denominator };
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0)
    {
        Logger::Log(
            LogLevel::Warning,
            "Non-finite variance in CalculateDataVarianceSquare for member -> "
            + m_member_info_list.at(static_cast<size_t>(member_id))
            + ", using fallback value");
        sigma_square = (std::isfinite(m_weight_data_min) && m_weight_data_min > 0.0)
            ? m_weight_data_min
            : std::numeric_limits<double>::max();
    }
    m_sigma_square_array(member_id) = sigma_square;
}

void HRLModelHelper::CalculateDataCovariance(int member_id)
{
    const auto & W{ m_W_list.at(static_cast<size_t>(member_id)) };
    const auto & sigma_square{ m_sigma_square_array(member_id) };
    VectorXd data_weight_array{ W.diagonal() };
    const auto data_size{ m_data_size_list.at(static_cast<size_t>(member_id)) };
    const auto W_inverse_trace{ EigenMatrixUtility::GetInverseDiagonalMatrix(W).diagonal().sum() };
    if (!std::isfinite(W_inverse_trace) || W_inverse_trace <= 0.0)
    {
        Logger::Log(LogLevel::Warning,
            "HRLModelHelper::CalculateDataCovariance : "
            "degenerate weights; using fallback covariance");
        auto & capital_sigma_matrix{ m_capital_sigma_list.at(static_cast<size_t>(member_id)) };
        if (capital_sigma_matrix.diagonal().size() != data_size)
        {
            capital_sigma_matrix.setIdentity(data_size);
        }
        return;
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
    m_capital_sigma_list.at(static_cast<size_t>(member_id)) = capital_sigma.asDiagonal();
}

void HRLModelHelper::CalculateMemberCovariance(double alpha_g)
{
    MatrixXd numerator{ MatrixXd::Zero(m_basis_size, m_basis_size) };
    double denominator{
        m_omega_sum - m_member_size * alpha_g * pow(1.0 + alpha_g, -1.0 - 0.5 * m_basis_size)
    };
    if (denominator <= 0.0)
    {
        Logger::Log(LogLevel::Warning,
            "HRLModelHelper::CalculateMemberCovariance : "
            "Member covariance denominator is non-positive; using previous value.");
        return; // Leave m_capital_lambda unchanged
    }
    MatrixXd residual_array{ m_beta_MDPDE_array.colwise() - m_mu_iter };
    for (int i = 0; i < m_member_size; i++)
    {
        VectorXd residual{ residual_array.col(i) };
        numerator += m_omega_array(i) * (residual * residual.transpose());
    }
    m_capital_lambda = numerator / denominator;
    if (!m_capital_lambda.array().allFinite())
    {
        Logger::Log(LogLevel::Warning,
            "HRLModelHelper::CalculateMemberCovariance : "
            "Resulting covariance has non-finite entries; resetting to identity.");
        m_capital_lambda = MatrixXd::Identity(m_basis_size, m_basis_size);
    }
}

void HRLModelHelper::CalculateBetaByMDPDE(int member_id)
{
    const auto & X{ m_X_list.at(static_cast<size_t>(member_id)) };
    const auto & y{ m_y_list.at(static_cast<size_t>(member_id)) };
    const auto & W{ m_W_list.at(static_cast<size_t>(member_id)) };
    MatrixXd gram_matrix{ (X.transpose() * W * X) };
    auto inverse_gram_matrix{ EigenMatrixUtility::GetInverseMatrix(gram_matrix) };
    m_beta_iter_array.col(member_id) = inverse_gram_matrix * (X.transpose() * W * y);
}

void HRLModelHelper::CalculateMuByMDPDE(void)
{
    MatrixXd numerator{ m_beta_MDPDE_array.array() / m_omega_sum };
    for (int i = 0; i < m_member_size; i++) numerator.col(i) *= m_omega_array(i);
    m_mu_iter = numerator.rowwise().sum();
}

void HRLModelHelper::CalculateMemberWeight(double alpha_g)
{
    auto inverse_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(m_capital_lambda) };
    MatrixXd residual_array{ m_beta_MDPDE_array.colwise() - m_mu_iter };
    auto omega_inverse_sum{ 0.0 };
    for (int i = 0; i < m_member_size; i++)
    {
        VectorXd residual{ residual_array.col(i) };
        auto exp_index{ static_cast<double>(residual.transpose() * inverse_capital_lambda * residual) };
        if (!std::isfinite(exp_index))
        {
            Logger::Log(LogLevel::Warning,
                        "HRLModelHelper::CalculateMemberWeight : non-finite exponent index");
            m_omega_array(i) = m_weight_member_min;
        }
        else
        {
            m_omega_array(i) = std::exp(-0.5 * alpha_g * exp_index);
            if (m_omega_array(i) < m_weight_member_min) m_omega_array(i) = m_weight_member_min;
        }
        omega_inverse_sum += (m_omega_array(i) == 0.0) ? 0.0 : 1.0 / m_omega_array(i);
    }
    m_omega_sum = m_omega_array.sum();
    m_omega_h = 1.0 / omega_inverse_sum;
}

void HRLModelHelper::CalculateDataWeight(int member_id, double alpha_r)
{
    const auto & X{ m_X_list.at(static_cast<size_t>(member_id)) };
    const auto & y{ m_y_list.at(static_cast<size_t>(member_id)) };
    auto sigma_square{ m_sigma_square_array(member_id) };
    if (!std::isfinite(sigma_square) || sigma_square <= 0.0)
    {
        if (!std::isfinite(m_weight_data_min) || m_weight_data_min <= 0.0)
        {
            throw std::runtime_error(
                "HRLModelHelper::CalculateDataWeight : invalid sigma square and weight min");
        }
        sigma_square = m_weight_data_min;
    }

    VectorXd beta{ m_beta_iter_array.col(member_id) };
    if (alpha_r == 0.0)
    {
        const auto data_size{ m_data_size_list.at(static_cast<size_t>(member_id)) };
        m_W_list.at(static_cast<size_t>(member_id)) = VectorXd::Ones(data_size).asDiagonal();
        return;
    }

    const double max_log{ std::log(std::numeric_limits<double>::max()) };
    const double fallback_weight{ m_weight_data_min };
    ArrayXd exponent{ -0.5 * alpha_r * (y - (X * beta)).array().square() / sigma_square };
    exponent = exponent.cwiseMin(max_log);
    ArrayXd W{ exponent.exp() };
    W = W.unaryExpr([&](double w) { return (std::isfinite(w)) ? w : fallback_weight; });
    W = W.cwiseMax(m_weight_data_min);
    m_W_list.at(static_cast<size_t>(member_id)) = W.matrix().asDiagonal();
}

void HRLModelHelper::CalculateStatisticalDistance(void)
{
    MatrixXd error_array{ m_beta_posterior_array.colwise() - m_mu_prior };
    if (m_member_size == 1)
    {
        m_statistical_distance_array.setZero();
        return;
    }
    
    auto inv_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(m_capital_lambda) };
    for (int i = 0; i < m_member_size; i++)
    {
        m_statistical_distance_array(i) =
            error_array.col(i).transpose() * inv_capital_lambda * error_array.col(i);
    }
}

void HRLModelHelper::LabelOutlierMember(void)
{
    auto chi_square_quantile = [](int df) {
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
        default:
            // Wilson-Hilferty approximation
            const double z{ 2.3263478740408408 }; // standard normal 0.99 quantile
            const double d{ static_cast<double>(df) };
            const double term{ 1.0 - 2.0 / (9.0 * d) + z * std::sqrt(2.0 / (9.0 * d)) };
            return d * term * term * term;
        }
    };
    const auto quantile{ chi_square_quantile(m_basis_size) };
    m_outlier_flag_array = (m_statistical_distance_array > quantile);
    for (int i = 0; i < m_member_size; ++i)
    {
        if (m_data_size_list.at(static_cast<size_t>(i)) < m_basis_size)
        {
            m_outlier_flag_array(i) = true;
        }
    }

    //m_outlier_flag_array = ((m_omega_array/m_omega_sum) < 0.05/static_cast<double>(m_member_size));
}

void HRLModelHelper::ValidateMemberId(int id) const
{
    if (id < 0 || id >= m_member_size)
    {
        throw std::out_of_range("member id out of range");
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

double HRLModelHelper::GetSigmaSquare(int id) const
{
    ValidateMemberId(id);
    return m_sigma_square_array(id);
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
    return m_W_list.at(static_cast<size_t>(id));
}

const Eigen::DiagonalMatrix<double, Eigen::Dynamic> &
HRLModelHelper::GetDataCovarianceMatrix(int id) const
{
    ValidateMemberId(id);
    return m_capital_sigma_list.at(static_cast<size_t>(id));
}

const Eigen::MatrixXd & HRLModelHelper::GetCapitalSigmaMatrixPosterior(int id) const
{
    ValidateMemberId(id);
    return m_capital_sigma_posterior_list.at(static_cast<size_t>(id));
}

Eigen::Ref<const Eigen::VectorXd> HRLModelHelper::GetBetaMatrixPosterior(int id) const
{
    ValidateMemberId(id);
    return Eigen::Ref<const Eigen::VectorXd>(m_beta_posterior_array.col(id));
}

Eigen::Ref<const Eigen::VectorXd> HRLModelHelper::GetBetaMatrixMDPDE(int id) const
{
    ValidateMemberId(id);
    return Eigen::Ref<const Eigen::VectorXd>(m_beta_MDPDE_array.col(id));
}

Eigen::Ref<const Eigen::VectorXd> HRLModelHelper::GetBetaMatrixOLS(int id) const
{
    ValidateMemberId(id);
    return Eigen::Ref<const Eigen::VectorXd>(m_beta_OLS_array.col(id));
}
