#include "HRLModelHelper.hpp"
#include "EigenMatrixUtility.hpp"
#include "Logger.hpp"

#include <cmath>
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

HRLModelHelper::HRLModelHelper(int basis_size, int member_size) :
    m_basis_size{ basis_size }, m_member_size{ member_size },
    m_maximum_iteration{ 100 }, m_tolerance{ 1.0e-5 },
    m_omega_sum{ 0.0 }, m_omega_h{ 0.0 },
    m_weight_data_min{ 1.0e-8 },
    m_weight_member_min{ (member_size <= 0) ? 0.0 : 1.0e-2 / member_size },
    m_sigma_square_array{ ArrayXd::Ones(member_size > 0 ? member_size : 1) },
    m_omega_array{ ArrayXd::Ones(member_size > 0 ? member_size : 1) },
    m_statistical_distance_array{ ArrayXd::Zero(member_size > 0 ? member_size : 1) },
    m_outlier_flag_array{ ArrayXb::Constant(member_size > 0 ? member_size : 1, false) },
    m_capital_lambda{ MatrixXd::Identity(basis_size > 0 ? basis_size : 1, basis_size > 0 ? basis_size : 1) },
    m_mu_iter{ VectorXd::Ones(basis_size > 0 ? basis_size : 1) },
    m_mu_MDPDE{ VectorXd::Ones(basis_size > 0 ? basis_size : 1) },
    m_mu_prior{ VectorXd::Ones(basis_size > 0 ? basis_size : 1) },
    m_mu_mean{ VectorXd::Ones(basis_size > 0 ? basis_size : 1) },
    m_beta_iter_array{ MatrixXd::Ones(basis_size > 0 ? basis_size : 1, member_size > 0 ? member_size : 1) },
    m_beta_OLS_array{ MatrixXd::Ones(basis_size > 0 ? basis_size : 1, member_size > 0 ? member_size : 1) },
    m_beta_MDPDE_array{ MatrixXd::Ones(basis_size > 0 ? basis_size : 1, member_size > 0 ? member_size : 1) },
    m_beta_posterior_array{ MatrixXd::Ones(basis_size > 0 ? basis_size : 1, member_size > 0 ? member_size : 1) }
{
    if (basis_size <= 0 || member_size <= 0)
    {
        throw std::invalid_argument("basis_size and member_size must be positive values");
    }

    m_data_size_list.reserve(static_cast<size_t>(m_member_size));
    m_member_info_list.reserve(static_cast<size_t>(m_member_size));
    m_X_list.reserve(static_cast<size_t>(m_member_size));
    m_y_list.reserve(static_cast<size_t>(m_member_size));
    m_W_list.reserve(static_cast<size_t>(m_member_size));
    m_capital_sigma_list.reserve(static_cast<size_t>(m_member_size));
    m_capital_sigma_posterior_list.reserve(static_cast<size_t>(m_member_size));
    m_capital_lambda_list.reserve(static_cast<size_t>(m_member_size));
}

void HRLModelHelper::SetDataArray(
    const std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> & data_array)
{
    int data_array_size{ static_cast<int>(data_array.size()) };
    if (data_array_size != m_member_size)
    {
        throw std::runtime_error("The input size of data list isn't consistent with member size.");
    }

    m_data_size_list.clear();
    m_member_info_list.clear();
    m_X_list.clear();
    m_y_list.clear();
    for (auto & [member_data, member_info] : data_array)
    {
        auto data_size{ static_cast<int>(member_data.size()) };
        MatrixXd x_data_matrix{ MatrixXd::Zero(data_size, m_basis_size) };
        VectorXd y_data_vector{ VectorXd::Zero(data_size) };
        for (int i = 0; i < data_size; i++)
        {
            if (member_data.at(static_cast<size_t>(i)).size() != m_basis_size + 1)
            {
                throw std::runtime_error("The input data size isn't consistent with basis size.");
            }
            for (int j = 0; j < m_basis_size; j++)
            {
                x_data_matrix(i, j) = member_data.at(static_cast<size_t>(i))(j);
            }
            y_data_vector(i) = member_data.at(static_cast<size_t>(i))(m_basis_size);
        }
        m_data_size_list.emplace_back(data_size);
        m_X_list.emplace_back(x_data_matrix);
        m_y_list.emplace_back(y_data_vector);
        m_member_info_list.emplace_back(member_info);
    }
}

void HRLModelHelper::RunEstimation(double alpha_r, double alpha_g)
{
    if (alpha_r < 0.0 || alpha_g < 0.0)
    {
        throw std::invalid_argument("Alpha parameters must be non-negative.");
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
    // Empty the containers of matrix objects
    m_W_list.clear();
    m_capital_sigma_list.clear();
    m_capital_sigma_posterior_list.clear();
    m_capital_lambda_list.clear();

    for (int i = 0; i < m_member_size; i++)
    { //=== Begin of member ID loop (0 ... I-1)
        const auto & data_size{ m_data_size_list.at(static_cast<size_t>(i)) };
        m_W_list.emplace_back(MatrixXd::Identity(data_size, data_size).diagonal());
        m_capital_sigma_list.emplace_back(MatrixXd::Identity(data_size, data_size).diagonal());
        m_capital_sigma_posterior_list.emplace_back(MatrixXd::Identity(m_basis_size, m_basis_size));
        m_capital_lambda_list.emplace_back(MatrixXd::Identity(m_basis_size, m_basis_size));
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
        m_sigma_square_array(i) = (n < 2) ? 1.0e+200 : (y - (X * m_beta_iter_array.col(i))).squaredNorm() / (n-1);
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
                            "Reach maximum iterations (Beta MDPDE) before achieving tolerance for member: "
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
        m_capital_lambda_list.at(static_cast<size_t>(i)) = m_member_size * m_omega_h / m_omega_array(i) * m_capital_lambda.array();
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
        auto inv_capital_sigma{ EigenMatrixUtility::GetInverseDiagonalMatrix(m_capital_sigma_list.at(i)) };
        auto inv_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(m_capital_lambda_list.at(i)) };
        if (m_member_size == 1) inv_capital_lambda = MatrixXd::Zero(m_basis_size, m_basis_size);
        MatrixXd gram_matrix{ X.transpose() * inv_capital_sigma * X };
        VectorXd moment_matrix{ X.transpose() * inv_capital_sigma * y };
        MatrixXd inv_capital_sigma_posterior{ gram_matrix + inv_capital_lambda };
        m_capital_sigma_posterior_list.at(i) = EigenMatrixUtility::GetInverseMatrix(inv_capital_sigma_posterior);
        const auto & capital_sigma_posterior{ m_capital_sigma_posterior_list.at(i) };
        m_beta_posterior_array.col(static_cast<int>(i)) = capital_sigma_posterior * (moment_matrix + inv_capital_lambda * m_mu_MDPDE);
        numerator += inv_capital_lambda * capital_sigma_posterior * moment_matrix;
        denominator += inv_capital_lambda * capital_sigma_posterior * gram_matrix;
    }
    MatrixXd inv_denominator{ EigenMatrixUtility::GetInverseMatrix(denominator) };
    m_mu_prior = inv_denominator * numerator;
    if (m_member_size == 1) m_mu_prior = m_beta_MDPDE_array.col(0);
    if (m_member_size == 2) m_mu_prior = m_mu_MDPDE;
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
    m_sigma_square_array(member_id) = (denominator == 0.0) ? 1.0e+200 : numerator / denominator;
}

void HRLModelHelper::CalculateDataCovariance(int member_id)
{
    const auto & W{ m_W_list.at(static_cast<size_t>(member_id)) };
    const auto & sigma_square{ m_sigma_square_array(member_id) };
    VectorXd data_weight_array{ W.diagonal() };
    const auto W_inverse_trace{ EigenMatrixUtility::GetInverseDiagonalMatrix(W).diagonal().sum() };
    const auto data_size{ m_data_size_list.at(static_cast<size_t>(member_id)) };
    VectorXd capital_sigma{ VectorXd::Zero(data_size) };
    for (int j = 0; j < data_size; j++)
    {
        if (data_weight_array(j) == 0.0 || W_inverse_trace == 0.0) continue;
        capital_sigma(j) = data_size * sigma_square / data_weight_array(j) / W_inverse_trace;
    }
    m_capital_sigma_list.at(static_cast<size_t>(member_id)) = capital_sigma.asDiagonal();
}

void HRLModelHelper::CalculateMemberCovariance(double alpha_g)
{
    MatrixXd numerator{ MatrixXd::Zero(m_basis_size, m_basis_size) };
    double denominator{ m_omega_sum - m_member_size * alpha_g * pow(1.0 + alpha_g, -1.0 - 0.5 * m_basis_size) };
    MatrixXd residual_array{ m_beta_MDPDE_array.colwise() - m_mu_iter };
    for (int i = 0; i < m_member_size; i++)
    {
        VectorXd residual{ residual_array.col(i) };
        numerator += m_omega_array(i) * (residual * residual.transpose());
    }
    m_capital_lambda = numerator/denominator;
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
        m_omega_array(i) = exp(-0.5 * alpha_g * exp_index);
        if (m_omega_array(i) < m_weight_member_min) m_omega_array(i) = m_weight_member_min;
        omega_inverse_sum += (m_omega_array(i) == 0.0) ? 0.0 : 1.0 / m_omega_array(i);
    }
    m_omega_sum = m_omega_array.sum();
    m_omega_h = 1.0 / omega_inverse_sum;
}

void HRLModelHelper::CalculateDataWeight(int member_id, double alpha_r)
{
    const auto & X{ m_X_list.at(static_cast<size_t>(member_id)) };
    const auto & y{ m_y_list.at(static_cast<size_t>(member_id)) };
    const auto & sigma_square{ m_sigma_square_array(member_id) };
    VectorXd beta{ m_beta_iter_array.col(member_id) };
    ArrayXd W{ (-0.5 * alpha_r * (y - (X * beta)).array().square() / sigma_square).exp() };
    W.max(m_weight_data_min);
    m_W_list.at(static_cast<size_t>(member_id)) = W.matrix().asDiagonal();
}

void HRLModelHelper::CalculateStatisticalDistance(void)
{
    MatrixXd error_array{ m_beta_posterior_array.colwise() - m_mu_prior };
    if (error_array.cols() <= 2)
    {
        for (int i = 0; i < m_member_size; i++) m_statistical_distance_array(i) = 1.0e-5;
    }
    else
    {
        auto inv_capital_lambda{ EigenMatrixUtility::GetInverseMatrix(m_capital_lambda) };
        for (int i = 0; i < m_member_size; i++)
        {
            m_statistical_distance_array(i) = error_array.col(i).transpose() * inv_capital_lambda * error_array.col(i);
        }
    }
}

void HRLModelHelper::LabelOutlierMember(void)
{
    auto quantile{ 9.21 };
    m_outlier_flag_array = (m_statistical_distance_array > quantile);

    //m_outlier_flag_array = ((m_omega_array/m_omega_sum) < 0.05/static_cast<double>(m_member_size));
}
