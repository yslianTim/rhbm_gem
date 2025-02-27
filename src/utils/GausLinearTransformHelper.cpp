#include "GausLinearTransformHelper.hpp"

#include <cmath>

Eigen::VectorXd GausLinearTransformHelper::BuildLinearModelDataVector(double gaus_x, double gaus_y)
{
    auto linear_model_data_vector_basis{ 2 + 1 };
    Eigen::VectorXd linear_model_data_vector{ Eigen::VectorXd::Zero(linear_model_data_vector_basis) };

    if (gaus_y <= 0.0)
    {
        throw std::runtime_error("The gaus y value should be positive value.");
    }

    linear_model_data_vector(0) = 1.0;
    linear_model_data_vector(1) = -0.5 * gaus_x * gaus_x;
    linear_model_data_vector(2) = std::log(gaus_y);
    
    return linear_model_data_vector;
}

Eigen::VectorXd GausLinearTransformHelper::BuildGausModel(const Eigen::VectorXd & linear_model)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
    if (linear_model(1) <= 0.0) return gaus_model;

    gaus_model(0) = std::exp(linear_model(0)) * std::pow(m_two_pi / linear_model(1), 1.5);
    gaus_model(1) = 1.0 / std::sqrt(linear_model(1));

    return gaus_model;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> GausLinearTransformHelper::BuildGausModelWithVariance(
    const Eigen::VectorXd & linear_model,
    const Eigen::MatrixXd & covariance_matrix)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
    Eigen::VectorXd gaus_model_variance{ Eigen::VectorXd::Zero(2) };

    gaus_model = BuildGausModel(linear_model);

    auto beta0{ linear_model(0) };
    auto beta1{ linear_model(1) };
    auto var_beta0{ covariance_matrix(0, 0) };
    auto var_beta1{ covariance_matrix(1, 1) };
    auto cov{ covariance_matrix(0, 1) };
    auto var_amplitude{
        std::pow(m_two_pi, 3) * std::exp(2.0*beta0) * std::pow(beta1, -5) *
        (std::pow(beta1, 2) * var_beta0 + 2.25 * var_beta1 - 3.0 * beta1 * cov)
    };
    auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
    gaus_model_variance(0) = std::sqrt(var_amplitude);
    gaus_model_variance(1) = std::sqrt(var_width);

    return std::make_tuple(gaus_model, gaus_model_variance);
}