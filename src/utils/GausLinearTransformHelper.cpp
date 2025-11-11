#include "GausLinearTransformHelper.hpp"
#include "Constants.hpp"
#include "Logger.hpp"

#include <stdexcept>

double GausLinearTransformHelper::GetGaussianResponseAtDistance(
    double r, double width, int dimension)
{
    if (width <= 0.0)
    {
        throw std::runtime_error("The gaus width should be positive value.");
    }

    double y{ 0.0 };
    double coeff{ 1.0 };
    double width_square{ width * width };
    coeff = 1.0 / std::pow(Constants::two_pi * width_square, 0.5 * dimension);
    y = coeff * std::exp(-0.5 * (r * r) / width_square);
    return y;
}

Eigen::VectorXd GausLinearTransformHelper::GetTylorSeriesBasisVector(
    double r, const Eigen::VectorXd & model_par)
{
    auto amplitude_0{ model_par(0) };
    auto width_0{ model_par(1) };
    auto G_0{ GetGaussianResponseAtDistance(r, width_0, 3) };
    Eigen::VectorXd basis_vector{ Eigen::VectorXd::Zero(3) };
    basis_vector(0) = G_0;
    basis_vector(1) = amplitude_0 * G_0 * (-3.0/width_0 + r*r/std::pow(width_0, 3));
    basis_vector(2) = 1.0;
    return basis_vector;
}

Eigen::VectorXd GausLinearTransformHelper::BuildLinearModelDataVector(
    double x, double y, const Eigen::VectorXd & model_par, int basis_dimension)
{
    auto data_vector_dimension{ basis_dimension + 1 };
    Eigen::VectorXd linear_model_data_vector{ Eigen::VectorXd::Zero(data_vector_dimension) };

    if (basis_dimension == 2)
    {
        if (y <= 0.0)
        {
            throw std::runtime_error("The gaus y value should be positive value.");
        }
        linear_model_data_vector(0) = 1.0;
        linear_model_data_vector(1) = -0.5 * x * x;
        linear_model_data_vector(2) = std::log(y);
        return linear_model_data_vector;
    }

    auto basis_vector{ GetTylorSeriesBasisVector(x, model_par) };
    auto amplitude_0{ model_par(0) };
    auto width_0{ model_par(1) };
    for (int i = 0; i < basis_dimension; i++)
    {
        linear_model_data_vector(i) = basis_vector(i);
    }
    auto intercept{ amplitude_0 * GetGaussianResponseAtDistance(x, width_0) + model_par(2) };
    linear_model_data_vector(basis_dimension) = y - intercept;
    return linear_model_data_vector;
}

Eigen::VectorXd GausLinearTransformHelper::BuildGaus2DModel(const Eigen::VectorXd & linear_model)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
    if (linear_model(1) <= 0.0) return gaus_model;

    auto model_dimension{ linear_model.rows() };
    if (model_dimension == 2)
    {
        gaus_model(0) = std::exp(linear_model(0)) * Constants::two_pi / linear_model(1);
        gaus_model(1) = 1.0 / std::sqrt(linear_model(1));
    }
    else
    {
        gaus_model(0) = linear_model(0);
        gaus_model(1) = linear_model(1);
    }

    return gaus_model;
}

Eigen::VectorXd GausLinearTransformHelper::BuildGaus3DModel(
    const Eigen::VectorXd & linear_model)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(3) };
    if (linear_model(1) <= 0.0) return gaus_model;

    auto model_dimension{ linear_model.rows() };
    if (model_dimension == 2)
    {
        gaus_model(0) = std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5);
        gaus_model(1) = 1.0 / std::sqrt(linear_model(1));
        gaus_model(2) = 0.0;
    }
    else
    {
        gaus_model(0) = linear_model(0);
        gaus_model(1) = linear_model(1);
        gaus_model(2) = linear_model(2);
    }

    return gaus_model;
}

Eigen::VectorXd GausLinearTransformHelper::BuildGaus3DModel(
    const Eigen::VectorXd & linear_model, const Eigen::VectorXd & model_par)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(3) };
    if (linear_model(1) <= 0.0) return gaus_model;

    auto model_dimension{ linear_model.rows() };
    if (model_dimension == 2)
    {
        gaus_model(0) = std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5);
        gaus_model(1) = 1.0 / std::sqrt(linear_model(1));
        gaus_model(2) = 0.0;
    }
    else
    {
        gaus_model(0) = std::exp(linear_model(0)) * model_par(0);
        gaus_model(1) = std::exp(linear_model(1)) + model_par(1);
        gaus_model(2) = linear_model(2) + model_par(2);
    }

    return gaus_model;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> GausLinearTransformHelper::BuildGaus2DModelWithVariance(
    const Eigen::VectorXd & linear_model,
    const Eigen::MatrixXd & covariance_matrix)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
    Eigen::VectorXd gaus_model_variance{ Eigen::VectorXd::Zero(2) };
    
    if (linear_model.rows() == 2)
    {
        if (covariance_matrix.rows() != 2 || covariance_matrix.cols() != 2)
        {
            throw std::invalid_argument("covariance_matrix must be 2x2");
        }
        if (!covariance_matrix.isApprox(covariance_matrix.transpose()))
        {
            throw std::invalid_argument("covariance_matrix must be symmetric");
        }

        gaus_model = BuildGaus2DModel(linear_model);
        auto beta0{ linear_model(0) };
        auto beta1{ linear_model(1) };
        if (beta1 <= 0.0)
        {
            return std::make_tuple(gaus_model, gaus_model_variance);
        }

        // TODO: modify with Gaus-2D error propagation
        auto var_beta0{ covariance_matrix(0, 0) };
        auto var_beta1{ covariance_matrix(1, 1) };
        auto cov{ covariance_matrix(0, 1) };
        auto var_amplitude{
            std::pow(Constants::two_pi, 3) * std::exp(2.0*beta0) * std::pow(beta1, -5) *
            (std::pow(beta1, 2) * var_beta0 + 2.25 * var_beta1 - 3.0 * beta1 * cov)
        };
        auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
        gaus_model_variance(0) = std::sqrt(var_amplitude);
        gaus_model_variance(1) = std::sqrt(var_width);
    }

    return std::make_tuple(gaus_model, gaus_model_variance);
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> GausLinearTransformHelper::BuildGaus3DModelWithVariance(
    const Eigen::VectorXd & linear_model,
    const Eigen::MatrixXd & covariance_matrix)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
    Eigen::VectorXd gaus_model_variance{ Eigen::VectorXd::Zero(2) };

    if (linear_model.rows() == 2)
    {
        if (covariance_matrix.rows() != 2 || covariance_matrix.cols() != 2)
        {
            throw std::invalid_argument("covariance_matrix must be 2x2");
        }
        if (!covariance_matrix.isApprox(covariance_matrix.transpose()))
        {
            throw std::invalid_argument("covariance_matrix must be symmetric");
        }

        gaus_model = BuildGaus3DModel(linear_model);
        auto beta0{ linear_model(0) };
        auto beta1{ linear_model(1) };
        if (beta1 <= 0.0)
        {
            return std::make_tuple(gaus_model, gaus_model_variance);
        }

        auto var_beta0{ covariance_matrix(0, 0) };
        auto var_beta1{ covariance_matrix(1, 1) };
        auto cov{ covariance_matrix(0, 1) };
        auto var_amplitude{
            std::pow(Constants::two_pi, 3) * std::exp(2.0*beta0) * std::pow(beta1, -5) *
            (std::pow(beta1, 2) * var_beta0 + 2.25 * var_beta1 - 3.0 * beta1 * cov)
        };
        auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
        gaus_model_variance(0) = std::sqrt(var_amplitude);
        gaus_model_variance(1) = std::sqrt(var_width);
    }

    return std::make_tuple(gaus_model, gaus_model_variance);
}
