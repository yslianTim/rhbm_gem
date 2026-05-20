#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::linearization_service
{

namespace {
constexpr int kLogQuadraticBasisSize{ 2 };
constexpr int kDifferentialBasisSize{ 3 };

SeriesPointList BuildLogQuadraticDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double range_min,
    double range_max)
{
    SeriesPointList basis_and_response_entry_list;
    basis_and_response_entry_list.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ sample.point.distance };
        const auto gaussian_response{ static_cast<double>(sample.response) };
        if (distance < static_cast<float>(range_min)) continue;
        if (distance > static_cast<float>(range_max)) continue;
        if (gaussian_response <= 0.0) continue;

        auto basis_vector{ std::vector<double>{ 1.0, -0.5 * distance * distance } };
        auto transformed_response{ std::log(gaussian_response) };
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{ std::move(basis_vector), transformed_response }
        );
    }

    if (basis_and_response_entry_list.empty())
    {
        Logger::Log(LogLevel::Warning,
            "linearization_service::BuildDatasetSeries : "
            "No valid gaus data entry in the specified range.");
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{ std::vector<double>(kLogQuadraticBasisSize, 0.0), 0.0 }
        );
    }
    basis_and_response_entry_list.shrink_to_fit();
    return basis_and_response_entry_list;
}

SeriesPointList BuildDifferentialDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double range_min,
    double range_max)
{
    std::map<double, std::pair<double, int>> response_by_distance;
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto response{ static_cast<double>(sample.response) };
        if (!numeric_validation::IsFinite(distance) || !numeric_validation::IsFinite(response))
        {
            continue;
        }
        if (distance < 0.0) continue;
        if (distance > range_max) continue;

        auto & sum_and_count{ response_by_distance[distance] };
        sum_and_count.first += response;
        sum_and_count.second++;
    }

    std::map<double, double> integral_by_distance;
    double previous_distance{ 0.0 };
    double previous_integrand{ 0.0 };
    double integral{ 0.0 };
    for (const auto & [distance, sum_and_count] : response_by_distance)
    {
        const auto mean_response{ sum_and_count.first / static_cast<double>(sum_and_count.second) };
        const auto integrand{ distance * mean_response };
        integral += 0.5 * (integrand + previous_integrand) * (distance - previous_distance);
        integral_by_distance.emplace(distance, integral);
        previous_distance = distance;
        previous_integrand = integrand;
    }

    SeriesPointList basis_and_response_entry_list;
    basis_and_response_entry_list.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto response{ static_cast<double>(sample.response) };
        if (!numeric_validation::IsFinite(distance) || !numeric_validation::IsFinite(response))
        {
            continue;
        }
        if (distance < range_min) continue;
        if (distance > range_max) continue;

        const auto integral_iter{ integral_by_distance.find(distance) };
        if (integral_iter == integral_by_distance.end()) continue;

        basis_and_response_entry_list.emplace_back(SeriesPoint{
            std::vector<double>{ 1.0, -integral_iter->second, distance * distance },
            response
        });
    }

    if (basis_and_response_entry_list.empty())
    {
        Logger::Log(LogLevel::Warning,
            "linearization_service::BuildDatasetSeries : "
            "No valid differential gaus data entry in the specified range.");
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{ std::vector<double>(kDifferentialBasisSize, 0.0), 0.0 }
        );
    }
    basis_and_response_entry_list.shrink_to_fit();
    return basis_and_response_entry_list;
}

RHBMParameterVector EncodeLogQuadraticGaussianToParameterVector(
    const GaussianModel3D & gaussian_model)
{
    RHBMParameterVector parameter_vector{ RHBMParameterVector::Zero(kLogQuadraticBasisSize) };
    const auto amplitude{ gaussian_model.GetAmplitude() };
    const auto width{ gaussian_model.GetWidth() };
    if (width == 0.0 || amplitude <= 0.0) return parameter_vector;

    const auto width_square{ width * width };
    parameter_vector(0) = std::log(amplitude) - 1.5 * std::log(Constants::two_pi * width_square);
    parameter_vector(1) = 1.0 / width_square;
    return parameter_vector;
}

RHBMParameterVector EncodeDifferentialGaussianToParameterVector(
    const GaussianModel3D & gaussian_model)
{
    RHBMParameterVector parameter_vector{ RHBMParameterVector::Zero(kDifferentialBasisSize) };
    const auto amplitude{ gaussian_model.GetAmplitude() };
    const auto width{ gaussian_model.GetWidth() };
    const auto intercept{ gaussian_model.GetIntercept() };
    if (width == 0.0) return parameter_vector;

    const auto width_square{ width * width };
    parameter_vector(0) =
        amplitude / std::pow(Constants::two_pi * width_square, 1.5) + intercept;
    parameter_vector(1) = 1.0 / width_square;
    parameter_vector(2) = intercept / (2.0 * width_square);
    return parameter_vector;
}

GaussianModel3D DecodeLogQuadraticParameterVector(const RHBMParameterVector & parameter_vector)
{
    eigen_validation::RequireVectorSize(parameter_vector, kLogQuadraticBasisSize, "parameter_vector");
    const auto beta1{ parameter_vector(1) };
    if (beta1 <= 0.0) return GaussianModel3D{ 0.0, 0.0, 0.0 };
    return GaussianModel3D{
        std::exp(parameter_vector(0)) * std::pow(Constants::two_pi / beta1, 1.5),
        1.0 / std::sqrt(beta1),
        0.0
    };
}

GaussianModel3D DecodeDifferentialParameterVector(const RHBMParameterVector & parameter_vector)
{
    eigen_validation::RequireVectorSize(parameter_vector, kDifferentialBasisSize, "parameter_vector");
    const auto beta0{ parameter_vector(0) };
    const auto beta1{ parameter_vector(1) };
    const auto beta2{ parameter_vector(2) };
    if (beta1 <= 0.0) return GaussianModel3D{ 0.0, 0.0, 0.0 };

    const auto amplitude{
        std::pow(Constants::two_pi, 1.5) *
        (beta0 / std::pow(beta1, 1.5) - 2.0 * beta2 / std::pow(beta1, 2.5))
    };
    return GaussianModel3D{
        amplitude,
        1.0 / std::sqrt(beta1),
        2.0 * beta2 / beta1
    };
}

double CalculateVariance(
    const Eigen::VectorXd & gradient,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    const auto variance{ static_cast<double>(gradient.transpose() * covariance_matrix * gradient) };
    return variance < 0.0 ? 0.0 : variance;
}

} // namespace

SeriesPointList BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double range_min,
    double range_max,
    LocalGaussianFitModel fit_model)
{
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "data range");

    switch (fit_model)
    {
    case LocalGaussianFitModel::LogQuadratic:
        return BuildLogQuadraticDatasetSeries(sampling_entries, range_min, range_max);
    case LocalGaussianFitModel::DifferentialMethod:
        return BuildDifferentialDatasetSeries(sampling_entries, range_min, range_max);
    }

    throw std::invalid_argument("Unsupported local Gaussian fit model for linearized dataset.");
}

RHBMParameterVector EncodeGaussianToParameterVector(const GaussianModel3D & gaussian_model)
{
    return EncodeGaussianToParameterVector(gaussian_model, LocalGaussianFitModel::LogQuadratic);
}

RHBMParameterVector EncodeGaussianToParameterVector(
    const GaussianModel3D & gaussian_model,
    LocalGaussianFitModel fit_model)
{
    switch (fit_model)
    {
    case LocalGaussianFitModel::LogQuadratic:
        return EncodeLogQuadraticGaussianToParameterVector(gaussian_model);
    case LocalGaussianFitModel::DifferentialMethod:
        return EncodeDifferentialGaussianToParameterVector(gaussian_model);
    }

    throw std::invalid_argument("Unsupported local Gaussian fit model for parameter encoding.");
}

GaussianModel3D DecodeParameterVector(const RHBMParameterVector & parameter_vector)
{
    return DecodeParameterVector(parameter_vector, LocalGaussianFitModel::LogQuadratic);
}

GaussianModel3D DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    LocalGaussianFitModel fit_model)
{
    switch (fit_model)
    {
    case LocalGaussianFitModel::LogQuadratic:
        return DecodeLogQuadraticParameterVector(parameter_vector);
    case LocalGaussianFitModel::DifferentialMethod:
        return DecodeDifferentialParameterVector(parameter_vector);
    }

    throw std::invalid_argument("Unsupported local Gaussian fit model for parameter decoding.");
}

GaussianModel3DWithUncertainty DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    return DecodeParameterVector(
        parameter_vector, covariance_matrix, LocalGaussianFitModel::LogQuadratic);
}

GaussianModel3DWithUncertainty DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix,
    LocalGaussianFitModel fit_model)
{
    eigen_validation::RequireShape(
        covariance_matrix,
        fit_model == LocalGaussianFitModel::DifferentialMethod ?
            kDifferentialBasisSize : kLogQuadraticBasisSize,
        fit_model == LocalGaussianFitModel::DifferentialMethod ?
            kDifferentialBasisSize : kLogQuadraticBasisSize,
        "covariance_matrix");
    const auto model{ DecodeParameterVector(parameter_vector, fit_model) };
    if (fit_model == LocalGaussianFitModel::DifferentialMethod)
    {
        const auto beta0{ parameter_vector(0) };
        const auto beta1{ parameter_vector(1) };
        const auto beta2{ parameter_vector(2) };
        if (beta1 <= 0.0)
        {
            return GaussianModel3DWithUncertainty{model, GaussianModel3DUncertainty{}};
        }

        const auto amplitude_factor{ std::pow(Constants::two_pi, 1.5) };
        Eigen::VectorXd amplitude_gradient{ Eigen::VectorXd::Zero(kDifferentialBasisSize) };
        amplitude_gradient(0) = amplitude_factor * std::pow(beta1, -1.5);
        amplitude_gradient(1) = amplitude_factor *
            (-1.5 * beta0 * std::pow(beta1, -2.5) +
                5.0 * beta2 * std::pow(beta1, -3.5));
        amplitude_gradient(2) = -2.0 * amplitude_factor * std::pow(beta1, -2.5);

        Eigen::VectorXd width_gradient{ Eigen::VectorXd::Zero(kDifferentialBasisSize) };
        width_gradient(1) = -0.5 * std::pow(beta1, -1.5);

        Eigen::VectorXd intercept_gradient{ Eigen::VectorXd::Zero(kDifferentialBasisSize) };
        intercept_gradient(1) = -2.0 * beta2 / (beta1 * beta1);
        intercept_gradient(2) = 2.0 / beta1;

        return GaussianModel3DWithUncertainty{
            model,
            GaussianModel3DUncertainty{
                std::sqrt(CalculateVariance(amplitude_gradient, covariance_matrix)),
                std::sqrt(CalculateVariance(width_gradient, covariance_matrix)),
                std::sqrt(CalculateVariance(intercept_gradient, covariance_matrix))
            }
        };
    }

    const auto beta1{ parameter_vector(1) };
    if (beta1 <= 0.0)
    {
        return GaussianModel3DWithUncertainty{model, GaussianModel3DUncertainty{}};
    }

    const auto amplitude{ model.GetAmplitude() };
    const auto var_beta0{ covariance_matrix(0, 0) };
    const auto var_beta1{ covariance_matrix(1, 1) };
    const auto cov{ covariance_matrix(0, 1) };
    const auto half_dimension{ 1.5 };
    const auto var_amplitude{
        amplitude * amplitude *
        (var_beta0 + half_dimension * half_dimension * var_beta1 / (beta1 * beta1) -
            3.0 * cov / beta1)
    };
    const auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
    return GaussianModel3DWithUncertainty{
        model,
        GaussianModel3DUncertainty{ std::sqrt(var_amplitude), std::sqrt(var_width), 0.0 }
    };
}

} // namespace rhbm_gem::linearization_service
