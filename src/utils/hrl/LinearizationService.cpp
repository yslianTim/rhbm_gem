#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace rhbm_gem::linearization_service
{

namespace {
constexpr int kLogQuadraticBasisSize{ 2 };

std::tuple<std::vector<double>, double> BuildLogQuadraticBasisVector(double x, double y)
{
    if (y <= 0.0)
    {
        throw std::runtime_error("The gaus y value should be positive value.");
    }
    return std::make_tuple(
        std::vector<double>{ 1.0, -0.5 * x * x },
        std::log(y));
}

RHBMParameterVector BuildParameterVector(const GaussianModel3D & model)
{
    RHBMParameterVector parameter_vector{ RHBMParameterVector::Zero(kLogQuadraticBasisSize) };
    const auto width{ model.GetWidth() };
    const auto amplitude{ model.GetAmplitude() };
    if (width == 0.0)
    {
        return parameter_vector;
    }

    const auto width_square{ width * width };
    parameter_vector(0) = (amplitude <= 0.0) ?
        0.0 : std::log(amplitude) - 1.5 * std::log(Constants::two_pi * width_square);
    parameter_vector(1) = 1.0 / width_square;
    return parameter_vector;
}

double ResolveLogQuadraticDimension(const LinearizationSpec & spec)
{
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return 2.0;
    case GaussianModelKind::MODEL_3D:
        return 3.0;
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

GaussianModel3D DecodeLogQuadratic(
    const RHBMParameterVector & linear_model,
    double dimension)
{
    const auto beta1{ linear_model(1) };
    if (beta1 <= 0.0)
    {
        return GaussianModel3D{ 0.0, 0.0, 0.0 };
    }
    return GaussianModel3D{
        std::exp(linear_model(0)) * std::pow(Constants::two_pi / beta1, 0.5 * dimension),
        1.0 / std::sqrt(beta1),
        0.0
    };
}

GaussianModel3DUncertainty CalculateLogQuadraticStandardDeviation(
    const RHBMParameterVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix,
    double dimension,
    const GaussianModel3D & model)
{
    const auto beta1{ linear_model(1) };
    if (beta1 <= 0.0)
    {
        return GaussianModel3DUncertainty{};
    }

    const auto amplitude{ model.GetAmplitude() };
    const auto var_beta0{ covariance_matrix(0, 0) };
    const auto var_beta1{ covariance_matrix(1, 1) };
    const auto cov{ covariance_matrix(0, 1) };
    const auto half_dimension{ 0.5 * dimension };
    const auto var_amplitude{
        amplitude * amplitude *
        (var_beta0 + half_dimension * half_dimension * var_beta1 / (beta1 * beta1) -
            dimension * cov / beta1)
    };
    const auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
    return GaussianModel3DUncertainty{ std::sqrt(var_amplitude), std::sqrt(var_width), 0.0 };
}

} // namespace

LinearizationSpec LinearizationSpec::BondDecode()
{
    auto spec{ LinearizationSpec{} };
    spec.model_kind = GaussianModelKind::MODEL_2D;
    return spec;
}

SeriesPointList BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    const LinearizationRange & fit_range)
{
    if (std::isnan(fit_range.min) || std::isnan(fit_range.max))
    {
        throw std::invalid_argument("linearization range values must not be NaN.");
    }
    if (fit_range.min > fit_range.max)
    {
        throw std::invalid_argument("linearization range min must not exceed max.");
    }

    SeriesPointList basis_and_response_entry_list;
    basis_and_response_entry_list.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ sample.distance };
        const auto gaussian_response{ static_cast<double>(sample.response) };
        if (distance < static_cast<float>(fit_range.min) ||
            distance > static_cast<float>(fit_range.max))
        {
            continue;
        }
        if (sample.score <= 0.0f || gaussian_response <= 0.0)
        {
            continue;
        }

        auto [basis_values, transformed_response]{
            BuildLogQuadraticBasisVector(static_cast<double>(distance), gaussian_response)
        };
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{
                std::move(basis_values),
                transformed_response,
                sample.score
            });
    }

    if (basis_and_response_entry_list.empty())
    {
        Logger::Log(LogLevel::Warning,
            "linearization_service::BuildDatasetSeries : "
            "No valid gaus data entry in the specified range.");
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{ std::vector<double>(kLogQuadraticBasisSize, 0.0), 0.0 });
    }
    basis_and_response_entry_list.shrink_to_fit();
    return basis_and_response_entry_list;
}

RHBMParameterVector EncodeGaussianToParameterVector(
    const LinearizationSpec & spec,
    const GaussianModel3D & gaussian_model)
{
    (void)spec;
    return BuildParameterVector(gaussian_model);
}

GaussianModel3D DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector)
{
    eigen_validation::RequireVectorSize(parameter_vector, kLogQuadraticBasisSize, "parameter_vector");
    const auto dimension{ ResolveLogQuadraticDimension(spec) };
    return DecodeLogQuadratic(parameter_vector, dimension);
}

GaussianModel3DWithUncertainty DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    eigen_validation::RequireVectorSize(parameter_vector, kLogQuadraticBasisSize, "parameter_vector");
    eigen_validation::RequireShape(
        covariance_matrix, kLogQuadraticBasisSize, kLogQuadraticBasisSize, "covariance_matrix");
    const auto dimension{ ResolveLogQuadraticDimension(spec) };
    const auto model{ DecodeLogQuadratic(parameter_vector, dimension) };
    return GaussianModel3DWithUncertainty{
        model,
        CalculateLogQuadraticStandardDeviation(parameter_vector, covariance_matrix, dimension, model)
    };
}

} // namespace rhbm_gem::linearization_service
