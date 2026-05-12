#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::linearization_service
{

namespace {
constexpr int kLogQuadraticBasisSize{ 2 };
} // namespace

SeriesPointList BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double range_min,
    double range_max)
{
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "data range");

    SeriesPointList basis_and_response_entry_list;
    basis_and_response_entry_list.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ sample.distance };
        const auto gaussian_response{ static_cast<double>(sample.response) };
        if (distance < static_cast<float>(range_min)) continue;
        if (distance > static_cast<float>(range_max)) continue;
        if (gaussian_response <= 0.0) continue;

        auto basis_vector{ std::vector<double>{ 1.0, -0.5 * distance * distance } };
        auto transformed_response{ std::log(gaussian_response) };
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{ std::move(basis_vector), transformed_response, sample.score }
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

RHBMParameterVector EncodeGaussianToParameterVector(const GaussianModel3D & gaussian_model)
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

GaussianModel3D DecodeParameterVector(const RHBMParameterVector & parameter_vector)
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

GaussianModel3DWithUncertainty DecodeParameterVector(
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    eigen_validation::RequireShape(
        covariance_matrix, kLogQuadraticBasisSize, kLogQuadraticBasisSize, "covariance_matrix");
    const auto model{ DecodeParameterVector(parameter_vector) };
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
