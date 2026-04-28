#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace rhbm_gem::linearization_service
{

namespace {

constexpr int kLogQuadraticBasisSize{ 2 };

void RequireParameterVectorSize(const RHBMParameterVector & parameter_vector)
{
    if (parameter_vector.rows() != kLogQuadraticBasisSize)
    {
        throw std::invalid_argument("parameter_vector size must match log-quadratic basis size.");
    }
}

void RequireLinearizationRange(const LinearizationRange & range)
{
    if (std::isnan(range.min) || std::isnan(range.max))
    {
        throw std::invalid_argument("linearization range values must not be NaN.");
    }
    if (range.min > range.max)
    {
        throw std::invalid_argument("linearization range min must not exceed max.");
    }
}

void RequireCovarianceMatrix(
    const RHBMPosteriorCovarianceMatrix & covariance_matrix,
    Eigen::Index expected_size)
{
    try
    {
        eigen_validation::RequireShape(
            covariance_matrix,
            expected_size,
            expected_size,
            "covariance_matrix");
    }
    catch (const std::invalid_argument &)
    {
        throw std::invalid_argument(
            "covariance_matrix must be " + std::to_string(expected_size) +
            "x" + std::to_string(expected_size));
    }
    if (!covariance_matrix.isApprox(covariance_matrix.transpose()))
    {
        throw std::invalid_argument("covariance_matrix must be symmetric");
    }
}

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

RHBMParameterVector BuildLinearModelCoefficientVector(const GaussianModel3D & model)
{
    RHBMParameterVector linear_model_coeff{ RHBMParameterVector::Zero(3) };
    const auto width{ model.GetWidth() };
    const auto amplitude{ model.GetAmplitude() };
    if (width == 0.0)
    {
        return linear_model_coeff;
    }

    const auto width_square{ width * width };
    linear_model_coeff(0) = (amplitude <= 0.0) ?
        0.0 : std::log(amplitude) - 1.5 * std::log(Constants::two_pi * width_square);
    linear_model_coeff(1) = 1.0 / width_square;
    linear_model_coeff(2) = model.GetIntercept();
    return linear_model_coeff;
}

GaussianModel3D DecodeLogQuadratic2D(const RHBMParameterVector & linear_model)
{
    if (linear_model(1) <= 0.0)
    {
        return GaussianModel3D{ 0.0, 0.0, 0.0 };
    }
    return GaussianModel3D{
        std::exp(linear_model(0)) * Constants::two_pi / linear_model(1),
        1.0 / std::sqrt(linear_model(1)),
        0.0
    };
}

GaussianModel3D DecodeLogQuadratic3D(const RHBMParameterVector & linear_model)
{
    const auto intercept{
        (linear_model.rows() > GaussianModel3D::InterceptIndex()) ?
            linear_model(GaussianModel3D::InterceptIndex()) : 0.0
    };
    if (linear_model(1) <= 0.0)
    {
        return GaussianModel3D{ 0.0, 0.0, intercept };
    }
    return GaussianModel3D{
        std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5),
        1.0 / std::sqrt(linear_model(1)),
        intercept
    };
}

GaussianModel3DUncertainty CalculateLogQuadratic2DStandardDeviation(
    const RHBMParameterVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    const auto beta1{ linear_model(1) };
    if (beta1 <= 0.0)
    {
        return GaussianModel3DUncertainty{};
    }

    const auto amplitude{ DecodeLogQuadratic2D(linear_model).GetAmplitude() };
    const auto var_beta0{ covariance_matrix(0, 0) };
    const auto var_beta1{ covariance_matrix(1, 1) };
    const auto cov{ covariance_matrix(0, 1) };
    const auto var_amplitude{
        amplitude * amplitude *
        (var_beta0 + var_beta1 / (beta1 * beta1) - 2.0 * cov / beta1)
    };
    const auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
    return GaussianModel3DUncertainty{
        std::sqrt(var_amplitude),
        std::sqrt(var_width),
        0.0
    };
}

GaussianModel3DUncertainty CalculateLogQuadratic3DStandardDeviation(
    const RHBMParameterVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    const auto beta1{ linear_model(1) };
    const auto intercept_standard_deviation{
        (linear_model.rows() > GaussianModel3D::InterceptIndex()) ?
            std::sqrt(covariance_matrix(
                GaussianModel3D::InterceptIndex(),
                GaussianModel3D::InterceptIndex())) : 0.0
    };
    if (beta1 <= 0.0)
    {
        return GaussianModel3DUncertainty{ 0.0, 0.0, intercept_standard_deviation };
    }

    const auto amplitude{ DecodeLogQuadratic3D(linear_model).GetAmplitude() };
    const auto var_beta0{ covariance_matrix(0, 0) };
    const auto var_beta1{ covariance_matrix(1, 1) };
    const auto cov{ covariance_matrix(0, 1) };
    const auto var_amplitude{
        amplitude * amplitude *
        (var_beta0 + 2.25 * var_beta1 / (beta1 * beta1) - 3.0 * cov / beta1)
    };
    const auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
    return GaussianModel3DUncertainty{
        std::sqrt(var_amplitude),
        std::sqrt(var_width),
        intercept_standard_deviation
    };
}

GaussianModel3D DecodeGaussianBySpec(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector)
{
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return DecodeLogQuadratic2D(parameter_vector);
    case GaussianModelKind::MODEL_3D:
        return DecodeLogQuadratic3D(parameter_vector);
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

GaussianModel3DWithUncertainty DecodeGaussianWithUncertaintyBySpec(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    RequireCovarianceMatrix(covariance_matrix, parameter_vector.rows());
    const auto model{ DecodeGaussianBySpec(spec, parameter_vector) };
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return GaussianModel3DWithUncertainty{
            model,
            CalculateLogQuadratic2DStandardDeviation(parameter_vector, covariance_matrix)
        };
    case GaussianModelKind::MODEL_3D:
        return GaussianModel3DWithUncertainty{
            model,
            CalculateLogQuadratic3DStandardDeviation(parameter_vector, covariance_matrix)
        };
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
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
    RequireLinearizationRange(fit_range);

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
    const auto encoded{ BuildLinearModelCoefficientVector(gaussian_model) };
    return encoded.head(kLogQuadraticBasisSize);
}

GaussianModel3D DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector)
{
    RequireParameterVectorSize(parameter_vector);
    return DecodeGaussianBySpec(spec, parameter_vector);
}

GaussianModel3DWithUncertainty DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    RequireParameterVectorSize(parameter_vector);
    return DecodeGaussianWithUncertaintyBySpec(spec, parameter_vector, covariance_matrix);
}

} // namespace rhbm_gem::linearization_service
