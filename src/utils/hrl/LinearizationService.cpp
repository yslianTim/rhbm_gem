#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace rhbm_gem::linearization_service
{

namespace {

GaussianModel3D BuildGaussianModel(const Eigen::VectorXd & gaussian_parameters)
{
    if (gaussian_parameters.rows() < 2)
    {
        throw std::invalid_argument("Gaussian parameter vector must have at least two entries.");
    }
    return GaussianModel3D{
        gaussian_parameters(GaussianModel3D::AmplitudeIndex()),
        gaussian_parameters(GaussianModel3D::WidthIndex()),
        (gaussian_parameters.rows() > GaussianModel3D::InterceptIndex()) ?
            gaussian_parameters(GaussianModel3D::InterceptIndex()) : 0.0
    };
}

GaussianModel3DUncertainty BuildGaussianModelUncertainty(
    const Eigen::VectorXd & gaussian_parameters)
{
    if (gaussian_parameters.rows() < 2)
    {
        throw std::invalid_argument("Gaussian parameter vector must have at least two entries.");
    }
    return GaussianModel3DUncertainty{
        gaussian_parameters(GaussianModel3D::AmplitudeIndex()),
        gaussian_parameters(GaussianModel3D::WidthIndex()),
        (gaussian_parameters.rows() > GaussianModel3D::InterceptIndex()) ?
            gaussian_parameters(GaussianModel3D::InterceptIndex()) : 0.0
    };
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

Eigen::VectorXd BuildGaus2DModel(const RHBMParameterVector & linear_model)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
    if (linear_model(1) <= 0.0)
    {
        return gaus_model;
    }
    gaus_model(0) = std::exp(linear_model(0)) * Constants::two_pi / linear_model(1);
    gaus_model(1) = 1.0 / std::sqrt(linear_model(1));
    return gaus_model;
}

GaussianModel3D BuildGaus3DModel(const RHBMParameterVector & linear_model)
{
    if (linear_model(1) <= 0.0)
    {
        return GaussianModel3D{ 0.0, 0.0, 0.0 };
    }
    return GaussianModel3D{
        std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5),
        1.0 / std::sqrt(linear_model(1)),
        0.0
    };
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus2DModelWithUncertainty(
    const RHBMParameterVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
    Eigen::VectorXd gaus_model_variance{ Eigen::VectorXd::Zero(2) };

    if (linear_model.rows() == 2)
    {
        try
        {
            eigen_validation::RequireShape(covariance_matrix, 2, 2, "covariance_matrix");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("covariance_matrix must be 2x2");
        }
        if (!covariance_matrix.isApprox(covariance_matrix.transpose()))
        {
            throw std::invalid_argument("covariance_matrix must be symmetric");
        }

        gaus_model = BuildGaus2DModel(linear_model);
        const auto beta0{ linear_model(0) };
        const auto beta1{ linear_model(1) };
        if (beta1 <= 0.0)
        {
            return std::make_tuple(gaus_model, gaus_model_variance);
        }

        const auto var_beta0{ covariance_matrix(0, 0) };
        const auto var_beta1{ covariance_matrix(1, 1) };
        const auto cov{ covariance_matrix(0, 1) };
        const auto var_amplitude{
            std::pow(Constants::two_pi, 3) * std::exp(2.0 * beta0) * std::pow(beta1, -5) *
            (std::pow(beta1, 2) * var_beta0 + 2.25 * var_beta1 - 3.0 * beta1 * cov)
        };
        const auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
        gaus_model_variance(0) = std::sqrt(var_amplitude);
        gaus_model_variance(1) = std::sqrt(var_width);
    }

    return std::make_tuple(gaus_model, gaus_model_variance);
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus3DModelWithUncertainty(
    const RHBMParameterVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()) };
    Eigen::VectorXd gaus_model_variance{
        Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize())
    };

    if (linear_model.rows() == 2)
    {
        try
        {
            eigen_validation::RequireShape(covariance_matrix, 2, 2, "covariance_matrix");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("covariance_matrix must be 2x2");
        }
        if (!covariance_matrix.isApprox(covariance_matrix.transpose()))
        {
            throw std::invalid_argument("covariance_matrix must be symmetric");
        }

        gaus_model = BuildGaus3DModel(linear_model).ToVector();
        const auto beta0{ linear_model(0) };
        const auto beta1{ linear_model(1) };
        if (beta1 <= 0.0)
        {
            return std::make_tuple(gaus_model, gaus_model_variance);
        }

        const auto var_beta0{ covariance_matrix(0, 0) };
        const auto var_beta1{ covariance_matrix(1, 1) };
        const auto cov{ covariance_matrix(0, 1) };
        const auto var_amplitude{
            std::pow(Constants::two_pi, 3) * std::exp(2.0 * beta0) * std::pow(beta1, -5) *
            (std::pow(beta1, 2) * var_beta0 + 2.25 * var_beta1 - 3.0 * beta1 * cov)
        };
        const auto var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
        gaus_model_variance(0) = std::sqrt(var_amplitude);
        gaus_model_variance(1) = std::sqrt(var_width);
        return std::make_tuple(gaus_model, gaus_model_variance);
    }

    if (linear_model.rows() == GaussianModel3D::ParameterSize())
    {
        try
        {
            eigen_validation::RequireShape(
                covariance_matrix,
                GaussianModel3D::ParameterSize(),
                GaussianModel3D::ParameterSize(),
                "covariance_matrix");
        }
        catch (const std::invalid_argument &)
        {
            throw std::invalid_argument("covariance_matrix must be 3x3");
        }
        if (!covariance_matrix.isApprox(covariance_matrix.transpose()))
        {
            throw std::invalid_argument("covariance_matrix must be symmetric");
        }

        gaus_model = BuildGaus3DModel(linear_model).ToVector();
        for (int i = 0; i < GaussianModel3D::ParameterSize(); i++)
        {
            gaus_model_variance(i) = std::sqrt(covariance_matrix(i, i));
        }
    }

    return std::make_tuple(gaus_model, gaus_model_variance);
}

void ValidateContextIfRequired(
    const LinearizationSpec & spec,
    const LinearizationContext & context)
{
    if (spec.requires_local_context && !context.HasModelParameters())
    {
        throw std::invalid_argument("Gaussian linearization context is required for this spec.");
    }
}

Eigen::VectorXd BuildGaussianVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & linear_model)
{
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return BuildGaus2DModel(linear_model);
    case GaussianModelKind::MODEL_3D:
        return BuildGaus3DModel(linear_model).ToVector();
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> DecodePosteriorParameters(
    const LinearizationSpec & spec,
    const RHBMParameterVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return BuildGaus2DModelWithUncertainty(linear_model, covariance_matrix);
    case GaussianModelKind::MODEL_3D:
        return BuildGaus3DModelWithUncertainty(linear_model, covariance_matrix);
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

} // namespace

LinearizationSpec LinearizationSpec::BondGroupDecode()
{
    auto spec{ LinearizationSpec{} };
    spec.model_kind = GaussianModelKind::MODEL_2D;
    return spec;
}

LinearizationContext LinearizationContext::FromModel(const GaussianModel3D & model)
{
    LinearizationContext context;
    context.model = model;
    return context;
}

SeriesPointList BuildDatasetSeries(
    const LinearizationSpec & spec,
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max,
    const LinearizationContext & context)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    if (spec.basis_size != 2)
    {
        throw std::invalid_argument("BuildDatasetSeries only supports log-quadratic basis_size == 2.");
    }
    ValidateContextIfRequired(spec, context);

    SeriesPointList basis_and_response_entry_list;
    basis_and_response_entry_list.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ sample.distance };
        const auto gaussian_response{ static_cast<double>(sample.response) };
        if (distance < static_cast<float>(x_min) || distance > static_cast<float>(x_max))
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
            SeriesPoint{ std::vector<double>(static_cast<std::size_t>(spec.basis_size), 0.0), 0.0 });
    }
    basis_and_response_entry_list.shrink_to_fit();
    return basis_and_response_entry_list;
}

RHBMParameterVector EncodeGaussianToParameterVector(
    const LinearizationSpec & spec,
    const GaussianModel3D & gaussian_model)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");

    const auto encoded{ BuildLinearModelCoefficientVector(gaussian_model) };
    if (spec.basis_size > encoded.rows())
    {
        throw std::invalid_argument("Requested basis size exceeds supported encoded Gaussian size.");
    }
    return encoded.head(spec.basis_size);
}

GaussianModel3D DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    return BuildGaussianModel(BuildGaussianVector(spec, parameter_vector));
}

GaussianModel3DWithUncertainty DecodeParameterVector(
    const LinearizationSpec & spec,
    const RHBMParameterVector & parameter_vector,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    const auto decoded{ DecodePosteriorParameters(spec, parameter_vector, covariance_matrix) };
    return GaussianModel3DWithUncertainty{
        BuildGaussianModel(std::get<0>(decoded)),
        BuildGaussianModelUncertainty(std::get<1>(decoded))
    };
}

} // namespace rhbm_gem::linearization_service
