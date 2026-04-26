#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::linearization_service
{

namespace {

GaussianEstimate BuildGaussianEstimate(const GaussianParameterVector & gaussian_parameters)
{
    if (gaussian_parameters.rows() < 2)
    {
        throw std::invalid_argument("Gaussian parameter vector must have at least two entries.");
    }
    return GaussianEstimate{ gaussian_parameters(0), gaussian_parameters(1) };
}

GaussianParameterVector GetTaylorSeriesBasisVector(
    double distance,
    const GaussianParameterVector & model_parameters)
{
    const auto amplitude_0{ model_parameters(0) };
    const auto width_0{ model_parameters(1) };
    const auto width_square{ width_0 * width_0 };
    const auto gaussian_0{
        1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * distance * distance / width_square)
    };
    GaussianParameterVector basis_vector{ GaussianParameterVector::Zero(3) };
    basis_vector(0) = gaussian_0;
    basis_vector(1) =
        amplitude_0 * gaussian_0 * (-3.0 / width_0 + distance * distance / std::pow(width_0, 3));
    basis_vector(2) = 1.0;
    return basis_vector;
}

GaussianParameterVector BuildLinearModelDataVector(
    double x,
    double y,
    const GaussianParameterVector & model_parameters,
    int basis_dimension)
{
    const auto data_vector_dimension{ basis_dimension + 1 };
    GaussianParameterVector linear_model_data_vector{
        GaussianParameterVector::Zero(data_vector_dimension)
    };

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

    const auto basis_vector{ GetTaylorSeriesBasisVector(x, model_parameters) };
    const auto amplitude_0{ model_parameters(0) };
    const auto width_0{ model_parameters(1) };
    for (int i = 0; i < basis_dimension; i++)
    {
        linear_model_data_vector(i) = basis_vector(i);
    }
    const auto width_square{ width_0 * width_0 };
    const auto gaussian_response{
        1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * x * x / width_square)
    };
    const auto intercept{ amplitude_0 * gaussian_response + model_parameters(2) };
    linear_model_data_vector(basis_dimension) = y - intercept;
    return linear_model_data_vector;
}

RHBMBetaVector BuildLinearModelCoefficientVector(double amplitude, double width)
{
    RHBMBetaVector linear_model_coeff{ RHBMBetaVector::Zero(3) };
    if (width == 0.0)
    {
        return linear_model_coeff;
    }

    const auto width_square{ width * width };
    linear_model_coeff(0) = (amplitude <= 0.0) ?
        0.0 : std::log(amplitude) - 1.5 * std::log(Constants::two_pi * width_square);
    linear_model_coeff(1) = 1.0 / width_square;
    linear_model_coeff(2) = 0.0;
    return linear_model_coeff;
}

GaussianParameterVector BuildGaus2DModel(const RHBMBetaVector & linear_model)
{
    GaussianParameterVector gaus_model{ GaussianParameterVector::Zero(2) };
    if (linear_model(1) <= 0.0)
    {
        return gaus_model;
    }

    const auto model_dimension{ linear_model.rows() };
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

GaussianParameterVector BuildGaus3DModel(const RHBMBetaVector & linear_model)
{
    GaussianParameterVector gaus_model{ GaussianParameterVector::Zero(3) };
    if (linear_model(1) <= 0.0)
    {
        return gaus_model;
    }

    const auto model_dimension{ linear_model.rows() };
    if (model_dimension == 2)
    {
        gaus_model(0) =
            std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5);
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

GaussianParameterVector BuildGaus3DModel(
    const RHBMBetaVector & linear_model,
    const GaussianParameterVector & model_parameters)
{
    GaussianParameterVector gaus_model{ GaussianParameterVector::Zero(3) };
    if (linear_model(1) <= 0.0)
    {
        return gaus_model;
    }

    const auto model_dimension{ linear_model.rows() };
    if (model_dimension == 2)
    {
        gaus_model(0) =
            std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5);
        gaus_model(1) = 1.0 / std::sqrt(linear_model(1));
        gaus_model(2) = 0.0;
    }
    else
    {
        gaus_model(0) = std::exp(linear_model(0)) * model_parameters(0);
        gaus_model(1) = std::exp(linear_model(1)) + model_parameters(1);
        gaus_model(2) = linear_model(2) + model_parameters(2);
    }

    return gaus_model;
}

std::tuple<GaussianParameterVector, GaussianParameterVector> BuildGaus2DModelWithVariance(
    const RHBMBetaVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    GaussianParameterVector gaus_model{ GaussianParameterVector::Zero(2) };
    GaussianParameterVector gaus_model_variance{ GaussianParameterVector::Zero(2) };

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

std::tuple<GaussianParameterVector, GaussianParameterVector> BuildGaus3DModelWithVariance(
    const RHBMBetaVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    GaussianParameterVector gaus_model{ GaussianParameterVector::Zero(2) };
    GaussianParameterVector gaus_model_variance{ GaussianParameterVector::Zero(2) };

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

        gaus_model = BuildGaus3DModel(linear_model);
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

void ValidateContextIfRequired(
    const LinearizationSpec & spec,
    const LinearizationContext & context)
{
    if (spec.requires_local_context && !context.HasModelParameters())
    {
        throw std::invalid_argument("Gaussian linearization context is required for this spec.");
    }
}

GaussianParameterVector BuildGaussianVector(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext * context)
{
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return BuildGaus2DModel(linear_model);
    case GaussianModelKind::MODEL_3D:
        if (context != nullptr && context->HasModelParameters())
        {
            return BuildGaus3DModel(linear_model, context->model_parameters);
        }
        return BuildGaus3DModel(linear_model);
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

} // namespace

LinearizationSpec LinearizationSpec::DefaultDataset()
{
    return LinearizationSpec{};
}

LinearizationSpec LinearizationSpec::DefaultMetricModel()
{
    return LinearizationSpec{};
}

LinearizationSpec LinearizationSpec::AtomLocalDecode()
{
    return LinearizationSpec{};
}

LinearizationSpec LinearizationSpec::AtomGroupDecode()
{
    return LinearizationSpec{};
}

LinearizationSpec LinearizationSpec::BondGroupDecode()
{
    auto spec{ LinearizationSpec{} };
    spec.model_kind = GaussianModelKind::MODEL_2D;
    return spec;
}

bool LinearizationContext::HasModelParameters() const
{
    return model_parameters.rows() > 0;
}

LinearizationContext LinearizationContext::FromModelParameters(
    const GaussianParameterVector & model_parameters)
{
    LinearizationContext context;
    context.model_parameters = model_parameters;
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
    ValidateContextIfRequired(spec, context);

    SeriesPointList basis_and_response_entry_list;
    basis_and_response_entry_list.reserve(sampling_entries.size());
    const auto model_parameters{
        context.HasModelParameters() ? context.model_parameters : Eigen::VectorXd::Zero(3)
    };
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ sample.distance };
        const auto gaussian_response{ static_cast<double>(sample.response) };
        if (distance < static_cast<float>(x_min) || distance > static_cast<float>(x_max))
        {
            continue;
        }
        if (gaussian_response <= 0.0)
        {
            continue;
        }

        const auto data_vector{
            BuildLinearModelDataVector(
                static_cast<double>(distance),
                gaussian_response,
                model_parameters,
                spec.basis_size)
        };
        std::vector<double> basis_values;
        basis_values.reserve(static_cast<std::size_t>(spec.basis_size));
        for (int i = 0; i < spec.basis_size; i++)
        {
            basis_values.emplace_back(data_vector(i));
        }
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{
                std::move(basis_values),
                data_vector(spec.basis_size),
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

SeriesPointList BuildLinearModelSeries(
    const LinearizationSpec & spec,
    const LocalPotentialSampleList & sampling_entries,
    const LinearizationContext & context)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    ValidateContextIfRequired(spec, context);
    if (spec.basis_size < 2)
    {
        throw std::invalid_argument("Linear model series requires basis_size >= 2.");
    }

    SeriesPointList linear_model_series;
    linear_model_series.reserve(sampling_entries.size());
    const auto model_parameters{
        context.HasModelParameters() ? context.model_parameters : Eigen::VectorXd::Zero(3)
    };
    for (const auto & sample : sampling_entries)
    {
        if (sample.score <= 0.0f || sample.response <= 0.0f)
        {
            continue;
        }

        const auto data_vector{
            BuildLinearModelDataVector(
                static_cast<double>(sample.distance),
                static_cast<double>(sample.response),
                model_parameters,
                spec.basis_size)
        };
        linear_model_series.emplace_back(
            SeriesPoint{ { data_vector(1) }, data_vector(spec.basis_size), sample.score });
    }
    return linear_model_series;
}

RHBMMemberDataset BuildDataset(
    const LinearizationSpec & spec,
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max,
    const LinearizationContext & context)
{
    return rhbm_helper::BuildMemberDataset(
        BuildDatasetSeries(spec, sampling_entries, x_min, x_max, context));
}

RHBMBetaVector EncodeGaussianToBeta(
    const LinearizationSpec & spec,
    const GaussianParameterVector & gaussian_parameters)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    if (spec.linearization_kind != LinearizationKind::LOG_QUADRATIC)
    {
        throw std::invalid_argument("EncodeGaussianToBeta only supports log-quadratic linearization.");
    }
    if (gaussian_parameters.rows() < 2)
    {
        throw std::invalid_argument("Gaussian parameter vector must have at least two entries.");
    }

    const auto encoded{
        BuildLinearModelCoefficientVector(
            gaussian_parameters(0),
            gaussian_parameters(1))
    };
    if (spec.basis_size > encoded.rows())
    {
        throw std::invalid_argument("Requested basis size exceeds supported encoded Gaussian size.");
    }
    return encoded.head(spec.basis_size);
}

RHBMBetaVector EncodeGaussianToBeta(
    const LinearizationSpec & spec,
    const GaussianEstimate & gaussian_estimate)
{
    GaussianParameterVector gaussian_parameters{ GaussianParameterVector::Zero(2) };
    gaussian_parameters(0) = gaussian_estimate.amplitude;
    gaussian_parameters(1) = gaussian_estimate.width;
    return EncodeGaussianToBeta(spec, gaussian_parameters);
}

GaussianParameterVector DecodeLocalBeta(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    ValidateContextIfRequired(spec, context);
    return BuildGaussianVector(spec, linear_model, &context);
}

GaussianParameterVector DecodeGroupBeta(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    return BuildGaussianVector(spec, linear_model, nullptr);
}

std::tuple<GaussianParameterVector, GaussianParameterVector> DecodePosterior(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return BuildGaus2DModelWithVariance(linear_model, covariance_matrix);
    case GaussianModelKind::MODEL_3D:
        return BuildGaus3DModelWithVariance(linear_model, covariance_matrix);
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

GaussianEstimate DecodeLocalEstimate(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context)
{
    return BuildGaussianEstimate(DecodeLocalBeta(spec, linear_model, context));
}

GaussianEstimate DecodeGroupEstimate(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model)
{
    return BuildGaussianEstimate(DecodeGroupBeta(spec, linear_model));
}

GaussianPosterior DecodePosteriorEstimate(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    const auto decoded{ DecodePosterior(spec, linear_model, covariance_matrix) };
    return GaussianPosterior{
        BuildGaussianEstimate(std::get<0>(decoded)),
        BuildGaussianEstimate(std::get<1>(decoded))
    };
}

} // namespace rhbm_gem::linearization_service
