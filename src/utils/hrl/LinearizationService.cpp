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

const GaussianModel3D & RequireContextModel(
    const LinearizationContext & context)
{
    if (!context.model.has_value())
    {
        throw std::invalid_argument("Gaussian linearization context is required for this spec.");
    }
    return context.model.value();
}

GaussianModel3D BuildGaussianModel(const Eigen::VectorXd & gaussian_parameters)
{
    if (gaussian_parameters.rows() < 2)
    {
        throw std::invalid_argument("Gaussian parameter vector must have at least two entries.");
    }
    return GaussianModel3D{ gaussian_parameters(0), gaussian_parameters(1) };
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

Eigen::VectorXd GetTaylorSeriesBasisVector(
    double distance,
    const GaussianModel3D & model)
{
    const auto amplitude_0{ model.GetAmplitude() };
    const auto width_0{ model.GetWidth() };
    const auto width_square{ width_0 * width_0 };
    const auto gaussian_0{
        1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * distance * distance / width_square)
    };
    Eigen::VectorXd basis_vector{
        Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize())
    };
    basis_vector(0) = gaussian_0;
    basis_vector(1) =
        amplitude_0 * gaussian_0 * (-3.0 / width_0 + distance * distance / std::pow(width_0, 3));
    basis_vector(2) = 1.0;
    return basis_vector;
}

Eigen::VectorXd BuildLinearModelDataVector(
    double x,
    double y,
    const GaussianModel3D & model,
    int basis_dimension)
{
    const auto data_vector_dimension{ basis_dimension + 1 };
    Eigen::VectorXd linear_model_data_vector{
        Eigen::VectorXd::Zero(data_vector_dimension)
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

    const auto basis_vector{ GetTaylorSeriesBasisVector(x, model) };
    const auto amplitude_0{ model.GetAmplitude() };
    const auto width_0{ model.GetWidth() };
    for (int i = 0; i < basis_dimension; i++)
    {
        linear_model_data_vector(i) = basis_vector(i);
    }
    const auto width_square{ width_0 * width_0 };
    const auto gaussian_response{
        1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * x * x / width_square)
    };
    const auto intercept{ amplitude_0 * gaussian_response + model.GetIntercept() };
    linear_model_data_vector(basis_dimension) = y - intercept;
    return linear_model_data_vector;
}

RHBMBetaVector BuildLinearModelCoefficientVector(const GaussianModel3D & model)
{
    RHBMBetaVector linear_model_coeff{ RHBMBetaVector::Zero(3) };
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

Eigen::VectorXd BuildGaus2DModel(const RHBMBetaVector & linear_model)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(2) };
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

GaussianModel3D BuildGaus3DModel(const RHBMBetaVector & linear_model)
{
    if (linear_model(1) <= 0.0)
    {
        return GaussianModel3D{ 0.0, 0.0, 0.0 };
    }

    const auto model_dimension{ linear_model.rows() };
    if (model_dimension == 2)
    {
        return GaussianModel3D{
            std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5),
            1.0 / std::sqrt(linear_model(1)),
            0.0
        };
    }

    return GaussianModel3D{ linear_model(0), linear_model(1), linear_model(2) };
}

GaussianModel3D BuildGaus3DModel(
    const RHBMBetaVector & linear_model,
    const GaussianModel3D & model)
{
    if (linear_model(1) <= 0.0)
    {
        return GaussianModel3D{ 0.0, 0.0, 0.0 };
    }

    const auto model_dimension{ linear_model.rows() };
    if (model_dimension == 2)
    {
        return GaussianModel3D{
            std::exp(linear_model(0)) * std::pow(Constants::two_pi / linear_model(1), 1.5),
            1.0 / std::sqrt(linear_model(1)),
            0.0
        };
    }

    return GaussianModel3D{
        std::exp(linear_model(0)) * model.GetAmplitude(),
        std::exp(linear_model(1)) + model.GetWidth(),
        linear_model(2) + model.GetIntercept()
    };
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus2DModelWithUncertainty(
    const RHBMBetaVector & linear_model,
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
    const RHBMBetaVector & linear_model,
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

void Validate3DModelSpec(const LinearizationSpec & spec)
{
    if (spec.model_kind != GaussianModelKind::MODEL_3D)
    {
        throw std::invalid_argument("GaussianModel3D decode requires MODEL_3D.");
    }
}

Eigen::VectorXd BuildGaussianVector(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext * context)
{
    switch (spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return BuildGaus2DModel(linear_model);
    case GaussianModelKind::MODEL_3D:
        if (context != nullptr && context->model.has_value())
        {
            return BuildGaus3DModel(linear_model, context->model.value()).ToVector();
        }
        return BuildGaus3DModel(linear_model).ToVector();
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
    return model.has_value();
}

LinearizationContext LinearizationContext::FromModelParameters(
    const Eigen::VectorXd & model_parameters)
{
    return FromModel(GaussianModel3D::FromVector(model_parameters));
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
    ValidateContextIfRequired(spec, context);

    SeriesPointList basis_and_response_entry_list;
    basis_and_response_entry_list.reserve(sampling_entries.size());
    const GaussianModel3D default_model{ 0.0, 0.0, 0.0 };
    const auto & model{ context.model.value_or(default_model) };
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
                model,
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
    const GaussianModel3D default_model{ 0.0, 0.0, 0.0 };
    const auto & model{ context.model.value_or(default_model) };
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
                model,
                spec.basis_size)
        };
        linear_model_series.emplace_back(
            SeriesPoint{ { data_vector(1) }, data_vector(spec.basis_size), sample.score });
    }
    return linear_model_series;
}

RHBMBetaVector EncodeGaussianToBeta(
    const LinearizationSpec & spec,
    const Eigen::VectorXd & gaussian_parameters)
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

    const auto gaussian_model{ (gaussian_parameters.rows() >= GaussianModel3D::ParameterSize()) ?
        GaussianModel3D::FromVector(gaussian_parameters) :
        GaussianModel3D{ gaussian_parameters(0), gaussian_parameters(1), 0.0 }
    };
    const auto encoded{ BuildLinearModelCoefficientVector(gaussian_model) };
    if (spec.basis_size > encoded.rows())
    {
        throw std::invalid_argument("Requested basis size exceeds supported encoded Gaussian size.");
    }
    return encoded.head(spec.basis_size);
}

RHBMBetaVector EncodeGaussianToBeta(
    const LinearizationSpec & spec,
    const GaussianModel3D & gaussian_model)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    if (spec.linearization_kind != LinearizationKind::LOG_QUADRATIC)
    {
        throw std::invalid_argument("EncodeGaussianToBeta only supports log-quadratic linearization.");
    }

    const auto encoded{ BuildLinearModelCoefficientVector(gaussian_model) };
    if (spec.basis_size > encoded.rows())
    {
        throw std::invalid_argument("Requested basis size exceeds supported encoded Gaussian size.");
    }
    return encoded.head(spec.basis_size);
}

Eigen::VectorXd DecodeLocalBeta(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    ValidateContextIfRequired(spec, context);
    return BuildGaussianVector(spec, linear_model, &context);
}

Eigen::VectorXd DecodeGroupBeta(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    return BuildGaussianVector(spec, linear_model, nullptr);
}

GaussianModel3D DecodeLocalModel3D(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    Validate3DModelSpec(spec);
    ValidateContextIfRequired(spec, context);
    if (context.HasModelParameters())
    {
        return BuildGaus3DModel(linear_model, RequireContextModel(context));
    }
    return BuildGaus3DModel(linear_model);
}

GaussianModel3D DecodeGroupModel3D(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model)
{
    numeric_validation::RequirePositive(spec.basis_size, "LinearizationSpec basis_size");
    Validate3DModelSpec(spec);
    return BuildGaus3DModel(linear_model);
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> DecodePosterior(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
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

GaussianModel3D DecodeLocalEstimate(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const LinearizationContext & context)
{
    return BuildGaussianModel(DecodeLocalBeta(spec, linear_model, context));
}

GaussianModel3D DecodeGroupEstimate(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model)
{
    return BuildGaussianModel(DecodeGroupBeta(spec, linear_model));
}

GaussianModel3DWithUncertainty DecodeGaussianModel3DWithUncertainty(
    const LinearizationSpec & spec,
    const RHBMBetaVector & linear_model,
    const RHBMPosteriorCovarianceMatrix & covariance_matrix)
{
    const auto decoded{ DecodePosterior(spec, linear_model, covariance_matrix) };
    return GaussianModel3DWithUncertainty{
        BuildGaussianModel(std::get<0>(decoded)),
        BuildGaussianModelUncertainty(std::get<1>(decoded))
    };
}

} // namespace rhbm_gem::linearization_service
