#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace rhbm_gem {

namespace {

GaussianEstimate BuildGaussianEstimate(const Eigen::VectorXd & gaussian_parameters)
{
    if (gaussian_parameters.rows() < 2)
    {
        throw std::invalid_argument("Gaussian parameter vector must have at least two entries.");
    }
    return GaussianEstimate{ gaussian_parameters(0), gaussian_parameters(1) };
}

Eigen::VectorXd GetTaylorSeriesBasisVector(
    double distance,
    const Eigen::VectorXd & model_parameters)
{
    const auto amplitude_0{ model_parameters(0) };
    const auto width_0{ model_parameters(1) };
    const auto width_square{ width_0 * width_0 };
    const auto gaussian_0{
        1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * distance * distance / width_square)
    };
    Eigen::VectorXd basis_vector{ Eigen::VectorXd::Zero(3) };
    basis_vector(0) = gaussian_0;
    basis_vector(1) =
        amplitude_0 * gaussian_0 * (-3.0 / width_0 + distance * distance / std::pow(width_0, 3));
    basis_vector(2) = 1.0;
    return basis_vector;
}

Eigen::VectorXd BuildLinearModelDataVector(
    double x,
    double y,
    const Eigen::VectorXd & model_parameters,
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

Eigen::VectorXd BuildLinearModelCoefficientVector(double amplitude, double width)
{
    Eigen::VectorXd linear_model_coeff{ Eigen::VectorXd::Zero(3) };
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

Eigen::VectorXd BuildGaus2DModel(const Eigen::VectorXd & linear_model)
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

Eigen::VectorXd BuildGaus3DModel(const Eigen::VectorXd & linear_model)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(3) };
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

Eigen::VectorXd BuildGaus3DModel(
    const Eigen::VectorXd & linear_model,
    const Eigen::VectorXd & model_parameters)
{
    Eigen::VectorXd gaus_model{ Eigen::VectorXd::Zero(3) };
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

std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus2DModelWithVariance(
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

std::tuple<Eigen::VectorXd, Eigen::VectorXd> BuildGaus3DModelWithVariance(
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

} // namespace

GaussianLinearizationSpec GaussianLinearizationSpec::DefaultDataset()
{
    return GaussianLinearizationSpec{};
}

GaussianLinearizationSpec GaussianLinearizationSpec::DefaultMetricModel()
{
    return GaussianLinearizationSpec{};
}

GaussianLinearizationSpec GaussianLinearizationSpec::AtomLocalDecode()
{
    return GaussianLinearizationSpec{};
}

GaussianLinearizationSpec GaussianLinearizationSpec::AtomGroupDecode()
{
    return GaussianLinearizationSpec{};
}

GaussianLinearizationSpec GaussianLinearizationSpec::BondGroupDecode()
{
    auto spec{ GaussianLinearizationSpec{} };
    spec.model_kind = GaussianModelKind::MODEL_2D;
    return spec;
}

bool GaussianLinearizationContext::HasModelParameters() const
{
    return model_parameters.rows() > 0;
}

GaussianLinearizationContext GaussianLinearizationContext::FromModelParameters(
    const Eigen::VectorXd & model_parameters)
{
    GaussianLinearizationContext context;
    context.model_parameters = model_parameters;
    return context;
}

GaussianLinearizationService::GaussianLinearizationService(GaussianLinearizationSpec spec) :
    m_spec{ std::move(spec) }
{
    ValidateSpec();
}

SeriesPointList GaussianLinearizationService::BuildDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max,
    const GaussianLinearizationContext & context) const
{
    ValidateContextIfRequired(context);

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
                m_spec.basis_size)
        };
        std::vector<double> basis_values;
        basis_values.reserve(static_cast<std::size_t>(m_spec.basis_size));
        for (int i = 0; i < m_spec.basis_size; i++)
        {
            basis_values.emplace_back(data_vector(i));
        }
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{
                std::move(basis_values),
                data_vector(m_spec.basis_size),
                sample.weight
            });
    }

    if (basis_and_response_entry_list.empty())
    {
        Logger::Log(LogLevel::Warning,
            "GaussianLinearizationService::BuildDatasetSeries : "
            "No valid gaus data entry in the specified range.");
        basis_and_response_entry_list.emplace_back(
            SeriesPoint{ std::vector<double>(static_cast<std::size_t>(m_spec.basis_size), 0.0), 0.0 });
    }
    basis_and_response_entry_list.shrink_to_fit();
    return basis_and_response_entry_list;
}

SeriesPointList GaussianLinearizationService::BuildLinearModelSeries(
    const LocalPotentialSampleList & sampling_entries,
    const GaussianLinearizationContext & context) const
{
    ValidateContextIfRequired(context);
    if (m_spec.basis_size < 2)
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
        if (sample.weight <= 0.0f || sample.response <= 0.0f)
        {
            continue;
        }

        const auto data_vector{
            BuildLinearModelDataVector(
                static_cast<double>(sample.distance),
                static_cast<double>(sample.response),
                model_parameters,
                m_spec.basis_size)
        };
        linear_model_series.emplace_back(
            SeriesPoint{ { data_vector(1) }, data_vector(m_spec.basis_size), sample.weight });
    }
    return linear_model_series;
}

HRLMemberDataset GaussianLinearizationService::BuildDataset(
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max,
    const GaussianLinearizationContext & context) const
{
    return HRLDataTransform::BuildMemberDataset(
        BuildDatasetSeries(sampling_entries, x_min, x_max, context));
}

Eigen::VectorXd GaussianLinearizationService::EncodeGaussianToBeta(
    const Eigen::VectorXd & gaussian_parameters) const
{
    if (m_spec.linearization_kind != GaussianLinearizationKind::LOG_QUADRATIC)
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
    if (m_spec.basis_size > encoded.rows())
    {
        throw std::invalid_argument("Requested basis size exceeds supported encoded Gaussian size.");
    }
    return encoded.head(m_spec.basis_size);
}

Eigen::VectorXd GaussianLinearizationService::EncodeGaussianToBeta(
    const GaussianEstimate & gaussian_estimate) const
{
    Eigen::VectorXd gaussian_parameters{ Eigen::VectorXd::Zero(2) };
    gaussian_parameters(0) = gaussian_estimate.amplitude;
    gaussian_parameters(1) = gaussian_estimate.width;
    return EncodeGaussianToBeta(gaussian_parameters);
}

Eigen::VectorXd GaussianLinearizationService::DecodeLocalBeta(
    const Eigen::VectorXd & linear_model,
    const GaussianLinearizationContext & context) const
{
    return BuildGaussianVector(linear_model, &context);
}

Eigen::VectorXd GaussianLinearizationService::DecodeGroupBeta(const Eigen::VectorXd & linear_model) const
{
    return BuildGaussianVector(linear_model, nullptr);
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> GaussianLinearizationService::DecodePosterior(
    const Eigen::VectorXd & linear_model,
    const Eigen::MatrixXd & covariance_matrix) const
{
    switch (m_spec.model_kind)
    {
    case GaussianModelKind::MODEL_2D:
        return BuildGaus2DModelWithVariance(linear_model, covariance_matrix);
    case GaussianModelKind::MODEL_3D:
        return BuildGaus3DModelWithVariance(linear_model, covariance_matrix);
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

GaussianEstimate GaussianLinearizationService::DecodeLocalEstimate(
    const Eigen::VectorXd & linear_model,
    const GaussianLinearizationContext & context) const
{
    return BuildGaussianEstimate(DecodeLocalBeta(linear_model, context));
}

GaussianEstimate GaussianLinearizationService::DecodeGroupEstimate(
    const Eigen::VectorXd & linear_model) const
{
    return BuildGaussianEstimate(DecodeGroupBeta(linear_model));
}

GaussianPosterior GaussianLinearizationService::DecodePosteriorEstimate(
    const Eigen::VectorXd & linear_model,
    const Eigen::MatrixXd & covariance_matrix) const
{
    const auto decoded{ DecodePosterior(linear_model, covariance_matrix) };
    return GaussianPosterior{
        BuildGaussianEstimate(std::get<0>(decoded)),
        BuildGaussianEstimate(std::get<1>(decoded))
    };
}

void GaussianLinearizationService::ValidateSpec() const
{
    if (m_spec.basis_size <= 0)
    {
        throw std::invalid_argument("GaussianLinearizationSpec basis_size must be positive.");
    }
}

void GaussianLinearizationService::ValidateContextIfRequired(
    const GaussianLinearizationContext & context) const
{
    if (m_spec.requires_local_context && !context.HasModelParameters())
    {
        throw std::invalid_argument("Gaussian linearization context is required for this spec.");
    }
}

Eigen::VectorXd GaussianLinearizationService::BuildGaussianVector(
    const Eigen::VectorXd & linear_model,
    const GaussianLinearizationContext * context) const
{
    switch (m_spec.model_kind)
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

} // namespace rhbm_gem
