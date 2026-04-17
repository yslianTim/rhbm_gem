#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>

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
            GausLinearTransformHelper::BuildLinearModelDataVector(
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
        GausLinearTransformHelper::BuildLinearModelCoefficentVector(
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
        return GausLinearTransformHelper::BuildGaus2DModelWithVariance(
            linear_model,
            covariance_matrix);
    case GaussianModelKind::MODEL_3D:
        return GausLinearTransformHelper::BuildGaus3DModelWithVariance(
            linear_model,
            covariance_matrix);
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
        return GausLinearTransformHelper::BuildGaus2DModel(linear_model);
    case GaussianModelKind::MODEL_3D:
        if (context != nullptr && context->HasModelParameters())
        {
            return GausLinearTransformHelper::BuildGaus3DModel(
                linear_model,
                context->model_parameters);
        }
        return GausLinearTransformHelper::BuildGaus3DModel(linear_model);
    }
    throw std::invalid_argument("Unsupported Gaussian model kind.");
}

} // namespace rhbm_gem
