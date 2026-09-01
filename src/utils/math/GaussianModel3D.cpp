#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace rhbm_gem
{

namespace
{

void RequireFiniteParameterValues(const Eigen::VectorXd & parameters, std::string_view value_name)
{
    const auto name{ std::string(value_name) };
    numeric_validation::RequireFinite(parameters(GaussianModel3D::AmplitudeIndex()), name + " amplitude");
    numeric_validation::RequireFinite(parameters(GaussianModel3D::WidthIndex()), name + " width");
    numeric_validation::RequireFinite(parameters(GaussianModel3D::OffsetIndex()), name + " offset");
}

} // namespace

GaussianModel3D::GaussianModel3D(double amplitude, double width, double offset) :
    m_amplitude{ amplitude }, m_width{ width }, m_offset{ offset }
{
}

void GaussianModel3D::RequireParameterVector(
    const Eigen::VectorXd & parameters,
    std::string_view value_name)
{
    eigen_validation::RequireVectorSize(parameters, kParameterSize, value_name);
    RequireFiniteParameterValues(parameters, value_name);
}

GaussianModel3D GaussianModel3D::FromVector(const Eigen::VectorXd & parameters)
{
    RequireParameterVector(parameters);
    return GaussianModel3D{
        parameters(kAmplitudeIndex),
        parameters(kWidthIndex),
        parameters(kOffsetIndex)
    };
}

GaussianModel3D GaussianModel3D::FromVectorPrefix(const Eigen::VectorXd & parameters)
{
    if (parameters.rows() < kParameterSize)
    {
        throw std::invalid_argument(
            "GaussianModel3D parameter vector must have at least three entries.");
    }
    RequireFiniteParameterValues(parameters, "GaussianModel3D parameter vector prefix");
    return GaussianModel3D{
        parameters(kAmplitudeIndex),
        parameters(kWidthIndex),
        parameters(kOffsetIndex)
    };
}

std::optional<GaussianModel3D> GaussianModel3D::FromTransformedCoordinates(
    const TransformedCoordinates & coordinates)
{
    if (!coordinates.allFinite()) return std::nullopt;

    const auto log_peak_height{
        coordinates(static_cast<Eigen::Index>(kLogPeakHeightCoordinateIndex))
    };
    const auto log_width{
        coordinates(static_cast<Eigen::Index>(kLogWidthCoordinateIndex))
    };
    const auto offset_to_peak_ratio{
        coordinates(static_cast<Eigen::Index>(kOffsetToPeakRatioCoordinateIndex))
    };
    const auto log_amplitude{
        log_peak_height + 1.5 * std::log(Constants::two_pi) + 3.0 * log_width
    };
    const auto amplitude{ std::exp(log_amplitude) };
    const auto width{ std::exp(log_width) };
    double offset{ 0.0 };
    if (offset_to_peak_ratio != 0.0)
    {
        const auto log_abs_offset{
            std::log(std::abs(offset_to_peak_ratio)) +
            log_peak_height + log_width - 0.5 * std::log(4.0 / Constants::two_pi)
        };
        offset = std::copysign(std::exp(log_abs_offset), offset_to_peak_ratio);
    }

    if (!std::isfinite(amplitude) || amplitude <= 0.0 ||
        !std::isfinite(width) || width <= 0.0 ||
        !std::isfinite(offset))
    {
        return std::nullopt;
    }
    return GaussianModel3D{ amplitude, width, offset };
}

void GaussianModel3D::RequireFiniteModel(
    const GaussianModel3D & model,
    std::string_view value_name)
{
    const auto name{ std::string(value_name) };
    numeric_validation::RequireFinite(model.GetAmplitude(), name + " amplitude");
    numeric_validation::RequireFinite(model.GetWidth(), name + " width");
    numeric_validation::RequireFinite(model.GetOffset(), name + " offset");
}

void GaussianModel3D::RequireFinitePositiveWidthModel(
    const GaussianModel3D & model,
    std::string_view value_name)
{
    const auto name{ std::string(value_name) };
    numeric_validation::RequireFinite(model.GetAmplitude(), name + " amplitude");
    numeric_validation::RequireFinitePositive(model.GetWidth(), name + " width");
    numeric_validation::RequireFinite(model.GetOffset(), name + " offset");
}

GaussianModel3D GaussianModel3D::WithAmplitude(double value) const
{
    return GaussianModel3D{ value, m_width, m_offset };
}

GaussianModel3D GaussianModel3D::WithWidth(double value) const
{
    return GaussianModel3D{ m_amplitude, value, m_offset };
}

GaussianModel3D GaussianModel3D::WithOffset(double value) const
{
    return GaussianModel3D{ m_amplitude, m_width, value };
}

Eigen::VectorXd GaussianModel3D::ToVector() const
{
    Eigen::VectorXd parameters{ Eigen::VectorXd::Zero(kParameterSize) };
    parameters(kAmplitudeIndex) = m_amplitude;
    parameters(kWidthIndex) = m_width;
    parameters(kOffsetIndex) = m_offset;
    return parameters;
}

std::optional<GaussianModel3D::TransformedCoordinates>
GaussianModel3D::ToTransformedCoordinates() const
{
    if (!std::isfinite(m_amplitude) || m_amplitude <= 0.0 ||
        !std::isfinite(m_width) || m_width <= 0.0 ||
        !std::isfinite(m_offset))
    {
        return std::nullopt;
    }

    const auto log_width{ std::log(m_width) };
    const auto log_peak_height{
        std::log(m_amplitude) - 1.5 * std::log(Constants::two_pi) - 3.0 * log_width
    };
    double offset_to_peak_ratio{ 0.0 };
    if (m_offset != 0.0)
    {
        const auto log_abs_offset_to_peak_ratio{
            std::log(std::abs(m_offset)) +
            0.5 * std::log(4.0 / Constants::two_pi) - log_width - log_peak_height
        };
        if (log_abs_offset_to_peak_ratio > std::log(std::numeric_limits<double>::max()))
        {
            return std::nullopt;
        }
        offset_to_peak_ratio = std::copysign(
            std::exp(log_abs_offset_to_peak_ratio),
            m_offset);
    }

    if (!std::isfinite(log_peak_height) ||
        !std::isfinite(log_width) ||
        !std::isfinite(offset_to_peak_ratio))
    {
        return std::nullopt;
    }

    return TransformedCoordinates{ log_peak_height, log_width, offset_to_peak_ratio };
}

double GaussianModel3D::GetModelParameter(int par_id) const
{
    switch (par_id)
    {
    case kAmplitudeIndex:
        return m_amplitude;
    case kWidthIndex:
        return m_width;
    case kOffsetIndex:
        return m_offset;
    default:
        throw std::out_of_range("GaussianModel3D parameter index is out of range.");
    }
}

double GaussianModel3D::GetDisplayParameter(int par_id) const
{
    switch (par_id)
    {
    case kAmplitudeIndex:
        return m_amplitude;
    case kWidthIndex:
        return m_width;
    case kOffsetIndex:
        return Intensity();
    default:
        throw std::out_of_range("GaussianModel3D display parameter index is out of range.");
    }
}

double GaussianModel3D::GetHeight() const
{
    return m_amplitude * std::pow(Constants::two_pi * m_width * m_width, -1.5);
}

double GaussianModel3D::Intensity() const
{
    return WithOffset(0.0).ResponseAtDistance(0.0);
}

double GaussianModel3D::SignalAtDistance(double distance) const
{
    return EvaluateAtDistance(distance).signal;
}

double GaussianModel3D::OffsetBasisAtDistance(double distance) const
{
    return EvaluateAtDistance(distance).offset_basis;
}

GaussianModel3DEvaluation GaussianModel3D::EvaluateAtDistance(double distance) const
{
    double offset_basis{ 0.0 };
    if (distance < 1.0e-5)
    {
        offset_basis = std::sqrt(2.0/M_PI) / m_width;
    }
    else
    {
        offset_basis = std::erf(distance/m_width/std::sqrt(2.0)) / distance;
    }

    if (m_width == 0.0)
    {
        return GaussianModel3DEvaluation{ 0.0, offset_basis, m_offset };
    }
    const auto width_square{ m_width * m_width };
    const auto signal{
        m_amplitude *
        std::pow(Constants::two_pi * width_square, -1.5) *
        std::exp(-0.5 * distance * distance / width_square)
    };
    return GaussianModel3DEvaluation{
        signal,
        offset_basis,
        signal + m_offset * offset_basis
    };
}

double GaussianModel3D::ResponseAtDistance(double distance) const
{
    return EvaluateAtDistance(distance).response;
}

GaussianModel3DUncertainty::GaussianModel3DUncertainty(
    double amplitude,
    double width,
    double offset) :
    m_amplitude{ amplitude },
    m_width{ width },
    m_offset{ offset }
{
}

void GaussianModel3DUncertainty::RequireFiniteNonNegativeUncertainty(
    const GaussianModel3DUncertainty & uncertainty,
    std::string_view value_name)
{
    const auto name{ std::string(value_name) };
    numeric_validation::RequireFiniteNonNegative(
        uncertainty.GetAmplitude(),
        name + " amplitude");
    numeric_validation::RequireFiniteNonNegative(
        uncertainty.GetWidth(),
        name + " width");
    numeric_validation::RequireFiniteNonNegative(
        uncertainty.GetOffset(),
        name + " offset");
}

Eigen::VectorXd GaussianModel3DUncertainty::ToVector() const
{
    Eigen::VectorXd parameters{ Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()) };
    parameters(GaussianModel3D::AmplitudeIndex()) = m_amplitude;
    parameters(GaussianModel3D::WidthIndex()) = m_width;
    parameters(GaussianModel3D::OffsetIndex()) = m_offset;
    return parameters;
}

double GaussianModel3DUncertainty::GetModelParameter(int par_id) const
{
    switch (par_id)
    {
    case GaussianModel3D::AmplitudeIndex():
        return m_amplitude;
    case GaussianModel3D::WidthIndex():
        return m_width;
    case GaussianModel3D::OffsetIndex():
        return m_offset;
    default:
        throw std::out_of_range(
            "GaussianModel3DUncertainty parameter index is out of range.");
    }
}

GaussianModel3DWithUncertainty::GaussianModel3DWithUncertainty(
    GaussianModel3D model,
    GaussianModel3DUncertainty standard_deviation) :
    m_model{ model },
    m_standard_deviation{ standard_deviation }
{
}

double GaussianModel3DWithUncertainty::GetDisplayParameter(int par_id) const
{
    return m_model.GetDisplayParameter(par_id);
}

double GaussianModel3DWithUncertainty::GetModelParameter(int par_id) const
{
    return m_model.GetModelParameter(par_id);
}

double GaussianModel3DWithUncertainty::GetModelStandardDeviation(int par_id) const
{
    return m_standard_deviation.GetModelParameter(par_id);
}

double GaussianModel3DWithUncertainty::GetDisplayStandardDeviation(int par_id) const
{
    switch (par_id)
    {
    case GaussianModel3D::AmplitudeIndex():
        return m_standard_deviation.GetAmplitude();
    case GaussianModel3D::WidthIndex():
        return m_standard_deviation.GetWidth();
    case GaussianModel3D::OffsetIndex():
        return IntensityStandardDeviation();
    default:
        throw std::out_of_range(
            "GaussianModel3DWithUncertainty display standard deviation index is out of range.");
    }
}

double GaussianModel3DWithUncertainty::IntensityStandardDeviation() const
{
    const auto sigma_amplitude{ m_standard_deviation.GetAmplitude() };
    const auto sigma_width{ m_standard_deviation.GetWidth() };
    const auto amplitude{ m_model.GetAmplitude() };
    const auto width{ m_model.GetWidth() };
    return std::sqrt(
        std::pow(std::pow(Constants::two_pi * width * width, -1.5) * sigma_amplitude, 2) +
        std::pow(
            -3.0 * amplitude * std::pow(Constants::two_pi, -1.5)
                * std::pow(width, -4) * sigma_width,
            2));
}

} // namespace rhbm_gem
