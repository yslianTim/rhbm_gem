#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace rhbm_gem
{

namespace
{

void RequireFiniteParameterValues(
    const Eigen::VectorXd & parameters,
    std::string_view value_name)
{
    const auto name{ std::string(value_name) };
    numeric_validation::RequireFinite(
        parameters(GaussianModel3D::AmplitudeIndex()),
        name + " amplitude");
    numeric_validation::RequireFinite(
        parameters(GaussianModel3D::WidthIndex()),
        name + " width");
    numeric_validation::RequireFinite(
        parameters(GaussianModel3D::InterceptIndex()),
        name + " intercept");
}

} // namespace

GaussianModel3D::GaussianModel3D(double amplitude, double width, double intercept) :
    m_amplitude{ amplitude },
    m_width{ width },
    m_intercept{ intercept }
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
        parameters(kInterceptIndex)
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
        parameters(kInterceptIndex)
    };
}

void GaussianModel3D::RequireFiniteModel(
    const GaussianModel3D & model,
    std::string_view value_name)
{
    const auto name{ std::string(value_name) };
    numeric_validation::RequireFinite(model.GetAmplitude(), name + " amplitude");
    numeric_validation::RequireFinite(model.GetWidth(), name + " width");
    numeric_validation::RequireFinite(model.GetIntercept(), name + " intercept");
}

void GaussianModel3D::RequireFinitePositiveWidthModel(
    const GaussianModel3D & model,
    std::string_view value_name)
{
    const auto name{ std::string(value_name) };
    numeric_validation::RequireFinite(model.GetAmplitude(), name + " amplitude");
    numeric_validation::RequireFinitePositive(model.GetWidth(), name + " width");
    numeric_validation::RequireFinite(model.GetIntercept(), name + " intercept");
}

GaussianModel3D GaussianModel3D::WithAmplitude(double value) const
{
    return GaussianModel3D{ value, m_width, m_intercept };
}

GaussianModel3D GaussianModel3D::WithWidth(double value) const
{
    return GaussianModel3D{ m_amplitude, value, m_intercept };
}

GaussianModel3D GaussianModel3D::WithIntercept(double value) const
{
    return GaussianModel3D{ m_amplitude, m_width, value };
}

double GaussianModel3D::GetAmplitude() const
{
    return m_amplitude;
}

double GaussianModel3D::GetWidth() const
{
    return m_width;
}

double GaussianModel3D::GetIntercept() const
{
    return m_intercept;
}

Eigen::VectorXd GaussianModel3D::ToVector() const
{
    Eigen::VectorXd parameters{ Eigen::VectorXd::Zero(kParameterSize) };
    parameters(kAmplitudeIndex) = m_amplitude;
    parameters(kWidthIndex) = m_width;
    parameters(kInterceptIndex) = m_intercept;
    return parameters;
}

double GaussianModel3D::GetModelParameter(int par_id) const
{
    switch (par_id)
    {
    case kAmplitudeIndex:
        return m_amplitude;
    case kWidthIndex:
        return m_width;
    case kInterceptIndex:
        return m_intercept;
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
    case kInterceptIndex:
        return Intensity();
    default:
        throw std::out_of_range("GaussianModel3D display parameter index is out of range.");
    }
}

double GaussianModel3D::Intensity() const
{
    if (m_width == 0.0)
    {
        return 0.0;
    }
    return m_amplitude * std::pow(Constants::two_pi * m_width * m_width, -1.5);
}

double GaussianModel3D::ResponseAtDistance(double distance) const
{
    if (m_width == 0.0)
    {
        return m_intercept;
    }
    return Intensity() * std::exp(-0.5 * distance * distance / (m_width * m_width))
        + m_intercept;
}

GaussianModel3DUncertainty::GaussianModel3DUncertainty(
    double amplitude,
    double width,
    double intercept) :
    m_amplitude{ amplitude },
    m_width{ width },
    m_intercept{ intercept }
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
        uncertainty.GetIntercept(),
        name + " intercept");
}

double GaussianModel3DUncertainty::GetAmplitude() const
{
    return m_amplitude;
}

double GaussianModel3DUncertainty::GetWidth() const
{
    return m_width;
}

double GaussianModel3DUncertainty::GetIntercept() const
{
    return m_intercept;
}

Eigen::VectorXd GaussianModel3DUncertainty::ToVector() const
{
    Eigen::VectorXd parameters{ Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()) };
    parameters(GaussianModel3D::AmplitudeIndex()) = m_amplitude;
    parameters(GaussianModel3D::WidthIndex()) = m_width;
    parameters(GaussianModel3D::InterceptIndex()) = m_intercept;
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
    case GaussianModel3D::InterceptIndex():
        return m_intercept;
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

const GaussianModel3D & GaussianModel3DWithUncertainty::GetModel() const
{
    return m_model;
}

const GaussianModel3DUncertainty &
GaussianModel3DWithUncertainty::GetStandardDeviationModel() const
{
    return m_standard_deviation;
}

double GaussianModel3DWithUncertainty::GetDisplayParameter(int par_id) const
{
    return m_model.GetDisplayParameter(par_id);
}

double GaussianModel3DWithUncertainty::GetDisplayStandardDeviation(int par_id) const
{
    switch (par_id)
    {
    case GaussianModel3D::AmplitudeIndex():
        return m_standard_deviation.GetAmplitude();
    case GaussianModel3D::WidthIndex():
        return m_standard_deviation.GetWidth();
    case GaussianModel3D::InterceptIndex():
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
