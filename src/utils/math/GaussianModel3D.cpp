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
    const GaussianParameterVector & parameters,
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
    const GaussianParameterVector & parameters,
    std::string_view value_name)
{
    eigen_validation::RequireVectorSize(parameters, kParameterSize, value_name);
    RequireFiniteParameterValues(parameters, value_name);
}

GaussianModel3D GaussianModel3D::FromVector(const GaussianParameterVector & parameters)
{
    RequireParameterVector(parameters);
    return GaussianModel3D{
        parameters(kAmplitudeIndex),
        parameters(kWidthIndex),
        parameters(kInterceptIndex)
    };
}

GaussianModel3D GaussianModel3D::FromVectorPrefix(const GaussianParameterVector & parameters)
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

GaussianParameterVector GaussianModel3D::ToVector() const
{
    GaussianParameterVector parameters{ GaussianParameterVector::Zero(kParameterSize) };
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

} // namespace rhbm_gem
