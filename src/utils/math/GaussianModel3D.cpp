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
        parameters(GaussianModel3D::kAmplitudeIndex),
        name + " amplitude");
    numeric_validation::RequireFinite(
        parameters(GaussianModel3D::kWidthIndex),
        name + " width");
    numeric_validation::RequireFinite(
        parameters(GaussianModel3D::kInterceptIndex),
        name + " intercept");
}

} // namespace

GaussianModel3D::GaussianModel3D(double amplitude, double width, double intercept) :
    amplitude{ amplitude },
    width{ width },
    intercept{ intercept }
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

GaussianParameterVector GaussianModel3D::ToVector() const
{
    GaussianParameterVector parameters{ GaussianParameterVector::Zero(kParameterSize) };
    parameters(kAmplitudeIndex) = amplitude;
    parameters(kWidthIndex) = width;
    parameters(kInterceptIndex) = intercept;
    return parameters;
}

double GaussianModel3D::GetModelParameter(int par_id) const
{
    switch (par_id)
    {
    case kAmplitudeIndex:
        return amplitude;
    case kWidthIndex:
        return width;
    case kInterceptIndex:
        return intercept;
    default:
        throw std::out_of_range("GaussianModel3D parameter index is out of range.");
    }
}

double GaussianModel3D::Intensity() const
{
    if (width == 0.0)
    {
        return 0.0;
    }
    return amplitude * std::pow(Constants::two_pi * width * width, -1.5);
}

double GaussianModel3D::ResponseAtDistance(double distance) const
{
    if (width == 0.0)
    {
        return intercept;
    }
    return Intensity() * std::exp(-0.5 * distance * distance / (width * width))
        + intercept;
}

} // namespace rhbm_gem
