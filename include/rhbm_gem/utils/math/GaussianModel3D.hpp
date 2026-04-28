#pragma once

#include <string_view>

#include <Eigen/Dense>

namespace rhbm_gem
{

// GaussianParameterVector stays in Gaussian-model space; its dimension depends on the model.
using GaussianParameterVector = Eigen::VectorXd;

class GaussianModel3D
{
    double m_amplitude{ 0.0 };
    double m_width{ 1.0 };
    double m_intercept{ 0.0 };

    static constexpr int kParameterSize{ 3 };
    static constexpr int kAmplitudeIndex{ 0 };
    static constexpr int kWidthIndex{ 1 };
    static constexpr int kInterceptIndex{ 2 };

public:
    GaussianModel3D() = default;
    GaussianModel3D(double amplitude, double width, double intercept);

    static constexpr int ParameterSize() { return kParameterSize; }
    static constexpr int AmplitudeIndex() { return kAmplitudeIndex; }
    static constexpr int WidthIndex() { return kWidthIndex; }
    static constexpr int InterceptIndex() { return kInterceptIndex; }

    static void RequireParameterVector(
        const GaussianParameterVector & parameters,
        std::string_view value_name = "GaussianModel3D parameter vector");

    static GaussianModel3D FromVector(const GaussianParameterVector & parameters);
    static GaussianModel3D FromVectorPrefix(const GaussianParameterVector & parameters);

    static void RequireFiniteModel(
        const GaussianModel3D & model,
        std::string_view value_name = "GaussianModel3D");
    static void RequireFinitePositiveWidthModel(
        const GaussianModel3D & model,
        std::string_view value_name = "GaussianModel3D");

    GaussianModel3D WithAmplitude(double value) const;
    GaussianModel3D WithWidth(double value) const;
    GaussianModel3D WithIntercept(double value) const;

    double GetAmplitude() const;
    double GetWidth() const;
    double GetIntercept() const;
    GaussianParameterVector ToVector() const;
    double GetModelParameter(int par_id) const;
    double Intensity() const;
    double ResponseAtDistance(double distance) const;
};

} // namespace rhbm_gem
