#pragma once

#include <string_view>

#include <Eigen/Dense>

namespace rhbm_gem
{

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
    GaussianModel3D(double amplitude, double width, double intercept = 0.0);

    static constexpr int ParameterSize() { return kParameterSize; }
    static constexpr int AmplitudeIndex() { return kAmplitudeIndex; }
    static constexpr int WidthIndex() { return kWidthIndex; }
    static constexpr int InterceptIndex() { return kInterceptIndex; }

    static void RequireParameterVector(
        const Eigen::VectorXd & parameters,
        std::string_view value_name = "GaussianModel3D parameter vector");

    static GaussianModel3D FromVector(const Eigen::VectorXd & parameters);
    static GaussianModel3D FromVectorPrefix(const Eigen::VectorXd & parameters);

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
    Eigen::VectorXd ToVector() const;
    double GetModelParameter(int par_id) const;
    double GetDisplayParameter(int par_id) const;
    double Intensity() const;
    double ResponseAtDistance(double distance) const;
};

class GaussianModel3DUncertainty
{
    double m_amplitude{ 0.0 };
    double m_width{ 0.0 };
    double m_intercept{ 0.0 };

public:
    GaussianModel3DUncertainty() = default;
    GaussianModel3DUncertainty(double amplitude, double width, double intercept = 0.0);

    static void RequireFiniteNonNegativeUncertainty(
        const GaussianModel3DUncertainty & uncertainty,
        std::string_view value_name = "GaussianModel3DUncertainty");

    double GetAmplitude() const;
    double GetWidth() const;
    double GetIntercept() const;
    Eigen::VectorXd ToVector() const;
    double GetModelParameter(int par_id) const;
};

class GaussianModel3DWithUncertainty
{
    GaussianModel3D m_model{ 0.0, 0.0 };
    GaussianModel3DUncertainty m_standard_deviation{};

public:
    GaussianModel3DWithUncertainty() = default;
    GaussianModel3DWithUncertainty(
        GaussianModel3D model,
        GaussianModel3DUncertainty standard_deviation);

    const GaussianModel3D & GetModel() const;
    const GaussianModel3DUncertainty & GetStandardDeviationModel() const;
    double GetDisplayParameter(int par_id) const;
    double GetDisplayStandardDeviation(int par_id) const;
    double IntensityStandardDeviation() const;
};

} // namespace rhbm_gem
