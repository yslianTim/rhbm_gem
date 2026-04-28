#pragma once

#include <string_view>

#include <Eigen/Dense>

namespace rhbm_gem
{

// GaussianParameterVector stays in Gaussian-model space; its dimension depends on the model.
using GaussianParameterVector = Eigen::VectorXd;

class GaussianModel3D
{
public:
    static constexpr int kParameterSize{ 3 };
    static constexpr int kAmplitudeIndex{ 0 };
    static constexpr int kWidthIndex{ 1 };
    static constexpr int kInterceptIndex{ 2 };

    double amplitude{ 0.0 };
    double width{ 1.0 };
    double intercept{ 0.0 };

    GaussianModel3D() = default;
    GaussianModel3D(double amplitude, double width, double intercept);

    static void RequireParameterVector(
        const GaussianParameterVector & parameters,
        std::string_view value_name = "GaussianModel3D parameter vector");

    static GaussianModel3D FromVector(const GaussianParameterVector & parameters);
    static GaussianModel3D FromVectorPrefix(const GaussianParameterVector & parameters);

    GaussianParameterVector ToVector() const;
    double GetModelParameter(int par_id) const;
    double Intensity() const;
    double ResponseAtDistance(double distance) const;
};

} // namespace rhbm_gem
