#pragma once

#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

namespace rhbm_gem {

// GaussianParameterVector stays in Gaussian-model space; its dimension depends on the model.
// HRLBetaVector stays in HRL linearized coefficient space.
using GaussianParameterVector = Eigen::VectorXd;

// Shared Gaussian value types used by HRL linearization and downstream consumers.
struct GaussianEstimate
{
    double amplitude{ 0.0 };
    double width{ 0.0 };

    double GetParameter(int par_id) const
    {
        switch (par_id)
        {
        case 0:
            return amplitude;
        case 1:
            return width;
        case 2:
            return Intensity();
        default:
            throw std::out_of_range("GaussianEstimate parameter index is out of range.");
        }
    }

    double Intensity() const
    {
        if (width == 0.0)
        {
            return 0.0;
        }
        return amplitude * std::pow(Constants::two_pi * width * width, -1.5);
    }

    HRLBetaVector ToBeta() const
    {
        HRLBetaVector coefficient_vector{ HRLBetaVector::Zero(3) };
        if (width == 0.0)
        {
            return coefficient_vector;
        }
        const auto width_square{ width * width };
        coefficient_vector(0) = (amplitude <= 0.0) ?
            0.0 : std::log(amplitude) - 1.5 * std::log(Constants::two_pi * width_square);
        coefficient_vector(1) = 1.0 / width_square;
        return coefficient_vector;
    }
};

struct GaussianPosterior
{
    GaussianEstimate estimate{};
    GaussianEstimate variance{};

    double GetEstimate(int par_id) const
    {
        return estimate.GetParameter(par_id);
    }

    double GetVariance(int par_id) const
    {
        switch (par_id)
        {
        case 0:
            return variance.amplitude;
        case 1:
            return variance.width;
        case 2:
            return IntensityVariance();
        default:
            throw std::out_of_range("GaussianPosterior variance index is out of range.");
        }
    }

    double IntensityVariance() const
    {
        const auto sigma_amplitude{ variance.amplitude };
        const auto sigma_width{ variance.width };
        const auto amplitude{ estimate.amplitude };
        const auto width{ estimate.width };
        return std::sqrt(
            std::pow(std::pow(Constants::two_pi * width * width, -1.5) * sigma_amplitude, 2) +
            std::pow(
                -3.0 * amplitude * std::pow(Constants::two_pi, -1.5)
                    * std::pow(width, -4) * sigma_width,
                2));
    }
};

} // namespace rhbm_gem
