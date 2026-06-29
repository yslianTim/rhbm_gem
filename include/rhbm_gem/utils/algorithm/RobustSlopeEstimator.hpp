#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <rhbm_gem/utils/algorithm/LinearRegressionSample.hpp>
#include <rhbm_gem/utils/algorithm/RobustSlopeOptions.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::algorithm {

class RobustSlopeEstimator
{
public:
    static bool EstimateOrdinarySlopeThroughOrigin(
        const std::vector<LinearRegressionSample> & sample_list,
        double & slope)
    {
        double numerator{ 0.0 };
        double denominator{ 0.0 };
        for (const auto & sample : sample_list)
        {
            numerator += sample.basis * sample.response;
            denominator += sample.basis * sample.basis;
        }
        if (!std::isfinite(numerator) ||
            !std::isfinite(denominator) ||
            denominator <= std::numeric_limits<double>::epsilon())
        {
            return false;
        }
        slope = numerator / denominator;
        return std::isfinite(slope);
    }

    static bool EstimateHuberSlopeThroughOrigin(
        const std::vector<LinearRegressionSample> & sample_list,
        const RobustSlopeOptions & options,
        double & slope)
    {
        if (sample_list.empty())
        {
            return false;
        }
        if (!EstimateOrdinarySlopeThroughOrigin(sample_list, slope))
        {
            return false;
        }

        for (int t = 0; t < options.maximum_iterations; t++)
        {
            const auto scale{ ComputeHuberResidualScale(sample_list, slope, options) };
            const auto cutoff{ options.cutoff_multiplier * scale };
            const auto lambda{ ComputeRegularizationLambda(scale, options) };
            double numerator{ 0.0 };
            double denominator{ 0.0 };
            for (const auto & sample : sample_list)
            {
                const auto error{ sample.response - slope * sample.basis };
                const auto abs_error{ std::abs(error) };
                const auto weight{ abs_error <= cutoff ? 1.0 : cutoff / abs_error };
                numerator += weight * sample.basis * sample.response;
                denominator += weight * sample.basis * sample.basis;
            }
            denominator += lambda;
            if (!std::isfinite(numerator) ||
                !std::isfinite(denominator) ||
                denominator <= std::numeric_limits<double>::epsilon())
            {
                return false;
            }

            const auto updated_slope{ numerator / denominator };
            if (!std::isfinite(updated_slope))
            {
                return false;
            }
            if (std::abs(updated_slope - slope) < options.tolerance)
            {
                slope = updated_slope;
                return true;
            }
            slope = updated_slope;
        }
        return true;
    }

private:
    static double ComputeHuberResidualScale(
        const std::vector<LinearRegressionSample> & sample_list,
        double slope,
        const RobustSlopeOptions & options)
    {
        std::vector<double> residual_list;
        residual_list.reserve(sample_list.size());
        for (const auto & sample : sample_list)
        {
            residual_list.emplace_back(sample.response - slope * sample.basis);
        }

        const auto median_residual{ array_helper::ComputeMedian(residual_list) };
        std::vector<double> deviation_list;
        deviation_list.reserve(residual_list.size());
        for (const auto residual : residual_list)
        {
            deviation_list.emplace_back(std::abs(residual - median_residual));
        }

        return std::max(
            options.scale_multiplier * array_helper::ComputeMedian(deviation_list),
            options.scale_min);
    }

    static double ComputeRegularizationLambda(
        double residual_scale,
        const RobustSlopeOptions & options)
    {
        if (!std::isfinite(options.regularization_prior_scale) ||
            options.regularization_prior_scale <= 0.0)
        {
            return 0.0;
        }
        const auto lambda{ residual_scale / options.regularization_prior_scale };
        return lambda * lambda;
    }
};

} // namespace rhbm_gem::algorithm
