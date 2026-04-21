#pragma once

#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>

namespace HRLModelMetrics
{

inline rhbm_gem::GaussianParameterVector CalculateNormalizedResidual(
    const HRLBetaVector & linear_estimate,
    const rhbm_gem::GaussianParameterVector & gaus_truth)
{
    static const rhbm_gem::GaussianLinearizationService service{
        rhbm_gem::GaussianLinearizationSpec::DefaultMetricModel()
    };
    const auto gaussian_estimate{ service.DecodeGroupBeta(linear_estimate) };
    rhbm_gem::EigenValidation::RequireVectorSize(gaussian_estimate, gaus_truth.rows(), "gaussian");
    return ((gaussian_estimate - gaus_truth).array() / gaus_truth.array()).matrix();
}

inline rhbm_gem::GaussianParameterVector CalculateAbsoluteGausDifference(
    const HRLBetaVector & linear_a,
    const HRLBetaVector & linear_b)
{
    static const rhbm_gem::GaussianLinearizationService service{
        rhbm_gem::GaussianLinearizationSpec::DefaultMetricModel()
    };
    const auto gaussian_a{ service.DecodeGroupBeta(linear_a) };
    const auto gaussian_b{ service.DecodeGroupBeta(linear_b) };
    rhbm_gem::EigenValidation::RequireVectorSize(gaussian_a, gaussian_b.rows(), "gaussian");
    return (gaussian_a - gaussian_b).array().abs().matrix();
}
} // namespace HRLModelMetrics
