#pragma once

#include <Eigen/Dense>

namespace HRLModelMetrics
{
Eigen::VectorXd CalculateNormalizedResidual(
    const Eigen::VectorXd & linear_estimate,
    const Eigen::VectorXd & gaus_truth
);

Eigen::VectorXd CalculateAbsoluteGausDifference(
    const Eigen::VectorXd & linear_a,
    const Eigen::VectorXd & linear_b
);
} // namespace HRLModelMetrics
