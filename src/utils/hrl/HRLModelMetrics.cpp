#include <rhbm_gem/utils/hrl/HRLModelMetrics.hpp>

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>

#include <stdexcept>
#include <string>

namespace
{
Eigen::VectorXd BuildGausModel(const Eigen::VectorXd & linear_model)
{
    static const rhbm_gem::GaussianLinearizationService service{
        rhbm_gem::GaussianLinearizationSpec::DefaultMetricModel()
    };
    return service.DecodeGroupBeta(linear_model);
}

void ValidateMatchingDimensions(
    const Eigen::VectorXd & estimate_gaus,
    const Eigen::VectorXd & gaus_truth)
{
    if (estimate_gaus.rows() != gaus_truth.rows())
    {
        Logger::Log(
            LogLevel::Error,
            "estimate size " + std::to_string(estimate_gaus.rows()) +
                " != model size " + std::to_string(gaus_truth.rows())
        );
        throw std::invalid_argument("model parameters size inconsistant.");
    }
}
} // namespace

namespace HRLModelMetrics
{
Eigen::VectorXd CalculateNormalizedResidual(
    const Eigen::VectorXd & linear_estimate,
    const Eigen::VectorXd & gaus_truth)
{
    const auto estimate_gaus{ BuildGausModel(linear_estimate) };
    ValidateMatchingDimensions(estimate_gaus, gaus_truth);
    return ((estimate_gaus - gaus_truth).array() / gaus_truth.array()).matrix();
}

Eigen::VectorXd CalculateAbsoluteGausDifference(
    const Eigen::VectorXd & linear_a,
    const Eigen::VectorXd & linear_b)
{
    const auto gaus_a{ BuildGausModel(linear_a) };
    const auto gaus_b{ BuildGausModel(linear_b) };
    ValidateMatchingDimensions(gaus_a, gaus_b);
    return (gaus_a - gaus_b).array().abs().matrix();
}
} // namespace HRLModelMetrics
