#pragma once

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

class HRLGroupEstimator final
{
    RHBMExecutionOptions m_options;

public:
    explicit HRLGroupEstimator(RHBMExecutionOptions options = {});
    RHBMGroupEstimationResult Estimate(const RHBMGroupEstimationInput & input, double alpha_g) const;

private:
    static RHBMGroupEstimationResult BuildFallbackResult(
        const RHBMGroupEstimationInput & input,
        const RHBMMuEstimateResult & mu_result
    );
};
