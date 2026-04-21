#pragma once

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

class HRLGroupEstimator final
{
    HRLExecutionOptions m_options;

public:
    explicit HRLGroupEstimator(HRLExecutionOptions options = {});
    HRLGroupEstimationResult Estimate(const HRLGroupEstimationInput & input, double alpha_g) const;

private:
    static HRLGroupEstimationResult BuildFallbackResult(
        const HRLGroupEstimationInput & input,
        const HRLMuEstimateResult & mu_result
    );
};
