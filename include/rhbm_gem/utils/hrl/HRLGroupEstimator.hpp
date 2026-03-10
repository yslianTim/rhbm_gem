#pragma once

#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>

class HRLGroupEstimator final
{
public:
    explicit HRLGroupEstimator(HRLExecutionOptions options = {});

    HRLGroupEstimationResult Estimate(
        const HRLGroupEstimationInput & input,
        double alpha_g
    ) const;

private:
    HRLExecutionOptions m_options;

    static void ValidateInput(const HRLGroupEstimationInput & input);
    static HRLGroupEstimationResult BuildFallbackResult(
        const HRLGroupEstimationInput & input,
        const HRLMuEstimateResult & mu_result
    );
};
