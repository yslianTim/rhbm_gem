#pragma once

#include "support/MDPDEExperiment.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include <memory>

namespace second_stage_test {

enum class EndpointPolicy { Legacy, FailedOnly, FreshResidual };
const char * EndpointPolicyName(EndpointPolicy);
bool ShouldRefineEndpoint(EndpointPolicy, rhbm_gem::RHBMEstimationStatus, bool fresh_pass);

// A process-scoped immutable policy is installed before workers start and removed
// after they join. Unlike a main-thread TLS switch, every OpenMP worker sees it.
class ScopedSecondStageEndpointExperiment
{
    struct Impl;
    std::unique_ptr<Impl> impl;
public:
    ScopedSecondStageEndpointExperiment();
    ~ScopedSecondStageEndpointExperiment();
};

bool IsEndpointOperatorProbe() noexcept;
bool IsEndpointExperimentActive() noexcept;
rhbm_gem::RHBMBetaEstimateResult ApplyEndpointExperiment(const ShapeFixture &);
void CaptureEndpointInitialState(const rhbm_gem::core::detail::SecondStageContext &,
    const rhbm_gem::core::detail::FitState &, const rhbm_gem::core::FitOptions &);
void CompareEndpointOperators(const char * label,
    const rhbm_gem::core::detail::SecondStageContext &,
    const std::vector<rhbm_gem::core::detail::ClusterKey> &,
    const rhbm_gem::core::detail::FitState &, const rhbm_gem::core::FitOptions &,
    const std::vector<double> & ridge,
    const rhbm_gem::core::detail::FixedPointOperatorEvidence & baseline);

} // namespace second_stage_test
