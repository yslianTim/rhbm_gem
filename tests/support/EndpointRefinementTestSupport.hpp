#pragma once

#include "support/MDPDETestSupport.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include <memory>

namespace second_stage_test {

enum class EndpointPolicy { Legacy, FailedOnly, FreshResidual };
const char * EndpointPolicyName(EndpointPolicy);
bool ShouldRefineEndpoint(EndpointPolicy, rhbm_gem::RHBMEstimationStatus, bool fresh_pass);

struct EndpointRefinementResult
{
    rhbm_gem::RHBMBetaEstimateResult result;
    bool accepted{};
    boost::json::object evidence;
};
EndpointRefinementResult RefineMDPDEEndpoint(const ShapeFixture &, int equation_budget,
    const MDPDEEquationEvidence * initial = nullptr);

// A process-scoped immutable policy is installed before workers start and removed
// after they join. Unlike a main-thread TLS switch, every OpenMP worker sees it.
class ScopedSecondStageEndpointTest
{
    struct Impl;
    std::unique_ptr<Impl> impl;
public:
    ScopedSecondStageEndpointTest();
    ~ScopedSecondStageEndpointTest();
};

bool IsEndpointOperatorProbe() noexcept;
bool IsEndpointTestActive() noexcept;
rhbm_gem::RHBMBetaEstimateResult ApplyEndpointTestPolicy(const ShapeFixture &);
void CaptureEndpointInitialState(const rhbm_gem::core::detail::SecondStageContext &,
    const rhbm_gem::core::detail::FitState &, const rhbm_gem::core::FitOptions &);
void CompareEndpointOperators(const char * label,
    const rhbm_gem::core::detail::SecondStageContext &,
    const std::vector<rhbm_gem::core::detail::ClusterKey> &,
    const rhbm_gem::core::detail::FitState &, const rhbm_gem::core::FitOptions &,
    const std::vector<double> & ridge,
    const rhbm_gem::core::detail::FixedPointOperatorEvidence & baseline);

} // namespace second_stage_test
