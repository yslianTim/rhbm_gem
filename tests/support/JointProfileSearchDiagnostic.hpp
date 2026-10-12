#pragma once
#include "core/detail/joint_component/Numerics.hpp"

namespace rhbm_gem::core::joint_component {
struct ProfileSearchDiagnosticResult
{
    SearchResult search;
    Assessment assessment;
    std::optional<Assessment> trusted_assessment;
    std::optional<Endpoint> trusted_state;
    std::optional<std::size_t> trusted_trial;
    double assessment_seconds{};
    std::optional<TrustEvidence> endpoint_trust;
    bool search_success{};
};

ProfileSearchDiagnosticResult AssessProfileSearchForDiagnostics(const Domain &,VectorRef,const EvaluationContext &,SearchResult,
    const JointProgressObserver & = {},const JointProgressComponent * = nullptr);
}
