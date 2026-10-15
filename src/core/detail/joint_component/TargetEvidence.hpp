#pragma once
#include "Problem.hpp"
namespace rhbm_gem::core::joint_component {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct TargetEvidenceWork
{
    std::size_t memberships_indexed{},atom_lookups{};
};
TargetEvidenceWork & TargetEvidenceWorkForTesting();
#endif
JointTargetEvidence AssessTargetPointState(const JointProblem &,const JointParameterLayout &,
    const std::optional<JointState> &,JointCheckStatus,const std::vector<JointBlockConditioning> &,
    const std::vector<JointBlockParameterStability> &,JointEvidenceScope,std::size_t);
JointCheckStatus TargetConvergence(const std::optional<JointTargetEvidence> &,bool);
void AddTargetEvidence(JointFitResult &);
}
