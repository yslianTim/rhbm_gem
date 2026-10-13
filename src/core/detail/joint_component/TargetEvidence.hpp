#pragma once
#include "Problem.hpp"
namespace rhbm_gem::core::joint_component {
MemoryEstimate TargetIdentifiabilityMemoryEstimate(const JointProblemInput &,const JointParameterLayout &,std::size_t);
std::string TargetRankContextIdentity(const JointProblemInput &,const JointParameterLayout &,const JointState &,std::size_t);
JointTargetEvidence DeriveTargetIdentifiability(const JointProblem &,const JointParameterLayout &,
    const std::optional<JointState> &,const std::vector<JointCheck> &,const std::vector<JointRankEvidence> &,
    JointEvidenceScope,std::size_t);
JointTargetEvidence AssessTargetState(const JointProblem &,const JointParameterLayout &,
    const std::optional<JointState> &,const std::vector<JointCheck> &,const std::vector<JointRankEvidence> &,
    JointEvidenceScope,std::size_t);
JointCheckStatus TargetConvergence(const std::optional<JointTargetEvidence> &,bool);
void AddTargetEvidence(JointFitResult &);
}
