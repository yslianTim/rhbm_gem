#pragma once
#include "Problem.hpp"
namespace rhbm_gem::core::joint_component {
bool ValidJointPointParameters(double amplitude,double tau,double charge,double log_tau);

JointTargetEvidence AssessTargetPointState(const JointProblem &,const JointParameterLayout &,
    const std::optional<JointState> &,JointCheckStatus,const std::vector<JointBlockConditioning> &,
    JointEvidenceScope,std::size_t);
JointCheckStatus TargetConvergence(const std::optional<JointTargetEvidence> &,bool);
void AddTargetEvidence(JointFitResult &);
}
