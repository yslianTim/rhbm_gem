#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include <cmath>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
using rhbm_gem::core::JointProblem;
void CheckSingleBlock(const std::string & topology,int atoms)
{
    JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector eta=n::Vector::Constant(atoms,std::log(.55));
    n::FixedNeighborPolicy policy; policy.maximum_sweeps=4;
    std::vector<n::FixedNeighborBlockSweep> observed;
    n::BlockCoordinateState observed_state; std::size_t state_callbacks{};
    policy.sweep_observer=[&](const auto & sweep){observed.push_back(sweep);};
    policy.state_observer=[&](std::size_t,const auto & state,const auto &,const auto &) {
        observed_state=state; ++state_callbacks;
    };
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
    auto unobserved_policy=policy; unobserved_policy.sweep_observer={}; unobserved_policy.state_observer={};
    const auto unobserved=n::SearchFixedNeighbor(problem,eta,unobserved_policy);
    ASSERT_TRUE(result.search_converged)<<result.reason;
    ASSERT_TRUE(result.endpoint_certified)<<result.endpoint_trust.reason;
    EXPECT_EQ(result.sweeps.size(),1u);
    ASSERT_EQ(observed.size(),result.sweeps.size());
    EXPECT_EQ(state_callbacks,result.sweeps.size());
    EXPECT_TRUE((observed_state.eta.array()==result.state.eta.array()).all());
    EXPECT_TRUE((observed_state.beta.array()==result.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(observed.back().global_ac_kkt,result.sweeps.back().global_ac_kkt);
    EXPECT_LE(result.sweeps.back().global_ac_kkt,1e-10);
    EXPECT_LE(result.sweeps.back().global_width_gradient_inf_norm,1e-12);
    EXPECT_LE(result.sweeps.back().cache_replay_error,2e-12);
    EXPECT_TRUE(result.assessment.primary.valid);
    EXPECT_TRUE(result.assessment.reference.valid);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(result.fit.search_completed);
    EXPECT_EQ(result.search_converged,unobserved.search_converged);
    EXPECT_EQ(result.endpoint_certified,unobserved.endpoint_certified);
    EXPECT_EQ(result.reason,unobserved.reason);
    EXPECT_EQ(result.sweeps.size(),unobserved.sweeps.size());
    EXPECT_EQ(result.state.eta.size(),unobserved.state.eta.size());
    EXPECT_TRUE((result.state.eta.array()==unobserved.state.eta.array()).all());
    EXPECT_TRUE((result.state.beta.array()==unobserved.state.beta.array()).all());
    EXPECT_DOUBLE_EQ(result.state.objective,unobserved.state.objective);
    auto stopped_search=result.fit; stopped_search.search_completed=false;
    EXPECT_EQ(stopped_search.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    for(const auto & sweep:result.sweeps)
        EXPECT_LE(sweep.objective_after,sweep.objective_before+n::BlockObjectiveReplayEnclosure(sweep.objective_before));
    EXPECT_EQ(data.partition.components.size(),1u);
}
void CheckMultiBlockBudgetFailure(const std::string & topology)
{
    JointProblem problem(second_stage_test::OperatorWorkload(topology,32));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto partition=n::BuildStructuralBlockPartition(*data.input,data.layout,16);
    ASSERT_EQ(partition.cores.size(),2u);
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=0;
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
    EXPECT_FALSE(result.search_converged);
    EXPECT_EQ(result.reason,"block-sweep-budget");
    EXPECT_FALSE(result.endpoint_certified);
    EXPECT_FALSE(result.fit.search_completed);
    EXPECT_EQ(result.fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Failed);
    EXPECT_EQ(result.fit.components.size(),1u);
}
}
TEST(JointFixedNeighborBlockCoordinateTest, SmallSingleBlockFixturesUseExistingEndpointCertification)
{
    CheckSingleBlock("chain",8); CheckSingleBlock("cube",8);
    CheckSingleBlock("chain",32); CheckSingleBlock("cube",32);
}
TEST(JointFixedNeighborBlockCoordinateTest, MultiBlockBudgetFailureRemainsDistinctFromRuntimeConvergence)
{
    CheckMultiBlockBudgetFailure("chain"); CheckMultiBlockBudgetFailure("cube");
}
TEST(JointFixedNeighborBlockCoordinateTest, InvalidCorePolicyKeepsFailureSemantics)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
    const n::Vector eta=n::Vector::Constant(8,std::log(.55)); n::FixedNeighborPolicy policy; policy.core_atoms=0;
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
    EXPECT_FALSE(result.search_converged); EXPECT_EQ(result.reason,"block-invalid-partition");
    EXPECT_FALSE(result.endpoint_certified); EXPECT_EQ(result.sweeps.size(),0u);
}
TEST(JointFixedNeighborBlockCoordinateTest, SameEtaRawPrimaryAndReferenceDiagnosticsPreserveSearchState)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const n::Vector initial_eta=n::Vector::Constant(8,std::log(.55));
    const auto search=n::SearchFixedNeighbor(problem,initial_eta);
    ASSERT_TRUE(search.search_converged);
    const n::Vector eta=search.state.eta,beta=search.state.beta;
    const n::Domain domain=n::ProfileDomain(data.domain,data.layout);
    const auto context=n::ProfileContext(data.context,data.layout,data.domain.rows);
    const n::Vector y=data.y;
    const auto raw=n::EvaluateState(domain,y,eta,beta,context);
    const auto primary=n::EvaluateProfile(domain,y,eta,false,&context);
    const auto reference=n::EvaluateProfile(domain,y,eta,true,&context);
    const auto raw_assessment=n::AssessEvaluated(domain,y,raw,reference,context,true);
    const auto primary_assessment=n::AssessEvaluated(domain,y,primary,reference,context,false);
    const auto reference_assessment=n::AssessEvaluated(domain,y,reference,primary,context,false);
    const auto raw_trust=n::CheckTrust(domain,y,raw,context,reference);
    const auto primary_trust=n::CheckTrust(domain,y,primary,context,reference);
    const auto reference_trust=n::CheckTrust(domain,y,reference,context,primary);

    EXPECT_TRUE((search.state.eta.array()==eta.array()).all());
    EXPECT_TRUE((raw.eta.array()==eta.array()).all());
    EXPECT_TRUE((primary.eta.array()==eta.array()).all());
    EXPECT_TRUE((reference.eta.array()==eta.array()).all());
    EXPECT_TRUE(raw.valid); EXPECT_TRUE(primary.valid); EXPECT_TRUE(reference.valid);
    EXPECT_TRUE(std::isfinite(raw.gradient.lpNorm<Eigen::Infinity>()));
    EXPECT_TRUE(std::isfinite(primary.gradient.lpNorm<Eigen::Infinity>()));
    EXPECT_TRUE(std::isfinite(reference.gradient.lpNorm<Eigen::Infinity>()));
    const auto difference=[](const n::Vector & lhs,const n::Vector & rhs) {
        return ((lhs-rhs).array().abs()/(1+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
    };
    EXPECT_TRUE(std::isfinite(difference(beta,primary.beta)));
    EXPECT_TRUE(std::isfinite(difference(primary.beta,reference.beta)));
    EXPECT_TRUE(std::isfinite(difference(beta,reference.beta)));
    for(const auto * assessment:{&raw_assessment,&primary_assessment,&reference_assessment})
    {
        EXPECT_TRUE(assessment->primary.valid);
        EXPECT_TRUE(assessment->reference.valid);
        EXPECT_TRUE(assessment->design.has_value());
        EXPECT_TRUE(assessment->widths.has_value());
        EXPECT_TRUE(assessment->jacobian.has_value());
        EXPECT_TRUE(assessment->normalized_widths.has_value());
        EXPECT_EQ(assessment->gradient,assessment->primary.gradient.lpNorm<Eigen::Infinity>()<=1e-12 &&
            assessment->reference.gradient.lpNorm<Eigen::Infinity>()<=1e-12);
    }
    EXPECT_TRUE(raw_trust.primary_valid); EXPECT_TRUE(primary_trust.primary_valid); EXPECT_TRUE(reference_trust.primary_valid);
}
TEST(JointFixedNeighborBlockCoordinateTest, CertifiedLocalCandidatesUseExistingChecksAndBoundedCores)
{
    JointProblem problem(second_stage_test::OperatorWorkload("chain",32));
    const n::Vector eta=n::Vector::Constant(32,std::log(.55));
    n::FixedNeighborPolicy policy; policy.core_atoms=16; policy.maximum_sweeps=3;
    policy.certify_local_candidates=true;
    const auto result=n::SearchFixedNeighbor(problem,eta,policy);
    ASSERT_FALSE(result.blocks.empty());
    for(const auto & block:result.blocks)
    {
        ASSERT_TRUE(block.local_assessment_attempted);
        EXPECT_LE(block.local_assessment_columns,2*static_cast<Eigen::Index>(policy.core_atoms));
        EXPECT_GT(block.local_assessment_rows,0);
        if(!block.local_assessment_passed)
        {
            EXPECT_FALSE(block.accepted);
            EXPECT_EQ(block.reason,"local-endpoint-uncertified");
        }
        if(block.local_assessment_passed)
        {
            EXPECT_TRUE(block.local_inner_passed);
            EXPECT_TRUE(block.local_gradient_passed);
            EXPECT_TRUE(block.local_correction_passed);
            EXPECT_TRUE(block.local_identified);
            EXPECT_TRUE(block.local_trust_passed);
        }
    }
    for(const auto & sweep:result.sweeps)
    {
        EXPECT_EQ(sweep.local_assessments,sweep.block_solves);
        EXPECT_LE(sweep.objective_after,sweep.objective_before+n::BlockObjectiveReplayEnclosure(sweep.objective_before));
        EXPECT_LE(sweep.maximum_local_assessment_columns,2*policy.core_atoms);
    }
}
TEST(JointFixedNeighborBlockCoordinateTest, LocalCertificationRequiresEveryExistingEndpointCheck)
{
    n::Assessment assessment; n::TrustEvidence trust;
    assessment.inner=assessment.gradient=assessment.local=assessment.identified=true; trust.passed=true;
    EXPECT_TRUE(n::IsCertifiedLocalEndpoint(assessment,trust));
    assessment.inner=false; EXPECT_FALSE(n::IsCertifiedLocalEndpoint(assessment,trust));
    assessment.inner=true; assessment.gradient=false; EXPECT_FALSE(n::IsCertifiedLocalEndpoint(assessment,trust));
    assessment.gradient=true; assessment.local=false; EXPECT_FALSE(n::IsCertifiedLocalEndpoint(assessment,trust));
    assessment.local=true; assessment.identified=false; EXPECT_FALSE(n::IsCertifiedLocalEndpoint(assessment,trust));
    assessment.identified=true; trust.passed=false; EXPECT_FALSE(n::IsCertifiedLocalEndpoint(assessment,trust));
}
