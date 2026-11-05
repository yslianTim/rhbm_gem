#include <gtest/gtest.h>
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/GaussSeidelBcd.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include "support/JointComponentChecks.hpp"
#include "support/JointRuntimeJson.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <cmath>
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>

namespace {
namespace core=rhbm_gem::core;
namespace n=core::joint_component;
namespace p=second_stage_test::matched::joint_component;
core::JointProblemInput Snapshot()
{
    core::JointProblemInput input; input.atom_ids={"a","b"}; input.support.resize(2);
    for(int a=0;a<2;++a) for(int k=0;k<40;++k)
    {
        const double square=.003*k*k; const auto basis=n::EvaluateKernel(square,.5,2.5);
        input.support[static_cast<std::size_t>(a)].push_back({input.observations.size(),square});
        input.row_ids.push_back(std::to_string(input.observations.size()));
        input.observations.push_back(2*basis.gaussian+.2*basis.charge);
    }
    return input;
}
;
}

TEST(JointComponentRuntimeTest, PublicObjectivesUseParentScaleAndIncludeConstantRows)
{
    auto input=Snapshot();
    input.row_ids.push_back("constant"); input.observations.push_back(7);
    const core::JointProblem problem(input); const double scale=problem.ObservationScale();
    ASSERT_GT(scale,1);
    const auto fit=p::FitBoundedJointComponentsForTesting(problem,{.55,.55},2.2);
    ASSERT_TRUE(fit.prediction && fit.objective && fit.assembled_state);
    double total{},components{};
    for(std::size_t r=0;r<input.observations.size();++r)
        total+=std::pow((*fit.prediction)[r]-input.observations[r],2)/(2*scale*scale);
    EXPECT_GT(total,1e-3);
    EXPECT_NEAR(*fit.objective,total,1e-15);
    EXPECT_DOUBLE_EQ(*fit.objective,fit.assembled_state->objective);
    for(const auto & component:fit.components)
    {
        ASSERT_TRUE(component.state); double local{};
        for(auto r:component.rows) local+=std::pow((*fit.prediction)[r]-input.observations[r],2)/(2*scale*scale);
        EXPECT_NEAR(component.state->objective,local,1e-16);
        components+=component.state->objective;
    }
    EXPECT_NEAR(components+49/(2*scale*scale),total,1e-15);
}

TEST(JointComponentRuntimeTest, AssemblyPreservesSuppliedCoefficientsWithoutGlobalReprofile)
{
    const core::JointProblem problem(Snapshot());
    const auto & data=core::JointProblemAccess::Get(problem);
    auto fit=p::FitBoundedJointComponentsForTesting(problem,{.55,.55},2.2);
    auto components=fit.components;
    ASSERT_TRUE(components[0].state && components[1].state);
    components[0].state->ac[0]+=.1;
    const double actual=components[0].state->ac[0];
    const auto assembled=n::AssembleComponents(data.domain,data.y,data.partition,data.context,components);
    ASSERT_TRUE(assembled.available);
    EXPECT_DOUBLE_EQ(assembled.beta(0),actual);
    EXPECT_FALSE(assembled.consistency_passed);
    const auto replay=n::EvaluateState(data.domain,data.y,assembled.eta,assembled.beta,data.context);
    EXPECT_LT((assembled.prediction-replay.x*assembled.beta).lpNorm<Eigen::Infinity>(),
        2e-12+2e-13*assembled.prediction.lpNorm<Eigen::Infinity>());
}

TEST(JointComponentRuntimeTest, AssemblyRejectsMappingMismatchAndMissingNumericalState)
{
    const core::JointProblem problem(Snapshot()); const auto & data=core::JointProblemAccess::Get(problem);
    const auto fit=p::FitBoundedJointComponentsForTesting(problem,{.55,.55},2.2);
    ASSERT_EQ(fit.components.size(),2u);
    auto mapping_mismatch=fit.components;
    mapping_mismatch.front().atoms.front()=mapping_mismatch.back().atoms.front();
    const auto bad_mapping=n::AssembleComponents(data.domain,data.y,data.partition,data.context,mapping_mismatch);
    EXPECT_FALSE(bad_mapping.available);
    EXPECT_FALSE(bad_mapping.mapping_valid);
    EXPECT_EQ(bad_mapping.mapping_reason,"component-identity-mapping-mismatch");

    auto missing_state=fit.components;
    missing_state.front().state.reset();
    const auto incomplete=n::AssembleComponents(data.domain,data.y,data.partition,data.context,missing_state);
    EXPECT_FALSE(incomplete.available);
    EXPECT_TRUE(incomplete.mapping_valid);
    EXPECT_EQ(incomplete.mapping_reason,"missing-component-state");

    auto unconverged=fit.components;
    std::erase_if(unconverged.front().evidence,[](const auto & check){return check.name=="endpoint-verification";});
    const auto rejected=n::AssembleComponents(data.domain,data.y,data.partition,data.context,unconverged);
    EXPECT_FALSE(rejected.available);
    EXPECT_EQ(rejected.mapping_reason,"component-not-converged");
}

TEST(JointComponentRuntimeTest, ReplayTrustRejectsIndependentDisagreements)
{
    const core::JointProblem problem(Snapshot()); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & view=data.partition.components[0]; const auto y=p::Select(data.y,view.rows);
    const auto context=p::ComponentContext(data.context,view,true);
    const auto e=n::EvaluateProfile(view.domain,y,n::Vector::Constant(1,std::log(.55)),&context);
    ASSERT_TRUE(e.valid); EXPECT_TRUE(n::CheckReplay(view.domain,y,e,context).passed);
    auto bad=e; bad.valid=false;
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).passed);
    bad=e; bad.certificate.available=false;
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).passed);
    bad=e; bad.certificate.kkt_passed=false;
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).passed);
    bad=e; bad.certificate.feasible=false;
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).passed);
    bad=e; bad.gradient(0)+=.01;
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).gradient_passed);
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).passed);
    bad=e; bad.x.coeffRef(0,0)+=.01;
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).prediction_passed);
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).passed);
    bad=e; bad.certificate.projected_kkt+=1e-11;
    EXPECT_FALSE(n::CheckReplay(view.domain,y,bad,context).passed);
    const auto independent=n::CheckIndependentReplay(view.domain,y,e,context);
    ASSERT_TRUE(independent.passed);
    EXPECT_TRUE(n::CheckTrust(e,independent).passed);
    auto bad_independent=independent;
    bad_independent.passed=false;
    bad_independent.prediction_passed=false;
    bad_independent.reason="independent-replay-disagreement";
    EXPECT_FALSE(n::CheckTrust(e,bad_independent).passed);
}

TEST(JointComponentRuntimeTest, RuntimeConvergenceIsDerivedFromActualEvidenceNotTermination)
{
    using Status=core::JointCheckStatus;
    auto fit=p::FitBoundedJointComponentsForTesting(core::JointProblem(Snapshot()),{.55,.55},2.2);
    ASSERT_EQ(fit.RuntimeConvergence(),Status::Passed);
    for(auto & component:fit.components)
    {
        component.search_completed=false; component.stop_reason="profile-budget";
        EXPECT_EQ(component.RuntimeConvergence(),Status::Passed);
        for(const char * name:{"two-step-derivative","richardson","precision-50-100","boundary-audit","regular-certificate"})
        {
            EXPECT_EQ(std::find_if(component.evidence.begin(),component.evidence.end(),[&](const auto & check) {
                return check.name==name;
            }),component.evidence.end());
            EXPECT_EQ(std::find_if(fit.evidence.begin(),fit.evidence.end(),[&](const auto & check) {
                return check.name==name;
            }),fit.evidence.end());
        }
        for(const char * name:{"ac-feasibility","conditional-ac-optimality","independent-ac-optimality",
            "independent-kkt-consistency","independent-replay","eta-gradient-consistency",
            "width-stationarity","endpoint-verification","numerical-qualification"})
            EXPECT_NE(std::find_if(component.evidence.begin(),component.evidence.end(),[&](const auto & check) {
                return check.name==name;
            }),component.evidence.end());
    }
    fit.search_completed=false;
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Passed);
    auto & evidence=fit.components[0].evidence;
    const auto set_status=[&](const std::string & name,Status status) {
        const auto found=std::find_if(evidence.begin(),evidence.end(),[&](const auto & check){return check.name==name;});
        ASSERT_NE(found,evidence.end()); found->status=status;
    };
    for(const auto status:{Status::Passed,Status::Failed,Status::Unavailable,Status::NotRun})
    {
        set_status("independent-replay",status);
        EXPECT_EQ(fit.RuntimeConvergence(),status);
    }
    set_status("independent-replay",Status::Passed);
    set_status("conditional-ac-optimality",Status::Unavailable);
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Unavailable);
    set_status("conditional-ac-optimality",Status::Passed);
    set_status("width-stationarity",Status::Failed);
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Failed);
    set_status("width-stationarity",Status::Passed);
    set_status("active-set-lm-subproblems",Status::Failed);
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Passed);
    set_status("active-set-lm-subproblems",Status::Passed);
    evidence.push_back({"local-correction",Status::Failed,rhbm_gem::JointEvidenceScope::ComponentLocal,.1,1e-10,
        "historical-correction-evidence"});
    set_status("local-correction",Status::Failed);
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Passed);
    fit.components[0].state.reset();
    EXPECT_EQ(fit.components[0].RuntimeConvergence(),Status::Unavailable);
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Unavailable);
}

TEST(JointComponentRuntimeTest, MissingOrWrongScopeEvidenceCannotVacuouslyConverge)
{
    using Status=core::JointCheckStatus;
    const auto fit=p::FitBoundedJointComponentsForTesting(core::JointProblem(Snapshot()),{.55,.55},2.2);
    auto component=fit.components[0];
    component.evidence.clear(); EXPECT_EQ(component.RuntimeConvergence(),Status::Unavailable);
    component.evidence=fit.evidence; EXPECT_EQ(component.RuntimeConvergence(),Status::Unavailable);
    auto global=fit;
    for(auto & check:global.evidence) if(check.name=="assembled-consistency") check.status=Status::Failed;
    EXPECT_EQ(global.RuntimeConvergence(),Status::Failed);
    global=fit;
    const auto mapping=std::find_if(global.evidence.begin(),global.evidence.end(),[](const auto & check) {
        return check.name=="assembly-mapping";
    });
    ASSERT_NE(mapping,global.evidence.end());
    mapping->status=Status::Unavailable;
    EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    global=fit; global.evidence=fit.components[0].evidence;
    EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    global=fit; global.prediction.reset(); EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    global=fit; global.available_row_mask[0]=false; EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    global=fit; global.components.clear(); EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    EXPECT_EQ(core::JointFitResult{}.RuntimeConvergence(),Status::Unavailable);
}

TEST(JointComponentRuntimeTest, GenericNumericalOutcomesAreNotPersistedAsFormalResults)
{
    namespace io=rhbm_gem::joint_result_io;
    for(int variant=0;variant<4;++variant)
    {
        auto input=Snapshot();
        input.row_ids.push_back("constant"); input.observations.push_back(3);
        if(variant==1) {input.atom_ids.push_back("unobserved"); input.support.emplace_back();}
        if(variant==2) for(auto & y:input.observations) y=0;
        std::vector<double> initial(input.atom_ids.size(),.55);
        if(variant==3) initial[0]=std::numeric_limits<double>::quiet_NaN();
        const auto fit=p::FitBoundedJointComponentsForTesting(core::JointProblem(input),initial,2.2);
        const auto saved=core::CaptureJointAnalysisResult(fit);
        ASSERT_TRUE(saved.metadata.software);
        EXPECT_FALSE(saved.metadata.model_sha256); EXPECT_FALSE(saved.metadata.map_sha256);
        EXPECT_FALSE(saved.metadata.map_normalization);
        EXPECT_EQ(saved.metadata.software->source_sha256.size(),64);
        EXPECT_EQ(saved.runtime_convergence,fit.RuntimeConvergence());
        EXPECT_EQ(saved.available_row_mask,fit.available_row_mask);
        EXPECT_EQ(saved.objective,fit.objective);
        EXPECT_EQ(saved.observation_scale,fit.observation_scale);
        ASSERT_EQ(saved.components.size(),fit.components.size());
        for(std::size_t k=0;k<fit.components.size();++k)
        {
            EXPECT_EQ(saved.components[k].runtime_convergence,fit.components[k].RuntimeConvergence());
            EXPECT_EQ(saved.components[k].state.has_value(),fit.components[k].state.has_value());
            if(fit.components[k].state) EXPECT_EQ(saved.components[k].state->ac,fit.components[k].state->ac);
        }
        EXPECT_THROW(io::Encode(saved),std::invalid_argument);
        if(variant==3) EXPECT_TRUE(std::isnan(saved.initialization.tau[0]));
    }
}

TEST(JointComponentRuntimeTest, PositiveMapScalingPreservesKernelCoefficientUnits)
{
    auto input=Snapshot();
    const auto original=p::FitBoundedJointComponentsForTesting(core::JointProblem(input),{.55,.55},2.2);
    constexpr double divisor=3.25;
    for(auto & y:input.observations) y/=divisor;
    const auto normalized=p::FitBoundedJointComponentsForTesting(core::JointProblem(input),{.55,.55},2.2);
    ASSERT_EQ(original.RuntimeConvergence(),core::JointCheckStatus::Passed);
    ASSERT_EQ(normalized.RuntimeConvergence(),core::JointCheckStatus::Passed);
    ASSERT_TRUE(original.assembled_state); ASSERT_TRUE(normalized.assembled_state);
    for(std::size_t i=0;i<original.assembled_state->ac.size();++i)
        EXPECT_NEAR(original.assembled_state->ac[i],divisor*normalized.assembled_state->ac[i],
            1e-10*(1+std::abs(original.assembled_state->ac[i])));
    for(std::size_t i=0;i<original.assembled_state->tau.size();++i)
        EXPECT_NEAR(original.assembled_state->tau[i],normalized.assembled_state->tau[i],1e-10);
    ASSERT_TRUE(original.prediction); ASSERT_TRUE(normalized.prediction);
    for(std::size_t i=0;i<original.prediction->size();++i)
        EXPECT_NEAR((*original.prediction)[i],divisor*(*normalized.prediction)[i],
            1e-10*(1+std::abs((*original.prediction)[i])));
}

TEST(JointComponentRuntimeTest, HistoricalLmFailureDoesNotVetoAVerifiedRecoveredEndpoint)
{
    using Status=core::JointCheckStatus;
    const auto fit=p::FitBoundedJointComponentsForTesting(core::JointProblem(Snapshot()),{.55,.55},2.2);
    ASSERT_EQ(fit.RuntimeConvergence(),Status::Passed);
    auto component=fit.components.front();
    const auto lm=std::find_if(component.evidence.begin(),component.evidence.end(),[](const auto & check) {
        return check.name=="active-set-lm-subproblems";
    });
    ASSERT_NE(lm,component.evidence.end());
    lm->status=Status::Failed; lm->reason="solves=3; failed=1; safely-rolled-back";
    EXPECT_EQ(component.RuntimeConvergence(),Status::Passed);
    const core::JointProblem problem(Snapshot());
    const auto & data=core::JointProblemAccess::Get(problem);
    auto components=fit.components; components.front()=component;
    EXPECT_TRUE(n::AssembleComponents(data.domain,data.y,data.partition,data.context,components).available);
    // Endpoint failures and incomplete recovery must continue to reject output.
    for(const char * name:{"width-stationarity","endpoint-verification","independent-replay","complete-sweep-stationarity","accepted-state-integrity"})
    {
        auto bad=component;
        for(auto & check:bad.evidence) if(check.name==name) check.status=Status::Failed;
        EXPECT_NE(bad.RuntimeConvergence(),Status::Passed)<<name;
        components.front()=bad;
        EXPECT_FALSE(n::AssembleComponents(data.domain,data.y,data.partition,data.context,components).available)<<name;
    }
    auto missing=component;
    std::erase_if(missing.evidence,[](const auto & check){return check.name=="accepted-state-integrity";});
    EXPECT_EQ(missing.RuntimeConvergence(),Status::Unavailable);
    component.state->tau[0]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_NE(component.RuntimeConvergence(),Status::Passed);
}

namespace {
struct RecoveryAssemblyRun
{
    core::JointFitResult fit;
    n::BcdSearchResult search;
    std::size_t injected{};
};
RecoveryAssemblyRun RunRecoveryAssembly(bool inject,bool persistent=false,std::size_t sweeps=200,
    bool corrupt_candidate=false)
{
    RecoveryAssemblyRun out;
    auto input=second_stage_test::SyntheticJointWorkload("chain",2);
    input.selection_domain.emplace(); input.selection_domain->target_indices={0,1};
    input.selection_domain->contributor_policy=std::string(rhbm_gem::PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->pure_halo_support_exclusion=rhbm_gem::JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    const core::JointProblem problem(std::move(input));
    n::BcdExecutionPolicy policy; policy.eta_upper_bound=std::log(2.2);
    policy.configure_for_testing=[&](n::BcdSearchPolicy & local) {
        local.core_atoms=1; local.maximum_sweeps=sweeps; local.collect_diagnostics=true;
        local.force_width_failure_for_testing=[&](std::size_t block,std::size_t sweep) {
            if(inject && block==0 && (persistent || sweep==0)) {++out.injected; return true;}
            return false;
        };
        if(corrupt_candidate) local.candidate_fault_for_testing=[](std::size_t sweep,std::size_t block,
            const n::BlockCoordinateState &) {
            return sweep==3 && block==1 ? n::BcdCandidateFaultForTesting::NaNParameter :
                n::BcdCandidateFaultForTesting::None;
        };
    };
    policy.observe_search_for_testing=[&](const auto & search) {out.search=search;};
    out.fit=n::FitBcdComponents(problem,{.55,.55},policy);
    return out;
}
}

TEST(JointComponentRuntimeTest, TransientWidthFailureRecoversAndAssemblesQualifiedComponent)
{
    using Status=core::JointCheckStatus;
    const auto baseline=RunRecoveryAssembly(false);
    ASSERT_EQ(baseline.fit.RuntimeConvergence(),Status::Passed);
    const auto recovered=RunRecoveryAssembly(true);
    const auto & search=recovered.search;
    EXPECT_EQ(recovered.injected,1u);
    ASSERT_FALSE(search.blocks.empty());
    const auto & first=search.blocks.front();
    EXPECT_EQ(first.status,"retry-pending"); EXPECT_FALSE(first.accepted);
    ASSERT_TRUE(first.failure_diagnostics);
    EXPECT_TRUE(first.failure_diagnostics->initial_profile_trusted);
    EXPECT_EQ(first.failure_diagnostics->local_search_stop_reason,
        "bounded-subproblem:injected-width-subproblem-failure");
    ASSERT_GT(search.blocks.size(),2u);
    const auto & retry=search.blocks[2];
    EXPECT_EQ(retry.block,1u); EXPECT_EQ(retry.sweep,2u);
    EXPECT_TRUE(retry.recovery_retry); EXPECT_TRUE(retry.recovery_success);
    EXPECT_TRUE(retry.accepted); EXPECT_GT(retry.accepted_updates,0);
    EXPECT_EQ(search.recovery_attempts,1u); EXPECT_EQ(search.recovery_successes,1u);
    EXPECT_EQ(search.unresolved_failures,0u); EXPECT_TRUE(search.accepted_state_integrity);
    EXPECT_EQ(search.work.bounded_lm.failed_subproblems,1u);
    ASSERT_TRUE(search.final_sweep);
    EXPECT_LE(search.final_sweep->global_ac_kkt,1e-7);
    EXPECT_LE(search.final_sweep->projected_eta_gradient_inf_norm,1e-9);
    EXPECT_EQ(search.final_sweep->maximum_width_bound_violation,0.0);
    EXPECT_GT(search.confirmed_stationarity_sweep,0u);
    ASSERT_TRUE(search.certification);
    EXPECT_TRUE(search.certification->qualified());
    EXPECT_TRUE(search.certification->assessment.endpoint_verified);
    EXPECT_TRUE(search.certification->assessment.independent_replay.passed);
    ASSERT_EQ(recovered.fit.components.size(),1u);
    const auto & component=recovered.fit.components.front();
    EXPECT_TRUE(component.state); EXPECT_EQ(component.RuntimeConvergence(),Status::Passed);
    const auto lm=std::find_if(component.evidence.begin(),component.evidence.end(),[](const auto & check) {
        return check.name=="active-set-lm-subproblems";
    });
    ASSERT_NE(lm,component.evidence.end());
    EXPECT_EQ(lm->status,Status::Failed);
    EXPECT_NE(lm->reason.find("; failed=1;"),std::string::npos);
    ASSERT_TRUE(recovered.fit.assembled_state); ASSERT_TRUE(recovered.fit.prediction);
    EXPECT_EQ(recovered.fit.RuntimeConvergence(),Status::Passed);
    EXPECT_EQ(recovered.fit.TargetRuntimeConvergence(),Status::Passed);
    for(const char * name:{"assembly-mapping","assembled-consistency","endpoint-verification","independent-replay"})
    {
        const auto check=std::find_if(recovered.fit.evidence.begin(),recovered.fit.evidence.end(),
            [&](const auto & item){return item.name==name;});
        ASSERT_NE(check,recovered.fit.evidence.end())<<name;
        EXPECT_EQ(check->status,Status::Passed)<<name;
    }
    ASSERT_TRUE(baseline.fit.assembled_state); ASSERT_TRUE(baseline.fit.prediction);
    EXPECT_NEAR(*recovered.fit.objective,*baseline.fit.objective,1e-12);
    for(std::size_t i=0;i<recovered.fit.prediction->size();++i)
        EXPECT_NEAR((*recovered.fit.prediction)[i],(*baseline.fit.prediction)[i],1e-6);
    for(std::size_t i=0;i<component.state->tau.size();++i)
        EXPECT_NEAR(component.state->tau[i],baseline.fit.assembled_state->tau[i],1e-6);
    for(std::size_t i=0;i<component.state->ac.size();++i)
        EXPECT_NEAR(component.state->ac[i],baseline.fit.assembled_state->ac[i],1e-5);
    const auto saved=core::CaptureJointAnalysisResult(recovered.fit);
    ASSERT_TRUE(saved.assembled_state); ASSERT_TRUE(saved.components.front().state);
    EXPECT_EQ(saved.runtime_convergence,Status::Passed);
    EXPECT_EQ(saved.target_runtime_convergence,Status::Passed);
    EXPECT_EQ(saved.components.front().runtime_convergence,Status::Passed);
    const auto historical=std::find_if(saved.components.front().evidence.begin(),saved.components.front().evidence.end(),
        [](const auto & check){return check.name=="active-set-lm-subproblems";});
    ASSERT_NE(historical,saved.components.front().evidence.end());
    EXPECT_EQ(historical->status,Status::Failed);
    EXPECT_EQ(historical->reason,lm->reason);
}

TEST(JointComponentRuntimeTest, WidthRecoveryExhaustionAndUnqualifiedTerminalStatesRejectAssembly)
{
    using Status=core::JointCheckStatus;
    const auto exhausted=RunRecoveryAssembly(true,true);
    EXPECT_EQ(exhausted.injected,2u);
    EXPECT_EQ(exhausted.search.reason,"block-recovery-exhausted");
    EXPECT_EQ(exhausted.search.recovery_attempts,1u);
    EXPECT_EQ(exhausted.search.recovery_successes,0u);
    EXPECT_EQ(exhausted.search.unresolved_failures,1u);
    EXPECT_FALSE(exhausted.fit.components.front().state);
    EXPECT_FALSE(exhausted.fit.assembled_state);
    EXPECT_NE(exhausted.fit.RuntimeConvergence(),Status::Passed);

    const auto early=RunRecoveryAssembly(true,false,2);
    EXPECT_EQ(early.search.recovery_successes,1u);
    ASSERT_TRUE(early.search.final_sweep);
    EXPECT_GT(early.search.final_sweep->projected_eta_gradient_inf_norm,1e-9);
    EXPECT_FALSE(early.fit.components.front().state); EXPECT_FALSE(early.fit.assembled_state);

    const auto corrupt=RunRecoveryAssembly(true,false,200,true);
    EXPECT_EQ(corrupt.search.recovery_successes,1u);
    EXPECT_FALSE(corrupt.search.accepted_state_integrity);
    EXPECT_FALSE(corrupt.fit.components.front().state); EXPECT_FALSE(corrupt.fit.assembled_state);

    // Start from genuinely recovered/certified evidence; only invalidate or remove it.
    const auto recovered=RunRecoveryAssembly(true);
    ASSERT_EQ(recovered.fit.RuntimeConvergence(),Status::Passed);
    const auto & data=core::JointProblemAccess::Get(*recovered.fit.problem);
    auto context=data.context; context.eta_upper_bound=std::log(2.2);
    auto altered=recovered.search.state;
    altered.beta(0)+=.1;
    const auto recertified=n::CertifyBcdEndpoint(data.domain,data.y,altered,context,
        0,n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(recertified.qualified());
    EXPECT_FALSE(recertified.assessment.endpoint_verified);
    EXPECT_FALSE(recertified.assessment.independent_replay.passed);
    auto rejected=recovered.fit.components;
    rejected.front().evidence=n::AssessmentEvidence(recertified.assessment,rhbm_gem::JointEvidenceScope::ComponentLocal);
    EXPECT_FALSE(n::AssembleComponents(data.domain,data.y,data.partition,context,rejected,
        n::PracticalJointConvergenceCriteria).available);
    for(const char * name:{"conditional-ac-optimality","width-stationarity","endpoint-verification",
        "independent-replay","accepted-state-integrity","complete-sweep-stationarity"})
    {
        for(bool missing:{false,true})
        {
            auto components=recovered.fit.components;
            auto & evidence=components.front().evidence;
            const auto check=std::find_if(evidence.begin(),evidence.end(),
                [&](const auto & item){return item.name==name;});
            ASSERT_NE(check,evidence.end());
            if(missing) evidence.erase(check); else check->status=Status::Failed;
            EXPECT_NE(components.front().RuntimeConvergence(),Status::Passed)<<name<<missing;
            EXPECT_FALSE(n::AssembleComponents(data.domain,data.y,data.partition,context,components,
                n::PracticalJointConvergenceCriteria).available)<<name<<missing;
        }
    }
}

TEST(JointComponentRuntimeTest, ComponentValidationProvidesReasonWithoutEnteringSearch)
{
    const core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=core::JointProblemAccess::Get(problem);
    EXPECT_TRUE(n::BcdComponentResult{}.search.reason.empty());
    for(int variant=0;variant<8;++variant)
    {
        SCOPED_TRACE(variant);
        auto prepared=n::MakeComponentProblemView(data.layout,data.partition.components.front());
        n::Vector widths=n::Vector::Constant(3,.55);
        if(variant==0) prepared.local_layout.full_atoms.clear();
        if(variant==1) prepared.local_layout.informative_rows.clear();
        if(variant==2) widths(0)=0;
        if(variant==3) widths(0)=-1;
        if(variant==4) widths(0)=std::numeric_limits<double>::quiet_NaN();
        if(variant==5) widths(0)=std::numeric_limits<double>::infinity();
        if(variant==6) prepared.local_to_parent_atoms[0]=-1;
        if(variant==7) prepared.local_to_parent_atoms[0]=widths.size();
        // Empty layouts retain precedence over invalid widths.
        if(variant<2) widths(0)=0;
        bool configured{},invoked{};
        n::BcdExecutionPolicy policy;
        policy.configure_for_testing=[&](auto &) {configured=true;};
        policy.observe_search_for_testing=[&](const auto &) {invoked=true;};
        std::vector<n::JointProgressEvent> events;
        const n::JointProgressComponent component{1,1,"component",3,prepared.local_to_parent_rows.size()};
        const auto result=n::SolveBcdComponent(prepared,widths,data.context,policy,
            [&](const auto & event) {events.push_back(event);},&component);
        EXPECT_EQ(result.search.reason,variant<2 ? "unobserved-full-parameters" : "invalid-initial-widths");
        EXPECT_FALSE(configured); EXPECT_FALSE(invoked); EXPECT_TRUE(events.empty());
        EXPECT_FALSE(result.succeeded); EXPECT_FALSE(result.search.search_converged);
        EXPECT_FALSE(result.search.accepted_state_integrity); EXPECT_FALSE(result.search.certification);
        EXPECT_EQ(result.search.sweep_count,0u); EXPECT_EQ(result.search.state.eta.size(),0);
        EXPECT_EQ(result.search.memory_estimate.status,n::MemoryStatus::Pass);
    }
}

TEST(JointComponentRuntimeTest, ComponentViewPreflightPreservesReasonAndMemoryDiagnostics)
{
    const core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=core::JointProblemAccess::Get(problem);
    for(int variant=0;variant<3;++variant)
    {
        SCOPED_TRACE(variant);
        auto view=data.partition.components.front(); auto layout=data.layout;
        auto mappings=std::make_shared<n::PartitionMappings>(*view.mappings); view.mappings=mappings;
        if(variant==0) mappings->row_to_local[0]=-1;
        if(variant==1) mappings->atom_to_local[0]=-1;
        // Repeated valid selectors exercise the estimate without allocating the rejected domain.
        if(variant==2) layout.full_atoms.assign(100000,0);
        bool configured{},invoked{};
        n::BcdExecutionPolicy policy;
        policy.configure_for_testing=[&](auto &) {configured=true;};
        policy.observe_search_for_testing=[&](const auto &) {invoked=true;};
        std::vector<n::JointProgressEvent> events;
        const n::JointProgressComponent component{1,1,view.id,view.atoms.size(),view.rows.size()};
        const auto result=n::SolveBcdComponentView(view,layout,n::Vector::Constant(3,.55),
            data.context,policy,[&](const auto & event) {events.push_back(event);},&component);
        EXPECT_EQ(result.search.reason,variant==0 ? "block-invalid-partition" : variant==1 ?
            "resource-limit:bcd-preparation:size-overflow" :
            "resource-limit:bcd-preparation:estimated-peak-exceeds-budget");
        EXPECT_FALSE(configured); EXPECT_FALSE(invoked); EXPECT_TRUE(events.empty());
        EXPECT_FALSE(result.succeeded); EXPECT_FALSE(result.search.search_converged);
        EXPECT_FALSE(result.search.accepted_state_integrity); EXPECT_FALSE(result.search.certification);
        EXPECT_EQ(result.search.sweep_count,0u); EXPECT_EQ(result.search.state.eta.size(),0);
        const auto & memory=result.search.memory_estimate;
        EXPECT_EQ(memory.status,variant==0 ? n::MemoryStatus::Pass : variant==1 ?
            n::MemoryStatus::Overflow : n::MemoryStatus::Reject);
        if(variant!=0)
        {
            EXPECT_EQ(memory.stage,n::MemoryStage::BcdPreparation);
            EXPECT_EQ(memory.reason,variant==1 ? "size-overflow" : "estimated-peak-exceeds-budget");
        }
        if(variant==2) EXPECT_GT(memory.estimated_peak_bytes,memory.budget_bytes);
    }
}

TEST(JointComponentRuntimeTest, ComponentSearchPreflightKeepsValidationOrderAndCertificationEvent)
{
    const core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=core::JointProblemAccess::Get(problem);
    for(int variant=0;variant<5;++variant)
    {
        SCOPED_TRACE(variant);
        auto prepared=n::MakeComponentProblemView(data.layout,data.partition.components.front());
        n::Vector widths=n::Vector::Constant(3,.55);
        n::BcdExecutionPolicy policy; policy.eta_upper_bound=std::log(2.2);
        if(variant==0) policy.eta_upper_bound.reset();
        if(variant==1) policy.eta_upper_bound=std::numeric_limits<double>::quiet_NaN();
        if(variant==2) policy.eta_upper_bound=std::numeric_limits<double>::infinity();
        if(variant==3) widths.setConstant(3);
        if(variant==4) prepared.row_ids=n::Identities{};
        int invocations{}; std::string observed_reason;
        policy.observe_search_for_testing=[&](const auto & search) {++invocations; observed_reason=search.reason;};
        std::vector<n::JointProgressEvent> events;
        const n::JointProgressComponent component{1,1,"component",3,prepared.local_to_parent_rows.size()};
        const auto result=n::SolveBcdComponent(prepared,widths,data.context,policy,
            [&](const auto & event) {events.push_back(event);},&component);
        EXPECT_EQ(result.search.reason,variant==0 ? "invalid-width-upper-bound" : variant<3 ?
            "invalid-global-width-bound" : variant==3 ? "bounded-initial-width-infeasible" : "block-invalid-partition");
        EXPECT_EQ(invocations,1); EXPECT_EQ(result.search.reason,observed_reason);
        EXPECT_FALSE(result.succeeded); EXPECT_FALSE(result.search.search_converged);
        EXPECT_FALSE(result.search.accepted_state_integrity); EXPECT_FALSE(result.search.certification);
        EXPECT_EQ(result.search.sweep_count,0u); EXPECT_EQ(result.search.state.eta.size(),0);
        EXPECT_EQ(result.search.memory_estimate.status,n::MemoryStatus::Pass);
        ASSERT_EQ(events.size(),1u);
        EXPECT_EQ(events[0].phase,n::JointProgressPhase::CertificationStarted);
        EXPECT_EQ(events[0].stop_reason,result.search.reason); EXPECT_FALSE(events[0].trusted_state);
        EXPECT_FALSE(events[0].bcd); EXPECT_FALSE(events[0].accepted_objective);
    }
}

TEST(JointComponentRuntimeTest, CanonicalReasonReachesComponentProgressAndCapturedResult)
{
    const core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto prepared=n::MakeComponentProblemView(data.layout,data.partition.components.front());
    const std::vector<std::string> reasons{"block-stationary","block-sweep-budget",
        "block-recovery-sweep-budget","block-recovery-exhausted","block-cache-replay-failed",
        "block-objective-replay-failed","resource-limit:global-state-evaluation:estimated-peak-exceeds-budget",
        "block-stationary"};
    for(std::size_t variant=0;variant<reasons.size();++variant)
    {
        SCOPED_TRACE(variant);
        n::BcdExecutionPolicy policy; policy.eta_upper_bound=std::log(2.2);
        policy.configure_for_testing=[&](auto & local) {
            if(variant==1 || variant==2) {local.maximum_sweeps=1; local.nominal_sweeps=0;}
            if(variant==2 || variant==3)
            {
                local.core_atoms=1;
                local.force_no_trusted_profile_for_testing=[](auto,auto) {return true;};
            }
            if(variant==4) local.perturb_sweep_replay_objective_for_testing=true;
            if(variant==5) local.candidate_fault_for_testing=[](auto,auto,const auto &) {
                return n::BcdCandidateFaultForTesting::NaNParameter;
            };
            if(variant==6) local.global_state_memory_budget_for_testing=1;
            if(variant==7) local.certify_endpoint=false;
        };
        std::string observed_reason; int invocations{};
        policy.observe_search_for_testing=[&](const auto & search) {++invocations; observed_reason=search.reason;};
        std::vector<n::JointProgressEvent> events;
        const n::JointProgressComponent component{1,1,data.partition.components.front().id,3,
            prepared.local_to_parent_rows.size()};
        const auto result=n::SolveBcdComponent(prepared,n::Vector::Constant(3,.55),data.context,policy,
            [&](const auto & event) {events.push_back(event);},&component);
        EXPECT_EQ(invocations,1); EXPECT_EQ(result.search.reason,reasons[variant]);
        EXPECT_EQ(result.search.reason,observed_reason);
        const bool qualified=result.search.accepted_state_integrity && result.search.certification &&
            result.search.certification->qualified();
        EXPECT_EQ(result.succeeded,result.search.search_converged && qualified);
        EXPECT_EQ(result.succeeded,variant==0);
        if(variant==7) {EXPECT_TRUE(result.search.search_converged); EXPECT_FALSE(result.search.certification);}
        ASSERT_FALSE(events.empty());
        const auto & certification=events.back();
        EXPECT_EQ(certification.phase,n::JointProgressPhase::CertificationStarted);
        EXPECT_EQ(certification.stop_reason,result.search.reason); EXPECT_EQ(certification.trusted_state,qualified);
        EXPECT_EQ(certification.component_id,component.id); EXPECT_EQ(certification.component_index,component.index);
        EXPECT_EQ(certification.profile_evaluations,static_cast<int>(result.search.total_profile_evaluations));
        EXPECT_EQ(certification.accepted_updates,static_cast<int>(result.search.total_accepted_local_updates));
        EXPECT_EQ(certification.accepted_objective,n::LastBcdObjective(result.search));
        EXPECT_EQ(certification.accepted_gradient_inf_norm,n::LastBcdGradient(result.search));
        events.clear();
        const auto fit=n::FitBcdComponents(problem,{.55,.55,.55},policy,
            [&](const auto & event) {events.push_back(event);});
        EXPECT_EQ(invocations,2);
        ASSERT_EQ(fit.components.size(),1u);
        EXPECT_EQ(fit.components[0].stop_reason,result.search.reason);
        EXPECT_EQ(fit.components[0].search_completed,result.succeeded);
        EXPECT_EQ(core::CaptureJointAnalysisResult(fit).components[0].stop_reason,result.search.reason);
        std::vector<n::JointProgressPhase> lifecycle;
        for(const auto & event:events)
        {
            if(event.phase==n::JointProgressPhase::SearchProgress) {EXPECT_TRUE(event.stop_reason.empty()); continue;}
            lifecycle.push_back(event.phase);
            if(event.phase==n::JointProgressPhase::CertificationStarted || event.phase==n::JointProgressPhase::ComponentCompleted)
                EXPECT_EQ(event.stop_reason,result.search.reason);
            else EXPECT_TRUE(event.stop_reason.empty());
        }
        EXPECT_EQ(lifecycle,(std::vector<n::JointProgressPhase>{n::JointProgressPhase::SolverConfigured,
            n::JointProgressPhase::ComponentStarted,n::JointProgressPhase::CertificationStarted,
            n::JointProgressPhase::ComponentCompleted,n::JointProgressPhase::AssemblyStarted,n::JointProgressPhase::AssemblyCompleted}));
    }
}
