#include <gtest/gtest.h>
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
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
namespace p=second_stage_test::matched::joint_abc;
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
    EXPECT_TRUE(n::CheckTrust(view.domain,y,e,context,independent).passed);
    auto bad_independent=independent;
    bad_independent.passed=false;
    bad_independent.prediction_passed=false;
    bad_independent.reason="independent-replay-disagreement";
    EXPECT_FALSE(n::CheckTrust(view.domain,y,e,context,bad_independent).passed);
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
            "independent-kkt-consistency","independent-replay","width-gradient-consistency",
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
        if(variant==3) EXPECT_TRUE(std::isnan(saved.initialization.b[0]));
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
    for(std::size_t i=0;i<original.assembled_state->b.size();++i)
        EXPECT_NEAR(original.assembled_state->b[i],normalized.assembled_state->b[i],1e-10);
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
    component.state->b[0]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_NE(component.RuntimeConvergence(),Status::Passed);
}

namespace {
struct RecoveryAssemblyRun
{
    core::JointFitResult fit;
    n::FixedNeighborSearchResult search;
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
    n::FixedNeighborSearchPolicy policy; policy.eta_upper_bound=std::log(2.2);
    policy.configure_for_testing=[&](n::FixedNeighborPolicy & local) {
        local.core_atoms=1; local.maximum_sweeps=sweeps; local.collect_diagnostics=true;
        local.force_width_failure_for_testing=[&](std::size_t block,std::size_t sweep) {
            if(inject && block==0 && (persistent || sweep==0)) {++out.injected; return true;}
            return false;
        };
        if(corrupt_candidate) local.candidate_fault_for_testing=[](std::size_t sweep,std::size_t block,
            const n::BlockCoordinateState &) {
            return sweep==3 && block==1 ? n::FixedNeighborCandidateFaultForTesting::NaNParameter :
                n::FixedNeighborCandidateFaultForTesting::None;
        };
    };
    policy.observe_search_for_testing=[&](const auto & search) {out.search=search;};
    out.fit=n::FitFixedNeighborComponents(problem,{.55,.55},policy);
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
    EXPECT_EQ(search.work.local_profile_work.bounded_lm.failed_subproblems,1u);
    ASSERT_TRUE(search.final_sweep);
    EXPECT_LE(search.final_sweep->global_ac_kkt,1e-7);
    EXPECT_LE(search.final_sweep->projected_width_gradient_inf_norm,1e-9);
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
    for(std::size_t i=0;i<component.state->b.size();++i)
        EXPECT_NEAR(component.state->b[i],baseline.fit.assembled_state->b[i],1e-6);
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
    EXPECT_GT(early.search.final_sweep->projected_width_gradient_inf_norm,1e-9);
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
    const auto recertified=n::CertifyFixedNeighborEndpoint(data.domain,data.y,altered,context,
        nullptr,0,n::PracticalJointConvergenceCriteria);
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
