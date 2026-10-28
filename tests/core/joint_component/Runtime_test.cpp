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
