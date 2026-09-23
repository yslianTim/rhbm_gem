#include <gtest/gtest.h>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/command/detail/MapSimulation.hpp"
#include "support/JointTestNumerics.hpp"
#include "support/JointRuntimeJson.hpp"
#include <cmath>
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/utils/math/EigenHelper.hpp>

namespace {
namespace core=rhbm_gem::core;
namespace n=core::joint_component;
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
}
TEST(JointComponentRuntimeTest, ImmutableSnapshotAndPartialFailure)
{
    auto input=Snapshot(); input.atom_ids.push_back("unobserved"); input.support.emplace_back();
    input.row_ids.push_back("constant"); input.observations.push_back(3);
    const core::JointProblem problem(input); input.observations[0]=999;
    EXPECT_NE(problem.Input().observations[0],999);
    const auto fit=core::FitJointComponents(problem,{.55,.55,.55});
    ASSERT_EQ(fit.components.size(),3);
    EXPECT_TRUE(fit.components[0].state.has_value()); EXPECT_TRUE(fit.components[1].state.has_value());
    EXPECT_FALSE(fit.components[2].state.has_value());
    EXPECT_FALSE(fit.prediction.has_value()); EXPECT_FALSE(fit.objective.has_value());
    EXPECT_EQ(fit.available_row_mask.size(),81); EXPECT_TRUE(fit.available_row_mask.back());
    EXPECT_EQ(fit.regular_certificate,core::JointCheckStatus::NotRun);
    const auto invalid=core::FitJointComponents(problem,{.55,0,.55});
    EXPECT_FALSE(invalid.initialization.valid); ASSERT_EQ(invalid.components.size(),3);
    EXPECT_TRUE(invalid.components[0].state); EXPECT_FALSE(invalid.components[1].state);
    EXPECT_EQ(invalid.components[1].stop_reason,"invalid-initial-widths");
    EXPECT_FALSE(invalid.objective.has_value());
}
TEST(JointComponentRuntimeTest, TypedSearchPreservesFrozenAdapterResults)
{
    const core::JointProblem problem(Snapshot()); const auto & data=core::JointProblemAccess::Get(problem);
    rhbm_gem::eigen_helper::ScopedEigenThreadCount caller_threads{2};
    const auto before_threads=Eigen::nbThreads();
    const auto fit=core::FitJointComponents(problem,{.55,.55}); ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(Eigen::nbThreads(),before_threads);
    auto original=second_stage_test::matched::joint_abc::Fit(data.domain,data.y,n::Vector::Constant(2,.55),&data.context);
    const auto & beta=original.at("primary").at("beta").as_array();
    for(std::size_t k=0;k<beta.size();++k) EXPECT_NEAR(fit.assembled_state->ac[k],boost::json::value_to<double>(beta[k]),1e-8);
    EXPECT_TRUE(fit.search_completed); EXPECT_EQ(fit.components.size(),2);
    ASSERT_EQ(fit.components[0].ranks.size(),4);
    EXPECT_EQ(fit.components[0].ranks[1].rank,1);
    EXPECT_GT(fit.components[0].ranks[1].threshold,0);
    for(const auto & check:fit.components[0].evidence)
        if(check.name=="richardson") EXPECT_EQ(check.status,core::JointCheckStatus::NotRun);
    EXPECT_EQ(fit.regular_certificate,core::JointCheckStatus::NotRun);
}
TEST(JointComponentRuntimeTest, GeometrySelectionAndNegativeObservations)
{
    auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(1); atom->SetElement(Element::CARBON); atom->SetPosition(0,0,0);
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms; atoms.push_back(std::move(atom));
    rhbm_gem::ModelObject model(std::move(atoms)); model.SelectAllAtoms();
    auto values=std::make_unique<double[]>(11); for(std::size_t k=0;k<11;++k) values[k]=k%2 ? -1 : 0;
    rhbm_gem::MapObject map({11,1,1},{.5,.5,.5},{-2.5,0,0},std::move(values));
    const auto problem=core::BuildJointProblem(map,model);
    ASSERT_EQ(problem.Input().observations.size(),11);
    EXPECT_EQ(problem.Input().support[0].front().squared_distance,6.25);
    EXPECT_EQ(problem.Input().support[0].back().squared_distance,6.25);
    EXPECT_EQ(problem.Input().observations[1],-1);
    model.SetAtomSelected(1,false);
    EXPECT_THROW(core::BuildJointProblem(map,model),std::invalid_argument);
}

TEST(JointComponentRuntimeTest, FreshInitializationPreservesSecondStageAndSelection)
{
    auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(1); atom->SetElement(Element::CARBON);
    atom->SetPosition(0,0,0); atom->SetChainID("A"); atom->SetComponentID("ALA"); atom->SetAtomID("C");
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms; atoms.push_back(std::move(atom));
    auto hydrogen=std::make_unique<rhbm_gem::AtomObject>(); hydrogen->SetSerialID(2); hydrogen->SetElement(Element::HYDROGEN);
    hydrogen->SetPosition(1,0,0); atoms.push_back(std::move(hydrogen));
    rhbm_gem::ModelObject model(std::move(atoms)); model.SelectAllAtoms();
    rhbm_gem::LocalGaussianResult sentinel; sentinel.alpha_r=.789;
    sentinel.mdpde=rhbm_gem::GaussianModel3DWithUncertainty{rhbm_gem::GaussianModel3D{7.,.91,.23},{}};
    model.EditAnalysis().SetAtomLocalGaussianResult(rhbm_gem::FittingStage::Second,*model.FindAtomPtr(1),sentinel);
    rhbm_gem::MapObject map({25,25,25},{.3,.3,.3},{-3.6,-3.6,-3.6});
    core::simulation::SimulationAtomPreparationResult generator;
    generator.atom_list.push_back(core::simulation::SimulationAtom{.serial_id=1,.element=Element::CARBON,.position={0,0,0},.charge_used=.2});
    core::MapSimulationRequest request; request.job_count=1; request.cutoff_distance=2.5;
    request.potential_model_choice=core::PotentialModel::SINGLE_GAUS;
    core::simulation::PopulateMapValueArray(map,generator,request,.5);
    const auto fit=core::EstimateJointComponents(map,model);
    ASSERT_TRUE(fit.initialization.valid) << fit.initialization.reason;
    EXPECT_EQ(model.GetSelectedAtoms().size(),2);
    EXPECT_EQ(fit.initialization.atoms.size(),1);
    const auto view=rhbm_gem::AtomLocalPotentialView::For(*model.FindAtomPtr(1));
    EXPECT_EQ(view.GetAlphaR(rhbm_gem::FittingStage::Second),.789);
    EXPECT_EQ(view.GetEstimateMDPDE(rhbm_gem::FittingStage::Second).GetWidth(),.91);
    EXPECT_EQ(view.GetEstimateMDPDE(rhbm_gem::FittingStage::First).GetWidth(),fit.initialization.b[0]);
    EXPECT_GT(view.GetRawSamplingEntries(false).size(),0);
}

TEST(JointComponentRuntimeTest, InnerEvidenceRemainsAvailableNearBoundary)
{
    auto input=Snapshot();
    for(const auto & support:input.support) for(const auto & point:support)
        input.observations[point.row]=.2*n::EvaluateKernel(point.squared_distance,.5,2.5).charge;
    const auto fit=core::FitJointComponents(core::JointProblem(std::move(input)),{.5,.5});
    ASSERT_EQ(fit.components.size(),2);
    for(const auto & component:fit.components)
    {
        ASSERT_TRUE(component.state);
        EXPECT_NEAR(component.state->ac[0],0,1e-10);
        bool found_inner=false;
        for(const auto & check:component.evidence) if(check.name=="inner")
        {
            found_inner=true;
            EXPECT_EQ(check.status,core::JointCheckStatus::Passed);
            ASSERT_TRUE(check.value); ASSERT_TRUE(check.threshold);
            EXPECT_LE(*check.value,*check.threshold);
        }
        EXPECT_TRUE(found_inner);
        EXPECT_EQ(component.regular_certificate,core::JointCheckStatus::NotRun);
    }
}

TEST(JointComponentRuntimeTest, PublicObjectivesUseParentScaleAndIncludeConstantRows)
{
    auto input=Snapshot();
    for(std::size_t k=0;k<input.observations.size();++k) input.observations[k]+=.03*std::sin(static_cast<double>(k));
    input.row_ids.push_back("constant"); input.observations.push_back(7);
    const core::JointProblem problem(input); const double scale=problem.ObservationScale();
    ASSERT_GT(scale,1);
    const auto fit=core::FitJointComponents(problem,{.55,.55});
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

TEST(JointComponentRuntimeTest, AssemblyPreservesSuppliedCoefficientsInsteadOfReprofiling)
{
    const core::JointProblem problem(Snapshot());
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto initial=n::Vector::Constant(2,.55);
    std::vector<n::ComponentResult> components;
    for(const auto & view:data.partition.components)
        components.push_back(n::SolveComponent(view,data.y,initial,data.context));
    ASSERT_TRUE(components[0].trusted_state && components[1].trusted_state);
    components[0].trusted_state->beta(0)+=.1;
    const double actual=components[0].trusted_state->beta(0);
    const auto assembled=n::AssembleComponents(data.domain,data.y,data.partition,data.context,components);
    ASSERT_TRUE(assembled.available);
    EXPECT_DOUBLE_EQ(assembled.beta(0),actual);
    EXPECT_FALSE(assembled.profile_agrees);
    const auto replay=n::EvaluateState(data.domain,data.y,assembled.eta,assembled.beta,data.context);
    EXPECT_EQ(assembled.prediction,replay.x*assembled.beta);
}

TEST(JointComponentRuntimeTest, AssessmentWorkIsSharedOnlyForIdenticalScopes)
{
    auto input=Snapshot(); input.atom_ids.resize(1); input.support.resize(1);
    input.row_ids.resize(40); input.observations.resize(40);
    n::AssessmentWorkForTesting()={};
    const auto fit=core::FitJointComponents(core::JointProblem(input),{.55});
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(n::AssessmentWorkForTesting().assessments,1);
    EXPECT_EQ(fit.components[0].reference_evaluations,0);
    EXPECT_EQ(fit.costs.search_reference_seconds,0);
    EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,fit.components[0].reference_evaluations+1);
    input.row_ids.push_back("constant"); input.observations.push_back(7);
    n::AssessmentWorkForTesting()={};
    const auto constant=core::FitJointComponents(core::JointProblem(input),{.55});
    ASSERT_TRUE(constant.assembled_state);
    EXPECT_EQ(n::AssessmentWorkForTesting().assessments,2);
    EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,2);
    n::AssessmentWorkForTesting()={};
    const auto multiple=core::FitJointComponents(core::JointProblem(Snapshot()),{.55,.55});
    ASSERT_TRUE(multiple.assembled_state);
    EXPECT_EQ(n::AssessmentWorkForTesting().assessments,3);
    EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,3);
}

TEST(JointComponentRuntimeTest, SearchReplaysWithoutReferenceIncludingInitialAndDetailedTrials)
{
    const core::JointProblem problem(Snapshot()); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & view=data.partition.components[0]; const auto y=n::SelectValues(data.y,view.rows);
    for(bool details:{false,true}) for(int budget:{0,100})
    {
        auto context=n::ChildContext(data.context,view,true);
        context.audit.trial_details=details; context.update_budget=budget;
        n::AssessmentWorkForTesting()={};
        const auto search=n::SearchProfile(view.domain,y,n::Vector::Constant(1,.55),context);
        EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,0);
        EXPECT_EQ(search.references,0); EXPECT_EQ(search.reference_seconds,0);
        ASSERT_TRUE(search.initial_accepted);
        if(budget==0) {EXPECT_EQ(search.accepted,0); EXPECT_EQ(search.stop_reason,"accepted-update-budget");}
        else EXPECT_GT(search.accepted,1);
        for(const auto & trial:search.trials)
        {
            if(trial.accepted || (details && trial.lm)) ASSERT_TRUE(trial.trust);
            if(trial.accepted) EXPECT_TRUE(trial.trust->passed);
            if(!trial.trust) continue;
            EXPECT_FALSE(trial.trust->reference);
            const auto json=second_stage_test::matched::runtime_json::Trust(*trial.trust);
            EXPECT_TRUE(json.at("reference").is_null());
            EXPECT_TRUE(json.at("scaled_coefficient_difference").is_null());
        }
        const auto fit=n::AssessComponentSearch(view.domain,y,context,search);
        ASSERT_TRUE(fit.trusted_state); ASSERT_TRUE(fit.endpoint_trust);
        EXPECT_TRUE(fit.endpoint_trust->passed); EXPECT_TRUE(fit.endpoint_trust->reference);
        EXPECT_EQ(n::AssessmentWorkForTesting().reference_evaluations,1);
        if(budget==0) {EXPECT_FALSE(fit.search_success); EXPECT_EQ(fit.search.stop_reason,"accepted-update-budget");}
    }
}

TEST(JointComponentRuntimeTest, ReplayAndReferenceRejectIndependentDisagreements)
{
    const core::JointProblem problem(Snapshot()); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & view=data.partition.components[0]; const auto y=n::SelectValues(data.y,view.rows);
    const auto context=n::ChildContext(data.context,view,true);
    const auto e=n::EvaluateProfile(view.domain,y,n::Vector::Constant(1,std::log(.55)),false,&context);
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
    const auto reference=n::EvaluateProfile(view.domain,y,e.eta,true,&context);
    EXPECT_TRUE(n::CheckTrust(view.domain,y,e,context,reference).passed);
    bad=reference; bad.valid=false;
    EXPECT_FALSE(n::CheckTrust(view.domain,y,e,context,bad).passed);
    bad=reference; bad.beta(0)+=.01;
    EXPECT_FALSE(n::CheckTrust(view.domain,y,e,context,bad).passed);
    bad=reference; bad.gradient(0)+=.01;
    EXPECT_FALSE(n::CheckTrust(view.domain,y,e,context,bad).passed);
    bad=reference; bad.certificate.kkt_passed=false;
    EXPECT_FALSE(n::CheckTrust(view.domain,y,e,context,bad).passed);
}

TEST(JointComponentRuntimeTest, FallbackCertifiesSavedCoefficientsNewestFirst)
{
    const core::JointProblem problem(Snapshot()); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & view=data.partition.components[0]; const auto y=n::SelectValues(data.y,view.rows);
    const auto context=n::ChildContext(data.context,view,true);
    auto search=n::SearchProfile(view.domain,y,n::Vector::Constant(1,.55),context);
    std::vector<std::size_t> accepted;
    for(std::size_t k=0;k<search.trials.size();++k) if(search.trials[k].accepted) accepted.push_back(k);
    ASSERT_GT(accepted.size(),1);
    const auto expected=search.trials[accepted[accepted.size()-2]].endpoint;
    search.trials[accepted.back()].endpoint.beta(0)+=1;
    search.eta=n::Vector::Constant(1,1000); search.stopped=false;
    const auto fit=n::AssessComponentSearch(view.domain,y,context,search);
    ASSERT_TRUE(fit.trusted_state && fit.trusted_assessment);
    EXPECT_EQ(fit.trusted_trial,accepted[accepted.size()-2]);
    EXPECT_EQ(fit.trusted_state->eta,expected.eta); EXPECT_EQ(fit.trusted_state->beta,expected.beta);
    EXPECT_EQ(fit.trusted_assessment->primary.beta,expected.beta);
    EXPECT_LE(fit.trusted_assessment->coefficient_difference,1e-10);
    EXPECT_FALSE(fit.search_success); EXPECT_TRUE(fit.search.stopped);
    EXPECT_EQ(fit.search.stop_reason,"endpoint-certification-failed");
    EXPECT_FALSE(fit.assessment.primary.valid);
    for(auto k:accepted) search.trials[k].endpoint.beta(0)=-1;
    const auto failed=n::AssessComponentSearch(view.domain,y,context,search);
    EXPECT_FALSE(failed.trusted_state); EXPECT_FALSE(failed.trusted_assessment);
    EXPECT_FALSE(failed.search_success);
    // Even the initial state is eligible, but it must pass fresh certification.
    search.trials[accepted.front()].endpoint=search.initial;
    search.stopped=true; search.stop_reason="profile-budget";
    const auto initial=n::AssessComponentSearch(view.domain,y,context,search);
    ASSERT_TRUE(initial.trusted_state && initial.trusted_assessment);
    EXPECT_EQ(initial.trusted_trial,accepted.front());
    EXPECT_EQ(initial.trusted_state->beta,search.initial.beta);
    EXPECT_EQ(initial.search.stop_reason,"profile-budget"); EXPECT_FALSE(initial.search_success);
}

TEST(JointComponentRuntimeTest, ReusedAssemblyMatchesFreshAssessmentAndRejectsChangedInputs)
{
    auto input=Snapshot(); input.atom_ids.resize(1); input.support.resize(1);
    input.row_ids.resize(40); input.observations.resize(40);
    const core::JointProblem problem(input); const auto & data=core::JointProblemAccess::Get(problem);
    auto context=n::ChildContext(data.context,data.partition.components[0],true);
    std::vector<n::ComponentResult> fits{n::SolveComponent(data.partition.components[0],data.y,n::Vector::Constant(1,.55),data.context)};
    ASSERT_TRUE(fits[0].trusted_assessment);
    const n::AssessmentReuse reuse{data.domain,data.y,context,*fits[0].trusted_assessment};
    const auto check=[&](const n::Domain & domain,const n::Vector & y,const n::EvaluationContext & policy,int assessments) {
        n::AssessmentWorkForTesting()={};
        const auto assembled=n::AssembleComponents(domain,y,data.partition,policy,fits,&reuse);
        EXPECT_EQ(n::AssessmentWorkForTesting().assessments,assessments);
        const auto fresh=n::AssessProfile(domain,y,assembled.eta,policy,&assembled.beta);
        EXPECT_EQ(second_stage_test::matched::runtime_json::Assessment(assembled.assessment),
            second_stage_test::matched::runtime_json::Assessment(fresh));
    };
    check(data.domain,data.y,data.context,0);
    auto policy=data.context; policy.rank.rows+=100; check(data.domain,data.y,policy,1);
    policy=data.context; policy.scale*=2; check(data.domain,data.y,policy,1);
    policy=data.context; policy.linear.release_factor*=2; check(data.domain,data.y,policy,1);
    policy=data.context; policy.audit.directions=n::Matrix::Ones(1,3); check(data.domain,data.y,policy,0);
    n::Vector y=data.y; y(0)+=.01; check(data.domain,y,data.context,1);
    auto support=data.domain.CopySupport(); support[0][1].square+=.001; n::Domain domain(data.domain.rows,std::move(support)); check(domain,data.y,data.context,1);
    fits[0].trusted_state->beta(0)+=.1; check(data.domain,data.y,data.context,1);
}

TEST(JointComponentRuntimeTest, ActualAssessmentFollowsFallbackAndLocalDirections)
{
    const core::JointProblem problem(Snapshot()); const auto & data=core::JointProblemAccess::Get(problem);
    const auto & view=data.partition.components[0];
    auto context=n::ChildContext(data.context,view,true); const auto y=n::SelectValues(data.y,view.rows);
    auto search=n::SearchProfile(view.domain,y,n::Vector::Constant(1,.55),context);
    ASSERT_TRUE(search.trials.back().accepted);
    const auto accepted=search.trials.back().endpoint;
    search.eta=n::Vector::Constant(1,1000); search.stopped=true; search.stop_reason="untrusted-trial";
    const auto fallback=n::AssessComponentSearch(view.domain,y,context,search);
    ASSERT_TRUE(fallback.trusted_state && fallback.trusted_assessment);
    EXPECT_FALSE(fallback.search_success); EXPECT_FALSE(fallback.assessment.primary.valid);
    EXPECT_EQ(fallback.trusted_state->eta,accepted.eta);
    EXPECT_EQ(fallback.trusted_assessment->primary.eta,accepted.eta);
    EXPECT_EQ(fallback.trusted_assessment->primary.beta,accepted.beta);
    context.audit.directions=n::Matrix::Zero(1,3);
    const auto local=n::SolveComponent(view,data.y,n::Vector::Constant(2,.55),data.context);
    const auto directed=n::AssessComponentSearch(view.domain,y,context,local.search);
    ASSERT_TRUE(directed.trusted_assessment);
    const auto fresh=n::AssessProfile(view.domain,y,directed.trusted_state->eta,
        n::ChildContext(data.context,view,true),&directed.trusted_state->beta);
    EXPECT_EQ(directed.trusted_assessment->primary.beta,fresh.primary.beta);
    EXPECT_EQ(directed.trusted_assessment->correction,fresh.correction);
    EXPECT_EQ(directed.trusted_assessment->weak_directions,fresh.weak_directions);
}

TEST(JointComponentRuntimeTest, RuntimeConvergenceIsDerivedFromActualEvidenceNotTermination)
{
    using Status=core::JointCheckStatus;
    auto fit=core::FitJointComponents(core::JointProblem(Snapshot()),{.55,.55});
    ASSERT_EQ(fit.RuntimeConvergence(),Status::Passed);
    for(auto & component:fit.components)
    {
        component.search_completed=false; component.stop_reason="profile-budget";
        EXPECT_EQ(component.RuntimeConvergence(),Status::Passed);
        for(const auto & check:component.evidence)
            if(check.name=="two-step-derivative" || check.name=="richardson" || check.name=="precision-50-100" ||
                check.name=="boundary-audit" || check.name=="regular-certificate")
                EXPECT_EQ(check.status,Status::NotRun);
    }
    fit.search_completed=false;
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Passed);
    EXPECT_EQ(fit.regular_certificate,Status::NotRun);
    auto & evidence=fit.components[0].evidence;
    evidence[0].status=Status::NotRun;
    EXPECT_EQ(fit.RuntimeConvergence(),Status::NotRun);
    evidence[1].status=Status::Unavailable;
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Unavailable);
    evidence[2].status=Status::Failed;
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Failed);
    fit.components[0].state.reset();
    EXPECT_EQ(fit.components[0].RuntimeConvergence(),Status::Unavailable);
    EXPECT_EQ(fit.RuntimeConvergence(),Status::Unavailable);
}

TEST(JointComponentRuntimeTest, MissingOrWrongScopeEvidenceCannotVacuouslyConverge)
{
    using Status=core::JointCheckStatus;
    const auto fit=core::FitJointComponents(core::JointProblem(Snapshot()),{.55,.55});
    auto component=fit.components[0];
    component.evidence.clear(); EXPECT_EQ(component.RuntimeConvergence(),Status::Unavailable);
    component.evidence=fit.evidence; EXPECT_EQ(component.RuntimeConvergence(),Status::Unavailable);
    auto global=fit;
    for(auto & check:global.evidence) if(check.name=="assembled-profile") check.status=Status::Failed;
    EXPECT_EQ(global.RuntimeConvergence(),Status::Failed);
    global=fit; global.evidence=fit.components[0].evidence;
    EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    global=fit; global.prediction.reset(); EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    global=fit; global.available_row_mask[0]=false; EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    global=fit; global.components.clear(); EXPECT_EQ(global.RuntimeConvergence(),Status::Unavailable);
    EXPECT_EQ(core::JointFitResult{}.RuntimeConvergence(),Status::Unavailable);
}

TEST(JointRuntimeTest, SnapshotViewsRetainStorageWithoutDuplicatingObservationsOrMemberships)
{
    std::optional<n::ComponentView> retained;
    std::optional<n::EvaluationContext> context;
    const double * observations=nullptr;
    {
        auto input=Snapshot(); core::JointProblem original(input);
        core::JointProblem copy=original; core::JointProblem moved=std::move(copy);
        const auto & data=core::JointProblemAccess::Get(moved);
        observations=moved.Input().observations.data();
        EXPECT_EQ(data.y.data(),observations);
        EXPECT_EQ(data.context.observations->data(),observations);
        EXPECT_EQ(data.domain.atoms.Storage(),&moved.Input());
        EXPECT_EQ(data.partition.components[0].domain.atoms.Storage(),&moved.Input());
        EXPECT_EQ(data.partition.components[0].mappings,data.partition.components[1].mappings);
        EXPECT_EQ(&data.context.atom_ids[0],&moved.Input().atom_ids[0]);
        EXPECT_EQ(data.partition.mappings->row_to_local.size(),input.observations.size());
        retained=data.partition.components[0]; context=n::ChildContext(data.context,*retained,true);
        input.observations[0]=999; input.support[0].clear();
    }
    ASSERT_TRUE(retained && context);
    EXPECT_EQ(context->observations->data(),observations);
    EXPECT_NE((*context->observations)(0),999);
    EXPECT_EQ(retained->domain.atoms[0].size(),40);
    const n::Vector y=n::SelectValues(*context->observations,retained->rows);
    EXPECT_TRUE(n::EvaluateProfile(retained->domain,y,n::Vector::Constant(1,std::log(.55)),false,&*context).valid);
}

TEST(JointRuntimeTest, ManyComponentsShareOneParentMapping)
{
    auto input=Snapshot();
    const auto base=input;
    for(std::size_t block=1;block<64;++block)
    {
        const auto offset=input.observations.size();
        for(std::size_t a=0;a<base.support.size();++a)
        {
            input.atom_ids.push_back(std::to_string(block)+"/"+base.atom_ids[a]);
            auto support=base.support[a]; for(auto & s:support) s.row+=offset;
            input.support.push_back(std::move(support));
        }
        for(std::size_t r=0;r<base.observations.size();++r)
        {input.observations.push_back(base.observations[r]); input.row_ids.push_back(std::to_string(block)+"/"+base.row_ids[r]);}
    }
    const core::JointProblem problem(std::move(input)); const auto & data=core::JointProblemAccess::Get(problem);
    ASSERT_EQ(data.partition.components.size(),128);
    for(const auto & c:data.partition.components)
    {
        EXPECT_EQ(c.mappings,data.partition.mappings);
        EXPECT_EQ(c.domain.atoms.Storage(),&problem.Input());
        for(std::size_t a=0;a<c.atoms.size();++a) EXPECT_EQ(c.LocalAtom(c.atoms[a]),static_cast<Eigen::Index>(a));
        for(std::size_t r=0;r<c.rows.size();++r) EXPECT_EQ(c.LocalRow(c.rows[r]),static_cast<Eigen::Index>(r));
    }
}

TEST(JointComponentRuntimeTest, SavedOutcomesPreserveActualResultsWithoutProblemOrPrediction)
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
        const auto fit=core::FitJointComponents(core::JointProblem(input),initial);
        const auto saved=core::CaptureJointAnalysisResult(fit);
        ASSERT_TRUE(saved.metadata.software);
        EXPECT_FALSE(saved.metadata.model_sha256); EXPECT_FALSE(saved.metadata.map_sha256);
        EXPECT_FALSE(saved.metadata.map_normalization);
        EXPECT_EQ(saved.metadata.software->source_sha256.size(),64);
        const auto encoded=io::Encode(saved);
        const auto decoded=io::Decode(encoded);
        EXPECT_EQ(io::Encode(decoded),encoded);
        EXPECT_EQ(decoded.runtime_convergence,fit.RuntimeConvergence());
        EXPECT_EQ(decoded.available_row_mask,fit.available_row_mask);
        EXPECT_EQ(decoded.objective,fit.objective);
        EXPECT_EQ(decoded.observation_scale,fit.observation_scale);
        ASSERT_EQ(decoded.components.size(),fit.components.size());
        for(std::size_t k=0;k<fit.components.size();++k)
        {
            EXPECT_EQ(decoded.components[k].runtime_convergence,fit.components[k].RuntimeConvergence());
            EXPECT_EQ(decoded.components[k].regular_certificate,core::JointCheckStatus::NotRun);
            EXPECT_EQ(decoded.components[k].state.has_value(),fit.components[k].state.has_value());
            if(fit.components[k].state) EXPECT_EQ(decoded.components[k].state->ac,fit.components[k].state->ac);
        }
        EXPECT_EQ(encoded.find("\"prediction\""),std::string::npos);
        EXPECT_EQ(encoded.find("\"observations\""),std::string::npos);
        if(variant==3) EXPECT_TRUE(std::isnan(decoded.initialization.b[0]));
    }
}

TEST(JointComponentRuntimeTest, PositiveMapScalingPreservesKernelCoefficientUnits)
{
    auto input=Snapshot();
    const auto original=core::FitJointComponents(core::JointProblem(input),{.55,.55});
    constexpr double divisor=3.25;
    for(auto & y:input.observations) y/=divisor;
    const auto normalized=core::FitJointComponents(core::JointProblem(input),{.55,.55});
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
