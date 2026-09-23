#include <gtest/gtest.h>
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include "support/JointPartialSelection.hpp"
#include "support/DataObjectTestSupport.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include <map>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include "support/CommandTestHelpers.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "core/detail/StageSummary.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/JointUncertainty.hpp"
#include "core/detail/MapInterpolation.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <set>

namespace {
namespace core=rhbm_gem::core;
namespace n=core::joint_component;
namespace sim=core::simulation;
// Exhaust the whole map/catalogue independently of the builder's bounding boxes.
void CheckReference(const rhbm_gem::MapObject & map,const rhbm_gem::ModelObject & model)
{
    const auto problem=core::BuildJointProblem(map,model); const auto & actual=problem.Input();
    std::set<const rhbm_gem::AtomObject *> selected(model.GetSelectedAtoms().begin(),model.GetSelectedAtoms().end());
    std::vector<std::string> rows,ids; std::vector<double> observations;
    std::vector<std::vector<core::JointSupport>> support(model.GetAtomList().size());
    const auto dims=map.GetGridSize();
    for(int z=0;z<dims[2];++z) for(int y=0;y<dims[1];++y) for(int x=0;x<dims[0];++x)
    {
        const auto p=sim::GridPosition({x,y,z},map.GetGridSpacing(),map.GetOrigin());
        bool target=false;
        for(const auto * atom:selected) if(atom->GetElement()!=Element::HYDROGEN && sim::SupportSquare(p,atom->GetPosition())<=6.25) target=true;
        if(!target) continue;
        const auto index=static_cast<std::size_t>(x+dims[0]*(y+dims[1]*z));
        for(std::size_t a=0;a<support.size();++a)
        {
            const auto & atom=model.GetAtomList()[a];
            const auto square=sim::SupportSquare(p,atom->GetPosition());
            if(atom->GetElement()!=Element::HYDROGEN && square<=6.25) support[a].push_back({rows.size(),square});
        }
        rows.push_back(std::to_string(index)); observations.push_back(map.GetMapValue(index));
    }
    std::vector<std::size_t> targets;
    for(std::size_t a=0;a<support.size();++a)
    {
        const auto & atom=model.GetAtomList()[a];
        if(atom->GetElement()==Element::HYDROGEN || (!selected.contains(atom.get()) && support[a].empty())) continue;
        const auto k=ids.size(); ids.push_back(std::to_string(atom->GetSerialID()));
        if(selected.contains(atom.get())) targets.push_back(k);
        ASSERT_LT(k,actual.support.size()); ASSERT_EQ(support[a].size(),actual.support[k].size());
        for(std::size_t i=0;i<support[a].size();++i)
        {
            EXPECT_EQ(support[a][i].row,actual.support[k][i].row);
            EXPECT_DOUBLE_EQ(support[a][i].squared_distance,actual.support[k][i].squared_distance);
        }
    }
    EXPECT_EQ(actual.row_ids,rows); EXPECT_EQ(actual.observations,observations); EXPECT_EQ(actual.atom_ids,ids);
    ASSERT_TRUE(actual.selection_domain); EXPECT_EQ(actual.selection_domain->target_indices,targets);
    double square=0; for(double value:observations) square+=value*value;
    EXPECT_NEAR(problem.ObservationScale(),std::max(1.,std::sqrt(square)),1e-12*problem.ObservationScale());
}
}
TEST(JointComponentPartialSelectionTest, FullCatalogueReferenceAndFixedRows)
{
    auto f=joint_partial_test::Make("partial"); CheckReference(*f.map,*f.model);
    const auto partial=core::BuildJointProblem(*f.map,*f.model);
    f.model->SelectAllAtoms(); CheckReference(*f.map,*f.model);
    const auto all=core::BuildJointProblem(*f.map,*f.model);
    EXPECT_GT(all.Input().row_ids.size(),partial.Input().row_ids.size());
    f.model->SelectAllAtoms(false); EXPECT_THROW(core::BuildJointProblem(*f.map,*f.model),std::invalid_argument);
}
TEST(JointComponentPartialSelectionTest, BoundariesClippingAndNonrecursiveHalo)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for(int i=0;i<6;++i)
    {
        auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(i+1);
        atom->SetElement(i==5 ? Element::HYDROGEN : Element::CARBON);
        atom->SetPosition(std::array<double,6>{0,4.8,7.2,30,-2.5}[static_cast<std::size_t>(i)],0,0);
        atoms.push_back(std::move(atom));
    }
    rhbm_gem::ModelObject model(std::move(atoms)); model.SelectAtoms([](const auto & a){return a.GetSerialID()==1 || a.GetSerialID()==4 || a.GetSerialID()==6;});
    auto values=std::make_unique<double[]>(11*3*3); for(std::size_t i=0;i<99;++i) values[i]=i%2 ? -1:0;
    rhbm_gem::MapObject map({11,3,3},{.25,.7,1.1},{0,-.7,-1.1},std::move(values));
    CheckReference(map,model);
    const auto p=core::BuildJointProblem(map,model);
    EXPECT_EQ(p.Input().atom_ids,(std::vector<std::string>{"1","2","4","5"}));
    EXPECT_TRUE(p.Input().support[2].empty());
    EXPECT_DOUBLE_EQ(p.Input().support[3][0].squared_distance,6.25);
    const auto fit=core::FitJointComponents(p,{.5,.5,.5,.5});
    EXPECT_FALSE(fit.assembled_state);
    const auto & data=core::JointProblemAccess::Get(p);
    EXPECT_EQ(data.partition.unobserved_atoms.size(),1);
}
TEST(JointComponentPartialSelectionTest, BridgeAndMonolithicEquivalence)
{
    auto f=joint_partial_test::Make("bridge"); CheckReference(*f.map,*f.model);
    const auto p=core::BuildJointProblem(*f.map,*f.model); const auto & data=core::JointProblemAccess::Get(p);
    ASSERT_EQ(data.partition.components.size(),1);
    auto without=p.Input(); without.selection_domain.reset(); without.atom_ids.pop_back(); without.support.pop_back();
    const core::JointProblem disconnected(without);
    EXPECT_EQ(core::JointProblemAccess::Get(disconnected).partition.components.size(),2);
    const auto fit=core::FitJointComponents(p,f.b);
    const n::Vector widths=Eigen::Map<n::Vector>(f.b.data(),static_cast<Eigen::Index>(f.b.size()));
    const auto mono=n::AssessComponentSearch(data.domain,data.y,data.context,n::SearchProfile(data.domain,data.y,widths,data.context));
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(mono.trusted_state);
    for(std::size_t i=0;i<fit.assembled_state->ac.size();++i)
        EXPECT_NEAR(fit.assembled_state->ac[i],mono.trusted_state->beta(static_cast<Eigen::Index>(i)),1e-10*(1+std::abs(fit.assembled_state->ac[i])));
    EXPECT_NEAR(*fit.objective,mono.trusted_state->certificate.objective/std::pow(p.ObservationScale(),2),1e-12);
}
TEST(JointComponentPartialSelectionTest, IdentifiableRecoveryOmissionAndOutsideDomain)
{
    auto f=joint_partial_test::Make("partial"); const auto p=core::BuildJointProblem(*f.map,*f.model);
    const auto fit=core::FitJointComponents(p,f.b);
    ASSERT_EQ(fit.RuntimeConvergence(),core::JointCheckStatus::Passed); ASSERT_TRUE(fit.assembled_state);
    for(std::size_t i=0;i<f.b.size();++i)
    {
        EXPECT_NEAR(fit.assembled_state->ac[2*i],f.a[i],1e-10*(1+f.a[i]));
        EXPECT_NEAR(fit.assembled_state->ac[2*i+1],f.c[i],1e-10*(1+f.c[i]));
        EXPECT_NEAR(fit.assembled_state->b[i],f.b[i],1e-10*(1+f.b[i]));
    }
    auto omitted=p.Input(); omitted.atom_ids.resize(1); omitted.support.resize(1);
    const auto wrong=core::FitJointComponents(core::JointProblem(omitted),{f.b[0]});
    ASSERT_TRUE(wrong.objective); ASSERT_TRUE(wrong.assembled_state);
    EXPECT_GT(*wrong.objective,*fit.objective+1e-8);
    EXPECT_GT(std::abs(wrong.assembled_state->ac[0]-f.a[0]),1e-5);
    std::set<std::string> rows(p.Input().row_ids.begin(),p.Input().row_ids.end());
    auto values=std::make_unique<double[]>(f.map->GetMapValueArraySize());
    for(std::size_t i=0;i<f.map->GetMapValueArraySize();++i) values[i]=rows.contains(std::to_string(i)) ? f.map->GetMapValue(i) : 12345;
    f.map->SetMapValueArray(std::move(values));
    const auto changed=core::BuildJointProblem(*f.map,*f.model);
    EXPECT_EQ(changed.Input().observations,p.Input().observations);
    const auto repeated=core::FitJointComponents(changed,f.b);
    EXPECT_EQ(repeated.assembled_state->ac,fit.assembled_state->ac); EXPECT_EQ(repeated.objective,fit.objective);
}
TEST(JointComponentPartialSelectionTest, WeakHaloAndInvalidWidthsAreNotSuccess)
{
    auto f=joint_partial_test::Make("weak"); const auto p=core::BuildJointProblem(*f.map,*f.model);
    ASSERT_EQ(p.Input().atom_ids.size(),2);
    for(const auto & point:p.Input().support[1]) EXPECT_GT(point.squared_distance,5.0);
    const auto weak=core::FitJointComponents(p,{.5,.05});
    EXPECT_NE(weak.RuntimeConvergence(),core::JointCheckStatus::Passed);
    EXPECT_FALSE(weak.components[0].state);
    for(double invalid:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()})
    {
        const auto fit=core::FitJointComponents(p,{.5,invalid});
        ASSERT_EQ(fit.components.size(),1); EXPECT_FALSE(fit.initialization.valid);
        EXPECT_EQ(fit.components[0].stop_reason,"invalid-initial-widths"); EXPECT_FALSE(fit.prediction);
    }
    EXPECT_TRUE(core::FitJointComponents(p,{.5}).components.empty());
}
TEST(JointComponentPartialSelectionTest, InitializationPreservesHaloHistoryAndRoles)
{
    auto f=joint_partial_test::Make("partial");
    f.model->EditAnalysis().InitializeFromSelection();
    f.model->EditAnalysis().InitializeGroupAlpha(.432);
    const auto groups=f.model->GetAnalysisView().CollectAtomGroupKeys();
    ASSERT_FALSE(groups.empty());
    const auto selected_bonds=f.model->GetSelectedBonds();
    rhbm_gem::LocalGaussianResult sentinel; sentinel.alpha_r=.789;
    sentinel.mdpde=rhbm_gem::GaussianModel3DWithUncertainty{rhbm_gem::GaussianModel3D{7.,.91,.23},{}};
    for(const auto & atom:f.model->GetAtomList())
        for(auto stage:{rhbm_gem::FittingStage::First,rhbm_gem::FittingStage::Second})
            f.model->EditAnalysis().SetAtomLocalGaussianResult(stage,*atom,sentinel);
    const auto fit=core::EstimateJointComponents(*f.map,*f.model);
    ASSERT_EQ(fit.initialization.atoms.size(),2); ASSERT_TRUE(fit.initialization.valid) << fit.initialization.reason;
    EXPECT_EQ(f.model->GetSelectedAtomCount(),1); EXPECT_EQ(f.model->GetSelectedAtoms()[0]->GetSerialID(),1);
    EXPECT_EQ(f.model->GetSelectedBonds(),selected_bonds);
    EXPECT_EQ(f.model->GetAnalysisView().CollectAtomGroupKeys(),groups);
    for(auto group:groups) EXPECT_EQ(f.model->GetAnalysisView().GetAtomAlphaG(group),.432);
    const auto halo=rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(2));
    EXPECT_EQ(halo.GetEstimateMDPDE(rhbm_gem::FittingStage::First).GetWidth(),.91);
    for(const auto & atom:f.model->GetAtomList())
        EXPECT_EQ(rhbm_gem::AtomLocalPotentialView::For(*atom).GetEstimateMDPDE(rhbm_gem::FittingStage::Second).GetWidth(),.91);
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto encoded=rhbm_gem::joint_result_io::Encode(saved);
    const auto decoded=rhbm_gem::joint_result_io::Decode(encoded);
    ASSERT_TRUE(decoded.selection_domain); EXPECT_EQ(decoded.selection_domain->target_indices,(std::vector<std::size_t>{0}));
    EXPECT_EQ(decoded.initialization.data_scope,"contributor-local-sampling-may-read-outside-target-domain");
    EXPECT_EQ(rhbm_gem::joint_result_io::Encode(decoded),encoded);
}

TEST(JointComponentPartialSelectionTest, InitializationExceptionsRemainLocal)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for(int i=0;i<2;++i)
    {
        auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(i+1); atom->SetElement(Element::CARBON);
        atom->SetPosition(9.*i,0,0); atom->SetChainID("A"); atom->SetComponentID("ALA"); atom->SetAtomID("CA");
        atoms.push_back(std::move(atom));
    }
    rhbm_gem::ModelObject model(std::move(atoms)); model.SelectAllAtoms();
    const std::array<int,3> dims{13,7,7}; const std::array<double,3> spacing{1.5,1.5,1.5},origin{-3,-4.5,-4.5};
    auto values=std::make_unique<double[]>(13*7*7);
    for(int z=0;z<7;++z) for(int y=0;y<7;++y) for(int x=0;x<13;++x)
    {
        const auto p=sim::GridPosition({x,y,z},spacing,origin);
        double value=0;
        for(const auto & atom:model.GetAtomList())
        {
            const auto kernel=n::EvaluateKernel(sim::SupportSquare(p,atom->GetPosition()),1.,2.5);
            value+=2*kernel.gaussian+.2*kernel.charge;
        }
        values[static_cast<std::size_t>(x+13*(y+7*z))]=value;
    }
    // Outside both target supports, but inside the first atom's cubic sampling stencil.
    values[4+13*(3+7*3)]=std::numeric_limits<double>::quiet_NaN();
    rhbm_gem::MapObject map(dims,spacing,origin,std::move(values));
    const auto fit=core::EstimateJointComponents(map,model);
    ASSERT_EQ(fit.initialization.atoms.size(),2);
    EXPECT_TRUE(fit.initialization.atoms[0].reason.starts_with("initialization-exception:"));
    EXPECT_EQ(fit.initialization.atoms[1].reason,"valid-width");
    ASSERT_EQ(fit.components.size(),2); EXPECT_FALSE(fit.components[0].state); EXPECT_TRUE(fit.components[1].state);
    EXPECT_FALSE(fit.initialization.valid); EXPECT_FALSE(fit.objective);
    EXPECT_EQ(model.GetSelectedAtomCount(),2);
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_EQ(decoded.initialization.atoms[0].reason,fit.initialization.atoms[0].reason);
    EXPECT_EQ(decoded.available_row_mask,fit.available_row_mask);
}

TEST(JointComponentPartialSelectionTest, ZeroSignalRankDeficiencyAndInsufficientRows)
{
    auto f=joint_partial_test::Make("partial"); const auto p=core::BuildJointProblem(*f.map,*f.model);
    for(int scenario=0;scenario<3;++scenario)
    {
        auto input=p.Input();
        if(scenario==0) std::fill(input.observations.begin(),input.observations.end(),0);
        if(scenario==1) input.support[1]=input.support[0];
        if(scenario==2)
        {
            input.row_ids.resize(2); input.observations.resize(2);
            for(auto & support:input.support) support={{0,0},{1,1}};
        }
        const auto fit=core::FitJointComponents(core::JointProblem(input),{.5,.5});
        EXPECT_NE(fit.RuntimeConvergence(),core::JointCheckStatus::Passed);
        EXPECT_EQ(fit.regular_certificate,core::JointCheckStatus::NotRun);
        if(scenario>0) EXPECT_FALSE(fit.assembled_state);
    }
}

TEST(JointComponentPartialSelectionTest, FullSelectionInitializerMatchesV1Workflow)
{
    auto f=joint_partial_test::Make("all");
    rhbm_gem::ModelObject original(*f.model);
    original.ApplyElementSelection(Element::HYDROGEN,true);
    original.EditAnalysis().InitializeFromSelection();
    core::RunPotentialSamplingWorkflow(*f.map,original,SphereSamplingMethod::FibonacciDeterministic,1);
    original.EditAnalysis().InitializeLocalFittingSeedModels();
    core::FitOptions options; options.thread_size=1; options.quiet_mode=true; options.exclude_hydrogen=true;
    core::RunLocalAlphaTraining(original,options,rhbm_gem::FittingStage::First);
    core::RunFixedOffsetLocalFitting(original,options,rhbm_gem::FittingStage::First);
    const auto fit=core::EstimateJointComponents(*f.map,*f.model);
    ASSERT_TRUE(fit.initialization.valid);
    for(std::size_t a=0;a<fit.initialization.atoms.size();++a)
    {
        const auto & atom=fit.initialization.atoms[a];
        const auto view=rhbm_gem::AtomLocalPotentialView::For(*original.FindAtomPtr(std::stoi(atom.id)));
        EXPECT_DOUBLE_EQ(fit.initialization.b[a],view.GetEstimateMDPDE(rhbm_gem::FittingStage::First).GetWidth());
        EXPECT_DOUBLE_EQ(atom.alpha,view.GetAlphaR(rhbm_gem::FittingStage::First));
        EXPECT_EQ(atom.sample_count,view.GetSamplingEntries(rhbm_gem::FittingStage::First).size());
    }
}

TEST(JointComponentPartialSelectionTest, OriginalBondSelectionIsPreserved)
{
    auto fixture=joint_partial_test::Make("partial");
    auto model=data_test::MakeModelWithBond();
    for(const auto & atom:model->GetAtomList())
    {
        atom->SetElement(Element::CARBON); atom->SetChainID("A"); atom->SetComponentID("ALA"); atom->SetAtomID("CA");
    }
    model->SelectAllAtoms(); model->SetAtomSelected(2,false);
    ASSERT_EQ(model->GetNumberOfBond(),1);
    for(bool selected:{true,false})
    {
        model->SelectAllBonds(selected);
        const auto before=model->GetSelectedBonds();
        const auto fit=core::EstimateJointComponents(*fixture.map,*model);
        EXPECT_EQ(model->GetSelectedBonds(),before);
        EXPECT_EQ(model->GetSelectedAtomCount(),1);
        EXPECT_EQ(fit.initialization.atoms.size(),2);
    }
}

TEST(JointComponentPartialSelectionTest, SharedWorkflowFitsEachContributorOnceWithoutChangingSelection)
{
    auto f = joint_partial_test::Make("partial");
    const auto selected = f.model->GetSelectedAtoms();
    const auto bonds = f.model->GetSelectedBonds();
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    std::map<std::pair<int, std::string>, int> calls;
    core::detail::FirstStageObserverForTesting() = [&](int id, std::string_view phase) {
        ++calls[{id, std::string(phase)}];
    };
    core::FitOptions options;
    options.estimator = core::PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode = true;
    core::RunPotentialFittingWorkflow(*f.map, *f.model, options);
    core::detail::FirstStageObserverForTesting() = {};
    ASSERT_TRUE(f.model->GetAnalysisView().GetJointResult());
    const auto & result = *f.model->GetAnalysisView().GetJointResult();
    EXPECT_EQ(f.model->GetSelectedAtoms(), selected);
    EXPECT_EQ(f.model->GetSelectedBonds(), bonds);
    EXPECT_EQ(result.row_ids, problem.Input().row_ids);
    ASSERT_EQ(result.initialization.atoms.size(), problem.Input().atom_ids.size());
    for (std::size_t i = 0; i < result.atom_ids.size(); ++i)
    {
        const auto id = std::stoi(result.atom_ids[i]);
        EXPECT_EQ((calls[{id, "raw"}]), 1);
        EXPECT_EQ((calls[{id, "first"}]), 1);
        const auto view = rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(id));
        EXPECT_DOUBLE_EQ(view.GetFinalModel(rhbm_gem::FittingStage::First).GetWidth(), result.initialization.b[i]);
        EXPECT_FALSE(view.GetRawSamplingEntries(false).empty());
        ASSERT_TRUE(view.GetPostFitPeeling());
        EXPECT_EQ(view.GetPostFitPeeling()->samples.size(), view.GetRawSamplingEntries(false).size());
        EXPECT_FALSE(view.GetGroupMemberResult());
    }
    const auto direct = core::FitJointComponents(problem, result.initialization.b);
    ASSERT_EQ(result.components.size(), direct.components.size());
    for (std::size_t c = 0; c < direct.components.size(); ++c)
    {
        ASSERT_EQ(result.components[c].state.has_value(), direct.components[c].state.has_value());
        if (!direct.components[c].state) continue;
        EXPECT_EQ(result.components[c].state->ac, direct.components[c].state->ac);
        EXPECT_EQ(result.components[c].state->b, direct.components[c].state->b);
        EXPECT_EQ(result.components[c].state->objective, direct.components[c].state->objective);
    }
    EXPECT_THROW(core::RunPotentialFittingWorkflow(*f.model, options), std::invalid_argument);
}

TEST(JointComponentPartialSelectionTest, StageAdapterUsesIdentityAndClearsMissingStates)
{
    auto f = joint_partial_test::Make("all");
    auto result = core::CaptureJointAnalysisResult(core::FitJointComponents(core::BuildJointProblem(*f.map, *f.model), f.b));
    ASSERT_EQ(result.atom_ids.size(), 2);
    std::swap(result.atom_ids[0], result.atom_ids[1]);
    auto & component = result.components.front();
    ASSERT_TRUE(component.state);
    component.atoms = {1, 0};
    result.selection_domain->target_indices = {1};
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, result, "test-run");
    const auto target = rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(1));
    const auto halo = rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(2));
    EXPECT_DOUBLE_EQ(target.GetFinalModel(rhbm_gem::FittingStage::Second).GetAmplitude(), component.state->ac[0]);
    EXPECT_DOUBLE_EQ(target.GetFinalModel(rhbm_gem::FittingStage::Second).GetOffset(), component.state->ac[1]);
    EXPECT_DOUBLE_EQ(target.GetFinalModel(rhbm_gem::FittingStage::Second).GetWidth(), component.state->b[0]);
    EXPECT_EQ(halo.GetStageEstimate(rhbm_gem::FittingStage::Second).source.role, rhbm_gem::FittingRole::Halo);
    auto summary = core::BuildSecondStageSpotSummary(*f.model);
    EXPECT_NE(summary.find("joint-components"), std::string::npos);
    EXPECT_NE(summary.find("| CA | 1 |"), std::string::npos);
    EXPECT_EQ(summary.find("| CB |"), std::string::npos);
    f.model->EditAnalysis().SetJointResult(result);
    f.model->EditAnalysis().ClearJointResult();
    EXPECT_TRUE(target.HasFinalModel(rhbm_gem::FittingStage::Second));
    component.state.reset(); component.stop_reason = "invalid-initial-widths";
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, result, "next-run");
    EXPECT_FALSE(target.HasFinalModel(rhbm_gem::FittingStage::Second));
    EXPECT_EQ(target.GetStageEstimate(rhbm_gem::FittingStage::Second).reason, "invalid-initial-widths");
    summary = core::BuildSecondStageSpotSummary(*f.model);
    EXPECT_NE(summary.find("| CA | 0 | 0 | 1"), std::string::npos);
}

TEST(JointComponentPartialSelectionTest, PostFitPeelingUsesGridOperatorAndRetainsMissingCoverage)
{
    auto f = joint_partial_test::Make("partial");
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    auto snapshot = core::CaptureJointAnalysisResult(core::FitJointComponents(problem, f.b));
    ASSERT_TRUE(snapshot.components.front().state);
    // Make the neighbor contribution strictly negative, including zero Gaussian A.
    snapshot.components.front().state->ac[2] = 0;
    snapshot.components.front().state->ac[3] = -0.4;
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, snapshot, "peeling-test");
    auto * atom = f.model->FindAtomPtr(1);
    auto editor = f.model->EditAnalysis();
    const SamplingPointList points{
        {0.2, {0.13, 0.07, 0.02}, false},
        {2.4, {2.4, 0, 0}, true},
        {2.9, {2.9, 0, 0}, true}};
    const auto raw = second_stage_test::SampleExperimentPoints(*f.map, points);
    editor.SetAtomLocalRawSamplingEntries(*atom, raw);
    const auto view = rhbm_gem::AtomLocalPotentialView::For(*atom);
    const auto before = view.GetFinalModel(rhbm_gem::FittingStage::Second).ToVector();
    const auto outputs = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem);
    const auto & output = outputs.at(1);
    ASSERT_EQ(output.samples.size(), raw.size());
    ASSERT_TRUE(output.samples[0].response);
    EXPECT_GT(*output.samples[0].response, raw[0].response);
    const auto & neighbor = rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(2))
        .GetFinalModel(rhbm_gem::FittingStage::Second);
    const auto stencil = core::detail::MakeTricubicStencil(*f.map, points[0].position);
    const double neighbor_prediction = core::detail::InterpolateTricubic(stencil, [&](const auto & node) {
        const auto position = core::simulation::GridPosition(node, f.map->GetGridSpacing(), f.map->GetOrigin());
        const double squared = core::simulation::SupportSquare(position, f.model->FindAtomPtr(2)->GetPosition());
        return squared <= 6.25 ? neighbor.ResponseAtDistance(std::sqrt(squared)) : 0.0;
    });
    EXPECT_NEAR(*output.samples[0].response, raw[0].response - neighbor_prediction, 1e-12);
    EXPECT_FALSE(output.samples[2].response);
    EXPECT_EQ(output.samples[2].reason, "outside-joint-domain");
    editor.SetAtomPostFitPeeling(*atom, output);
    ASSERT_TRUE(view.GetLocalFittingPeelingRatio(0, 1));
    EXPECT_LT(*view.GetLocalFittingPeelingRatio(0, 1), 0);
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(0, 3));
    EXPECT_EQ(view.GetFinalModel(rhbm_gem::FittingStage::Second).ToVector(), before);
    const auto after = view.GetRawSamplingEntries(false);
    for (std::size_t i = 0; i < raw.size(); ++i)
    {
        EXPECT_EQ(after[i].response, raw[i].response);
        EXPECT_EQ(after[i].point.position, raw[i].point.position);
        EXPECT_EQ(after[i].point.is_selected, raw[i].point.is_selected);
    }
    auto missing = rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(2)).GetStageEstimate(rhbm_gem::FittingStage::Second);
    missing.point.reset(); missing.reason = "missing-state";
    editor.SetAtomStageEstimate(rhbm_gem::FittingStage::Second, *f.model->FindAtomPtr(2), missing);
    const auto incomplete = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem);
    EXPECT_FALSE(incomplete.at(1).samples[0].response);
    EXPECT_EQ(incomplete.at(1).samples[0].reason, "missing-contributor-state");
}

TEST(JointComponentPartialSelectionTest, PeelingIncludesDistantHaloThroughInterpolationStencil)
{
    auto f = joint_partial_test::Make("weak");
    // Binary-exact spacing makes the y/z zero-weight stencil nodes exact.
    f.map = std::make_unique<rhbm_gem::MapObject>(std::array<int, 3>{49, 25, 25},
        std::array<double, 3>{.25, .25, .25}, std::array<double, 3>{-6, -3, -3},
        std::make_unique<double[]>(49 * 25 * 25));
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    auto editor = f.model->EditAnalysis();
    for (std::size_t i = 0; i < f.a.size(); ++i)
    {
        rhbm_gem::LocalStageEstimate estimate;
        estimate.point = rhbm_gem::GaussianModel3D{f.a[i], f.b[i], f.c[i]};
        estimate.source.method = rhbm_gem::EstimateMethod::JointComponents;
        editor.SetAtomStageEstimate(rhbm_gem::FittingStage::Second,
            *f.model->FindAtomPtr(static_cast<int>(i + 1)), estimate);
    }
    auto * target = f.model->FindAtomPtr(1);
    const SamplingPointList points{{2.11, {2.11, 0, 0}, true}};
    auto raw = second_stage_test::SampleExperimentPoints(*f.map, points);
    editor.SetAtomLocalRawSamplingEntries(*target, raw);
    auto result = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem).at(1);
    ASSERT_TRUE(result.samples[0].response) << result.samples[0].reason;
    // Halo center 4.8 A away, with zero analytic support at the sample position.
    // Its contribution at the 2.5 A stencil node must still be subtracted.
    EXPECT_GT(std::abs(*result.samples[0].response - raw[0].response), 1e-6);
    EXPECT_EQ(result.neighbor_count, 1);
    raw[0].response = 0;
    editor.SetAtomLocalRawSamplingEntries(*target, raw);
    result = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem).at(1);
    editor.SetAtomPostFitPeeling(*target, result);
    EXPECT_FALSE(rhbm_gem::AtomLocalPotentialView::For(*target).GetLocalFittingPeelingRatio(0, 3));
    raw[0].response = -100;
    editor.SetAtomLocalRawSamplingEntries(*target, raw);
    result = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem).at(1);
    ASSERT_TRUE(result.samples[0].response);
    EXPECT_LT(*result.samples[0].response, 0);
}

TEST(JointComponentPartialSelectionTest, FullComponentUncertaintyMatchesDenseReferenceAndReportsGates)
{
    auto f = joint_partial_test::Make("all");
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    auto snapshot = core::CaptureJointAnalysisResult(core::FitJointComponents(problem, f.b));
    auto & component = snapshot.components.front();
    ASSERT_TRUE(component.state);
    component.runtime_convergence = core::JointCheckStatus::Passed;
    // A fixed, slightly perturbed endpoint gives a known nonzero residual variance.
    component.state->ac[0] += 0.01;
    const auto & input = problem.Input();
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(input.observations.size()), 6);
    Eigen::VectorXd residual = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(input.observations.size()));
    for (std::size_t a = 0; a < 2; ++a)
        for (const auto & support : input.support[a])
        {
            const auto k = n::EvaluateKernel(support.squared_distance, component.state->b[a], 2.5);
            const auto A = component.state->ac[2*a], C = component.state->ac[2*a+1];
            jacobian.row(static_cast<Eigen::Index>(support.row)).segment<3>(static_cast<Eigen::Index>(3*a)) << k.gaussian, k.charge,
                A*k.gaussian_log_width + C*k.charge_log_width;
            residual(static_cast<Eigen::Index>(support.row)) += A*k.gaussian + C*k.charge;
        }
    for (std::size_t r = 0; r < input.observations.size(); ++r) residual(static_cast<Eigen::Index>(r)) -= input.observations[r];
    const Eigen::MatrixXd dense = (jacobian.transpose()*jacobian).inverse() *
        (residual.squaredNorm() / static_cast<double>(input.observations.size()-6));
    auto uncertainty = core::detail::ComputeJointUncertainty(problem, snapshot);
    for (int a = 0; a < 2; ++a)
    {
        ASSERT_EQ(uncertainty.at(a+1).status, rhbm_gem::EvidenceStatus::Available);
        ASSERT_TRUE(uncertainty.at(a+1).covariance);
        EXPECT_TRUE(uncertainty.at(a+1).covariance->isApprox(dense.block<3,3>(3*a,3*a), 1e-8));
    }
    auto exact = input;
    // Use identical arithmetic to the service when constructing a zero-residual fixture.
    std::fill(exact.observations.begin(), exact.observations.end(), 0.0);
    for (std::size_t a = 0; a < 2; ++a)
        for (const auto & support : exact.support[a])
        {
            const auto k = n::EvaluateKernel(support.squared_distance, component.state->b[a], 2.5);
            exact.observations[support.row] += component.state->ac[2*a]*k.gaussian + component.state->ac[2*a+1]*k.charge;
        }
    EXPECT_EQ(core::detail::ComputeJointUncertainty(core::JointProblem(exact), snapshot).at(1).reason,
        "residual-variance-unavailable");
    auto deficient = input; deficient.support[1] = deficient.support[0];
    auto deficient_snapshot = snapshot;
    deficient_snapshot.components[0].state->b[1] = deficient_snapshot.components[0].state->b[0];
    EXPECT_EQ(core::detail::ComputeJointUncertainty(core::JointProblem(deficient), deficient_snapshot).at(1).reason,
        "rank-deficient-jacobian");
    auto short_snapshot = snapshot; short_snapshot.components[0].rows.resize(6);
    EXPECT_EQ(core::detail::ComputeJointUncertainty(problem, short_snapshot).at(1).reason,
        "insufficient-residual-degrees-of-freedom");
    component.runtime_convergence = core::JointCheckStatus::Failed;
    EXPECT_EQ(core::detail::ComputeJointUncertainty(problem, snapshot).at(1).reason, "component-not-converged");
    component.runtime_convergence = core::JointCheckStatus::Passed;
    component.state->ac[0] = 0;
    EXPECT_EQ(core::detail::ComputeJointUncertainty(problem, snapshot).at(2).reason, "active-amplitude-boundary");
    component.state.reset();
    EXPECT_EQ(core::detail::ComputeJointUncertainty(problem, snapshot).at(1).reason, "missing-component-state");
    rhbm_gem::LocalStageEstimate stage;
    stage.source.role = rhbm_gem::FittingRole::Target;
    stage.point = rhbm_gem::GaussianModel3D(0, 0.5);
    const auto evidence = core::detail::BuildJointParameterEvidence(stage);
    EXPECT_EQ(evidence.reason, "nonpositive-amplitude");
    EXPECT_FALSE(evidence.estimate);
}

TEST(JointComponentPartialSelectionTest, GroupPosteriorUsesEvidenceWithoutSamplesAndPreservesSecond)
{
    using namespace rhbm_gem;
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int i = 1; i <= 5; ++i)
    {
        auto atom = std::make_unique<AtomObject>();
        atom->SetSerialID(i); atom->SetElement(Element::CARBON);
        atom->SetChainID("A"); atom->SetComponentID("ALA"); atom->SetAtomID("CA");
        atoms.push_back(std::move(atom));
    }
    ModelObject model(std::move(atoms)); model.SelectAllAtoms();
    auto editor = model.EditAnalysis(); editor.InitializeFromSelection();
    std::map<int, LocalStageEstimate> stages;
    for (int i = 1; i <= 5; ++i)
    {
        LocalStageEstimate stage;
        stage.source = {EstimateMethod::JointComponents, std::to_string(i), "component", "evidence-test", FittingRole::Target};
        stage.point = GaussianModel3D(2 + i*0.2, 0.5 + 0.02*(i*i % 7), -0.1*i);
        stage.uncertainty.status = EvidenceStatus::Available;
        stage.uncertainty.method = "iid-ls-linearized";
        stage.uncertainty.covariance = Eigen::Matrix3d::Identity()*0.001;
        if (i == 3) stage.uncertainty.status = EvidenceStatus::Unavailable;
        stages[i] = std::move(stage);
    }
    editor.ApplySecondStageEstimates(stages);
    for (const auto & [id, stage] : stages)
        editor.SetAtomGroupEvidence(*model.FindAtomPtr(id), core::detail::BuildJointParameterEvidence(stage));
    core::FitOptions options; options.estimator = core::PotentialEstimator::JOINT_COMPONENTS; options.quiet_mode = true;
    core::RunGroupPotentialFitting(model, options);
    const auto keys = model.GetAnalysisView().CollectAtomGroupKeys();
    ASSERT_EQ(keys.size(), 1);
    const auto summary = *model.GetAnalysisView().GetGroupParameterSummary(keys[0]);
    ASSERT_EQ(summary.status, EvidenceStatus::Available) << summary.reason;
    EXPECT_EQ(summary.point_count, 5); EXPECT_EQ(summary.eligible_count, 4); EXPECT_EQ(summary.excluded_count, 1);
    EXPECT_EQ(summary.member_ids, (std::vector<int>{1,2,4,5}));
    EXPECT_FALSE(AtomLocalPotentialView::For(*model.FindAtomPtr(3)).GetGroupMemberResult());
    const auto first = AtomLocalPotentialView::For(*model.FindAtomPtr(1));
    const auto before = first.GetGroupMemberResult()->posterior.GetModel().ToVector();
    const auto second = first.GetFinalModel(FittingStage::Second).ToVector();
    auto detached_summary = model.GetAnalysisView().GetGroupParameterSummary(keys[0]);
    detached_summary->inference->member_results.front().posterior = {};
    EXPECT_EQ(first.GetGroupMemberResult()->posterior.GetModel().ToVector(), before);
    EXPECT_EQ(model.GetAnalysisView().GetGroupParameterSummary(keys[0])->inference->member_results.front().posterior.GetModel().ToVector(), before);
    EXPECT_FALSE(first.GetGroupMemberResult()->charge_inferred);
    EXPECT_EQ(first.GetGroupMemberResult()->evidence_source_id, "evidence-test");
    {
        const command_test::ScopedTempDir directory{"joint_posterior_roundtrip"};
        DataRepository repository(directory.path()/"posterior.sqlite"); repository.SaveModel(model,"model");
        const auto loaded=repository.LoadModel("model");
        const auto saved=AtomLocalPotentialView::For(*loaded->FindAtomPtr(1));
        ASSERT_TRUE(saved.GetGroupMemberResult());
        EXPECT_EQ(saved.GetGroupMemberResult()->posterior.GetModel().ToVector(),before);
        EXPECT_EQ(*saved.GetGroupMemberResult()->parameter_covariance,*first.GetGroupMemberResult()->parameter_covariance);
        EXPECT_FALSE(saved.GetGroupMemberResult()->charge_inferred);
        EXPECT_TRUE(std::isnan(saved.GetGroupMemberResult()->posterior.GetStandardDeviationModel().GetOffset()));
        ASSERT_TRUE(loaded->GetAnalysisView().GetGroupParameterSummary(keys[0])->inference);
        EXPECT_EQ(loaded->GetAnalysisView().GetGroupParameterSummary(keys[0])->member_ids,summary.member_ids);
    }
    {
        ModelObject changed_model(model);
        auto changed_editor = changed_model.EditAnalysis();
        const auto changed_view = AtomLocalPotentialView::For(*changed_model.FindAtomPtr(1));
        auto stage = changed_view.GetStageEstimate(FittingStage::Second);
        stage.uncertainty.covariance = *stage.uncertainty.covariance * 2;
        changed_editor.SetAtomStageEstimate(FittingStage::Second, *changed_model.FindAtomPtr(1), stage);
        EXPECT_FALSE(changed_view.GetGroupEvidence());
        EXPECT_FALSE(changed_view.GetGroupMemberResult());
        EXPECT_FALSE(changed_model.GetAnalysisView().GetGroupParameterSummary(keys[0]));
        EXPECT_FALSE(AtomLocalPotentialView::For(*changed_model.FindAtomPtr(2)).GetGroupMemberResult());
    }
    editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(1), {{-12345, {0.3, {0,0,0}, true}}});
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(first.GetGroupMemberResult()->posterior.GetModel().ToVector(), before);
    EXPECT_EQ(first.GetFinalModel(FittingStage::Second).ToVector(), second);
    auto changed = *first.GetGroupEvidence(); changed.covariance = *changed.covariance * 10;
    editor.SetAtomGroupEvidence(*model.FindAtomPtr(1), changed);
    core::RunGroupPotentialFitting(model, options);
    EXPECT_FALSE(first.GetGroupMemberResult()->posterior.GetModel().ToVector().isApprox(before, 1e-9));
    for (int i = 2; i <= 5; ++i)
    {
        auto evidence = *AtomLocalPotentialView::For(*model.FindAtomPtr(i)).GetGroupEvidence();
        evidence.status = EvidenceStatus::Unavailable;
        editor.SetAtomGroupEvidence(*model.FindAtomPtr(i), evidence);
    }
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(model.GetAnalysisView().GetGroupParameterSummary(keys[0])->reason, "single-member");
    EXPECT_FALSE(first.GetGroupMemberResult());
    EXPECT_FALSE(model.GetAnalysisView().HasAtomGroupPrior(keys[0]));
    for (int i = 1; i <= 5; ++i)
    {
        auto evidence = *AtomLocalPotentialView::For(*model.FindAtomPtr(i)).GetGroupEvidence();
        evidence.status = EvidenceStatus::Available;
        evidence.estimate = Eigen::Vector2d(1, 4);
        evidence.covariance = Eigen::Matrix2d::Identity() * 0.01;
        editor.SetAtomGroupEvidence(*model.FindAtomPtr(i), evidence);
    }
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(model.GetAnalysisView().GetGroupParameterSummary(keys[0])->reason, "singular-group-covariance");
    EXPECT_FALSE(first.GetGroupMemberResult());
}

TEST(JointComponentPartialSelectionTest, SummaryReportsMixedSourcesWithoutSnapshot)
{
    auto f = joint_partial_test::Make("all");
    const auto result = core::CaptureJointAnalysisResult(core::FitJointComponents(core::BuildJointProblem(*f.map, *f.model), f.b));
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, result, "mixed-test");
    auto editor = f.model->EditAnalysis();
    auto local = rhbm_gem::LocalGaussianResult{};
    local.mdpde = {rhbm_gem::GaussianModel3D{2, .5, 0}, {}};
    editor.SetAtomLocalGaussianResult(rhbm_gem::FittingStage::Second, *f.model->FindAtomPtr(2), local);
    EXPECT_FALSE(f.model->GetAnalysisView().GetJointResult());
    EXPECT_NE(core::BuildSecondStageSpotSummary(*f.model).find("Estimator: mixed"), std::string::npos);
}
