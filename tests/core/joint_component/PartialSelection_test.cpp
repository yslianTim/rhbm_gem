#include <gtest/gtest.h>
#include <boost/json.hpp>
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <algorithm>
#include "support/JointPartialSelection.hpp"
#include "support/JointHaloObservationExclusion.hpp"
#include "support/JointProfileSearchDiagnostic.hpp"
#include "support/DataObjectTestSupport.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/HaloEligibility.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include <map>
#include <numeric>
#include <utility>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include "support/CommandTestHelpers.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/detail/PostJointLocalValidation.hpp"
#include "core/detail/StageSummary.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "core/detail/MapInterpolation.hpp"
#include "core/detail/MapSampler.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include "data/io/detail/StageResultJson.hpp"
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <set>

namespace {
namespace core=rhbm_gem::core;
namespace n=core::joint_component;
namespace sim=core::simulation;
void ResetProgressPercentForTest()
{
    testing::internal::CaptureStdout();
    Logger::ProgressPercent(0, 1);
    testing::internal::GetCapturedStdout();
}
core::JointProblem FullContributorProblem(const core::JointProblem & problem)
{
    auto input=problem.Input();
    if(input.selection_domain)
    {
        input.selection_domain->contributor_policy="all-non-hydrogen";
        input.selection_domain->fixed_halo_selection.reset();
    }
    return core::JointProblem(std::move(input));
}
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
    std::size_t candidate_halos{},retained_halos{},zero_support_halos{};
    for(std::size_t a=0;a<support.size();++a)
    {
        const auto & atom=model.GetAtomList()[a];
        if(atom->GetElement()==Element::HYDROGEN) continue;
        if(!selected.contains(atom.get()))
        {
            ++candidate_halos;
            if(support[a].empty()) {++zero_support_halos; continue;}
            ++retained_halos;
        }
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
    ASSERT_TRUE(actual.selection_domain->fixed_halo_selection);
    const auto & provenance=*actual.selection_domain->fixed_halo_selection;
    EXPECT_EQ(actual.selection_domain->contributor_policy,rhbm_gem::FixedHaloContributorPolicy);
    EXPECT_EQ(provenance.minimum_observations,1u);
    EXPECT_EQ(provenance.candidate_halos,candidate_halos);
    EXPECT_EQ(provenance.retained_halos,retained_halos);
    EXPECT_EQ(provenance.zero_support_halos,zero_support_halos);
    double square=0; for(double value:observations) square+=value*value;
    EXPECT_NEAR(problem.ObservationScale(),std::max(1.,std::sqrt(square)),1e-12*problem.ObservationScale());
}
void SetFixedHaloExactMap(joint_partial_test::Fixture & f)
{
    const auto dims=f.map->GetGridSize();
    const auto spacing=f.map->GetGridSpacing(),origin=f.map->GetOrigin();
    auto values=std::make_unique<double[]>(f.map->GetMapValueArraySize());
    for(int z=0;z<dims[2];++z) for(int y=0;y<dims[1];++y) for(int x=0;x<dims[0];++x)
    {
        const auto position=sim::GridPosition({x,y,z},spacing,origin);
        double value{};
        for(std::size_t atom=0;atom<f.a.size();++atom)
        {
            const auto squared=sim::SupportSquare(position,f.model->GetAtomList()[atom]->GetPosition());
            if(squared>6.25) continue;
            const auto basis=n::EvaluateKernel(squared,f.b[atom],2.5);
            value+=f.a[atom]*basis.gaussian+(atom==0 ? f.c[atom]*basis.charge : 0.0);
        }
        values[static_cast<std::size_t>(x+dims[0]*(y+dims[1]*z))]=value;
    }
    f.map->SetMapValueArray(std::move(values));
}
rhbm_gem::JointFixedHaloParameter FixedHaloFor(const joint_partial_test::Fixture & f,std::size_t index)
{
    return {index,std::to_string(f.model->GetAtomList()[index]->GetSerialID()),f.a[index],f.b[index],0.0,17,0,
        std::string(rhbm_gem::FixedHaloParameterSourceContractVersion),true};
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
    ASSERT_TRUE(p.Input().selection_domain->fixed_halo_selection);
    EXPECT_EQ(p.Input().selection_domain->fixed_halo_selection->zero_support_halos,1u);
    const auto fit=core::FitJointComponents(p,{.5,.5});
    EXPECT_FALSE(fit.assembled_state);
    const auto & data=core::JointProblemAccess::Get(p);
    EXPECT_EQ(data.partition.unobserved_atoms.size(),1);
}
TEST(JointComponentPartialSelectionTest, RetainsEverySupportedHaloAndKeepsFixedRows)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for(const auto & [id,x]:std::vector<std::pair<int,double>>{{10,0.0},{20,4.5},{30,3.5},{40,30.0}})
    {
        auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(id);
        atom->SetElement(Element::CARBON); atom->SetPosition(x,0,0); atoms.push_back(std::move(atom));
    }
    rhbm_gem::ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==10;});
    auto values=std::make_unique<double[]>(11);
    for(std::size_t row=0;row<11;++row) values[row]=static_cast<double>(row)-3.0;
    rhbm_gem::MapObject map({11,1,1},{.5,.5,.5},{-2.5,0,0},std::move(values));

    const auto problem=core::BuildJointProblem(map,model);
    const auto & input=problem.Input();
    EXPECT_EQ(input.row_ids,(std::vector<std::string>{"0","1","2","3","4","5","6","7","8","9","10"}));
    EXPECT_EQ(input.observations,(std::vector<double>{-3,-2,-1,0,1,2,3,4,5,6,7}));
    EXPECT_EQ(input.atom_ids,(std::vector<std::string>{"10","20","30"}));
    ASSERT_TRUE(input.selection_domain);
    EXPECT_EQ(input.selection_domain->target_indices,(std::vector<std::size_t>{0}));
    ASSERT_TRUE(input.selection_domain->fixed_halo_selection);
    const auto & provenance=*input.selection_domain->fixed_halo_selection;
    EXPECT_EQ(input.selection_domain->contributor_policy,rhbm_gem::FixedHaloContributorPolicy);
    EXPECT_EQ(provenance.candidate_halos,3u);
    EXPECT_EQ(provenance.retained_halos,2u);
    EXPECT_EQ(provenance.zero_support_halos,1u);
    ASSERT_EQ(input.support.size(),3u);
    EXPECT_EQ(input.support[1].size(),2u);
    EXPECT_EQ(input.support[2].size(),4u);

    const auto workset=core::detail::MakeJointFittingWorkset(model,problem);
    EXPECT_EQ(workset.contributors.size(),input.atom_ids.size());
    EXPECT_EQ(core::detail::CollectFirstStageAtoms(workset).size(),input.atom_ids.size());
    const auto & data=core::JointProblemAccess::Get(problem);
    EXPECT_EQ(data.construction_memory.shape.component_atoms,1u);
    EXPECT_EQ(data.construction_memory.shape.component_support_entries,17u);

    auto invalid=input; invalid.support[2].clear();
    EXPECT_THROW(core::JointProblem(std::move(invalid)),std::invalid_argument);
}

TEST(JointComponentPartialSelectionTest, SelectedTargetWithOneObservationIsAlwaysRetained)
{
    auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(1);
    atom->SetElement(Element::CARBON); atom->SetPosition(0,0,0);
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms; atoms.push_back(std::move(atom));
    rhbm_gem::ModelObject model(std::move(atoms)); model.SelectAllAtoms();
    rhbm_gem::MapObject map({1,1,1},{1,1,1},{0,0,0});
    const auto problem=core::BuildJointProblem(map,model);
    EXPECT_EQ(problem.Input().atom_ids,(std::vector<std::string>{"1"}));
    ASSERT_TRUE(problem.Input().selection_domain);
    EXPECT_EQ(problem.Input().selection_domain->target_indices,(std::vector<std::size_t>{0}));
    EXPECT_EQ(problem.Input().support[0].size(),1u);
}

TEST(JointComponentPartialSelectionTest, SelectedTargetWithZeroObservationsKeepsIdentity)
{
    auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(1);
    atom->SetElement(Element::CARBON); atom->SetPosition(100,0,0);
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms; atoms.push_back(std::move(atom));
    rhbm_gem::ModelObject model(std::move(atoms)); model.SelectAllAtoms();
    rhbm_gem::MapObject map({1,1,1},{1,1,1},{0,0,0});
    const auto problem=core::BuildJointProblem(map,model);
    EXPECT_EQ(problem.Input().atom_ids,(std::vector<std::string>{"1"}));
    EXPECT_TRUE(problem.Input().row_ids.empty());
    EXPECT_TRUE(problem.Input().observations.empty());
    EXPECT_TRUE(problem.Input().support[0].empty());
    ASSERT_TRUE(problem.Input().selection_domain);
    EXPECT_EQ(problem.Input().selection_domain->target_indices,(std::vector<std::size_t>{0}));
}

TEST(JointComponentPartialSelectionTest, FixedHaloDoesNotConnectTargetComponents)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for(const auto & [id,x]:std::vector<std::pair<int,double>>{{1,0.0},{2,7.5},{3,3.75}})
    {
        auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(id);
        atom->SetElement(Element::CARBON); atom->SetPosition(x,0,0); atoms.push_back(std::move(atom));
    }
    rhbm_gem::ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()!=3;});
    rhbm_gem::MapObject map({6,1,1},{2.5,.1,.1},{-2.5,0,0});
    const auto problem=core::BuildJointProblem(map,model);
    EXPECT_EQ(problem.Input().atom_ids,(std::vector<std::string>{"1","2","3"}));
    const auto & data=core::JointProblemAccess::Get(problem);
    ASSERT_EQ(data.domain.atoms.size(),2u);
    ASSERT_EQ(data.partition.components.size(),2u);
    EXPECT_EQ(data.partition.components[0].atoms.size(),1u);
    EXPECT_EQ(data.partition.components[1].atoms.size(),1u);
}

TEST(JointComponentPartialSelectionTest, SharedFixedHaloPersistsOutsideTargetComponents)
{
    auto f=joint_partial_test::Make("bridge");
    const auto dims=f.map->GetGridSize();
    const auto spacing=f.map->GetGridSpacing(),origin=f.map->GetOrigin();
    auto values=std::make_unique<double[]>(f.map->GetMapValueArraySize());
    for(int z=0;z<dims[2];++z) for(int y=0;y<dims[1];++y) for(int x=0;x<dims[0];++x)
    {
        const auto position=sim::GridPosition({x,y,z},spacing,origin);
        double value{};
        for(std::size_t atom=0;atom<f.a.size();++atom)
        {
            const auto squared=sim::SupportSquare(position,f.model->GetAtomList()[atom]->GetPosition());
            if(squared>6.25) continue;
            const auto basis=n::EvaluateKernel(squared,f.b[atom],2.5);
            value+=f.a[atom]*basis.gaussian+(atom<2 ? f.c[atom]*basis.charge : 0.0);
        }
        values[static_cast<std::size_t>(x+dims[0]*(y+dims[1]*z))]=value;
    }
    f.map->SetMapValueArray(std::move(values));
    const auto problem=core::BuildJointProblem(*f.map,*f.model);
    ASSERT_EQ(problem.Input().atom_ids,(std::vector<std::string>{"1","2","3"}));
    const auto & input=problem.Input();
    auto overlaps=[&](std::size_t a,std::size_t b) {
        return std::any_of(input.support[a].begin(),input.support[a].end(),[&](const auto & left) {
            return std::any_of(input.support[b].begin(),input.support[b].end(),[&](const auto & right) {
                return left.row==right.row;
            });
        });
    };
    EXPECT_FALSE(overlaps(0,1));
    EXPECT_TRUE(overlaps(0,2));
    EXPECT_TRUE(overlaps(1,2));
    rhbm_gem::JointFixedHaloParameter recovered_fixed{
        2,"3",f.a[2],f.b[2],0.0,17,
        static_cast<int>(rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED),
        std::string(rhbm_gem::FixedHaloParameterSourceContractVersion),true};
    recovered_fixed.first_stage_qualification=rhbm_gem::RHBMSolveQualification::RefinedSuccess;
    rhbm_gem::RHBMFailedOnlyRecoveryDiagnostics recovery;
    recovery.accepted=true;
    recovery.reason="endpoint-refinement-accepted";
    recovery.endpoint_refinement.emplace();
    recovery.endpoint_refinement->accepted=true;
    recovery.endpoint_refinement->reason="accepted";
    recovery.endpoint_refinement->candidate_residual=1e-14;
    recovery.endpoint_refinement->reference_residual=1e-11;
    recovery.endpoint_refinement->floor_masks_equal=true;
    recovered_fixed.first_stage_recovery=recovery;
    const std::vector<rhbm_gem::JointFixedHaloParameter> fixed{recovered_fixed};
    const auto initialized=problem.WithFixedHalos(fixed);
    const auto & data=core::JointProblemAccess::Get(initialized);
    ASSERT_EQ(data.partition.components.size(),2u);
    EXPECT_EQ(data.domain.atoms.size(),2u);
    ASSERT_EQ(data.partition.components[0].atoms.size(),1u);
    ASSERT_EQ(data.partition.components[1].atoms.size(),1u);
    EXPECT_EQ(data.partition.components[0].atoms.front(),0);
    EXPECT_EQ(data.partition.components[1].atoms.front(),1);
    EXPECT_EQ(data.fixed_halos,fixed);

    const auto fit=core::FitJointComponents(initialized,{f.b[0],f.b[1]});
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_EQ(fit.components.size(),2u);
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(fit.assembled_state->b.size(),2u);
    const auto saved=core::CaptureJointAnalysisResult(fit);
    EXPECT_EQ(saved.fixed_halos,fixed);
    for(const auto & component:saved.components)
    {
        ASSERT_EQ(component.atoms.size(),1u);
        EXPECT_LT(component.atoms.front(),2u);
    }
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_EQ(decoded.fixed_halos,fixed);
    EXPECT_EQ(decoded.fixed_halos.front().first_stage_status,
        static_cast<int>(rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED));
    EXPECT_EQ(decoded.fixed_halos.front().first_stage_qualification,
        rhbm_gem::RHBMSolveQualification::RefinedSuccess);
    ASSERT_TRUE(decoded.fixed_halos.front().first_stage_recovery);
    EXPECT_EQ(decoded.fixed_halos.front().first_stage_recovery->reason,"endpoint-refinement-accepted");
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model,decoded,"shared-fixed-halo");
    const auto halo_stage=rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(3))
        .GetStageEstimate(rhbm_gem::FittingStage::Second);
    ASSERT_TRUE(halo_stage.point);
    EXPECT_EQ(halo_stage.source.role,rhbm_gem::FittingRole::Halo);
    EXPECT_TRUE(halo_stage.source.fixed_from_first_stage);
    EXPECT_EQ(halo_stage.convergence,rhbm_gem::JointCheckStatus::NotRun);

    const std::array<std::pair<int,double>,2> samples{{{1,-1.2},{2,1.2}}};
    for(const auto & [target_id,x]:samples)
    {
        const SamplingPointList point{{0.0,{x,0.0,0.0},true}};
        const auto raw=core::detail::BuildLocalPotentialSampleList(*f.map,point);
        f.model->EditAnalysis().SetAtomLocalRawSamplingEntries(*f.model->FindAtomPtr(target_id),raw);
        const auto peelings=core::detail::BuildPostFitPeelingSamples(
            *f.map,*f.model,problem,std::nullopt,&decoded);
        const auto & peeled=peelings.at(target_id);
        ASSERT_EQ(peeled.neighbor_count,1u);
        ASSERT_TRUE(peeled.samples.front().response) << peeled.samples.front().reason;
        const auto squared=sim::SupportSquare(point.front().position,f.model->FindAtomPtr(3)->GetPosition());
        const auto basis=n::EvaluateKernel(squared,fixed.front().width,2.5);
        EXPECT_NEAR(*peeled.samples.front().response,raw.front().response-fixed.front().amplitude*basis.gaussian,1e-11);
        const auto prepared=core::detail::BuildPostJointLocalSamples(
            *f.map,*f.model,problem,target_id,peeled,decoded);
        EXPECT_EQ(prepared.samples.size(),1u);
        EXPECT_THROW(core::detail::BuildPostJointLocalSamples(
            *f.map,*f.model,problem,3,peeled,decoded),std::invalid_argument);
    }
    f.model->SelectAllAtoms();
    f.model->EditAnalysis().InitializeFromSelection();
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model,decoded,"shared-fixed-halo");
    core::FitOptions group_options; group_options.estimator=core::PotentialEstimator::JOINT_COMPONENTS;
    group_options.quiet_mode=true;
    core::RunGroupPotentialFitting(*f.model,group_options);
    std::size_t target_exclusions{};
    for(const auto & key:f.model->GetAnalysisView().CollectAtomGroupKeys())
    {
        const auto summary=f.model->GetAnalysisView().GetGroupParameterSummary(key);
        ASSERT_TRUE(summary);
        target_exclusions+=summary->excluded_count;
        EXPECT_EQ(std::find(summary->member_ids.begin(),summary->member_ids.end(),3),summary->member_ids.end());
    }
    EXPECT_EQ(target_exclusions,2u);
}

TEST(JointComponentPartialSelectionTest, FixedHaloResponseUsesRawScaleAndReplaysFullResidual)
{
    auto f=joint_partial_test::Make("partial");
    const auto dims=f.map->GetGridSize();
    const auto spacing=f.map->GetGridSpacing(),origin=f.map->GetOrigin();
    auto values=std::make_unique<double[]>(f.map->GetMapValueArraySize());
    for(int z=0;z<dims[2];++z) for(int y=0;y<dims[1];++y) for(int x=0;x<dims[0];++x)
    {
        const auto position=sim::GridPosition({x,y,z},spacing,origin);
        double value{};
        for(std::size_t atom=0;atom<f.a.size();++atom)
        {
            const auto squared=sim::SupportSquare(position,f.model->GetAtomList()[atom]->GetPosition());
            if(squared>6.25) continue;
            const auto basis=n::EvaluateKernel(squared,f.b[atom],2.5);
            value+=f.a[atom]*basis.gaussian+(atom==0 ? f.c[atom]*basis.charge : 0.0);
        }
        values[static_cast<std::size_t>(x+dims[0]*(y+dims[1]*z))]=value;
    }
    f.map->SetMapValueArray(std::move(values));
    const auto problem=core::BuildJointProblem(*f.map,*f.model);
    ASSERT_EQ(problem.Input().atom_ids,(std::vector<std::string>{"1","2"}));
    const std::vector<rhbm_gem::JointFixedHaloParameter> fixed{{1,"2",f.a[1],f.b[1],0.,17,0,
        std::string(rhbm_gem::FixedHaloParameterSourceContractVersion),true}};
    const auto initialized=problem.WithFixedHalos(fixed);
    const auto & data=core::JointProblemAccess::Get(initialized);
    ASSERT_EQ(data.domain.atoms.size(),1u);
    ASSERT_EQ(data.layout.full_atoms,(std::vector<std::size_t>{0}));
    ASSERT_EQ(data.partition.components.size(),1u);
    EXPECT_NEAR(problem.ObservationScale(),std::max(1.0,std::sqrt(std::inner_product(
        problem.Input().observations.begin(),problem.Input().observations.end(),
        problem.Input().observations.begin(),0.0))),1e-12*problem.ObservationScale());
    ASSERT_EQ(data.fixed_halos,fixed);
    for(std::size_t row=0;row<problem.Input().observations.size();++row)
        EXPECT_DOUBLE_EQ(data.y(static_cast<Eigen::Index>(row)),problem.Input().observations[row]-
            data.fixed_prediction(static_cast<Eigen::Index>(row)));

    const auto fit=core::FitJointComponents(initialized,{f.b[0]});
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(fit.assembled_state);
    ASSERT_TRUE(fit.prediction);
    ASSERT_TRUE(fit.objective);
    ASSERT_EQ(fit.assembled_state->b.size(),1u);
    ASSERT_EQ(fit.components.size(),1u);
    EXPECT_EQ(fit.components[0].atoms,(std::vector<std::size_t>{0}));
    EXPECT_EQ(initialized.FixedHalos(),fixed);
    const auto replay=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="full-observation-replay";
    });
    ASSERT_NE(replay,fit.evidence.end());
    EXPECT_EQ(replay->status,rhbm_gem::JointCheckStatus::Passed);
    for(const auto * name:{"fixed-halo-replay","full-objective-replay","target-ac-feasibility",
        "target-ac-optimality","target-width-stationarity","target-contributor-mapping",
        "original-observation-identities","fixed-halo-provenance"})
    {
        const auto check=std::find_if(fit.evidence.begin(),fit.evidence.end(),[&](const auto & item) {
            return item.name==name;
        });
        ASSERT_NE(check,fit.evidence.end()) << name;
        EXPECT_EQ(check->status,rhbm_gem::JointCheckStatus::Passed) << name;
    }

    double rss{};
    for(std::size_t row=0;row<problem.Input().observations.size();++row)
    {
        double target_prediction{};
        for(const auto & support:problem.Input().support[0]) if(support.row==row)
        {
            const auto basis=n::EvaluateKernel(support.squared_distance,fit.assembled_state->b[0],2.5);
            target_prediction=fit.assembled_state->ac[0]*basis.gaussian+
                fit.assembled_state->ac[1]*basis.charge;
        }
        const double effective_residual=target_prediction-data.y(static_cast<Eigen::Index>(row));
        const double full_residual=fit.prediction->at(row)-problem.Input().observations[row];
        EXPECT_NEAR(full_residual,effective_residual,2e-12*(1+std::abs(full_residual)));
        rss+=full_residual*full_residual;
    }
    EXPECT_NEAR(*fit.objective,rss/(2*problem.ObservationScale()*problem.ObservationScale()),1e-12);

    auto noisy_input=problem.Input();
    const auto noisy_row=noisy_input.support[1].front().row;
    noisy_input.observations[noisy_row]+=.05;
    const core::JointProblem noisy_problem(std::move(noisy_input));
    const auto noisy_fixed=noisy_problem.WithFixedHalos(fixed);
    const auto noisy_fit=core::FitJointComponents(noisy_fixed,{f.b[0]});
    EXPECT_EQ(noisy_fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(noisy_fit.assembled_state);
    ASSERT_TRUE(noisy_fit.prediction);
    double fixed_halo_amplitude_gradient{},fixed_halo_charge_gradient{},fixed_halo_width_gradient{};
    for(const auto & support:noisy_problem.Input().support[1])
    {
        const auto basis=n::EvaluateKernel(support.squared_distance,fixed[0].width,2.5);
        const double residual=noisy_fit.prediction->at(support.row)-noisy_problem.Input().observations[support.row];
        fixed_halo_amplitude_gradient+=basis.gaussian*residual;
        fixed_halo_charge_gradient+=basis.charge*residual;
        fixed_halo_width_gradient+=fixed[0].amplitude*basis.gaussian_log_width*residual;
    }
    fixed_halo_width_gradient/=std::pow(noisy_problem.ObservationScale(),2);
    const double fixed_halo_gradient=std::max({std::abs(fixed_halo_amplitude_gradient),
        std::abs(fixed_halo_charge_gradient),std::abs(fixed_halo_width_gradient)});
    EXPECT_GT(fixed_halo_gradient,1e-12);
    EXPECT_EQ(noisy_fixed.FixedHalos(),fixed);

    auto snapshot=core::CaptureJointAnalysisResult(fit);
    EXPECT_EQ(snapshot.parameterization_contract,rhbm_gem::JointTargetOnlyFixedHaloParameterizationContract);
    EXPECT_EQ(snapshot.fixed_halos,fixed);
    ASSERT_TRUE(snapshot.assembled_state);
    ASSERT_EQ(snapshot.assembled_state->b.size(),1u);
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(snapshot));
    EXPECT_EQ(decoded.parameterization_contract,rhbm_gem::JointTargetOnlyFixedHaloParameterizationContract);
    EXPECT_EQ(decoded.fixed_halos,fixed);
    ASSERT_TRUE(decoded.assembled_state);
    ASSERT_EQ(decoded.assembled_state->b.size(),1u);

    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model,decoded,"fixed-halo-roundtrip");
    const auto target_stage=rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(1))
        .GetStageEstimate(rhbm_gem::FittingStage::Second);
    ASSERT_TRUE(target_stage.point);
    EXPECT_EQ(target_stage.source.role,rhbm_gem::FittingRole::Target);
    EXPECT_FALSE(target_stage.source.fixed_from_first_stage);
    const auto halo_stage=rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(2))
        .GetStageEstimate(rhbm_gem::FittingStage::Second);
    ASSERT_TRUE(halo_stage.point);
    EXPECT_DOUBLE_EQ(halo_stage.point->GetAmplitude(),fixed[0].amplitude);
    EXPECT_DOUBLE_EQ(halo_stage.point->GetWidth(),fixed[0].width);
    EXPECT_DOUBLE_EQ(halo_stage.point->GetOffset(),0.0);
    EXPECT_EQ(halo_stage.reason,"fixed-from-first-stage");
    EXPECT_EQ(halo_stage.convergence,rhbm_gem::JointCheckStatus::NotRun);
    EXPECT_EQ(halo_stage.source.role,rhbm_gem::FittingRole::Halo);
    EXPECT_TRUE(halo_stage.source.fixed_from_first_stage);

    const SamplingPointList center{{0.0,{0.0,0.0,0.0},true}};
    auto raw=core::detail::BuildLocalPotentialSampleList(*f.map,center);
    f.model->EditAnalysis().SetAtomLocalRawSamplingEntries(*f.model->FindAtomPtr(1),raw);
    const auto peeled=core::detail::BuildPostFitPeelingSamples(
        *f.map,*f.model,problem,std::nullopt,&decoded).at(1);
    ASSERT_EQ(peeled.samples.size(),1u);
    ASSERT_TRUE(peeled.samples[0].response) << peeled.samples[0].reason;
    EXPECT_EQ(peeled.neighbor_count,1u);
    const auto halo_square=sim::SupportSquare(center[0].position,f.model->FindAtomPtr(2)->GetPosition());
    const auto halo_basis=n::EvaluateKernel(halo_square,fixed[0].width,2.5);
    EXPECT_NEAR(*peeled.samples[0].response,raw[0].response-fixed[0].amplitude*halo_basis.gaussian,1e-12);
}

TEST(JointComponentPartialSelectionTest, NonleadingTargetMapsBackToFullContributorIdentity)
{
    auto f=joint_partial_test::Make("partial");
    f.model->SelectAtoms([](const auto & atom){return atom.GetSerialID()==2;});
    const auto dims=f.map->GetGridSize();
    const auto spacing=f.map->GetGridSpacing(),origin=f.map->GetOrigin();
    auto values=std::make_unique<double[]>(f.map->GetMapValueArraySize());
    for(int z=0;z<dims[2];++z) for(int y=0;y<dims[1];++y) for(int x=0;x<dims[0];++x)
    {
        const auto position=sim::GridPosition({x,y,z},spacing,origin);
        double value{};
        for(std::size_t atom=0;atom<f.a.size();++atom)
        {
            const auto squared=sim::SupportSquare(position,f.model->GetAtomList()[atom]->GetPosition());
            if(squared>6.25) continue;
            const auto basis=n::EvaluateKernel(squared,f.b[atom],2.5);
            value+=f.a[atom]*basis.gaussian+(atom==1 ? f.c[atom]*basis.charge : 0.0);
        }
        values[static_cast<std::size_t>(x+dims[0]*(y+dims[1]*z))]=value;
    }
    f.map->SetMapValueArray(std::move(values));
    const auto problem=core::BuildJointProblem(*f.map,*f.model);
    EXPECT_EQ(problem.Input().atom_ids,(std::vector<std::string>{"1","2"}));
    ASSERT_EQ(problem.Input().selection_domain->target_indices,(std::vector<std::size_t>{1}));
    const auto initialized=problem.WithFixedHalos({{0,"1",f.a[0],f.b[0],0.,17,0,
        std::string(rhbm_gem::FixedHaloParameterSourceContractVersion),true}});
    const auto & data=core::JointProblemAccess::Get(initialized);
    EXPECT_EQ(data.target_to_contributor,(std::vector<std::size_t>{1}));
    EXPECT_EQ(data.context.atom_ids,(std::vector<std::string>{"2"}));
    const auto fit=core::FitJointComponents(initialized,{f.b[1]});
    EXPECT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(fit.assembled_state);
    ASSERT_EQ(fit.components.size(),1u);
    EXPECT_EQ(fit.components[0].atoms,(std::vector<std::size_t>{1}));
    ASSERT_TRUE(fit.target_evidence);
    ASSERT_EQ(fit.target_evidence->point_atoms.size(),1u);
    EXPECT_EQ(fit.target_evidence->point_atoms[0].atom,1u);
    for(const auto & block:fit.components[0].local_conditioning)
        EXPECT_EQ(block.core_atoms,(std::vector<std::size_t>{1}));
    ASSERT_TRUE(fit.problem);
    EXPECT_EQ(fit.problem->FixedHalos(),initialized.FixedHalos());
    EXPECT_NEAR(fit.assembled_state->ac[0],f.a[1],1e-9*(1+f.a[1]));
    EXPECT_NEAR(fit.assembled_state->ac[1],f.c[1],1e-9*(1+std::abs(f.c[1])));
    EXPECT_NEAR(fit.assembled_state->b[0],f.b[1],1e-9*(1+f.b[1]));
}

TEST(JointComponentPartialSelectionTest, ContributorOrderIsIndependentOfModelTraversal)
{
    auto make_model=[](std::vector<int> order) {
        std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
        for(const int id:order)
        {
            auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(id);
            atom->SetElement(Element::CARBON); atom->SetPosition(id==1 ? 0.0 : id==2 ? 1.0 : 3.5,0,0);
            atoms.push_back(std::move(atom));
        }
        auto model=std::make_unique<rhbm_gem::ModelObject>(std::move(atoms));
        model->SelectAtoms([](const auto & atom){return atom.GetSerialID()==1;});
        return model;
    };
    auto values=std::make_unique<double[]>(11);
    rhbm_gem::MapObject map({11,1,1},{.5,.5,.5},{-2.5,0,0},std::move(values));
    const auto forward=core::BuildJointProblem(map,*make_model({1,2,3}));
    const auto reverse=core::BuildJointProblem(map,*make_model({3,2,1}));
    EXPECT_EQ(forward.Input().atom_ids,reverse.Input().atom_ids);
    EXPECT_EQ(forward.Input().row_ids,reverse.Input().row_ids);
    EXPECT_EQ(forward.Input().observations,reverse.Input().observations);
    ASSERT_EQ(forward.Input().support.size(),reverse.Input().support.size());
    for(std::size_t atom=0;atom<forward.Input().support.size();++atom)
    {
        ASSERT_EQ(forward.Input().support[atom].size(),reverse.Input().support[atom].size());
        for(std::size_t row=0;row<forward.Input().support[atom].size();++row)
        {
            EXPECT_EQ(forward.Input().support[atom][row].row,reverse.Input().support[atom][row].row);
            EXPECT_DOUBLE_EQ(forward.Input().support[atom][row].squared_distance,
                reverse.Input().support[atom][row].squared_distance);
        }
    }
}
TEST(JointComponentPartialSelectionTest, BridgeAndMonolithicEquivalence)
{
    auto f=joint_partial_test::Make("bridge"); CheckReference(*f.map,*f.model);
    const auto p=FullContributorProblem(core::BuildJointProblem(*f.map,*f.model)); const auto & data=core::JointProblemAccess::Get(p);
    ASSERT_EQ(data.partition.components.size(),1);
    auto without=p.Input(); without.selection_domain.reset(); without.atom_ids.pop_back(); without.support.pop_back();
    const core::JointProblem disconnected(without);
    EXPECT_EQ(core::JointProblemAccess::Get(disconnected).partition.components.size(),2);
    const auto fit=core::FitJointComponents(p,f.b);
    const n::Vector widths=Eigen::Map<n::Vector>(f.b.data(),static_cast<Eigen::Index>(f.b.size()));
    const auto mono=n::AssessProfileSearchForDiagnostics(data.domain,data.y,data.context,n::SearchProfile(data.domain,data.y,widths,data.context));
    ASSERT_TRUE(fit.assembled_state); ASSERT_TRUE(mono.trusted_state);
    for(std::size_t i=0;i<fit.assembled_state->ac.size();++i)
        EXPECT_NEAR(fit.assembled_state->ac[i],mono.trusted_state->beta(static_cast<Eigen::Index>(i)),1e-10*(1+std::abs(fit.assembled_state->ac[i])));
    EXPECT_NEAR(*fit.objective,mono.trusted_state->certificate.objective/std::pow(p.ObservationScale(),2),1e-12);
}
TEST(JointComponentPartialSelectionTest, IdentifiableRecoveryOmissionAndOutsideDomain)
{
    auto f=joint_partial_test::Make("partial"); const auto p=FullContributorProblem(core::BuildJointProblem(*f.map,*f.model));
    const auto fit=core::FitJointComponents(p,f.b);
    ASSERT_EQ(fit.RuntimeConvergence(),core::JointCheckStatus::Passed); ASSERT_TRUE(fit.assembled_state);
    for(std::size_t i=0;i<f.b.size();++i)
    {
        EXPECT_NEAR(fit.assembled_state->ac[2*i],f.a[i],1e-10*(1+f.a[i]));
        EXPECT_NEAR(fit.assembled_state->ac[2*i+1],f.c[i],1e-10*(1+f.c[i]));
        EXPECT_NEAR(fit.assembled_state->b[i],f.b[i],1e-10*(1+f.b[i]));
    }
    auto omitted=p.Input(); omitted.atom_ids.resize(1); omitted.support.resize(1);
    omitted.selection_domain->contributor_policy="all-non-hydrogen";
    omitted.selection_domain->halo_eligibility.reset();
    const auto wrong=core::FitJointComponents(core::JointProblem(omitted),{f.b[0]});
    ASSERT_FALSE(wrong.objective); ASSERT_FALSE(wrong.assembled_state);
    ASSERT_TRUE(wrong.components[0].candidate_diagnostic);
    EXPECT_TRUE(wrong.components[0].candidate_diagnostic->candidate_available);
    ASSERT_GE(wrong.components[0].candidate_diagnostic->candidate_ac.size(),2);
    EXPECT_GT(std::abs(wrong.components[0].candidate_diagnostic->candidate_ac[0]-f.a[0]),1e-5);
    EXPECT_NE(wrong.RuntimeConvergence(),core::JointCheckStatus::Passed);
    std::set<std::string> rows(p.Input().row_ids.begin(),p.Input().row_ids.end());
    auto values=std::make_unique<double[]>(f.map->GetMapValueArraySize());
    for(std::size_t i=0;i<f.map->GetMapValueArraySize();++i) values[i]=rows.contains(std::to_string(i)) ? f.map->GetMapValue(i) : 12345;
    f.map->SetMapValueArray(std::move(values));
    const auto changed=FullContributorProblem(core::BuildJointProblem(*f.map,*f.model));
    EXPECT_EQ(changed.Input().observations,p.Input().observations);
    const auto repeated=core::FitJointComponents(changed,f.b);
    EXPECT_EQ(repeated.assembled_state->ac,fit.assembled_state->ac); EXPECT_EQ(repeated.objective,fit.objective);
}
TEST(JointComponentPartialSelectionTest, InvalidWidthsAreNotSuccessForRetainedContributors)
{
    auto f=joint_partial_test::Make("weak"); const auto p=FullContributorProblem(core::BuildJointProblem(*f.map,*f.model));
    ASSERT_EQ(p.Input().atom_ids,(std::vector<std::string>{"1","2"}));
    ASSERT_TRUE(p.Input().selection_domain);
    ASSERT_EQ(p.Input().selection_domain->target_indices,(std::vector<std::size_t>{0}));
    for(double invalid:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()})
    {
        const auto fit=core::FitJointComponents(p,{.5,invalid});
        ASSERT_EQ(fit.components.size(),1); EXPECT_FALSE(fit.initialization.valid);
        EXPECT_EQ(fit.components[0].stop_reason,"invalid-initial-widths"); EXPECT_FALSE(fit.prediction);
    }
    EXPECT_TRUE(core::FitJointComponents(p,{}).components.empty());
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
    EXPECT_EQ(boost::json::parse(encoded).at("schema_version").as_int64(),16);
    const auto decoded=rhbm_gem::joint_result_io::Decode(encoded);
    ASSERT_TRUE(decoded.selection_domain); EXPECT_EQ(decoded.selection_domain->target_indices,(std::vector<std::size_t>{0}));
    ASSERT_TRUE(decoded.selection_domain->fixed_halo_selection);
    EXPECT_EQ(decoded.selection_domain->fixed_halo_selection,saved.selection_domain->fixed_halo_selection);
    EXPECT_EQ(decoded.initialization.data_scope,"contributor-local-sampling-may-read-outside-target-domain");
    EXPECT_EQ(rhbm_gem::joint_result_io::Encode(decoded),encoded);
    auto missing_provenance=boost::json::parse(encoded).as_object();
    missing_provenance.at("selection_domain").as_object()["fixed_halo_selection"]=nullptr;
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(missing_provenance)),std::invalid_argument);
    missing_provenance=boost::json::parse(encoded).as_object();
    missing_provenance["ranks"]=missing_provenance.at("legacy_global_rank_evidence");
    missing_provenance.erase("legacy_global_rank_evidence");
    for(auto & item:missing_provenance.at("components").as_array())
    {
        auto & component=item.as_object();
        component["ranks"]=component.at("legacy_global_rank_evidence");
        component.erase("legacy_global_rank_evidence");
    }
    missing_provenance["schema_version"]=10;
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(missing_provenance)),std::invalid_argument);
    f.model->EditAnalysis().SetJointResult(saved);
}

TEST(JointComponentPartialSelectionTest, FixedHaloSelectionPersistsThroughJsonAndSqlite)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for(const auto & [id,x]:std::vector<std::pair<int,double>>{{10,0.0},{20,4.5},{30,3.5},{40,30.0}})
    {
        auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(id);
        atom->SetElement(Element::CARBON); atom->SetPosition(x,0,0); atoms.push_back(std::move(atom));
    }
    rhbm_gem::ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==10;});
    auto values=std::make_unique<double[]>(11);
    for(std::size_t row=0;row<11;++row) values[row]=static_cast<double>(row)-3.0;
    rhbm_gem::MapObject map({11,1,1},{.5,.5,.5},{-2.5,0,0},std::move(values));
    const auto problem=core::BuildJointProblem(map,model);
    const auto initialized=problem.WithFixedHalos({
        {1,"20",2.,.5,0.,17,0,std::string(rhbm_gem::FixedHaloParameterSourceContractVersion),true},
        {2,"30",2.,.5,0.,17,0,std::string(rhbm_gem::FixedHaloParameterSourceContractVersion),true}});
    const auto fit=core::FitJointComponents(initialized,{.5,.5,.5});
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_EQ(decoded.fixed_halos,initialized.FixedHalos());
    auto historical=boost::json::parse(rhbm_gem::joint_result_io::Encode(saved)).as_object();
    historical["schema_version"]=14;
    for(auto & item:historical.at("fixed_halos").as_array())
    {
        auto & fixed=item.as_object();
        fixed.erase("first_stage_qualification");
        fixed.erase("first_stage_recovery");
    }
    for(auto & item:historical.at("initialization").as_object().at("atoms").as_array())
    {
        auto & atom=item.as_object();
        atom.erase("first_stage_fit_qualified");
        atom.erase("first_stage_qualification");
        atom.erase("first_stage_recovery");
    }
    auto malformed_current=historical;
    malformed_current["schema_version"]=15;
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(malformed_current)),std::invalid_argument);
    const auto legacy_decoded=rhbm_gem::joint_result_io::Decode(boost::json::serialize(historical));
    ASSERT_EQ(legacy_decoded.fixed_halos.size(),initialized.FixedHalos().size());
    EXPECT_TRUE(rhbm_gem::IsQualifiedFirstStageFixedHalo(legacy_decoded.fixed_halos.front()));
    ASSERT_TRUE(decoded.selection_domain->fixed_halo_selection);
    EXPECT_EQ(decoded.selection_domain->fixed_halo_selection->candidate_halos,3u);
    EXPECT_EQ(decoded.selection_domain->fixed_halo_selection->retained_halos,2u);
    model.EditAnalysis().ApplyJointResult(decoded,"provenance");

    const command_test::ScopedTempDir directory{"joint_halo_provenance_roundtrip"};
    rhbm_gem::DataRepository repository(directory.path()/"provenance.sqlite");
    repository.SaveModel(model,"joint");
    const auto loaded=repository.LoadModel("joint");
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(loaded->GetAnalysisView().GetJointResult()->fixed_halos,decoded.fixed_halos);
    EXPECT_EQ(loaded->GetAnalysisView().GetJointResult()->selection_domain->fixed_halo_selection,
        decoded.selection_domain->fixed_halo_selection);
    const auto loaded_target=rhbm_gem::AtomLocalPotentialView::For(*loaded->FindAtomPtr(10))
        .GetStageEstimate(rhbm_gem::FittingStage::Second);
    EXPECT_EQ(loaded_target.source.role,rhbm_gem::FittingRole::Target);
    EXPECT_FALSE(loaded_target.source.fixed_from_first_stage);
    const auto loaded_halo=rhbm_gem::AtomLocalPotentialView::For(*loaded->FindAtomPtr(20))
        .GetStageEstimate(rhbm_gem::FittingStage::Second);
    EXPECT_EQ(loaded_halo.source.role,rhbm_gem::FittingRole::Halo);
    EXPECT_TRUE(loaded_halo.source.fixed_from_first_stage);
    ASSERT_TRUE(loaded_halo.point);
    EXPECT_EQ(loaded_halo.reason,"fixed-from-first-stage");

    auto radial=decoded;
    radial.selection_domain->fixed_halo_selection->zero_support_halos=1;
    const auto radial_decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(radial));
    EXPECT_EQ(radial_decoded.selection_domain->fixed_halo_selection,radial.selection_domain->fixed_halo_selection);
    model.EditAnalysis().ApplyJointResult(radial,"radial-provenance");
    repository.SaveModel(model,"radial-provenance");
    const auto radial_loaded=repository.LoadModel("radial-provenance");
    ASSERT_TRUE(radial_loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(radial_loaded->GetAnalysisView().GetJointResult()->selection_domain->fixed_halo_selection,
        radial.selection_domain->fixed_halo_selection);
}

TEST(JointComponentPartialSelectionTest, InitializationExceptionsKeepProvenanceWithMedianSeed)
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
    const auto problem = core::BuildJointProblem(map, model);
    const auto workset = core::detail::MakeJointFittingWorkset(model, problem);
    core::FitOptions options; options.thread_size = 1;
    const auto first_stage_atoms = core::detail::CollectFirstStageAtoms(workset);
    core::RunPotentialSamplingWorkflow(map, model, first_stage_atoms, options.sampling_method, 1);
    ResetProgressPercentForTest();
    testing::internal::CaptureStdout();
    const auto initialization = core::detail::RunJointFirstStageInitializationFromPreparedSamples(
        model, workset, options);
    const auto progress = testing::internal::GetCapturedStdout();
    EXPECT_NE(progress.find("50%"), std::string::npos);
    EXPECT_NE(progress.find("100%"), std::string::npos);
    EXPECT_EQ(progress.find("Run local alpha training for 1 atoms."), std::string::npos);
    EXPECT_EQ(std::count(progress.begin(), progress.end(), '\r'), 2);
    const auto fit = core::FitJointComponents(problem, initialization.b);
    ASSERT_EQ(initialization.atoms.size(),2);
    EXPECT_TRUE(initialization.atoms[0].reason.starts_with("initialization-exception:"));
    EXPECT_EQ(initialization.atoms[1].reason,"valid-width");
    ASSERT_EQ(fit.components.size(),2);
    for(const auto & component:fit.components)
    {
        ASSERT_TRUE(component.state);
        ASSERT_TRUE(component.candidate_diagnostic);
        EXPECT_TRUE(component.candidate_diagnostic->candidate_available);
        EXPECT_EQ(component.candidate_diagnostic->endpoint_verification_status,core::JointCheckStatus::Passed);
        EXPECT_EQ(component.candidate_diagnostic->numerical_qualification_status,core::JointCheckStatus::Passed);
    }
    EXPECT_TRUE(initialization.valid); EXPECT_TRUE(fit.objective);
    EXPECT_EQ(fit.RuntimeConvergence(),core::JointCheckStatus::Passed);
    EXPECT_EQ(initialization.atoms[0].seed_source,"median-fallback");
    EXPECT_EQ(initialization.atoms[0].donor_count,1);
    EXPECT_FALSE(initialization.atoms[0].original_b);
    EXPECT_DOUBLE_EQ(initialization.b[0],initialization.b[1]);
    EXPECT_EQ(model.GetSelectedAtomCount(),2);
    auto saved=core::CaptureJointAnalysisResult(fit);
    saved.initialization = initialization;
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    EXPECT_EQ(decoded.initialization.atoms[0].reason,initialization.atoms[0].reason);
    EXPECT_EQ(decoded.available_row_mask,fit.available_row_mask);
}

TEST(JointComponentPartialSelectionTest, LocalAlphaSingleAtomAndBatchPreserveNumericalResults)
{
    auto f = joint_partial_test::Make("all");
    const auto contributors = f.model->GetSelectedAtoms();
    ASSERT_EQ(contributors.size(), 2);
    core::FitOptions options; options.thread_size = 1;
    std::vector<double> expected_alpha;
    {
        auto editor = f.model->EditAnalysis();
        for (auto * atom : contributors)
        {
            editor.SetAtomLocalRawSamplingEntries(*atom,
                core::SampleAtomMapValues(*f.map, *atom, options.sampling_method));
            const auto entries = rhbm_gem::AtomLocalPotentialView::For(*atom)
                .GetSamplingEntries(rhbm_gem::FittingStage::First);
            ASSERT_GE(entries.size(), 10);
            expected_alpha.push_back(core::TrainAlphaR(entries, options));
        }
        core::detail::TrainLocalAlphaForAtom(editor, options,
            rhbm_gem::FittingStage::First, *contributors.front());
    }
    EXPECT_DOUBLE_EQ(rhbm_gem::AtomLocalPotentialView::For(*contributors.front())
        .GetAlphaR(rhbm_gem::FittingStage::First), expected_alpha.front());

    ResetProgressPercentForTest();
    testing::internal::CaptureStdout();
    core::RunLocalAlphaTraining(*f.model, options, rhbm_gem::FittingStage::First, contributors);
    const auto progress = testing::internal::GetCapturedStdout();
    EXPECT_NE(progress.find("Run local alpha training for 2 atoms."), std::string::npos);
    EXPECT_NE(progress.find("50%"), std::string::npos);
    EXPECT_NE(progress.find("100%"), std::string::npos);
    EXPECT_EQ(std::count(progress.begin(), progress.end(), '\r'), 2);
    for (std::size_t i = 0; i < contributors.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(rhbm_gem::AtomLocalPotentialView::For(*contributors[i])
            .GetAlphaR(rhbm_gem::FittingStage::First), expected_alpha[i]);
    }
}

TEST(JointComponentPartialSelectionTest, IsolatedProgressCountsSkippedContributorsOnce)
{
    auto f = joint_partial_test::Make("partial");
    const core::detail::FittingWorkset workset{
        {f.model->FindAtomPtr(1), f.model->FindAtomPtr(2)}, {true, false}};
    core::FitOptions options; options.thread_size = 1;
    std::map<std::pair<int, std::string>, int> calls;
    std::map<int, int> sampling_calls;
    core::detail::FirstStageObserverForTesting() = [&](int id, std::string_view phase) {
        ++calls[{id, std::string(phase)}];
    };
    core::detail::PotentialSamplingObserverForTesting() = [&](int id) { ++sampling_calls[id]; };
    const auto first_stage_atoms = core::detail::CollectFirstStageAtoms(workset);
    core::RunPotentialSamplingWorkflow(*f.map, *f.model, first_stage_atoms, options.sampling_method, 1);
    core::detail::PotentialSamplingObserverForTesting() = {};
    ResetProgressPercentForTest();
    testing::internal::CaptureStdout();
    const auto initialization = core::detail::RunJointFirstStageInitializationFromPreparedSamples(
        *f.model, workset, options);
    const auto progress = testing::internal::GetCapturedStdout();
    core::detail::FirstStageObserverForTesting() = {};

    ASSERT_EQ(initialization.atoms.size(), 2);
    EXPECT_NE(initialization.atoms[0].reason, "not-required-observable-contribution");
    EXPECT_EQ((calls[{1, "first"}]), 1);
    EXPECT_EQ((calls[{2, "first"}]), 1);
    EXPECT_EQ(sampling_calls[1], 1);
    EXPECT_EQ(sampling_calls[2], 1);
    EXPECT_NE(progress.find("50%"), std::string::npos);
    EXPECT_NE(progress.find("100%"), std::string::npos);
    EXPECT_EQ(std::count(progress.begin(), progress.end(), '\r'), 2);

    options.quiet_mode = true;
    testing::internal::CaptureStdout();
    core::detail::RunJointFirstStageInitializationFromPreparedSamples(*f.model, workset, options);
    EXPECT_TRUE(testing::internal::GetCapturedStdout().empty());
}

TEST(JointComponentPartialSelectionTest, WorksetSamplesAndInitializesEveryFullAbcContributor)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for (int id = 1; id <= 3; ++id)
    {
        auto atom = std::make_unique<rhbm_gem::AtomObject>();
        atom->SetSerialID(id);
        atom->SetElement(Element::CARBON);
        atom->SetPosition(0.7 * (id - 1), 0, 0);
        atoms.push_back(std::move(atom));
    }
    rhbm_gem::ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom) { return atom.GetSerialID() == 1; });
    rhbm_gem::MapObject map({17, 17, 17}, {.5, .5, .5}, {-4, -4, -4});

    core::JointProblemInput input;
    input.atom_ids = {"1", "2", "3"};
    input.row_ids = {"0", "1", "2", "3"};
    input.observations = {.3, .5, .7, .9};
    input.support = {
        {{0, .1}, {1, .4}, {2, .9}, {3, 1.6}},
        {{2, .25}, {3, 1.0}},
        {{3, 2.0}}
    };
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0};
    const core::JointProblem problem(std::move(input));
    const auto workset = core::detail::MakeJointFittingWorkset(model, problem);
    ASSERT_EQ(workset.contributors.size(), 3);
    EXPECT_EQ(workset.target_mask, (std::vector<bool>{true, false, false}));
    EXPECT_EQ(problem.ParameterLayout().full_atoms, (std::vector<std::size_t>{0, 1, 2}));
    const auto selected = model.GetSelectedAtoms();
    const auto bonds = model.GetSelectedBonds();

    model.EditAnalysis().InitializeFromSelection();
    core::FitOptions options;
    options.thread_size = 1;
    options.quiet_mode = true;
    std::map<int, int> sampling_calls;
    std::map<int, int> first_fit_calls;
    core::detail::PotentialSamplingObserverForTesting() = [&](int id) { ++sampling_calls[id]; };
    core::detail::FirstStageObserverForTesting() = [&](int id, std::string_view phase) {
        if (phase == "first") ++first_fit_calls[id];
    };
    const auto initialization = core::detail::PrepareJointWidthInitialization(map, model, workset, options);
    core::detail::PotentialSamplingObserverForTesting() = {};
    core::detail::FirstStageObserverForTesting() = {};

    EXPECT_EQ(sampling_calls[1], 1);
    EXPECT_EQ(sampling_calls[2], 1);
    EXPECT_EQ(sampling_calls[3], 1);
    EXPECT_EQ(first_fit_calls[1], 1);
    EXPECT_EQ(first_fit_calls[2], 1);
    EXPECT_EQ(first_fit_calls[3], 1);
    ASSERT_EQ(initialization.atoms.size(), 3);
    EXPECT_NE(initialization.atoms[2].reason, "not-required-observable-contribution");
    EXPECT_EQ(model.GetSelectedAtoms(), selected);
    EXPECT_EQ(model.GetSelectedBonds(), bonds);
}

TEST(JointComponentPartialSelectionTest, ZeroSignalRankDeficiencyAndInsufficientRows)
{
    auto f=joint_partial_test::Make("partial"); const auto p=core::BuildJointProblem(*f.map,*f.model);
    for(int scenario=0;scenario<3;++scenario)
    {
        auto input=p.Input();
        input.selection_domain->contributor_policy="all-non-hydrogen";
        input.selection_domain->halo_eligibility.reset();
        input.selection_domain->fixed_halo_selection.reset();
        if(scenario==0) std::fill(input.observations.begin(),input.observations.end(),0);
        if(scenario==1) input.support[1]=input.support[0];
        if(scenario==2)
        {
            input.row_ids.resize(2); input.observations.resize(2);
            for(auto & support:input.support) support={{0,0},{1,1}};
        }
        const auto fit=core::FitJointComponents(core::JointProblem(input),{.5,.5});
        if(scenario==0) EXPECT_EQ(fit.RuntimeConvergence(),core::JointCheckStatus::Passed);
        else EXPECT_NE(fit.RuntimeConvergence(),core::JointCheckStatus::Passed);
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
    auto f = joint_partial_test::Make("weak");
    f.map=std::make_unique<rhbm_gem::MapObject>(std::array<int,3>{101,41,41},
        std::array<double,3>{.3,.3,.3},std::array<double,3>{-9,-6,-6},
        std::make_unique<double[]>(101*41*41));
    SetFixedHaloExactMap(f);
    const auto selected = f.model->GetSelectedAtoms();
    const auto bonds = f.model->GetSelectedBonds();
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    std::map<std::pair<int, std::string>, int> calls;
    std::map<int, int> sampling_calls;
    core::detail::FirstStageObserverForTesting() = [&](int id, std::string_view phase) {
        ++calls[{id, std::string(phase)}];
    };
    core::detail::PotentialSamplingObserverForTesting() = [&](int id) { ++sampling_calls[id]; };
    core::FitOptions options;
    options.estimator = core::PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode = true;
    core::detail::JointPostprocessingMetrics postprocessing_metrics;
    core::detail::RunJointComponentWorkflow(*f.map, *f.model, options, &postprocessing_metrics);
    core::detail::FirstStageObserverForTesting() = {};
    core::detail::PotentialSamplingObserverForTesting() = {};
    ASSERT_TRUE(f.model->GetAnalysisView().GetJointResult());
    const auto & result = *f.model->GetAnalysisView().GetJointResult();
    ASSERT_TRUE(result.assembled_state) << result.initialization.reason << " components=" << result.components.size()
        << (result.components.empty() ? std::string{} : ", stop="+result.components.front().stop_reason);
    ASSERT_EQ(result.fixed_halos.size(),1u);
    const auto expected_stages = rhbm_gem::data_internal::BuildJointStageEstimates(result, "expected-run");
    EXPECT_EQ(f.model->GetSelectedAtoms(), selected);
    EXPECT_EQ(f.model->GetSelectedBonds(), bonds);
    EXPECT_EQ(postprocessing_metrics.snapshot_identity_validation_count, 1u);
    EXPECT_EQ(postprocessing_metrics.full_row_identity_scan_count, 1u);
    EXPECT_EQ(postprocessing_metrics.observation_row_mapping_builds, 1u);
    EXPECT_EQ(postprocessing_metrics.contributor_mapping_builds, 1u);
    EXPECT_EQ(postprocessing_metrics.joint_stage_mapping_builds, 1u);
    EXPECT_EQ(postprocessing_metrics.target_preparation_count, selected.size());
    EXPECT_EQ(result.row_ids, problem.Input().row_ids);
    ASSERT_EQ(result.initialization.atoms.size(), problem.Input().atom_ids.size());
    std::size_t local_regression_count{};
    for (std::size_t i = 0; i < result.atom_ids.size(); ++i)
    {
        const auto id = std::stoi(result.atom_ids[i]);
        EXPECT_EQ((calls[{id, "first"}]), 1);
        EXPECT_EQ(sampling_calls[id], 1);
        const auto view = rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(id));
        EXPECT_DOUBLE_EQ(view.GetFinalModel(rhbm_gem::FittingStage::First).GetWidth(),
            result.initialization.atoms[i].mdpde[1]);
        EXPECT_FALSE(view.GetRawSamplingEntries(false).empty());
        const auto & stage = view.GetStageEstimate(rhbm_gem::FittingStage::Second);
        const auto & expected = expected_stages.at(id);
        EXPECT_EQ(stage.point.has_value(), expected.point.has_value());
        if (stage.point && expected.point) EXPECT_EQ(stage.point->ToVector(), expected.point->ToVector());
        EXPECT_EQ(stage.uncertainty.status, rhbm_gem::EvidenceStatus::NotRun);
        EXPECT_FALSE(stage.uncertainty.covariance);
        const bool target = stage.source.role == rhbm_gem::FittingRole::Target;
        EXPECT_EQ(view.GetPostFitPeeling().has_value(), target);
        EXPECT_EQ(view.GetPostJointLocalResult().has_value(), target);
        if (target)
        {
            EXPECT_EQ(view.GetPostFitPeeling()->samples.size(), view.GetRawSamplingEntries(false).size());
            const auto & local = *view.GetPostJointLocalResult();
            EXPECT_EQ(local.source.run_id, stage.source.run_id);
            EXPECT_EQ(local.source.component_id, stage.source.component_id);
            if (stage.point)
            {
                EXPECT_DOUBLE_EQ(local.fixed_charge.coefficient, stage.point->GetOffset());
                EXPECT_DOUBLE_EQ(local.fixed_charge.width, stage.point->GetWidth());
            }
            if (local.status == rhbm_gem::EvidenceStatus::Available)
            {
                const auto prepared = core::detail::BuildPostJointLocalSamples(
                    *f.map, *f.model, problem, id, *view.GetPostFitPeeling(), result);
                const auto gaussian_only = prepared.BuildGaussianSamples(view.GetRawSamplingEntriesView());
                const double alpha = core::TrainAlphaR(gaussian_only, options);
                const auto direct = core::EstimateLocalGaussian(
                    gaussian_only, alpha, options, rhbm_gem::GaussianModel3D{0.0, 1.0, 0.0});
                ASSERT_TRUE(direct.fit_result);
                EXPECT_DOUBLE_EQ(local.fit.alpha_r, direct.alpha_r);
                EXPECT_EQ(local.fit.ols.GetModel().ToVector(), direct.ols.GetModel().ToVector());
                EXPECT_EQ(local.fit.mdpde.GetModel().ToVector(), direct.mdpde.GetModel().ToVector());
                ASSERT_TRUE(local.fit.diagnostics);
                EXPECT_EQ(local.fit.diagnostics->status, direct.fit_result->status);
                EXPECT_DOUBLE_EQ(local.fit.diagnostics->sigma_square, direct.fit_result->sigma_square);
                EXPECT_EQ(local.fit.diagnostics->iterations.iterations, direct.fit_result->diagnostics.iterations);
                EXPECT_EQ(local.fit.diagnostics->refinement.has_value(), direct.fit_result->refinement.has_value());
                ++local_regression_count;
            }
        }
        EXPECT_FALSE(view.GetGroupMemberResult());
    }
    EXPECT_GT(local_regression_count, 0u);
    const auto initialized_problem=problem.WithFixedHalos(result.fixed_halos);
    const auto direct = core::FitJointComponents(initialized_problem, result.initialization.b);
    EXPECT_EQ(result.objective, direct.objective);
    ASSERT_EQ(result.components.size(), direct.components.size());
    for (std::size_t c = 0; c < direct.components.size(); ++c)
    {
        ASSERT_EQ(result.components[c].state.has_value(), direct.components[c].state.has_value());
        if (!direct.components[c].state) continue;
        EXPECT_EQ(result.components[c].state->ac, direct.components[c].state->ac);
        EXPECT_EQ(result.components[c].state->b, direct.components[c].state->b);
        EXPECT_EQ(result.components[c].state->objective, direct.components[c].state->objective);
    }
}

TEST(JointComponentPartialSelectionTest, UnassignedTargetDoesNotEnterPostFitPeelingBatch)
{
    auto f = joint_partial_test::Make("weak");
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for (const auto & atom : f.model->GetAtomList())
        atoms.push_back(std::make_unique<rhbm_gem::AtomObject>(*atom));
    auto unobserved = std::make_unique<rhbm_gem::AtomObject>();
    unobserved->SetSerialID(3); unobserved->SetElement(Element::CARBON);
    unobserved->SetPosition(100.0, 0.0, 0.0); unobserved->SetChainID("A");
    unobserved->SetComponentID("ALA"); unobserved->SetAtomID("CG");
    atoms.push_back(std::move(unobserved));
    f.model = std::make_unique<rhbm_gem::ModelObject>(std::move(atoms));
    f.model->SelectAtoms([](const auto & atom) { return atom.GetSerialID() != 2; });

    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    ASSERT_EQ(problem.Input().atom_ids, (std::vector<std::string>{"1", "2", "3"}));
    ASSERT_TRUE(problem.Input().selection_domain);
    EXPECT_EQ(problem.Input().selection_domain->target_indices, (std::vector<std::size_t>{0, 2}));

    rhbm_gem::JointAnalysisResult result;
    result.parameterization_contract = std::string(rhbm_gem::JointTargetOnlyFixedHaloParameterizationContract);
    result.atom_ids = problem.Input().atom_ids;
    result.row_ids = problem.Input().row_ids;
    result.selection_domain = problem.Input().selection_domain;
    result.fixed_halos = {FixedHaloFor(f, 1)};
    rhbm_gem::JointAnalysisComponent component;
    component.id = "target-component"; component.atoms = {0};
    component.state = rhbm_gem::JointState{{f.a[0], f.c[0]}, {f.b[0]}, {std::log(f.b[0])}, {0.0}, 0.0};
    result.components.push_back(component);

    f.model->EditAnalysis().ApplySecondStageEstimates(
        rhbm_gem::data_internal::BuildJointStageEstimates(result, "unassigned-target-test"));
    f.model->EditAnalysis().SetJointResult(result);
    const core::detail::JointPostprocessingContext postprocessing(*f.map, *f.model, problem, &result);
    auto peelings = postprocessing.BuildPostFitPeelingSamples();
    ASSERT_EQ(peelings.size(), 1u);
    EXPECT_TRUE(peelings.contains(1));
    EXPECT_FALSE(peelings.contains(3));
    std::vector<std::pair<int, rhbm_gem::PostFitPeelingResult>> updates;
    for (auto & [id, peeling] : peelings) updates.emplace_back(id, std::move(peeling));
    EXPECT_NO_THROW(f.model->EditAnalysis().ApplyPostFitPeelingsBatch(std::move(updates)));
    EXPECT_TRUE(rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(1)).GetPostFitPeeling());
    EXPECT_FALSE(rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(3)).GetPostFitPeeling());
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
    EXPECT_NE(summary.find("\n| Spot | valid | not-converged | unavailable | A mean / s.d. | B mean / s.d. | C(J) mean / s.d. |"), std::string::npos);
    EXPECT_NE(summary.find("\n| CA   |"), std::string::npos);
    EXPECT_EQ(summary.find("\n| CB   |"), std::string::npos);
    const auto expectSummaryRowWidths = [](const std::string & output, bool unavailable) {
        const auto row_start = output.find("\n| CA   |");
        ASSERT_NE(row_start, std::string::npos);
        const auto row_end = output.find('\n', row_start + 1);
        const auto row = output.substr(row_start + 1, row_end - row_start - 1);
        std::array<std::size_t, 8> pipe_positions{};
        std::size_t pipe_count = 0;
        for (auto position = row.find('|'); position != std::string::npos; position = row.find('|', position + 1))
        {
            ASSERT_LT(pipe_count, pipe_positions.size());
            pipe_positions[pipe_count++] = position;
        }
        ASSERT_EQ(pipe_count, pipe_positions.size());
        EXPECT_EQ(pipe_positions[1] - pipe_positions[0], 7);
        EXPECT_EQ(pipe_positions[2] - pipe_positions[1], 8);
        EXPECT_EQ(pipe_positions[3] - pipe_positions[2], 16);
        EXPECT_EQ(pipe_positions[4] - pipe_positions[3], 14);
        for (std::size_t column = 4; column < 7; ++column)
        {
            EXPECT_EQ(pipe_positions[column + 1] - pipe_positions[column], 16);
            if (unavailable)
                EXPECT_EQ(row.substr(pipe_positions[column] + 2, 11), "unavailable");
            else
                EXPECT_EQ(row[pipe_positions[column] + 9], '/');
        }
    };
    expectSummaryRowWidths(summary, true);
    f.model->EditAnalysis().SetJointResult(result);
    f.model->EditAnalysis().ClearJointResult();
    EXPECT_TRUE(target.HasFinalModel(rhbm_gem::FittingStage::Second));
    component.state.reset(); component.stop_reason = "invalid-initial-widths";
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, result, "next-run");
    EXPECT_FALSE(target.HasFinalModel(rhbm_gem::FittingStage::Second));
    EXPECT_EQ(target.GetStageEstimate(rhbm_gem::FittingStage::Second).reason, "invalid-initial-widths");
    summary = core::BuildSecondStageSpotSummary(*f.model);
    EXPECT_NE(summary.find("\n| CA   |     0 |             0 |           1 |"), std::string::npos);
    expectSummaryRowWidths(summary, true);
}

TEST(JointComponentPartialSelectionTest, PostFitPeelingUsesGridOperatorAndRetainsMissingCoverage)
{
    auto f = joint_partial_test::Make("partial");
    SetFixedHaloExactMap(f);
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    const auto initialized=problem.WithFixedHalos({FixedHaloFor(f,1)});
    const auto fit=core::FitJointComponents(initialized,{f.b[0]});
    ASSERT_TRUE(fit.assembled_state) << fit.initialization.reason << " search=" << fit.search_completed
        << " components=" << fit.components.size()
        << (fit.components.empty() ? std::string{} : ", stop="+fit.components.front().stop_reason);
    const auto snapshot=core::CaptureJointAnalysisResult(fit);
    ASSERT_TRUE(snapshot.components.front().state);
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, snapshot, "peeling-test");
    auto * atom = f.model->FindAtomPtr(1);
    auto editor = f.model->EditAnalysis();
    const SamplingPointList points{
        {0.2, {0.13, 0.07, 0.02}, false},
        {2.4, {2.4, 0, 0}, true},
        {2.9, {2.9, 0, 0}, true}};
    const auto raw = core::detail::BuildLocalPotentialSampleList(*f.map, points);
    editor.SetAtomLocalRawSamplingEntries(*atom, raw);
    const auto view = rhbm_gem::AtomLocalPotentialView::For(*atom);
    const auto before = view.GetFinalModel(rhbm_gem::FittingStage::Second).ToVector();
    const auto outputs = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem,std::nullopt,&snapshot);
    const auto & output = outputs.at(1);
    ASSERT_EQ(output.samples.size(), raw.size());
    ASSERT_TRUE(output.samples[0].response);
    EXPECT_LT(*output.samples[0].response, raw[0].response);
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
    EXPECT_GT(*view.GetLocalFittingPeelingRatio(0, 1), 0);
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

TEST(JointComponentPartialSelectionTest, ExcludedFitRowsAreUnavailableToPostFitPeeling)
{
    auto f=joint_partial_test::Make("partial");
    SetFixedHaloExactMap(f);
    const auto original=core::BuildJointProblem(*f.map,*f.model);
    const auto retained=joint_test::halo_exclusion::HaloSupportExclusionRows(original.Input());
    auto masked=joint_test::halo_exclusion::BuildMaskedInput(original.Input(),retained,false);
    auto fit_problem=core::JointProblem(std::move(masked.input));
    fit_problem=core::JointProblemAccess::WithReferenceNormalization(fit_problem,original);
    ASSERT_LT(fit_problem.Input().row_ids.size(),original.Input().row_ids.size());

    const auto fit=core::FitJointComponents(fit_problem,{f.b[0]});
    ASSERT_TRUE(fit.assembled_state) << fit.initialization.reason;
    ASSERT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_EQ(fit.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    const auto snapshot=core::CaptureJointAnalysisResult(fit);
    f.model->EditAnalysis().ApplyJointResult(snapshot,"masked-post-fit-domain-test");
    const SamplingPointList points{{0.0,{0.0,0.0,0.0},true}};
    f.model->EditAnalysis().SetAtomLocalRawSamplingEntries(*f.model->FindAtomPtr(1),
        core::detail::BuildLocalPotentialSampleList(*f.map,points));

    const core::detail::JointPostprocessingContext context(*f.map,*f.model,fit_problem,&snapshot);
    const auto peelings=context.BuildPostFitPeelingSamples();
    ASSERT_TRUE(peelings.contains(1));
    ASSERT_EQ(peelings.at(1).samples.size(),1u);
    EXPECT_FALSE(peelings.at(1).samples.front().response);
    EXPECT_EQ(peelings.at(1).samples.front().reason,"outside-joint-domain");
    EXPECT_THROW((core::detail::JointPostprocessingContext(*f.map,*f.model,original,&snapshot)),
        std::invalid_argument);
}

TEST(JointComponentPartialSelectionTest, PeelingIncludesDistantHaloThroughInterpolationStencil)
{
    auto f = joint_partial_test::Make("weak");
    // Binary-exact spacing makes the y/z zero-weight stencil nodes exact.
    f.map = std::make_unique<rhbm_gem::MapObject>(std::array<int, 3>{111, 61, 61},
        std::array<double, 3>{.1, .1, .1}, std::array<double, 3>{-5.5, -3, -3},
        std::make_unique<double[]>(111 * 61 * 61));
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
    auto raw = core::detail::BuildLocalPotentialSampleList(*f.map, points);
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

TEST(JointComponentPartialSelectionTest, SharedPostprocessingContextValidatesIdentityAndBuildsMappingsOnce)
{
    auto f = joint_partial_test::Make("all");
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    const auto snapshot = core::CaptureJointAnalysisResult(core::FitJointComponents(problem, f.b));
    f.model->EditAnalysis().ApplyJointResult(snapshot, "shared-context-test");
    for (const auto * atom : f.model->GetSelectedAtoms())
        f.model->EditAnalysis().SetAtomLocalRawSamplingEntries(*atom,
            core::SampleAtomMapValues(*f.map, *atom, SphereSamplingMethod::FibonacciDeterministic));

    core::detail::JointPostprocessingMetrics metrics;
    const core::detail::JointPostprocessingContext context(
        *f.map, *f.model, problem, &snapshot, &metrics);
    EXPECT_EQ(metrics.snapshot_identity_validation_count, 1u);
    EXPECT_EQ(metrics.full_row_identity_scan_count, 1u);
    EXPECT_EQ(metrics.observation_row_mapping_builds, 1u);
    EXPECT_EQ(metrics.contributor_mapping_builds, 1u);
    const auto peelings = context.BuildPostFitPeelingSamples(
        problem.Input().selection_domain->target_indices);
    ASSERT_EQ(peelings.size(), 2u);
    EXPECT_EQ(metrics.snapshot_identity_validation_count, 1u);
    EXPECT_EQ(metrics.full_row_identity_scan_count, 1u);
    EXPECT_EQ(metrics.observation_row_mapping_builds, 1u);
    EXPECT_EQ(metrics.contributor_mapping_builds, 1u);
    EXPECT_EQ(metrics.target_preparation_count, 0u);

    const auto standalone_peelings = core::detail::BuildPostFitPeelingSamples(
        *f.map, *f.model, problem, problem.Input().selection_domain->target_indices, &snapshot);
    std::size_t expected_prepared_samples{};
    for (const auto * atom : f.model->GetSelectedAtoms())
    {
        const int id = atom->GetSerialID();
        const auto & shared_peeling = peelings.at(id);
        const auto & standalone_peeling = standalone_peelings.at(id);
        const auto local = rhbm_gem::AtomLocalPotentialView::For(*atom);
        const auto & raw_samples = local.GetRawSamplingEntriesView();
        ASSERT_EQ(shared_peeling.samples.size(), standalone_peeling.samples.size());
        for (std::size_t i = 0; i < shared_peeling.samples.size(); ++i)
        {
            EXPECT_EQ(shared_peeling.samples[i].response, standalone_peeling.samples[i].response);
            EXPECT_EQ(shared_peeling.samples[i].reason, standalone_peeling.samples[i].reason);
        }

        const auto shared = context.BuildPostJointLocalSamples(id, shared_peeling);
        const auto standalone = core::detail::BuildPostJointLocalSamples(
            *f.map, *f.model, problem, id, standalone_peeling, snapshot);
        EXPECT_EQ(shared.samples.size(), raw_samples.size());
        ASSERT_EQ(shared.samples.size(), standalone.samples.size());
        EXPECT_DOUBLE_EQ(shared.fixed_charge.coefficient, standalone.fixed_charge.coefficient);
        EXPECT_DOUBLE_EQ(shared.fixed_charge.width, standalone.fixed_charge.width);
        for (std::size_t i = 0; i < shared.samples.size(); ++i)
        {
            EXPECT_EQ(shared.samples[i].fixed_charge_response, standalone.samples[i].fixed_charge_response);
            EXPECT_EQ(shared.samples[i].gaussian_response, standalone.samples[i].gaussian_response);
            EXPECT_EQ(shared.samples[i].reason, standalone.samples[i].reason);
        }
        expected_prepared_samples += shared.samples.size();
        EXPECT_EQ(metrics.snapshot_identity_validation_count, 1u);
        EXPECT_EQ(metrics.full_row_identity_scan_count, 1u);
        EXPECT_EQ(metrics.observation_row_mapping_builds, 1u);
        EXPECT_EQ(metrics.contributor_mapping_builds, 1u);
    }
    EXPECT_EQ(metrics.snapshot_identity_validation_count, 1u);
    EXPECT_EQ(metrics.full_row_identity_scan_count, 1u);
    EXPECT_EQ(metrics.observation_row_mapping_builds, 1u);
    EXPECT_EQ(metrics.contributor_mapping_builds, 1u);
    EXPECT_EQ(metrics.target_preparation_count, 2u);
    EXPECT_EQ(metrics.total_prepared_sample_count, expected_prepared_samples);
    EXPECT_EQ(metrics.total_prepared_sample_count,
        metrics.valid_local_sample_count + metrics.invalid_local_sample_count);
}

TEST(JointComponentPartialSelectionTest, PostJointEvidenceAssessmentEnforcesAvailableContract)
{
    using namespace rhbm_gem;
    LocalStageEstimate joint;
    joint.source = {EstimateMethod::JointComponents, "1", "component", "validation-run", FittingRole::Target};
    joint.point = GaussianModel3D{2.0, 0.5, -0.1};
    joint.convergence = JointCheckStatus::Passed;
    LocalPotentialSampleList raw(4);
    for (std::size_t i = 0; i < raw.size(); ++i)
        raw[i].point.distance = 0.1 + 0.2 * static_cast<double>(i);
    PostFitPeelingResult peeling;
    peeling.source = joint.source;
    peeling.samples.resize(raw.size());
    for (auto & sample : peeling.samples) sample.response = 1.0;

    PostJointLocalResult result;
    result.source = joint.source;
    result.status = EvidenceStatus::Available;
    result.reason.clear();
    result.fixed_charge = {joint.point->GetOffset(), joint.point->GetWidth()};
    result.sample_count = raw.size();
    result.effective_sample_count = 3;
    result.fit.alpha_r = 0.2;
    result.fit.ols = {GaussianModel3D{2.1, 0.47, 0.0}, {}};
    result.fit.mdpde = {GaussianModel3D{2.0, 0.48, 0.0}, {}};
    result.fit.diagnostics = LocalFitDiagnostics{RHBMEstimationStatus::SUCCESS, 0.1, {8, 1e-8, 1e-8}, {}, std::nullopt};
    RHBMInformation information;
    information.precision = 20.0 * Eigen::Matrix2d::Identity();
    information.moment = information.precision *
        linearization_service::EncodeGaussianToParameterVector(result.fit.mdpde.GetModel());
    result.information = information;

    const auto assess = [&](const PostJointLocalResult & candidate, const LocalStageEstimate & endpoint,
                            const PostFitPeelingResult * paired) {
        return data_internal::AssessPostJointLocalResult(candidate, 1, endpoint, paired, raw);
    };
    const auto assessCurrent = [&](const PostJointLocalResult & candidate) {
        return assess(candidate, joint, &peeling);
    };
    const auto valid = assessCurrent(result);
    EXPECT_TRUE(valid.valid) << valid.reason;
    EXPECT_TRUE(valid.group_eligible);

    auto bad = result;
    EXPECT_FALSE(assess(bad, joint, nullptr).valid);
    auto unqualified = joint;
    unqualified.convergence = JointCheckStatus::NotRun;
    EXPECT_FALSE(assess(result, unqualified, &peeling).valid);
    bad = result; bad.source.run_id = "different-run";
    EXPECT_FALSE(assessCurrent(bad).valid);
    auto halo = joint; halo.source.role = FittingRole::Halo;
    EXPECT_FALSE(assess(result, halo, &peeling).valid);
    bad = result; bad.fit.mdpde = {GaussianModel3D{
        std::numeric_limits<double>::quiet_NaN(), 0.48, 0.0}, {}};
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.fit.alpha_r = -0.1;
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.fit.diagnostics->sigma_square = 0.0;
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.effective_sample_count = 2;
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.fit.diagnostics.reset();
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.information->precision = Eigen::Matrix2d::Zero();
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.information->moment(0) = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result;
    bad.sample_failures = {{0, raw[0].point.distance, "invalid"}, {0, raw[0].point.distance, "duplicate"}};
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.fixed_charge.width += 0.01;
    EXPECT_FALSE(assessCurrent(bad).valid);
    bad = result; bad.reason = "available-with-reason";
    EXPECT_FALSE(assessCurrent(bad).valid);

    bad = result;
    bad.status = EvidenceStatus::Unavailable;
    bad.reason = "insufficient-positive-responses";
    bad.information.reset();
    const auto unavailable = assessCurrent(bad);
    EXPECT_TRUE(unavailable.valid) << unavailable.reason;
    EXPECT_FALSE(unavailable.group_eligible);
}

TEST(JointComponentPartialSelectionTest, GroupPosteriorUsesPostJointLocalInformationAndPreservesJoint)
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
        stage.source = {EstimateMethod::JointComponents, std::to_string(i), "component", "local-test-run", FittingRole::Target};
        stage.point = GaussianModel3D(2 + i*0.2, 0.5 + 0.02*(i*i % 7), -0.1*i);
        stage.convergence = JointCheckStatus::Passed;
        stages[i] = std::move(stage);
    }
    editor.ApplySecondStageEstimates(stages);

    std::vector<std::pair<int, PostFitPeelingResult>> peelings;
    std::vector<std::pair<int, PostJointLocalResult>> local_results;
    for (const auto & [id, unused] : stages)
    {
        (void)unused;
        const auto & stage = AtomLocalPotentialView::For(*model.FindAtomPtr(id))
            .GetStageEstimate(FittingStage::Second);
        LocalPotentialSampleList raw;
        PostFitPeelingResult peeling;
        peeling.source = stage.source;
        for (int sample_index = 0; sample_index < 6; ++sample_index)
        {
            LocalPotentialSample sample;
            sample.point.distance = 0.12 + 0.43 * sample_index;
            sample.response = 2.0 + 0.1 * id + 0.02 * sample_index;
            raw.push_back(sample);
            peeling.samples.push_back({sample.response, {}});
        }
        editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(id), std::move(raw));
        peelings.emplace_back(id, std::move(peeling));

        PostJointLocalResult local;
        local.source = stage.source;
        local.fixed_charge = {stage.point->GetOffset(), stage.point->GetWidth()};
        local.sample_count = 6;
        local.effective_sample_count = 6;
        local.fit.alpha_r = 0.2;
        local.fit.ols = {GaussianModel3D(1.9 + 0.1 * id, 0.40 + 0.01 * id, 0.0), {}};
        local.fit.mdpde = {GaussianModel3D(2.0 + 0.1 * id, 0.42 + 0.01 * id, 0.0), {}};
        local.fit.diagnostics = LocalFitDiagnostics{RHBMEstimationStatus::SUCCESS, 0.1, {8, 1e-8, 1e-8}, {}, std::nullopt};
        if (id == 3)
        {
            local.status = EvidenceStatus::Unavailable;
            local.reason = "local-mdpde-unqualified";
        }
        else
        {
            local.status = EvidenceStatus::Available;
            local.reason.clear();
            RHBMInformation information;
            information.precision = 20.0 * Eigen::Matrix2d::Identity();
            const auto beta = linearization_service::EncodeGaussianToParameterVector(local.fit.mdpde.GetModel());
            information.moment = information.precision * beta;
            local.information = std::move(information);
        }
        local_results.emplace_back(id, std::move(local));
    }

    ResetGroupInvalidationStatsForTesting();
    editor.ApplyPostFitPeelingsBatch(std::move(peelings));
    const auto peeling_stats = GetGroupInvalidationStatsForTesting();
    EXPECT_EQ(peeling_stats.membership_traversals, 1u);
    EXPECT_EQ(peeling_stats.groups_invalidated, 1u);
    editor.ApplyPostJointLocalResultsBatch(local_results);

    core::FitOptions options; options.estimator = core::PotentialEstimator::JOINT_COMPONENTS; options.quiet_mode = true;
    core::RunGroupPotentialFitting(model, options);
    const auto keys = model.GetAnalysisView().CollectAtomGroupKeys();
    ASSERT_EQ(keys.size(), 1);
    auto summary = *model.GetAnalysisView().GetGroupParameterSummary(keys[0]);
    ASSERT_EQ(summary.status, EvidenceStatus::Available) << summary.reason;
    EXPECT_EQ(summary.point_count, 4); EXPECT_EQ(summary.eligible_count, 4); EXPECT_EQ(summary.excluded_count, 1);
    EXPECT_EQ(summary.member_ids, (std::vector<int>{1,2,4,5}));
    EXPECT_FALSE(AtomLocalPotentialView::For(*model.FindAtomPtr(3)).GetGroupMemberResult());
    ASSERT_TRUE(summary.descriptive_charge_median);
    EXPECT_DOUBLE_EQ(*summary.descriptive_charge_median, -0.3);
    ASSERT_TRUE(summary.inference);
    EXPECT_DOUBLE_EQ(summary.inference->mean.GetOffset(), 0.0);
    EXPECT_DOUBLE_EQ(summary.inference->mdpde.GetOffset(), 0.0);
    EXPECT_DOUBLE_EQ(summary.inference->prior.GetModel().GetOffset(), 0.0);

    const auto first = AtomLocalPotentialView::For(*model.FindAtomPtr(1));
    ASSERT_TRUE(first.GetPostJointLocalResult());
    const auto joint_point = first.GetFinalModel(FittingStage::Second).ToVector();
    const auto local_point = first.GetPostJointLocalResult()->fit.mdpde.GetModel().ToVector();
    EXPECT_FALSE(joint_point.isApprox(local_point, 1.0e-10));
    const auto posterior = first.GetGroupMemberResult()->posterior.GetModel().ToVector();
    const auto evaluateMember = [](const GroupGaussianMemberResult & member, double distance) {
        EXPECT_TRUE(member.fixed_charge);
        if (!member.fixed_charge) return;
        const auto & gaussian = member.posterior.GetModel();
        const double charge = core::joint_component::EvaluateKernel(
            distance * distance, member.fixed_charge->width, 2.5).charge;
        EXPECT_DOUBLE_EQ(core::EvaluateGroupMemberResponseAtDistance(member, distance),
            gaussian.SignalAtDistance(distance) + member.fixed_charge->coefficient * charge);
    };
    for (const int id : summary.member_ids)
    {
        const auto member = AtomLocalPotentialView::For(*model.FindAtomPtr(id)).GetGroupMemberResult();
        ASSERT_TRUE(member);
        ASSERT_TRUE(member->fixed_charge);
        const auto local = AtomLocalPotentialView::For(*model.FindAtomPtr(id)).GetPostJointLocalResult();
        ASSERT_TRUE(local);
        EXPECT_DOUBLE_EQ(member->fixed_charge->coefficient, local->fixed_charge.coefficient);
        EXPECT_DOUBLE_EQ(member->fixed_charge->width, local->fixed_charge.width);
        EXPECT_DOUBLE_EQ(member->posterior.GetModel().GetOffset(), 0.0);
        evaluateMember(*member, 0.0);
        evaluateMember(*member, 1.2);
        evaluateMember(*member, 2.5);
        evaluateMember(*member, 2.6);
    }
    EXPECT_NE(first.GetGroupMemberResult()->fixed_charge->width,
        first.GetGroupMemberResult()->posterior.GetModel().GetWidth());

    const auto before_invalid_batch = rhbm_gem::stage_result_io::Encode(model);
    auto invalid_batch = local_results;
    invalid_batch[0].second.fit.mdpde = {GaussianModel3D(9.0, 0.7, 0.0), {}};
    invalid_batch[1].second.information->precision << -1.0, 0.0, 0.0, 1.0;
    EXPECT_THROW(editor.ApplyPostJointLocalResultsBatch(std::move(invalid_batch)), std::invalid_argument);
    EXPECT_EQ(rhbm_gem::stage_result_io::Encode(model), before_invalid_batch);
    EXPECT_THROW(editor.ApplyPostJointLocalResultsBatch({local_results[0], local_results[0]}), std::invalid_argument);
    EXPECT_EQ(rhbm_gem::stage_result_io::Encode(model), before_invalid_batch);

    {
        auto stale = *first.GetPostJointLocalResult();
        stale.source.run_id = "stale-run";
        EXPECT_THROW(editor.SetAtomPostJointLocalResult(*model.FindAtomPtr(1), std::move(stale)), std::invalid_argument);
        EXPECT_TRUE(model.GetAnalysisView().GetGroupParameterSummary(keys[0]));
        EXPECT_EQ(first.GetGroupMemberResult()->posterior.GetModel().ToVector(), posterior);
    }
    {
        auto stale = *first.GetPostFitPeeling();
        stale.source.run_id = "stale-run";
        EXPECT_THROW(editor.ApplyPostFitPeelingsBatch({{1, std::move(stale)}}), std::invalid_argument);
        EXPECT_TRUE(first.GetPostJointLocalResult());
    }

    std::vector<std::pair<int, PostFitPeelingResult>> republished_peelings;
    for (int id = 1; id <= 5; ++id)
        republished_peelings.emplace_back(id, *AtomLocalPotentialView::For(*model.FindAtomPtr(id)).GetPostFitPeeling());
    editor.ApplyPostFitPeelingsBatch(std::move(republished_peelings));
    EXPECT_FALSE(first.GetPostJointLocalResult());
    EXPECT_FALSE(first.GetGroupMemberResult());
    EXPECT_FALSE(model.GetAnalysisView().GetGroupParameterSummary(keys[0]));
    editor.ApplyPostJointLocalResultsBatch(local_results);
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(first.GetGroupMemberResult()->posterior.GetModel().ToVector(), posterior);

    {
        const command_test::ScopedTempDir directory{"post_joint_local_group_roundtrip"};
        DataRepository repository(directory.path() / "model.sqlite");
        repository.SaveModel(model, "model");
        const auto loaded = repository.LoadModel("model");
        EXPECT_EQ(rhbm_gem::stage_result_io::Encode(*loaded), rhbm_gem::stage_result_io::Encode(model));
        const auto saved = AtomLocalPotentialView::For(*loaded->FindAtomPtr(1));
        ASSERT_TRUE(saved.GetPostJointLocalResult());
        EXPECT_EQ(saved.GetPostJointLocalResult()->source.run_id, "local-test-run");
        EXPECT_EQ(saved.GetPostJointLocalResult()->status, EvidenceStatus::Available);
        ASSERT_TRUE(saved.GetPostJointLocalResult()->information);
        EXPECT_EQ(saved.GetPostJointLocalResult()->information->precision,
            first.GetPostJointLocalResult()->information->precision);
        ASSERT_TRUE(saved.GetGroupMemberResult());
        EXPECT_EQ(saved.GetGroupMemberResult()->posterior.GetModel().ToVector(), posterior);
        ASSERT_TRUE(saved.GetGroupMemberResult()->fixed_charge);
        EXPECT_DOUBLE_EQ(saved.GetGroupMemberResult()->fixed_charge->coefficient,
            first.GetGroupMemberResult()->fixed_charge->coefficient);
        EXPECT_DOUBLE_EQ(saved.GetGroupMemberResult()->fixed_charge->width,
            first.GetGroupMemberResult()->fixed_charge->width);
        EXPECT_DOUBLE_EQ(core::EvaluateGroupMemberResponseAtDistance(*saved.GetGroupMemberResult(), 1.2),
            core::EvaluateGroupMemberResponseAtDistance(*first.GetGroupMemberResult(), 1.2));
    }

    auto version_three = boost::json::parse(rhbm_gem::stage_result_io::Encode(model)).as_object();
    version_three["version"] = 3;
    for (auto & value : version_three["atoms"].as_array())
    {
        auto & atom_record = value.as_object();
        if (atom_record.at("posterior").is_null()) continue;
        auto & posterior_record = atom_record.at("posterior").as_object();
        const auto & charge = atom_record.at("post_joint_local").as_object().at("fixed_charge").as_array();
        posterior_record.at("posterior").as_object().at("point").as_array()[2] = charge[0];
        posterior_record.erase("fixed_charge");
    }
    for (auto & value : version_three["groups"].as_array())
    {
        auto & group_record = value.as_object();
        if (group_record.at("summary").is_null()) continue;
        auto & summary_record = group_record.at("summary").as_object();
        if (!summary_record.at("charge_median").is_null())
        {
            const auto median = summary_record.at("charge_median");
            summary_record.at("mean").as_array()[2] = median;
            auto & inference = summary_record.at("inference").as_object();
            inference.at("mean").as_array()[2] = median;
            inference.at("mdpde").as_array()[2] = median;
            inference.at("prior").as_object().at("point").as_array()[2] = median;
        }
        summary_record.erase("charge_median");
    }
    ModelObject migrated_v3{model};
    rhbm_gem::stage_result_io::Decode(
        migrated_v3, boost::json::serialize(version_three), 4, {});
    const auto migrated_summary = migrated_v3.GetAnalysisView().GetGroupParameterSummary(keys[0]);
    ASSERT_TRUE(migrated_summary);
    ASSERT_TRUE(migrated_summary->inference);
    EXPECT_EQ(migrated_summary->member_ids, summary.member_ids);
    EXPECT_DOUBLE_EQ(migrated_summary->descriptive_charge_median.value(),
        summary.descriptive_charge_median.value());
    EXPECT_DOUBLE_EQ(migrated_summary->inference->mean.GetOffset(), 0.0);
    EXPECT_DOUBLE_EQ(migrated_summary->inference->mdpde.GetOffset(), 0.0);
    for (const int id : summary.member_ids)
    {
        const auto original_member = AtomLocalPotentialView::For(*model.FindAtomPtr(id)).GetGroupMemberResult();
        const auto migrated_member = AtomLocalPotentialView::For(*migrated_v3.FindAtomPtr(id)).GetGroupMemberResult();
        ASSERT_TRUE(original_member);
        ASSERT_TRUE(migrated_member);
        ASSERT_TRUE(original_member->fixed_charge);
        ASSERT_TRUE(migrated_member->fixed_charge);
        EXPECT_DOUBLE_EQ(migrated_member->fixed_charge->coefficient,
            original_member->fixed_charge->coefficient);
        EXPECT_DOUBLE_EQ(migrated_member->fixed_charge->width,
            original_member->fixed_charge->width);
        EXPECT_DOUBLE_EQ(migrated_member->posterior.GetModel().GetAmplitude(),
            original_member->posterior.GetModel().GetAmplitude());
        EXPECT_DOUBLE_EQ(migrated_member->posterior.GetModel().GetWidth(),
            original_member->posterior.GetModel().GetWidth());
        for (const double distance : {0.0, 1.2, 2.5, 2.6})
            EXPECT_DOUBLE_EQ(core::EvaluateGroupMemberResponseAtDistance(*migrated_member, distance),
                core::EvaluateGroupMemberResponseAtDistance(*original_member, distance));
    }

    auto invalid_v4 = boost::json::parse(rhbm_gem::stage_result_io::Encode(model)).as_object();
    for (auto & value : invalid_v4.at("atoms").as_array())
    {
        auto & atom_record = value.as_object();
        if (atom_record.at("id").as_int64() != 1) continue;
        auto & precision = atom_record.at("post_joint_local").as_object()
            .at("information").as_object().at("precision").as_array();
        precision[0] = -1.0; precision[1] = 0.0; precision[2] = 0.0; precision[3] = 1.0;
    }
    ModelObject invalid_persistence{model};
    EXPECT_THROW(rhbm_gem::stage_result_io::Decode(
        invalid_persistence, boost::json::serialize(invalid_v4), 4, {}), std::invalid_argument);

    editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(1), {{-12345, {0.3, {0,0,0}, true}}});
    EXPECT_FALSE(first.GetPostFitPeeling());
    EXPECT_FALSE(first.GetPostJointLocalResult());
    EXPECT_FALSE(first.GetGroupMemberResult());
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(model.GetAnalysisView().GetGroupParameterSummary(keys[0])->eligible_count, 3);
    EXPECT_EQ(first.GetFinalModel(FittingStage::Second).ToVector(), joint_point);
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

TEST(JointComponentPartialSelectionTest, TargetOnlyPostprocessingRetainsFullContributorCoupling)
{
    auto f = joint_partial_test::Make("partial");
    SetFixedHaloExactMap(f);
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    const auto initialized=problem.WithFixedHalos({FixedHaloFor(f,1)});
    const auto fit=core::FitJointComponents(initialized,{f.b[0]});
    ASSERT_TRUE(fit.assembled_state) << fit.initialization.reason << " search=" << fit.search_completed
        << " components=" << fit.components.size()
        << (fit.components.empty() ? std::string{} : ", stop="+fit.components.front().stop_reason);
    const auto snapshot = core::CaptureJointAnalysisResult(fit);
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, snapshot, "output-scope");
    const auto target_stage=rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(1))
        .GetStageEstimate(rhbm_gem::FittingStage::Second);
    ASSERT_TRUE(target_stage.point) << target_stage.reason;
    for (const auto & atom : f.model->GetAtomList())
        f.model->EditAnalysis().SetAtomLocalRawSamplingEntries(*atom,
            core::SampleAtomMapValues(*f.map, *atom, SphereSamplingMethod::FibonacciDeterministic));
    const std::vector<std::size_t> all{0, 1};
    const auto target = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem,std::nullopt,&snapshot);
    const auto complete = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem,all,&snapshot);
    ASSERT_EQ(target.size(), 1); ASSERT_EQ(complete.size(), 1);
    EXPECT_EQ(target.at(1).neighbor_count, complete.at(1).neighbor_count);
    for (std::size_t i=0; i<target.at(1).samples.size(); ++i)
    {
        EXPECT_EQ(target.at(1).samples[i].response, complete.at(1).samples[i].response);
        EXPECT_EQ(target.at(1).samples[i].reason, complete.at(1).samples[i].reason);
    }
    f.model->EditAnalysis().SetJointResult(snapshot);
    auto metadata = snapshot.metadata; metadata.model_path = "metadata-only";
    f.model->EditAnalysis().UpdateJointMetadata(metadata);
    EXPECT_EQ(f.model->GetAnalysisView().GetJointResult()->objective, snapshot.objective);
    EXPECT_EQ(f.model->GetAnalysisView().GetJointResult()->metadata.model_path, "metadata-only");
    EXPECT_TRUE(rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(1)).HasFinalModel(rhbm_gem::FittingStage::Second));
}

TEST(JointComponentPartialSelectionTest, FixedChargePosteriorEvaluatorUsesJointWidthAndKernelSupport)
{
    const std::array<rhbm_gem::JointFixedChargeProfile, 4> profiles{{
        {0.0, 0.6}, {0.4, 0.6}, {-0.3, 0.42}, {0.4, 0.9}}};
    const std::array<double, 4> posterior_widths{0.6, 0.48, 0.51, 0.48};
    for (std::size_t i = 0; i < profiles.size(); ++i)
    {
        rhbm_gem::GroupGaussianMemberResult member;
        member.posterior = {rhbm_gem::GaussianModel3D{2.3, posterior_widths[i], 0.0}, {}};
        member.fixed_charge = profiles[i];
        member.charge_inferred = false;
        for (const double distance : {0.0, 1.2, 2.5, 2.6})
        {
            const double joint_charge = n::EvaluateKernel(
                distance * distance, profiles[i].width, 2.5).charge;
            const double gaussian_only = member.posterior.GetModel().SignalAtDistance(distance);
            const double reconstructed = core::EvaluateGroupMemberResponseAtDistance(member, distance);
            const double expected = member.posterior.GetModel().SignalAtDistance(distance) +
                profiles[i].coefficient * joint_charge;
            EXPECT_DOUBLE_EQ(reconstructed, expected);
            if (joint_charge == 0.0 || profiles[i].coefficient == 0.0)
                EXPECT_DOUBLE_EQ(reconstructed, gaussian_only);
        }
        const double gaussian_only = member.posterior.GetModel().SignalAtDistance(0.8);
        const double reconstructed = core::EvaluateGroupMemberResponseAtDistance(member, 0.8);
        if (profiles[i].coefficient == 0.0)
            EXPECT_DOUBLE_EQ(reconstructed, gaussian_only);
        else
            EXPECT_NE(reconstructed, gaussian_only);
    }
    EXPECT_NE(profiles[2].width, posterior_widths[2]);

    rhbm_gem::GroupGaussianMemberResult narrow, broad;
    narrow.posterior = broad.posterior = {rhbm_gem::GaussianModel3D{2.3, 0.48, 0.0}, {}};
    narrow.fixed_charge = rhbm_gem::JointFixedChargeProfile{0.4, 0.42};
    broad.fixed_charge = rhbm_gem::JointFixedChargeProfile{0.4, 0.9};
    EXPECT_NE(core::EvaluateGroupMemberResponseAtDistance(narrow, 1.2),
        core::EvaluateGroupMemberResponseAtDistance(broad, 1.2));
    EXPECT_NE(core::EvaluateGroupMemberResponseAtDistance(narrow, 0.8),
        narrow.posterior.GetModel().SignalAtDistance(0.8));
    EXPECT_NE(core::EvaluateGroupMemberResponseAtDistance(broad, 0.8),
        broad.posterior.GetModel().SignalAtDistance(0.8));
}

TEST(JointComponentPartialSelectionTest, ExistingSampleFirstExecutorUsesOnlyTheExplicitWorkset)
{
    auto f = joint_partial_test::Make("partial");
    auto editor = f.model->EditAnalysis();
    editor.InitializeFromSelection();
    rhbm_gem::LocalGaussianResult seed;
    seed.mdpde = {rhbm_gem::GaussianModel3D{7, .73, 0}, {}};
    editor.SetAtomLocalGaussianResult(rhbm_gem::FittingStage::First, *f.model->FindAtomPtr(1), seed);
    seed.mdpde = {rhbm_gem::GaussianModel3D{0, 1, 0}, {}};
    const auto & halo = *f.model->FindAtomPtr(2);
    editor.SetAtomLocalGaussianResult(rhbm_gem::FittingStage::First, halo, seed);
    editor.SetAtomLocalRawSamplingEntries(halo, core::SampleAtomMapValues(*f.map, halo, SphereSamplingMethod::FibonacciDeterministic));
    const auto selected = f.model->GetSelectedAtoms();
    std::map<std::pair<int, std::string>, int> calls;
    core::detail::FirstStageObserverForTesting() = [&](int id, std::string_view phase) { ++calls[{id, std::string(phase)}]; };
    const std::vector<rhbm_gem::AtomObject *> atoms{f.model->FindAtomPtr(2)};
    core::FitOptions options; options.quiet_mode = true; options.thread_size = 1;
    core::detail::RunBatchFirstStageFromPreparedSamples(*f.model, atoms, options);
    core::detail::FirstStageObserverForTesting() = {};
    EXPECT_EQ(f.model->GetSelectedAtoms(), selected);
    EXPECT_DOUBLE_EQ(rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(1)).GetFinalModel(rhbm_gem::FittingStage::First).GetWidth(), .73);
    EXPECT_EQ((calls[{2, "first"}]), 1);
    EXPECT_EQ((calls[{1, "first"}]), 0);
}
