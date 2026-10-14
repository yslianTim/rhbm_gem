#include <gtest/gtest.h>
#include <boost/json.hpp>
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <algorithm>
#include "support/JointPartialSelection.hpp"
#include "support/JointProfileSearchDiagnostic.hpp"
#include "support/DataObjectTestSupport.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/HaloEligibility.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include <map>
#include <utility>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include "support/CommandTestHelpers.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "core/detail/StageSummary.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/LocalConditionalCovariance.hpp"
#include "core/detail/MapInterpolation.hpp"
#include "core/detail/MapSampler.hpp"
#include "data/io/detail/JointResultJson.hpp"
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
    std::vector<rhbm_gem::ExcludedHaloRecord> excluded;
    std::size_t candidate_halos{};
    for(std::size_t a=0;a<support.size();++a)
    {
        const auto & atom=model.GetAtomList()[a];
        if(atom->GetElement()==Element::HYDROGEN) continue;
        if(!selected.contains(atom.get()))
        {
            ++candidate_halos;
            auto candidate_support=support[a];
            const auto eligibility=n::EvaluateFullABCHaloEligibility(candidate_support);
            if(!eligibility.eligible)
            {
                excluded.push_back({std::to_string(atom->GetSerialID()),eligibility.support_rows,
                    eligibility.distinct_squared_distances,eligibility.reason});
                continue;
            }
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
    ASSERT_TRUE(actual.selection_domain->halo_eligibility);
    const auto & provenance=*actual.selection_domain->halo_eligibility;
    std::sort(excluded.begin(),excluded.end(),[](const auto & a,const auto & b){return a.atom_id<b.atom_id;});
    EXPECT_EQ(actual.selection_domain->contributor_policy,rhbm_gem::FullABCHaloContributorPolicy);
    EXPECT_EQ(provenance.minimum_observations,3u);
    EXPECT_EQ(provenance.minimum_distinct_radial_information,3u);
    EXPECT_EQ(provenance.candidate_halos,candidate_halos);
    EXPECT_EQ(provenance.retained_halos,ids.size()-targets.size());
    ASSERT_EQ(provenance.excluded_halos.size(),excluded.size());
    for(std::size_t i=0;i<excluded.size();++i)
    {
        EXPECT_EQ(provenance.excluded_halos[i].atom_id,excluded[i].atom_id);
        EXPECT_EQ(provenance.excluded_halos[i].support_rows,excluded[i].support_rows);
        EXPECT_EQ(provenance.excluded_halos[i].distinct_radial_distances,excluded[i].distinct_radial_distances);
        EXPECT_EQ(provenance.excluded_halos[i].reason,excluded[i].reason);
    }
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
    EXPECT_EQ(p.Input().atom_ids,(std::vector<std::string>{"1","4"}));
    ASSERT_TRUE(p.Input().selection_domain->halo_eligibility);
    ASSERT_EQ(p.Input().selection_domain->halo_eligibility->excluded_halos.size(),3u);
    const auto fit=core::FitJointComponents(p,{.5,.5});
    EXPECT_FALSE(fit.assembled_state);
    const auto & data=core::JointProblemAccess::Get(p);
    EXPECT_EQ(data.partition.unobserved_atoms.size(),1);
}
TEST(JointComponentPartialSelectionTest, PrunesIneligibleHalosBeforeProblemConstructionAndKeepsFixedRows)
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
    EXPECT_EQ(input.atom_ids,(std::vector<std::string>{"10","30"}));
    ASSERT_TRUE(input.selection_domain);
    EXPECT_EQ(input.selection_domain->target_indices,(std::vector<std::size_t>{0}));
    ASSERT_TRUE(input.selection_domain->halo_eligibility);
    const auto & provenance=*input.selection_domain->halo_eligibility;
    EXPECT_EQ(input.selection_domain->contributor_policy,rhbm_gem::FullABCHaloContributorPolicy);
    EXPECT_EQ(provenance.candidate_halos,3u);
    EXPECT_EQ(provenance.retained_halos,1u);
    ASSERT_EQ(provenance.excluded_halos.size(),2u);
    EXPECT_EQ(provenance.excluded_halos[0].atom_id,"20");
    EXPECT_EQ(provenance.excluded_halos[0].support_rows,2u);
    EXPECT_EQ(provenance.excluded_halos[0].reason,"insufficient-observations");
    EXPECT_EQ(provenance.excluded_halos[1].atom_id,"40");
    EXPECT_EQ(provenance.excluded_halos[1].support_rows,0u);
    EXPECT_EQ(provenance.excluded_halos[1].reason,"insufficient-observations");
    ASSERT_EQ(input.support.size(),2u);
    auto retained= input.support[1];
    const auto eligibility=n::EvaluateFullABCHaloEligibility(retained);
    EXPECT_TRUE(eligibility.eligible);
    EXPECT_EQ(eligibility.support_rows,4u);
    EXPECT_EQ(eligibility.distinct_squared_distances,4u);

    const auto workset=core::detail::MakeJointFittingWorkset(model,problem);
    EXPECT_EQ(workset.contributors.size(),input.atom_ids.size());
    EXPECT_EQ(core::detail::CollectFirstStageAtoms(workset).size(),input.atom_ids.size());
    const auto & data=core::JointProblemAccess::Get(problem);
    EXPECT_EQ(data.construction_memory.shape.component_atoms,2u);
    EXPECT_EQ(data.construction_memory.shape.component_support_entries,15u);

    auto invalid=input; invalid.support[1].resize(1);
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

TEST(JointComponentPartialSelectionTest, PrunedBridgeRebuildsStructuralComponents)
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
    const auto pruned=core::BuildJointProblem(map,model);
    EXPECT_EQ(pruned.Input().atom_ids,(std::vector<std::string>{"1","2"}));
    const auto & after=core::JointProblemAccess::Get(pruned);
    ASSERT_EQ(after.partition.components.size(),2u);

    core::JointProblemInput unpruned=pruned.Input();
    unpruned.selection_domain->contributor_policy="all-non-hydrogen";
    unpruned.selection_domain->halo_eligibility.reset();
    unpruned.atom_ids.push_back("3");
    unpruned.support.emplace_back();
    for(std::size_t row=0;row<unpruned.row_ids.size();++row)
    {
        const double x=-2.5+2.5*static_cast<double>(row);
        const double squared=(x-3.75)*(x-3.75);
        if(squared<=6.25) unpruned.support.back().push_back({row,squared});
    }
    unpruned.selection_domain->target_indices={0,1};
    const core::JointProblem unpruned_problem(std::move(unpruned));
    const auto & before=core::JointProblemAccess::Get(unpruned_problem);
    EXPECT_EQ(before.partition.components.size(),1u);
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
    const auto p=core::BuildJointProblem(*f.map,*f.model); const auto & data=core::JointProblemAccess::Get(p);
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
    const auto changed=core::BuildJointProblem(*f.map,*f.model);
    EXPECT_EQ(changed.Input().observations,p.Input().observations);
    const auto repeated=core::FitJointComponents(changed,f.b);
    EXPECT_EQ(repeated.assembled_state->ac,fit.assembled_state->ac); EXPECT_EQ(repeated.objective,fit.objective);
}
TEST(JointComponentPartialSelectionTest, InvalidWidthsAreNotSuccessForRetainedContributors)
{
    auto f=joint_partial_test::Make("weak"); const auto p=core::BuildJointProblem(*f.map,*f.model);
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
    EXPECT_EQ(boost::json::parse(encoded).at("schema_version").as_int64(),11);
    const auto decoded=rhbm_gem::joint_result_io::Decode(encoded);
    ASSERT_TRUE(decoded.selection_domain); EXPECT_EQ(decoded.selection_domain->target_indices,(std::vector<std::size_t>{0}));
    ASSERT_TRUE(decoded.selection_domain->halo_eligibility);
    EXPECT_EQ(decoded.selection_domain->halo_eligibility->candidate_halos,
        saved.selection_domain->halo_eligibility->candidate_halos);
    EXPECT_EQ(decoded.selection_domain->halo_eligibility->excluded_halos.size(),
        saved.selection_domain->halo_eligibility->excluded_halos.size());
    EXPECT_EQ(decoded.initialization.data_scope,"contributor-local-sampling-may-read-outside-target-domain");
    EXPECT_EQ(rhbm_gem::joint_result_io::Encode(decoded),encoded);
    auto missing_provenance=boost::json::parse(encoded).as_object();
    missing_provenance.at("selection_domain").as_object()["halo_eligibility"]=nullptr;
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(missing_provenance)),std::invalid_argument);
    missing_provenance=boost::json::parse(encoded).as_object();
    missing_provenance["schema_version"]=10;
    EXPECT_THROW(rhbm_gem::joint_result_io::Decode(boost::json::serialize(missing_provenance)),std::invalid_argument);
    f.model->EditAnalysis().SetJointResult(saved);
}

TEST(JointComponentPartialSelectionTest, ExcludedHaloProvenancePersistsThroughJsonAndSqlite)
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
    const auto fit=core::FitJointComponents(problem,{.5,.5});
    const auto saved=core::CaptureJointAnalysisResult(fit);
    const auto decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(saved));
    ASSERT_TRUE(decoded.selection_domain->halo_eligibility);
    ASSERT_EQ(decoded.selection_domain->halo_eligibility->excluded_halos.size(),2u);
    EXPECT_EQ(decoded.selection_domain->halo_eligibility->excluded_halos[0].atom_id,"20");
    model.EditAnalysis().ApplyJointResult(decoded,"provenance");

    const command_test::ScopedTempDir directory{"joint_halo_provenance_roundtrip"};
    rhbm_gem::DataRepository repository(directory.path()/"provenance.sqlite");
    repository.SaveModel(model,"joint");
    const auto loaded=repository.LoadModel("joint");
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(loaded->GetAnalysisView().GetJointResult()->selection_domain->halo_eligibility->excluded_halos,
        decoded.selection_domain->halo_eligibility->excluded_halos);

    auto radial=decoded;
    radial.selection_domain->halo_eligibility->excluded_halos[0]={"20",3,2,"insufficient-radial-information"};
    const auto radial_decoded=rhbm_gem::joint_result_io::Decode(rhbm_gem::joint_result_io::Encode(radial));
    EXPECT_EQ(radial_decoded.selection_domain->halo_eligibility->excluded_halos[0],
        radial.selection_domain->halo_eligibility->excluded_halos[0]);
    model.EditAnalysis().ApplyJointResult(radial,"radial-provenance");
    repository.SaveModel(model,"radial-provenance");
    const auto radial_loaded=repository.LoadModel("radial-provenance");
    ASSERT_TRUE(radial_loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(radial_loaded->GetAnalysisView().GetJointResult()->selection_domain->halo_eligibility->excluded_halos[0],
        radial.selection_domain->halo_eligibility->excluded_halos[0]);
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
    std::map<int, int> sampling_calls;
    core::detail::FirstStageObserverForTesting() = [&](int id, std::string_view phase) {
        ++calls[{id, std::string(phase)}];
    };
    core::detail::PotentialSamplingObserverForTesting() = [&](int id) { ++sampling_calls[id]; };
    core::FitOptions options;
    options.estimator = core::PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode = true;
    core::RunPotentialFittingWorkflow(*f.map, *f.model, options);
    core::detail::FirstStageObserverForTesting() = {};
    core::detail::PotentialSamplingObserverForTesting() = {};
    ASSERT_TRUE(f.model->GetAnalysisView().GetJointResult());
    const auto & result = *f.model->GetAnalysisView().GetJointResult();
    const auto expected_stages = rhbm_gem::data_internal::BuildJointStageEstimates(result, "expected-run");
    EXPECT_EQ(f.model->GetSelectedAtoms(), selected);
    EXPECT_EQ(f.model->GetSelectedBonds(), bonds);
    EXPECT_EQ(result.row_ids, problem.Input().row_ids);
    ASSERT_EQ(result.initialization.atoms.size(), problem.Input().atom_ids.size());
    for (std::size_t i = 0; i < result.atom_ids.size(); ++i)
    {
        const auto id = std::stoi(result.atom_ids[i]);
        EXPECT_EQ((calls[{id, "first"}]), 1);
        EXPECT_EQ(sampling_calls[id], 1);
        const auto view = rhbm_gem::AtomLocalPotentialView::For(*f.model->FindAtomPtr(id));
        EXPECT_DOUBLE_EQ(view.GetFinalModel(rhbm_gem::FittingStage::First).GetWidth(), result.initialization.b[i]);
        EXPECT_FALSE(view.GetRawSamplingEntries(false).empty());
        const auto & stage = view.GetStageEstimate(rhbm_gem::FittingStage::Second);
        const auto & expected = expected_stages.at(id);
        EXPECT_EQ(stage.point.has_value(), expected.point.has_value());
        if (stage.point && expected.point) EXPECT_EQ(stage.point->ToVector(), expected.point->ToVector());
        EXPECT_EQ(stage.uncertainty.status, rhbm_gem::EvidenceStatus::NotRun);
        EXPECT_FALSE(stage.uncertainty.covariance);
        const bool target = stage.source.role == rhbm_gem::FittingRole::Target;
        EXPECT_EQ(view.GetPostFitPeeling().has_value(), target);
        EXPECT_EQ(view.GetGroupEvidence().has_value(), target);
        if (target)
        {
            EXPECT_EQ(view.GetPostFitPeeling()->samples.size(), view.GetRawSamplingEntries(false).size());
            const auto & evidence = *view.GetGroupEvidence();
            EXPECT_EQ(evidence.uncertainty_method, "local-conditional-iid-ls");
            EXPECT_EQ(evidence.atom_id, stage.source.atom_id);
            EXPECT_EQ(evidence.component_id, stage.source.component_id);
            EXPECT_EQ(evidence.source_id, stage.source.run_id);
            if (stage.point && evidence.estimate)
                EXPECT_TRUE(evidence.estimate->isApprox(
                    rhbm_gem::linearization_service::EncodeGaussianToParameterVector(*stage.point), 1.0e-12));
        }
        EXPECT_FALSE(view.GetGroupMemberResult());
    }
    const auto direct = n::FitFixedNeighborComponents(problem, result.initialization.b, {});
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
    EXPECT_NE(summary.find("\n| Spot | valid | not-converged | unavailable | A mean / s.d. | B mean / s.d. | C mean / s.d. |"), std::string::npos);
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
    expectSummaryRowWidths(summary, false);
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
    const auto raw = core::detail::BuildLocalPotentialSampleList(*f.map, points);
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

TEST(JointComponentPartialSelectionTest, GroupPosteriorUsesLocalEvidenceAndPreservesSecond)
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
        stage.convergence = JointCheckStatus::Passed;
        stages[i] = std::move(stage);
    }
    editor.ApplySecondStageEstimates(stages);
    std::vector<std::pair<int, PostFitPeelingResult>> peeling_updates;
    std::vector<std::pair<int, GroupParameterEvidence>> evidence_updates;
    peeling_updates.reserve(stages.size());
    evidence_updates.reserve(stages.size());
    for (const auto & [id, unused] : stages)
    {
        (void)unused;
        const auto & stage = AtomLocalPotentialView::For(*model.FindAtomPtr(id)).GetStageEstimate(FittingStage::Second);
        LocalPotentialSampleList raw;
        PostFitPeelingResult peeling;
        peeling.source = stage.source;
        for (int sample_index = 0; sample_index < 6; ++sample_index)
        {
            const double distance = 0.12 + 0.43 * sample_index;
            const auto basis = core::joint_component::EvaluateKernel(distance * distance, stage.point->GetWidth(), 2.5);
            LocalPotentialSample sample;
            sample.point.distance = distance;
            sample.response = stage.point->GetAmplitude() * basis.gaussian +
                stage.point->GetOffset() * basis.charge + 0.01 * std::sin(sample_index + id);
            raw.push_back(sample);
            peeling.samples.push_back({sample.response, {}});
        }
        if (id == 3)
        {
            for (std::size_t i = 2; i < peeling.samples.size(); ++i)
            {
                peeling.samples[i].response.reset();
                peeling.samples[i].reason = "test-missing-coverage";
            }
        }
        editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(id), raw);
        const auto local = AtomLocalPotentialView::For(*model.FindAtomPtr(id));
        const auto covariance = core::detail::ComputeLocalConditionalCovariance(
            stage, raw, peeling, local.HasSampleGeometry());
        evidence_updates.emplace_back(id, core::detail::BuildLocalGroupEvidence(stage, &peeling, covariance));
        peeling_updates.emplace_back(id, std::move(peeling));
    }
    ResetGroupInvalidationStatsForTesting();
    editor.ApplyJointPostFitEvidenceBatch(std::move(peeling_updates), std::move(evidence_updates));
    const auto publication_stats = GetGroupInvalidationStatsForTesting();
    EXPECT_EQ(publication_stats.membership_traversals, 1u);
    EXPECT_EQ(publication_stats.groups_invalidated, 1u);
    core::FitOptions options; options.estimator = core::PotentialEstimator::JOINT_COMPONENTS; options.quiet_mode = true;
    core::RunGroupPotentialFitting(model, options);
    const auto keys = model.GetAnalysisView().CollectAtomGroupKeys();
    ASSERT_EQ(keys.size(), 1);
    auto summary = *model.GetAnalysisView().GetGroupParameterSummary(keys[0]);
    ASSERT_EQ(summary.status, EvidenceStatus::Available) << summary.reason;
    EXPECT_EQ(summary.point_count, 5); EXPECT_EQ(summary.eligible_count, 4); EXPECT_EQ(summary.excluded_count, 1);
    EXPECT_EQ(summary.member_ids, (std::vector<int>{1,2,4,5}));
    EXPECT_FALSE(AtomLocalPotentialView::For(*model.FindAtomPtr(3)).GetGroupMemberResult());

    const auto first = AtomLocalPotentialView::For(*model.FindAtomPtr(1));
    const auto original_peeling = *first.GetPostFitPeeling();
    const auto original_evidence = *first.GetGroupEvidence();
    auto assert_unchanged = [&] {
        const auto current = AtomLocalPotentialView::For(*model.FindAtomPtr(1));
        ASSERT_TRUE(current.GetPostFitPeeling());
        ASSERT_TRUE(current.GetGroupEvidence());
        EXPECT_EQ(current.GetPostFitPeeling()->source.run_id, original_peeling.source.run_id);
        ASSERT_EQ(current.GetPostFitPeeling()->samples.size(), original_peeling.samples.size());
        for (std::size_t i = 0; i < original_peeling.samples.size(); ++i)
        {
            EXPECT_EQ(current.GetPostFitPeeling()->samples[i].response, original_peeling.samples[i].response);
            EXPECT_EQ(current.GetPostFitPeeling()->samples[i].reason, original_peeling.samples[i].reason);
        }
        EXPECT_EQ(current.GetGroupEvidence()->source_id, original_evidence.source_id);
        EXPECT_EQ(current.GetGroupEvidence()->status, original_evidence.status);
        EXPECT_EQ(current.GetGroupEvidence()->covariance, original_evidence.covariance);
        EXPECT_TRUE(model.GetAnalysisView().GetGroupParameterSummary(keys[0]));
    };
    ResetGroupInvalidationStatsForTesting();
    {
        auto changed = original_peeling;
        *changed.samples.front().response += 0.5;
        auto invalid = changed;
        invalid.source.atom_id = "999";
        auto invalid_evidence = original_evidence;
        invalid_evidence.atom_id = "999";
        std::vector<std::pair<int, PostFitPeelingResult>> peelings{{1, std::move(changed)}, {999, std::move(invalid)}};
        std::vector<std::pair<int, GroupParameterEvidence>> evidence{{1, original_evidence}, {999, std::move(invalid_evidence)}};
        EXPECT_THROW(editor.ApplyJointPostFitEvidenceBatch(std::move(peelings), std::move(evidence)), std::invalid_argument);
        assert_unchanged();
    }
    {
        std::vector<std::pair<int, PostFitPeelingResult>> peelings{{1, original_peeling}, {1, original_peeling}};
        std::vector<std::pair<int, GroupParameterEvidence>> evidence{{1, original_evidence}, {1, original_evidence}};
        EXPECT_THROW(editor.ApplyJointPostFitEvidenceBatch(std::move(peelings), std::move(evidence)), std::invalid_argument);
        assert_unchanged();
    }
    {
        auto stale = original_peeling;
        stale.source.run_id = "stale-run";
        EXPECT_THROW(editor.ApplyJointPostFitEvidenceBatch({{1, std::move(stale)}}, {{1, original_evidence}}),
            std::invalid_argument);
        assert_unchanged();
    }
    {
        auto mismatch = original_peeling;
        mismatch.samples.pop_back();
        EXPECT_THROW(editor.ApplyJointPostFitEvidenceBatch({{1, std::move(mismatch)}}, {{1, original_evidence}}),
            std::invalid_argument);
        assert_unchanged();
    }
    const auto failed_publication_stats = GetGroupInvalidationStatsForTesting();
    EXPECT_EQ(failed_publication_stats.membership_traversals, 0u);
    EXPECT_EQ(failed_publication_stats.groups_invalidated, 0u);

    const auto before = first.GetGroupMemberResult()->posterior.GetModel().ToVector();
    const auto second = first.GetFinalModel(FittingStage::Second).ToVector();
    const auto posterior_before_republish = *first.GetGroupMemberResult();
    std::vector<std::pair<int, PostFitPeelingResult>> refreshed_peelings;
    std::vector<std::pair<int, GroupParameterEvidence>> refreshed_evidence;
    for (const auto & [id, unused] : stages)
    {
        (void)unused;
        const auto view = AtomLocalPotentialView::For(*model.FindAtomPtr(id));
        refreshed_peelings.emplace_back(id, *view.GetPostFitPeeling());
        refreshed_evidence.emplace_back(id, *view.GetGroupEvidence());
    }
    ResetGroupInvalidationStatsForTesting();
    editor.ApplyJointPostFitEvidenceBatch(std::move(refreshed_peelings), std::move(refreshed_evidence));
    const auto summary_invalidation_stats = GetGroupInvalidationStatsForTesting();
    EXPECT_EQ(summary_invalidation_stats.membership_traversals, 1u);
    EXPECT_EQ(summary_invalidation_stats.groups_invalidated, 1u);
    EXPECT_FALSE(model.GetAnalysisView().GetGroupParameterSummary(keys[0]));
    for (const auto & atom : model.GetAtomList())
        EXPECT_FALSE(AtomLocalPotentialView::For(*atom).GetGroupMemberResult());
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(first.GetGroupMemberResult()->posterior.GetModel().ToVector(), before);
    EXPECT_EQ(*first.GetGroupMemberResult()->parameter_covariance, *posterior_before_republish.parameter_covariance);
    summary = *model.GetAnalysisView().GetGroupParameterSummary(keys[0]);

    auto legacy_evidence = *AtomLocalPotentialView::For(*model.FindAtomPtr(5)).GetGroupEvidence();
    legacy_evidence.uncertainty_method = "iid-ls-linearized";
    editor.SetAtomGroupEvidence(*model.FindAtomPtr(5), legacy_evidence);
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(model.GetAnalysisView().GetGroupParameterSummary(keys[0])->eligible_count, 3);
    legacy_evidence.uncertainty_method = "local-conditional-iid-ls";
    editor.SetAtomGroupEvidence(*model.FindAtomPtr(5), legacy_evidence);
    core::RunGroupPotentialFitting(model, options);

    {
        const command_test::ScopedTempDir directory{"joint_posterior_roundtrip"};
        DataRepository repository(directory.path()/"posterior.sqlite"); repository.SaveModel(model,"model");
        const auto loaded=repository.LoadModel("model");
        const auto saved=AtomLocalPotentialView::For(*loaded->FindAtomPtr(1));
        ASSERT_TRUE(saved.GetGroupMemberResult());
        ASSERT_TRUE(saved.GetGroupEvidence());
        EXPECT_EQ(saved.GetGroupEvidence()->uncertainty_method, "local-conditional-iid-ls");
        EXPECT_EQ(*saved.GetGroupEvidence()->covariance, *first.GetGroupEvidence()->covariance);
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
        stage.point = stage.point->WithAmplitude(stage.point->GetAmplitude() + 0.1);
        changed_editor.SetAtomStageEstimate(FittingStage::Second, *changed_model.FindAtomPtr(1), stage);
        EXPECT_FALSE(changed_view.GetGroupEvidence());
        EXPECT_FALSE(changed_view.GetGroupMemberResult());
        EXPECT_FALSE(changed_model.GetAnalysisView().GetGroupParameterSummary(keys[0]));
        EXPECT_FALSE(AtomLocalPotentialView::For(*changed_model.FindAtomPtr(2)).GetGroupMemberResult());
    }
    auto changed_evidence = *first.GetGroupEvidence();
    changed_evidence.covariance = *changed_evidence.covariance * 10.0;
    editor.SetAtomGroupEvidence(*model.FindAtomPtr(1), changed_evidence);
    core::RunGroupPotentialFitting(model, options);
    EXPECT_FALSE(first.GetGroupMemberResult()->posterior.GetModel().ToVector().isApprox(before, 1.0e-9));
    EXPECT_EQ(first.GetFinalModel(FittingStage::Second).ToVector(), second);

    auto changed_peeling = *first.GetPostFitPeeling();
    *changed_peeling.samples[0].response += 0.02;
    editor.SetAtomPostFitPeeling(*model.FindAtomPtr(1), changed_peeling);
    EXPECT_FALSE(first.GetGroupEvidence());
    EXPECT_FALSE(first.GetGroupMemberResult());
    const auto raw = first.GetRawSamplingEntries(false);
    const auto & stage = first.GetStageEstimate(FittingStage::Second);
    const auto covariance = core::detail::ComputeLocalConditionalCovariance(stage, raw, changed_peeling, true);
    editor.SetAtomGroupEvidence(*model.FindAtomPtr(1),
        core::detail::BuildLocalGroupEvidence(stage, &changed_peeling, covariance));
    core::RunGroupPotentialFitting(model, options);

    editor.SetAtomLocalRawSamplingEntries(*model.FindAtomPtr(1), {{-12345, {0.3, {0,0,0}, true}}});
    EXPECT_FALSE(first.GetPostFitPeeling());
    EXPECT_FALSE(first.GetGroupEvidence());
    EXPECT_FALSE(first.GetGroupMemberResult());
    core::RunGroupPotentialFitting(model, options);
    EXPECT_EQ(model.GetAnalysisView().GetGroupParameterSummary(keys[0])->eligible_count, 3);
    EXPECT_EQ(first.GetFinalModel(FittingStage::Second).ToVector(), second);
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
    const auto problem = core::BuildJointProblem(*f.map, *f.model);
    const auto snapshot = core::CaptureJointAnalysisResult(core::FitJointComponents(problem, f.b));
    rhbm_gem::data_internal::ApplyJointStageEstimates(*f.model, snapshot, "output-scope");
    for (const auto & atom : f.model->GetAtomList())
        f.model->EditAnalysis().SetAtomLocalRawSamplingEntries(*atom,
            core::SampleAtomMapValues(*f.map, *atom, SphereSamplingMethod::FibonacciDeterministic));
    const std::vector<std::size_t> all{0, 1};
    const auto target = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem);
    const auto complete = core::detail::BuildPostFitPeelingSamples(*f.map, *f.model, problem, all);
    ASSERT_EQ(target.size(), 1); ASSERT_EQ(complete.size(), 2);
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
