#include <gtest/gtest.h>
#include <boost/json.hpp>

#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/MapSampler.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include <cmath>
#include <memory>
#include <vector>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace core=rhbm_gem::core;

JointProblemInput PureInput()
{
    JointProblemInput input;
    input.atom_ids={"1"};
    input.row_ids={"0","1","2","3","4","5","6","7","8"};
    input.support.resize(1);
    for (std::size_t row=0;row<input.row_ids.size();++row)
    {
        const auto basis=core::joint_component::EvaluateKernel(.04*static_cast<double>(row*row),.65,2.5);
        input.support[0].push_back({row,.04*static_cast<double>(row*row)});
        input.observations.push_back(2.0*basis.gaussian-.2*basis.charge);
    }
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->target_indices={0};
    input.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    return input;
}

JointInitialization TargetInitialization()
{
    JointInitialization initialization;
    initialization.valid=true;
    initialization.b={.65};
    JointInitializationAtom atom;
    atom.id="1";
    atom.target=true;
    atom.reason="valid-width";
    atom.used_b=.65;
    atom.original_b=.65;
    initialization.atoms.push_back(std::move(atom));
    initialization.reason="valid-widths";
    return initialization;
}

JointAnalysisResult PurePeelingSnapshot(const JointProblem & problem,const ModelObject & model)
{
    JointAnalysisResult result;
    result.parameterization_contract=std::string(JointPureHaloSupportExclusionParameterizationContract);
    result.atom_ids=problem.OriginalInput().atom_ids;
    result.row_ids=problem.OriginalInput().row_ids;
    result.fitting_row_ids=problem.Input().row_ids;
    result.fitting_to_original_contributors=problem.FittingToOriginalContributors();
    result.fitting_to_original_rows=problem.FittingToOriginalRows();
    result.selection_domain=problem.OriginalInput().selection_domain;
    result.available_row_mask.assign(result.row_ids.size(),false);
    for (const auto row:problem.FittingToOriginalRows()) result.available_row_mask[row]=true;
    result.runtime_convergence=JointCheckStatus::Passed;
    result.target_runtime_convergence=JointCheckStatus::Passed;
    JointAnalysisComponent component;
    component.id="pure-component";
    component.atoms={0};
    component.rows=problem.FittingToOriginalRows();
    component.runtime_convergence=JointCheckStatus::Passed;
    component.target_runtime_convergence=JointCheckStatus::Passed;
    component.state=JointState{{1.0,0.0},{1.0},{0.0},{0.0},0.0,{}};
    result.components.push_back(component);
    result.assembled_state=component.state;
    JointHaloMedianSnapshot medians;
    medians.joint_run_id="pure-peeling-run";
    medians.valid_target_count=1;
    medians.global_median={.5,.8,0.0};
    const auto * halo=model.FindAtomPtr(2);
    medians.halos.push_back({1,"2",data_internal::GetGroupKey(halo),.5,.8,0.0,
        HaloMedianSource::GlobalTargetMedian,1,"pure-peeling-run",
        std::string(HaloMedianReconstructionContractVersion)});
    result.halo_median_reconstruction=std::move(medians);
    return result;
}

struct MapFixture
{
    ModelObject model;
    MapObject map;
};

MapFixture BoundaryFixture()
{
    auto target=std::make_unique<AtomObject>();
    target->SetSerialID(1);
    target->SetElement(Element::CARBON);
    target->SetPosition(0,0,0);
    auto halo=std::make_unique<AtomObject>();
    halo->SetSerialID(2);
    halo->SetElement(Element::CARBON);
    halo->SetPosition(5,0,0);
    std::vector<std::unique_ptr<AtomObject>> atoms;
    atoms.push_back(std::move(target));
    atoms.push_back(std::move(halo));
    ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==1;});
    auto values=std::make_unique<double[]>(11);
    for (std::size_t row=0;row<11;++row) values[row]=static_cast<double>(row+1);
    MapObject map({11,1,1},{.5,.5,.5},{-2.5,0,0},std::move(values));
    return {std::move(model),std::move(map)};
}
}

TEST(PureHaloSupportExclusionTest, ExcludesExactBoundarySupportAndPreservesOriginalDomain)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto & fitting=problem.Input();
    const auto & original=problem.OriginalInput();
    ASSERT_TRUE(fitting.selection_domain);
    EXPECT_EQ(fitting.selection_domain->contributor_policy,PureHaloSupportExclusionContributorPolicy);
    EXPECT_EQ(original.row_ids.size(),11u);
    EXPECT_EQ(fitting.row_ids.size(),10u);
    EXPECT_EQ(fitting.row_ids.back(),"9");
    EXPECT_EQ(original.row_ids.back(),"10");
    EXPECT_EQ(problem.FittingToOriginalRows(),(std::vector<std::size_t>{0,1,2,3,4,5,6,7,8,9}));
    EXPECT_EQ(problem.FittingToOriginalContributors(),(std::vector<std::size_t>{0}));
    ASSERT_TRUE(original.selection_domain->pure_halo_support_exclusion);
    const auto & provenance=*original.selection_domain->pure_halo_support_exclusion;
    EXPECT_EQ(provenance.original_rows,11u);
    EXPECT_EQ(provenance.fit_rows,10u);
    EXPECT_EQ(provenance.excluded_rows,1u);
    EXPECT_EQ(provenance.candidate_halos,1u);
    EXPECT_EQ(provenance.retained_halos,1u);
    EXPECT_EQ(provenance.zero_support_halos,0u);
    double norm{};
    for (const double value:original.observations) norm+=value*value;
    EXPECT_DOUBLE_EQ(problem.ObservationScale(),std::max(1.0,std::sqrt(norm)));
    for (const auto & support:fitting.support)
        for (const auto & entry:support) EXPECT_LT(entry.row,fitting.row_ids.size());
}

TEST(PureHaloSupportExclusionTest, FitsTargetOnlyWithoutFixedHaloSnapshot)
{
    const JointProblem problem(PureInput());
    EXPECT_TRUE(problem.FixedHalos().empty());
    const auto fit=FitJointComponentsInitialized(problem,TargetInitialization());
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(fit.RuntimeConvergence(),JointCheckStatus::Passed);
    EXPECT_EQ(fit.TargetRuntimeConvergence(),JointCheckStatus::Passed);
    EXPECT_TRUE(fit.available_row_mask==std::vector<bool>(problem.Input().row_ids.size(),true));
    EXPECT_TRUE(fit.evidence.end()!=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="independent-replay" && check.status==JointCheckStatus::Passed;
    }));
    EXPECT_FALSE(std::any_of(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="fixed-halo-replay";
    }));
}

TEST(PureHaloSupportExclusionTest, CapturesOriginalAndFittingIdentitiesSeparately)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto snapshot=CaptureJointAnalysisResult(FitJointComponents(problem,{.65}));
    EXPECT_EQ(snapshot.parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    EXPECT_EQ(snapshot.atom_ids,problem.OriginalInput().atom_ids);
    EXPECT_EQ(snapshot.row_ids,problem.OriginalInput().row_ids);
    EXPECT_EQ(snapshot.fitting_row_ids,problem.Input().row_ids);
    EXPECT_EQ(snapshot.fitting_to_original_contributors,problem.FittingToOriginalContributors());
    EXPECT_EQ(snapshot.fitting_to_original_rows,problem.FittingToOriginalRows());
    ASSERT_EQ(snapshot.available_row_mask.size(),snapshot.row_ids.size());
    EXPECT_FALSE(snapshot.available_row_mask.back());
    for (const auto & component:snapshot.components)
    {
        for (const auto atom:component.atoms) EXPECT_LT(atom,snapshot.atom_ids.size());
        for (const auto row:component.rows) EXPECT_LT(row,snapshot.row_ids.size());
    }
}

TEST(PureHaloSupportExclusionTest, PeelingUsesOriginalRowsAndReconstructedHalo)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    data_internal::ApplyJointStageEstimates(fixture.model,snapshot,"pure-peeling-run");
    const SamplingPointList points{{2.5,{2.5,0,0},true}};
    fixture.model.EditAnalysis().SetAtomLocalRawSamplingEntries(*fixture.model.FindAtomPtr(1),
        core::detail::BuildLocalPotentialSampleList(fixture.map,points));
    const core::detail::JointPostprocessingContext context(
        fixture.map,fixture.model,problem,&snapshot);
    const auto peelings=context.BuildPostFitPeelingSamples(std::vector<std::size_t>{0});
    ASSERT_TRUE(peelings.contains(1));
    ASSERT_EQ(peelings.at(1).samples.size(),1u);
    EXPECT_TRUE(peelings.at(1).samples.front().response);
    EXPECT_TRUE(peelings.at(1).samples.front().reason.empty());
    const auto prepared=context.BuildPostJointLocalSamples(1,peelings.at(1));
    ASSERT_EQ(prepared.samples.size(),1u);
    EXPECT_TRUE(prepared.samples.front().fixed_charge_response);
    EXPECT_TRUE(prepared.samples.front().gaussian_response);
}

TEST(PureHaloSupportExclusionTest, PeelingRejectsFitReconstructionIdentityMismatch)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    snapshot.fitting_row_ids.front()="stale";
    EXPECT_THROW((core::detail::JointPostprocessingContext(fixture.map,fixture.model,problem,&snapshot)),
        std::invalid_argument);
}

TEST(PureHaloSupportExclusionTest, PersistsVersionedDualDomainAndMedianProvenance)
{
    auto fixture=BoundaryFixture();
    const auto problem=BuildPureHaloSupportExclusionProblem(fixture.map,fixture.model);
    const auto snapshot=PurePeelingSnapshot(problem,fixture.model);
    const auto encoded=joint_result_io::Encode(snapshot);
    EXPECT_EQ(boost::json::parse(encoded).at("schema_version").as_int64(),17);
    const auto decoded=joint_result_io::Decode(encoded);
    EXPECT_EQ(decoded.parameterization_contract,JointPureHaloSupportExclusionParameterizationContract);
    EXPECT_EQ(decoded.fitting_to_original_contributors,problem.FittingToOriginalContributors());
    ASSERT_TRUE(decoded.halo_median_reconstruction);
    EXPECT_EQ(decoded.halo_median_reconstruction->halos,snapshot.halo_median_reconstruction->halos);
}
