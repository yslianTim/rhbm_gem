#include <gtest/gtest.h>

#include "core/detail/HaloMedianReconstruction.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <memory>
#include <limits>
#include <numeric>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <string>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
using namespace rhbm_gem::core::detail;

struct Fixture
{
    ModelObject model;
    JointProblem problem;
    JointAnalysisResult result;
};

Fixture MakeFixture(std::size_t target_count=3)
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int id=1;id<=5;++id)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(id);
        atom->SetElement(Element::CARBON);
        atom->SetComponentKey(id<=4 ? 1 : 3);
        atom->SetAtomKey(id<=3 ? 7 : id==4 ? 7 : 8);
        atom->SetPosition(static_cast<double>(id),0,0);
        atoms.push_back(std::move(atom));
    }
    ModelObject model(std::move(atoms));
    JointProblemInput input;
    input.atom_ids={"1","2","3","4","5"};
    input.row_ids={"0","1","2","3","4"};
    input.observations={0,0,0,0,0};
    input.support={{{0,0}},{{1,0}},{{2,0}},{{3,0}},{{4,0}}};
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->target_indices={0,1,2};
    input.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        5,5,0,2,2,0};
    JointProblem problem(std::move(input));

    JointAnalysisResult result;
    result.parameterization_contract=std::string(JointPureHaloSupportExclusionParameterizationContract);
    result.atom_ids=problem.Input().atom_ids;
    result.row_ids=problem.Input().row_ids;
    result.selection_domain=problem.Input().selection_domain;
    result.runtime_convergence=JointCheckStatus::Passed;
    result.target_runtime_convergence=JointCheckStatus::Passed;
    result.initialization.valid=true;
    result.initialization.b={1,1,1};
    result.components.emplace_back();
    auto & component=result.components.back();
    component.id="component";
    component.atoms={0,1,2};
    component.rows={0,1,2,3,4};
    component.runtime_convergence=JointCheckStatus::Passed;
    component.target_runtime_convergence=JointCheckStatus::Passed;
    component.state=JointState{{1.0,2.0,3.0,5.0,5.0,8.0},{1,2,3},{0,0,0},{0,0,0},0,{}};
    result.assembled_state=component.state;
    result.target_evidence.emplace();
    result.target_evidence->contract=std::string(JointTargetPointValidationContractVersion);
    result.target_evidence->original_rows=5;
    for (std::size_t atom=0;atom<target_count;++atom)
    {
        JointAtomPointValidation point;
        point.atom=atom;
        point.component_convergence=JointCheckStatus::Passed;
        point.parameter_validity=JointCheckStatus::Passed;
        point.local_conditioning="Adequate";
        point.parameter_stability="NotAssessed";
        point.point_quality="NotAssessed";
        point.reason="test";
        result.target_evidence->point_atoms.push_back(std::move(point));
    }
    return {std::move(model),std::move(problem),std::move(result)};
}

Fixture MakeObservationCountFixture(bool all_fixed=false,bool all_free=false,bool include_zero_halo=true,
    bool reverse_model=false)
{
    auto base=MakeFixture();
    JointProblemInput input;
    input.atom_ids={"1","2","3","4","5","6","7"};
    input.row_ids={"0","1","2","3","4","5","6","7","8"};
    input.observations.assign(input.row_ids.size(),0.0);
    input.support={{{0,0}},{{1,0}},{{2,0}},{{3,.1}},{{3,.3},{4,.4}},
        {{5,.5},{6,.6},{7,.7}},{{5,.8},{6,.9},{7,1.0},{8,1.1}}};
    if(all_fixed)
    {
        input.support[5]={{5,.5},{6,.6}};
        input.support[6]={{7,.8}};
    }
    if(all_free)
    {
        input.support[3]={{3,.1},{4,.2},{5,.3}};
        input.support[4]={{3,.4},{4,.5},{5,.6},{6,.7}};
    }
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->target_indices={0,1,2};
    const std::size_t candidate_halos=include_zero_halo?5:4;
    input.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size()-3,3,candidate_halos,4,
        candidate_halos-4};
    JointProblem problem(input);

    std::vector<std::unique_ptr<AtomObject>> atoms;
    const int maximum_id=include_zero_halo?8:7;
    for(int step=0;step<maximum_id;++step)
    {
        const int id=reverse_model ? maximum_id-step : step+1;
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(id);
        atom->SetElement(Element::CARBON);
        atom->SetComponentKey(id%2?1:2);
        atom->SetAtomKey(static_cast<AtomKey>(id));
        atom->SetPosition(static_cast<double>(id),0,0);
        atoms.push_back(std::move(atom));
    }
    ModelObject model(std::move(atoms));
    auto result=std::move(base.result);
    result.atom_ids=input.atom_ids;
    result.row_ids=input.row_ids;
    result.selection_domain=input.selection_domain;
    result.target_evidence->original_rows=input.observations.size();
    return {std::move(model),std::move(problem),std::move(result)};
}

;

}

TEST(HaloMedianReconstructionTest, ComputesSeparateComponentAndWidthAndChargeMedians)
{
    auto fixture=MakeFixture();
    const auto built=BuildHaloMedianReconstruction(fixture.model,fixture.problem,fixture.result,"run-1");
    ASSERT_EQ(built.snapshot.halos.size(),2u);
    EXPECT_EQ(built.audit.valid_target_count,3u);
    EXPECT_EQ(built.audit.geometry_local_halos,2u);
    EXPECT_EQ(built.snapshot.neighbor_count,4u);
    EXPECT_EQ(built.snapshot.donor_policy,GeometryLocalTargetMedianDonorPolicy);
    EXPECT_EQ(built.snapshot.halos[0].source,HaloMedianSource::GeometryLocalTargetMedian);
    EXPECT_DOUBLE_EQ(built.snapshot.halos[0].amplitude,3.0);
    EXPECT_DOUBLE_EQ(built.snapshot.halos[0].width,2.0);
    EXPECT_DOUBLE_EQ(built.snapshot.halos[0].charge,5.0);
    EXPECT_EQ(built.snapshot.halos[0].donor_count,3u);
    EXPECT_EQ(built.snapshot.halos[1].source,HaloMedianSource::GeometryLocalTargetMedian);
    EXPECT_DOUBLE_EQ(built.snapshot.halos[1].amplitude,3.0);
    EXPECT_DOUBLE_EQ(built.snapshot.halos[1].width,2.0);
    EXPECT_DOUBLE_EQ(built.snapshot.halos[1].charge,5.0);
}

TEST(HaloMedianReconstructionTest, AcceptsNotAssessedStabilityButRejectsUnqualifiedDonor)
{
    auto fixture=MakeFixture();
    fixture.result.target_evidence->point_atoms[0].parameter_validity=JointCheckStatus::Failed;
    const auto built=BuildHaloMedianReconstruction(fixture.model,fixture.problem,fixture.result,"run-2");
    EXPECT_EQ(built.snapshot.valid_target_count,2u);
    EXPECT_EQ(built.snapshot.global_median[0],4.0);
    EXPECT_EQ(built.snapshot.global_median[1],2.5);
    EXPECT_EQ(built.snapshot.global_median[2],6.5);
}

TEST(HaloMedianReconstructionTest, ExcludesInvalidWidthAndNonfiniteDonors)
{
    auto invalid_width=MakeFixture();
    invalid_width.result.components.front().state->b[0]=0.0;
    const auto width_built=BuildHaloMedianReconstruction(
        invalid_width.model,invalid_width.problem,invalid_width.result,"invalid-width");
    EXPECT_EQ(width_built.snapshot.valid_target_count,2u);
    EXPECT_EQ(width_built.snapshot.global_median[0],4.0);

    auto nonfinite=MakeFixture();
    nonfinite.result.components.front().state->ac[3]=std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite_built=BuildHaloMedianReconstruction(
        nonfinite.model,nonfinite.problem,nonfinite.result,"nonfinite-charge");
    EXPECT_EQ(nonfinite_built.snapshot.valid_target_count,2u);
    EXPECT_EQ(nonfinite_built.snapshot.global_median[0],3.0);
    EXPECT_EQ(nonfinite_built.snapshot.global_median[2],5.0);
}

TEST(HaloMedianReconstructionTest, PreservesMixedSignChargeMedian)
{
    auto fixture=MakeFixture();
    fixture.result.components.front().state->ac={1.0,-2.0,3.0,4.0,5.0,6.0};
    const auto built=BuildHaloMedianReconstruction(fixture.model,fixture.problem,fixture.result,"mixed-sign");
    EXPECT_EQ(built.snapshot.global_median[2],4.0);
    EXPECT_EQ(built.snapshot.halos.front().charge,4.0);
}

TEST(HaloMedianReconstructionTest, EmptyQualifiedTargetPoolIsExplicitlyRejected)
{
    auto fixture=MakeFixture();
    for (auto & point:fixture.result.target_evidence->point_atoms)
        point.local_conditioning="Warning";
    try
    {
        (void)BuildHaloMedianReconstruction(fixture.model,fixture.problem,fixture.result,"run-3");
        FAIL() << "expected explicit no-donor rejection";
    }
    catch (const std::invalid_argument & error)
    { EXPECT_STREQ(error.what(),"halo-median-no-valid-targets"); }
}

TEST(HaloMedianReconstructionTest, ResultIsIndependentOfModelAtomTraversalOrder)
{
    auto first=MakeFixture();
    auto second=MakeFixture();
    auto atoms=std::vector<std::unique_ptr<AtomObject>>{};
    for (int id=5;id>=1;--id)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(id);
        atom->SetElement(Element::CARBON);
        atom->SetComponentKey(id<=4 ? 1 : 3);
        atom->SetAtomKey(id<=3 ? 7 : id==4 ? 7 : 8);
        atom->SetPosition(static_cast<double>(id),0,0);
        atoms.push_back(std::move(atom));
    }
    ModelObject reversed(std::move(atoms));
    const auto left=BuildHaloMedianReconstruction(first.model,first.problem,first.result,"run-4");
    const auto right=BuildHaloMedianReconstruction(reversed,second.problem,second.result,"run-4");
    EXPECT_EQ(left.snapshot,right.snapshot);
}

TEST(HaloMedianReconstructionTest, ChemicalLabelsDoNotChangeGeometryLocalHaloEstimates)
{
    auto baseline=MakeFixture();
    auto relabeled=MakeFixture();
    relabeled.model.FindAtomPtr(3)->SetComponentKey(2);
    relabeled.model.FindAtomPtr(3)->SetAtomKey(19);
    relabeled.model.FindAtomPtr(3)->SetElement(Element::OXYGEN);
    relabeled.model.FindAtomPtr(3)->SetResidue(Residue::GLY);
    relabeled.model.FindAtomPtr(3)->SetSpot(Spot::O);
    relabeled.model.FindAtomPtr(4)->SetComponentKey(7);
    relabeled.model.FindAtomPtr(4)->SetAtomKey(28);

    const auto original=BuildHaloMedianReconstruction(
        baseline.model,baseline.problem,baseline.result,"label-baseline");
    const auto changed=BuildHaloMedianReconstruction(
        relabeled.model,relabeled.problem,relabeled.result,"label-baseline");

    EXPECT_EQ(original.snapshot,changed.snapshot);
}

TEST(HaloMedianReconstructionTest, MissingQualifiedTargetDonorsRemainUnavailableWithoutDefaultParameters)
{
    auto fixture=MakeObservationCountFixture();
    for(auto & point:fixture.result.target_evidence->point_atoms)
        point.local_conditioning="Warning";
    const auto classified=BuildHaloFullABCPreclassification(
        fixture.model,fixture.problem,fixture.result);
    EXPECT_FALSE(classified.available);
    EXPECT_EQ(classified.unavailable_reason,"no-qualified-joint-target-donors");
    EXPECT_FALSE(classified.global_target_median);
    EXPECT_EQ(classified.halos.size(),5u);
    for(const auto & halo:classified.halos) EXPECT_FALSE(halo.parameter_profile);
}
