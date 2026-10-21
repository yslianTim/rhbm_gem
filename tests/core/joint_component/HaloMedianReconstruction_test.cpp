#include <gtest/gtest.h>

#include "core/detail/HaloMedianReconstruction.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <memory>
#include <limits>
#include <vector>
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

Fixture MakeFixture()
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
    component.state=JointState{{1.0,2.0,3.0,5.0,5.0,8.0},{1,2,3},{0,0,0},{0,0,0},0};
    result.assembled_state=component.state;
    result.target_evidence.emplace();
    result.target_evidence->contract=std::string(JointTargetPointValidationContractVersion);
    result.target_evidence->original_rows=5;
    for (std::size_t atom=0;atom<3;++atom)
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


}

TEST(HaloMedianReconstructionTest, ComputesSeparateComponentAndWidthAndChargeMedians)
{
    auto fixture=MakeFixture();
    const auto built=BuildHaloMedianReconstruction(fixture.model,fixture.problem,fixture.result,"run-1");
    ASSERT_TRUE(built.snapshot);
    const auto & snapshot=*built.snapshot;
    ASSERT_EQ(snapshot.halos.size(),2u);
    EXPECT_EQ(built.audit.valid_target_count,3u);
    EXPECT_EQ(built.audit.geometry_local_halos,2u);
    EXPECT_EQ(snapshot.neighbor_count,4u);
    EXPECT_EQ(snapshot.donor_policy,GeometryLocalTargetMedianDonorPolicy);
    EXPECT_EQ(snapshot.halos[0].source,HaloMedianSource::GeometryLocalTargetMedian);
    EXPECT_DOUBLE_EQ(snapshot.halos[0].amplitude,3.0);
    EXPECT_DOUBLE_EQ(snapshot.halos[0].width,2.0);
    EXPECT_DOUBLE_EQ(snapshot.halos[0].charge,5.0);
    EXPECT_EQ(snapshot.halos[0].donor_count,3u);
    EXPECT_EQ(snapshot.halos[1].source,HaloMedianSource::GeometryLocalTargetMedian);
    EXPECT_DOUBLE_EQ(snapshot.halos[1].amplitude,3.0);
    EXPECT_DOUBLE_EQ(snapshot.halos[1].width,2.0);
    EXPECT_DOUBLE_EQ(snapshot.halos[1].charge,5.0);
}

TEST(HaloMedianReconstructionTest, AcceptsNotAssessedStabilityButRejectsUnqualifiedDonor)
{
    auto fixture=MakeFixture();
    fixture.result.target_evidence->point_atoms[0].parameter_validity=JointCheckStatus::Failed;
    const auto built=BuildHaloMedianReconstruction(fixture.model,fixture.problem,fixture.result,"run-2");
    ASSERT_TRUE(built.snapshot);
    EXPECT_EQ(built.snapshot->valid_target_count,2u);
    EXPECT_EQ(built.snapshot->global_median[0],4.0);
    EXPECT_EQ(built.snapshot->global_median[1],2.5);
    EXPECT_EQ(built.snapshot->global_median[2],6.5);
}

TEST(HaloMedianReconstructionTest, ExcludesInvalidWidthAndNonfiniteDonors)
{
    auto invalid_width=MakeFixture();
    invalid_width.result.components.front().state->b[0]=0.0;
    const auto width_built=BuildHaloMedianReconstruction(
        invalid_width.model,invalid_width.problem,invalid_width.result,"invalid-width");
    ASSERT_TRUE(width_built.snapshot);
    EXPECT_EQ(width_built.snapshot->valid_target_count,2u);
    EXPECT_EQ(width_built.snapshot->global_median[0],4.0);

    auto nonfinite=MakeFixture();
    nonfinite.result.components.front().state->ac[3]=std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite_built=BuildHaloMedianReconstruction(
        nonfinite.model,nonfinite.problem,nonfinite.result,"nonfinite-charge");
    ASSERT_TRUE(nonfinite_built.snapshot);
    EXPECT_EQ(nonfinite_built.snapshot->valid_target_count,2u);
    EXPECT_EQ(nonfinite_built.snapshot->global_median[0],3.0);
    EXPECT_EQ(nonfinite_built.snapshot->global_median[2],5.0);
}

TEST(HaloMedianReconstructionTest, PreservesMixedSignChargeMedian)
{
    auto fixture=MakeFixture();
    fixture.result.components.front().state->ac={1.0,-2.0,3.0,4.0,5.0,6.0};
    const auto built=BuildHaloMedianReconstruction(fixture.model,fixture.problem,fixture.result,"mixed-sign");
    ASSERT_TRUE(built.snapshot);
    EXPECT_EQ(built.snapshot->global_median[2],4.0);
    EXPECT_EQ(built.snapshot->halos.front().charge,4.0);
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
    ASSERT_TRUE(left.snapshot);
    ASSERT_TRUE(right.snapshot);
    EXPECT_EQ(*left.snapshot,*right.snapshot);
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

    ASSERT_TRUE(original.snapshot);
    ASSERT_TRUE(changed.snapshot);
    EXPECT_EQ(*original.snapshot,*changed.snapshot);
}
