#include <gtest/gtest.h>

#include "core/detail/joint_component/HaloMedianReconstruction.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/LocalConditioning.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "support/JointTestContext.hpp"
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

TEST(HaloMedianReconstructionTest, AcceptsNotAssessedQualityButRejectsUnqualifiedDonor)
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

TEST(HaloMedianReconstructionTest, FullModelConditioningGuardsZeroContributionWidthMedians)
{
    namespace n=rhbm_gem::core::joint_component;
    namespace tc=second_stage_test::matched::joint_abc;
    for(double charge:{0.0,.2})
    {
        JointProblemInput input; input.atom_ids={"1"}; input.support.resize(1);
        for(double square:{0.0,.25,1.0,2.25,4.0,6.25})
        {
            input.support[0].push_back({input.observations.size(),square});
            input.row_ids.push_back(std::to_string(input.observations.size()));
            input.observations.push_back(charge*n::EvaluateKernel(square,.7,2.5).charge);
        }
        input.selection_domain.emplace(); input.selection_domain->target_indices={0};
        JointProblem target(std::move(input));
        const auto & data=JointProblemAccess::Get(target);
        n::Vector eta=n::Vector::Constant(1,std::log(.7)),beta(2); beta<<0.,charge;
        auto block=n::AssessFixedNeighborLocalConditioning(data.domain,data.y,eta,beta,
            data.context,"test/block-1").evidence;
        block.core_atoms={0};
        if(charge==0) EXPECT_EQ(block.conditioning_classification,"Singular");
        else
        {
            EXPECT_EQ(block.projected_width_rank,1u);
            EXPECT_GT(block.projected_width_column_norms.front(),0.0);
            EXPECT_EQ(block.conditioning_classification,"ActiveFace");
        }
        JointState state{{0.,charge},{.7},{std::log(.7)},{0.},0.};
        const auto evidence=n::AssessTargetPointState(target,target.ParameterLayout(),state,
            JointCheckStatus::Passed,{block},JointEvidenceScope::AssembledGlobal,6);
        ASSERT_EQ(evidence.point_atoms.size(),1u);
        EXPECT_EQ(evidence.point_atoms.front().local_conditioning,"Warning");

        auto f=MakeFixture();
        auto & point=*f.result.components.front().state;
        point.ac[0]=0.; point.ac[1]=charge; point.b[0]=.7;
        f.result.target_evidence->point_atoms[0].local_conditioning=
            evidence.point_atoms.front().local_conditioning;
        const auto guarded=BuildHaloMedianReconstruction(f.model,f.problem,f.result,"guarded");
        ASSERT_TRUE(guarded.snapshot);
        EXPECT_EQ(guarded.snapshot->valid_target_count,2u);
        EXPECT_EQ(guarded.snapshot->halos.front().donor_atom_ids,(std::vector<std::string>{"3","2"}));
        EXPECT_DOUBLE_EQ(guarded.snapshot->global_median[1],2.5);
        EXPECT_DOUBLE_EQ(guarded.snapshot->halos.front().width,2.5);

        // Counterfactual gate bypass: a zero-contribution target's B can change
        // another atom's median, even though its own model remains zero.
        f.result.target_evidence->point_atoms[0].local_conditioning="Adequate";
        const auto admitted=QualifiedJointTargetDonors(f.model,f.problem,f.result);
        EXPECT_EQ(admitted.size(),3u); // No amplitude-only exclusion in donor policy.
        if(charge==0)
        {
            point.b[0]=4.;
            const auto bypassed=BuildHaloMedianReconstruction(f.model,f.problem,f.result,"counterfactual");
            EXPECT_DOUBLE_EQ(bypassed.snapshot->halos.front().width,3.0);
        }
    }
}

TEST(HaloMedianReconstructionTest, ConditioningExclusionPreservesSmallPoolProvenanceAndRejectsEmptyPool)
{
    auto f=MakeFixture();
    f.result.target_evidence->point_atoms[0].local_conditioning="Warning";
    f.result.target_evidence->point_atoms[1].local_conditioning="Warning";
    const auto single=BuildHaloMedianReconstruction(f.model,f.problem,f.result,"one-donor");
    ASSERT_TRUE(single.snapshot);
    EXPECT_EQ(single.snapshot->valid_target_count,1u);
    EXPECT_EQ(single.snapshot->global_median,(std::array<double,3>{5.,3.,8.}));
    for(const auto & halo:single.snapshot->halos)
    {
        EXPECT_EQ(halo.donor_count,1u); EXPECT_EQ(halo.donor_atom_ids,(std::vector<std::string>{"3"}));
        EXPECT_EQ(halo.joint_run_id,"one-donor");
        EXPECT_EQ(halo.source,HaloMedianSource::GeometryLocalTargetMedian);
        EXPECT_DOUBLE_EQ(halo.width,3.0);
    }
    f.result.target_evidence->point_atoms[2].local_conditioning="Unavailable";
    EXPECT_THROW(BuildHaloMedianReconstruction(f.model,f.problem,f.result,"empty"),std::invalid_argument);
}
