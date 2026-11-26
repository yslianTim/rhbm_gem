#include <gtest/gtest.h>

#include "core/estimator/joint_component/ContributorReconstruction.hpp"
#include "core/estimator/joint_component/Problem.hpp"
#include "core/estimator/joint_component/LocalConditioning.hpp"
#include "core/estimator/joint_component/TargetEvidence.hpp"
#include "support/JointTestContext.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <memory>
#include <cmath>
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

Fixture MakeFixture(bool buffer=false)
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
    input.selection_domain->contributor_policy=std::string(BufferExpansionContributorPolicy);
    input.selection_domain->output_indices=buffer ? std::vector<std::size_t>{0,1} :
        std::vector<std::size_t>{0,1,2};
    if (buffer) input.selection_domain->buffer_indices={2};
    input.selection_domain->buffer_expansion=JointBufferExpansionProvenance{
        5,5,0,buffer ? 3u : 2u,2,0,buffer ? 1u : 0u};
    JointProblem problem(std::move(input));

    JointAnalysisResult result;
    result.parameterization_contract=std::string(JointBufferExpansionParameterizationContract);
    result.atom_ids=problem.Input().atom_ids;
    result.row_ids=problem.Input().row_ids;
    result.selection_domain=problem.Input().selection_domain;
    result.runtime_convergence=JointCheckStatus::Passed;
    result.target_runtime_convergence=JointCheckStatus::Passed;
    result.initialization.valid=true;
    result.search_completed=true;
    result.evidence.push_back({"global-width-calibration",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,1.25,5.0,"global-initial-stage-q50-k4"});
    result.initialization.tau={1,1,1};
    result.components.emplace_back();
    auto & component=result.components.back();
    component.id="component";
    component.atoms={0,1,2};
    component.rows={0,1,2,3,4};
    component.runtime_convergence=JointCheckStatus::Passed;
    component.target_runtime_convergence=JointCheckStatus::Passed;
    component.state=JointState{{1.0,2.0,3.0,5.0,5.0,8.0},{1,2,3},{0,std::log(2.),std::log(3.)},{0,0,0},0};
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

TEST(ContributorReconstructionTest, ComputesSeparateComponentAndWidthAndChargeMedians)
{
    auto fixture=MakeFixture();
    const auto built=BuildContributorReconstruction(fixture.model,fixture.problem,fixture.result,"run-1");
    ASSERT_TRUE(built.snapshot);
    const auto & snapshot=*built.snapshot;
    ASSERT_EQ(snapshot.replacements.size(),2u);
    EXPECT_EQ(built.audit.qualified_donor_count,3u);
    EXPECT_EQ(built.audit.reconstructed_halos,2u);
    EXPECT_EQ(snapshot.neighbor_count,4u);
    EXPECT_EQ(snapshot.donor_policy,ContributorReconstructionDonorPolicy);
    EXPECT_EQ(snapshot.replacements[0].source,ReconstructionSource::GeometryLocalFittingMedian);
    EXPECT_DOUBLE_EQ(snapshot.replacements[0].amplitude,3.0);
    EXPECT_DOUBLE_EQ(snapshot.replacements[0].width,2.0);
    EXPECT_DOUBLE_EQ(snapshot.replacements[0].charge,5.0);
    EXPECT_EQ(snapshot.replacements[0].donor_count,3u);
    EXPECT_EQ(snapshot.replacements[1].source,ReconstructionSource::GeometryLocalFittingMedian);
    EXPECT_DOUBLE_EQ(snapshot.replacements[1].amplitude,3.0);
    EXPECT_DOUBLE_EQ(snapshot.replacements[1].width,2.0);
    EXPECT_DOUBLE_EQ(snapshot.replacements[1].charge,5.0);
}

TEST(ContributorReconstructionTest, DistinguishesQualityWarningsFromNumericalFailures)
{
    auto f=MakeFixture();
    for (const std::string quality:{"Warning","Unavailable","NotAssessed"})
    {
        f.result.target_evidence->point_atoms[0].local_conditioning=quality;
        const auto built=BuildContributorReconstruction(f.model,f.problem,f.result,"warning");
        ASSERT_EQ(built.snapshot->qualified_donor_count,2u);
        const auto & replacement=built.snapshot->replacements.front();
        EXPECT_EQ(replacement.role,FittingRole::OutputTarget);
        EXPECT_EQ(replacement.rejection_reasons,ConditioningNotAdequateReconstruction);
        EXPECT_EQ(replacement.donor_atom_ids,(std::vector<std::string>{"2","3"}));
        EXPECT_EQ((std::array{replacement.amplitude,replacement.width,replacement.charge}),
            (std::array{4.,2.5,6.5}));
        EXPECT_EQ(*built.joint_points[0],(std::array{1.,1.,2.}));
    }
    auto invalid=f.result;
    invalid.target_evidence->point_atoms[0].parameter_validity=JointCheckStatus::Failed;
    EXPECT_THROW(BuildContributorReconstruction(f.model,f.problem,invalid,"bad-validity"),std::invalid_argument);
    invalid=f.result;
    invalid.components.front().state->tau[0]=0.0;
    EXPECT_THROW(BuildContributorReconstruction(f.model,f.problem,invalid,"bad-width"),std::invalid_argument);
    invalid=f.result;
    invalid.components.front().state->ac[3]=std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(BuildContributorReconstruction(f.model,f.problem,invalid,"nan-charge"),std::invalid_argument);
    invalid=f.result;
    invalid.components.front().runtime_convergence=JointCheckStatus::Failed;
    EXPECT_THROW(BuildContributorReconstruction(f.model,f.problem,invalid,"bad-component"),std::invalid_argument);
    invalid=f.result;
    invalid.runtime_convergence=JointCheckStatus::Failed;
    EXPECT_THROW(BuildContributorReconstruction(f.model,f.problem,invalid,"bad-joint"),std::invalid_argument);
}

TEST(ContributorReconstructionTest, WidthBoundsArePerAtomAndRetainOverlappingRejectionReasons)
{
    auto f=MakeFixture();
    f.result.evidence.push_back({"global-width-lower-bound",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,4.0,1.0,"isotropic-voxel-spacing-k0.25"});
    auto & state=*f.result.components.front().state;
    state.log_tau[0]=std::nextafter(0.,1.);
    state.tau[2]=5.;
    state.log_tau[2]=std::log(5.);
    f.result.target_evidence->point_atoms[0].local_conditioning="Warning";
    const auto built=BuildContributorReconstruction(f.model,f.problem,f.result,"boundary");
    ASSERT_EQ(built.snapshot->qualified_donor_count,1u);
    ASSERT_EQ(built.snapshot->replacements.size(),4u);
    const auto & lower=built.snapshot->replacements[0];
    const auto & upper=built.snapshot->replacements[1];
    EXPECT_EQ(lower.contributor_index,0u);
    EXPECT_EQ(lower.rejection_reasons,WidthLowerBoundReconstruction|ConditioningNotAdequateReconstruction);
    EXPECT_EQ(upper.contributor_index,2u);
    EXPECT_EQ(upper.rejection_reasons,WidthUpperBoundReconstruction);
    for (const auto & replacement:built.snapshot->replacements)
        EXPECT_EQ(replacement.donor_atom_ids,(std::vector<std::string>{"2"}));
    EXPECT_EQ(built.audit.overlapping_rejections,1u);
}

TEST(ContributorReconstructionTest, QualifiedBufferDonatesAndBothWidthFacesRetainBufferRole)
{
    auto f=MakeFixture(true);
    f.result.evidence.push_back({"global-width-lower-bound",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,2.0,.5,"isotropic-voxel-spacing-k0.25"});
    auto interior=BuildContributorReconstruction(f.model,f.problem,f.result,"buffer-interior");
    ASSERT_EQ(interior.donor_index->Donors().size(),3u);
    EXPECT_EQ(interior.donor_index->Donors().back().contributor_index,2u);

    for (const auto & [width,reason]:{std::pair{.5,WidthLowerBoundReconstruction},
        std::pair{5.,WidthUpperBoundReconstruction}})
    {
        auto & state=*f.result.components.front().state;
        state.tau[2]=width;
        state.log_tau[2]=std::log(width);
        const auto built=BuildContributorReconstruction(f.model,f.problem,f.result,"buffer-face");
        ASSERT_EQ(built.audit.replaced_buffers,1u);
        ASSERT_EQ(built.donor_index->Donors().size(),2u);
        const auto & buffer=built.snapshot->replacements.front();
        EXPECT_EQ(buffer.role,FittingRole::BufferTarget);
        EXPECT_EQ(buffer.rejection_reasons,reason);
        EXPECT_EQ(buffer.donor_atom_ids,(std::vector<std::string>{"2","1"}));
        EXPECT_EQ((*built.joint_points[2])[1],width);
    }
}

TEST(ContributorReconstructionTest, EquidistantAndDuplicatePositionsUseCanonicalAtomId)
{
    const ReconstructionDonorIndex index({
        {0,"8",{1,0,0},{1,1,1}}, {1,"7",{-1,0,0},{1,1,1}}, {2,"12",{1,0,0},{1,1,1}},
        {3,"6",{0,1,0},{1,1,1}}, {4,"9",{0,-1,0},{1,1,1}}});
    std::vector<std::string> ids;
    for (const auto * donor:index.Query({0,0,0})) ids.push_back(donor->atom_id);
    EXPECT_EQ(ids,(std::vector<std::string>{"12","6","7","8"}));
}

TEST(ContributorReconstructionTest, PreservesMixedSignChargeMedian)
{
    auto fixture=MakeFixture();
    fixture.result.components.front().state->ac={1.0,-2.0,3.0,4.0,5.0,6.0};
    const auto built=BuildContributorReconstruction(fixture.model,fixture.problem,fixture.result,"mixed-sign");
    ASSERT_TRUE(built.snapshot);
    EXPECT_EQ(built.snapshot->replacements.front().charge,4.0);
}

TEST(ContributorReconstructionTest, EmptyQualifiedTargetPoolIsExplicitlyRejected)
{
    auto fixture=MakeFixture();
    for (auto & point:fixture.result.target_evidence->point_atoms)
        point.local_conditioning="Warning";
    const auto built=BuildContributorReconstruction(fixture.model,fixture.problem,fixture.result,"run-3");
    EXPECT_EQ(built.audit.unavailable,fixture.result.atom_ids.size());
    for (const auto & replacement:built.snapshot->replacements)
    {
        EXPECT_EQ(replacement.source,ReconstructionSource::Unavailable);
        EXPECT_TRUE(replacement.donor_atom_ids.empty());
    }
}

TEST(ContributorReconstructionTest, ResultIsIndependentOfModelAtomTraversalOrder)
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
    const auto left=BuildContributorReconstruction(first.model,first.problem,first.result,"run-4");
    const auto right=BuildContributorReconstruction(reversed,second.problem,second.result,"run-4");
    ASSERT_TRUE(left.snapshot);
    ASSERT_TRUE(right.snapshot);
    EXPECT_EQ(*left.snapshot,*right.snapshot);
}

TEST(ContributorReconstructionTest, ChemicalLabelsDoNotChangeGeometryLocalHaloEstimates)
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

    const auto original=BuildContributorReconstruction(
        baseline.model,baseline.problem,baseline.result,"label-baseline");
    const auto changed=BuildContributorReconstruction(
        relabeled.model,relabeled.problem,relabeled.result,"label-baseline");

    ASSERT_TRUE(original.snapshot);
    ASSERT_TRUE(changed.snapshot);
    EXPECT_EQ(*original.snapshot,*changed.snapshot);
}

TEST(ContributorReconstructionTest, FullModelConditioningGuardsZeroContributionWidthMedians)
{
    namespace n=rhbm_gem::core::joint_component;
    namespace tc=joint_component_test::matched::joint_component;
    for(double charge:{0.0,.2})
    {
        JointProblemInput input; input.atom_ids={"1"}; input.support.resize(1);
        for(double square:{0.0,.25,1.0,2.25,4.0,6.25})
        {
            input.support[0].push_back({input.observations.size(),square});
            input.row_ids.push_back(std::to_string(input.observations.size()));
            input.observations.push_back(charge*n::EvaluateKernel(square,.7,2.5).charge);
        }
        input.selection_domain.emplace(); input.selection_domain->output_indices={0};
        JointProblem target(std::move(input));
        const auto & data=JointProblemAccess::Get(target);
        n::Vector eta=n::Vector::Constant(1,std::log(.7)),beta(2); beta<<0.,charge;
        auto block=n::AssessBcdLocalConditioning(data.domain,data.y,eta,beta,
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
        point.ac[0]=0.; point.ac[1]=charge; point.tau[0]=.7; point.log_tau[0]=std::log(.7);
        f.result.target_evidence->point_atoms[0].local_conditioning=
            evidence.point_atoms.front().local_conditioning;
        const auto guarded=BuildContributorReconstruction(f.model,f.problem,f.result,"guarded");
        ASSERT_TRUE(guarded.snapshot);
        EXPECT_EQ(guarded.snapshot->qualified_donor_count,2u);
        EXPECT_EQ(guarded.snapshot->replacements.back().donor_atom_ids,(std::vector<std::string>{"3","2"}));
        EXPECT_DOUBLE_EQ(guarded.snapshot->replacements.back().width,2.5);

        // Counterfactual gate bypass: a zero-contribution target's tau can change
        // another atom's median, even though its own model remains zero.
        f.result.target_evidence->point_atoms[0].local_conditioning="Adequate";
        const auto admitted=CollectQualifiedReconstructionDonors(f.model,f.result);
        EXPECT_EQ(admitted.size(),3u); // No amplitude-only exclusion in donor policy.
        if(charge==0)
        {
            point.tau[0]=4.; point.log_tau[0]=std::log(4.);
            const auto bypassed=BuildContributorReconstruction(f.model,f.problem,f.result,"counterfactual");
            EXPECT_DOUBLE_EQ(bypassed.snapshot->replacements.front().width,3.0);
        }
    }
}

TEST(ContributorReconstructionTest, ConditioningExclusionPreservesSmallPoolProvenanceAndRejectsEmptyPool)
{
    auto f=MakeFixture();
    f.result.target_evidence->point_atoms[0].local_conditioning="Warning";
    f.result.target_evidence->point_atoms[1].local_conditioning="Warning";
    const auto single=BuildContributorReconstruction(f.model,f.problem,f.result,"one-donor");
    ASSERT_TRUE(single.snapshot);
    EXPECT_EQ(single.snapshot->qualified_donor_count,1u);
    for(const auto & halo:single.snapshot->replacements)
    {
        EXPECT_EQ(halo.donor_count,1u); EXPECT_EQ(halo.donor_atom_ids,(std::vector<std::string>{"3"}));
        EXPECT_EQ(halo.joint_run_id,"one-donor");
        EXPECT_EQ(halo.source,ReconstructionSource::GeometryLocalFittingMedian);
        EXPECT_DOUBLE_EQ(halo.width,3.0);
    }
    f.result.target_evidence->point_atoms[2].local_conditioning="Unavailable";
    EXPECT_EQ(BuildContributorReconstruction(f.model,f.problem,f.result,"empty").audit.unavailable,5u);
}
