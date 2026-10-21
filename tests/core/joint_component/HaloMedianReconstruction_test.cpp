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

struct FullABCFixture
{
    ModelObject model;
    JointProblem problem;
    JointAnalysisResult result;
    std::array<double,3> free_truth{};
    double zero_row_residual_square{};
};

FullABCFixture MakeFullABCFixture(std::size_t free_support_count=30)
{
    constexpr std::size_t target_count=20,row_count=32;
    constexpr std::array<double,3> target_point{1.0,.8,.1};
    constexpr std::array<double,3> free_truth{1.4,.95,-.12};
    const std::array<double,3> fixed_point=target_point;
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int id=1;id<=23;++id)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(id);
        atom->SetElement(Element::CARBON);
        atom->SetPosition(4.0*static_cast<double>(id),0,0);
        atoms.push_back(std::move(atom));
    }
    ModelObject model(std::move(atoms));
    JointProblemInput input;
    input.row_ids.reserve(row_count);
    input.observations.assign(row_count,0.0);
    input.atom_ids.reserve(target_count+2);
    input.support.resize(target_count+2);
    for (std::size_t row=0;row<row_count;++row)
    {
        input.row_ids.push_back(std::to_string(row));
        const double target_square=.03*static_cast<double>(row);
        const auto target_basis=joint_component::EvaluateKernel(target_square,target_point[1],2.5);
        input.observations[row]=static_cast<double>(target_count)*
            (target_point[0]*target_basis.gaussian+target_point[2]*target_basis.charge);
        for (std::size_t target=0;target<target_count;++target)
            input.support[target].push_back({row,target_square});
    }
    for (std::size_t target=0;target<target_count;++target)
        input.atom_ids.push_back(std::to_string(target+1));
    input.atom_ids.push_back("21");
    input.atom_ids.push_back("22");
    for (std::size_t row=0;row<free_support_count;++row)
    {
        const double square=free_support_count==1 ? .2 :
            .02+6.2*static_cast<double>(row)/static_cast<double>(free_support_count-1);
        input.support[target_count].push_back({row,square});
        const auto basis=joint_component::EvaluateKernel(square,free_truth[1],2.5);
        input.observations[row]+=free_truth[0]*basis.gaussian+free_truth[2]*basis.charge;
    }
    for (const auto [row,square]:std::array<std::pair<std::size_t,double>,2>{{{3,.35},{7,.8}}})
    {
        input.support[target_count+1].push_back({row,square});
        const auto basis=joint_component::EvaluateKernel(square,fixed_point[1],2.5);
        input.observations[row]+=fixed_point[0]*basis.gaussian+fixed_point[2]*basis.charge;
    }
    input.observations[30]+=.25;
    input.observations[31]-=.2;
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->target_indices.resize(target_count);
    std::iota(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),0);
    input.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        row_count,row_count,0,3,2,1};
    JointProblem problem(input);

    JointAnalysisResult result;
    result.parameterization_contract=std::string(JointPureHaloSupportExclusionParameterizationContract);
    result.atom_ids=input.atom_ids;
    result.row_ids=input.row_ids;
    result.selection_domain=input.selection_domain;
    result.initialization.valid=true;
    result.initialization.reason="test-qualified-first-stage-widths";
    result.initialization.data_scope="first-stage-target-widths";
    result.initialization.b.assign(target_count,.7);
    for (std::size_t target=0;target<target_count;++target)
    {
        JointInitializationAtom initialized;
        initialized.id=input.atom_ids[target];
        initialized.target=true;
        initialized.first_stage_fit_qualified=true;
        initialized.native_status=static_cast<int>(RHBMEstimationStatus::SUCCESS);
        initialized.first_stage_qualification=RHBMSolveQualification::NativeSuccess;
        initialized.mdpde={target_point[0],.7,target_point[2]};
        result.initialization.atoms.push_back(std::move(initialized));
    }
    result.search_completed=true;
    result.runtime_convergence=JointCheckStatus::Passed;
    result.target_runtime_convergence=JointCheckStatus::Passed;
    JointAnalysisComponent component;
    component.id="frozen-targets";
    component.atoms=input.selection_domain->target_indices;
    component.state.emplace();
    component.state->b.assign(target_count,target_point[1]);
    component.state->ac.resize(2*target_count);
    component.state->log_b.assign(target_count,std::log(target_point[1]));
    component.state->width_gradient.assign(target_count,0.0);
    result.target_evidence.emplace();
    for (std::size_t target=0;target<target_count;++target)
    {
        component.state->ac[2*target]=target_point[0];
        component.state->ac[2*target+1]=target_point[2];
        JointAtomPointValidation point;
        point.atom=target;
        point.component_convergence=JointCheckStatus::Passed;
        point.parameter_validity=JointCheckStatus::Passed;
        point.local_conditioning="Adequate";
        point.point_quality="Available";
        result.target_evidence->point_atoms.push_back(std::move(point));
    }
    component.runtime_convergence=JointCheckStatus::Passed;
    component.target_runtime_convergence=JointCheckStatus::Passed;
    result.components.push_back(std::move(component));
    result.assembled_state=result.components.front().state;
    result.target_evidence->contract=std::string(JointTargetPointValidationContractVersion);
    result.target_evidence->original_rows=row_count;
    const double raw_norm=std::sqrt(std::inner_product(input.observations.begin(),input.observations.end(),
        input.observations.begin(),0.0));
    return {std::move(model),std::move(problem),std::move(result),free_truth,
        (.25*.25+.2*.2)/(2.0*std::max(1.0,raw_norm)*std::max(1.0,raw_norm))};
}
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

TEST(HaloFullABCJointEstimationTest, FitsFreeHaloOnOriginalRowsAfterFrozenTargetAndFixedHaloSubtraction)
{
    auto fixture=MakeFullABCFixture();
    const auto target_before=fixture.result.components.front().state;
    const auto estimated=EstimateHaloFullABCJoint(
        fixture.model,fixture.problem,fixture.result,"fullabc-test-run");
    ASSERT_TRUE(estimated.available) << estimated.unavailable_reason;
    ASSERT_TRUE(estimated.diagnostics.solver_invoked);
    ASSERT_TRUE(estimated.diagnostics.solver_qualified);
    ASSERT_EQ(estimated.diagnostics.solver_status,JointCheckStatus::Passed);
    ASSERT_EQ(estimated.diagnostics.original_observation_rows,32u);
    EXPECT_EQ(estimated.diagnostics.zero_sensitivity_rows,2u);
    ASSERT_EQ(estimated.fixed_halo_count,2u);
    ASSERT_EQ(estimated.free_halo_count,1u);
    ASSERT_EQ(estimated.diagnostics.qualified_first_stage_width_donors,20u);
    ASSERT_DOUBLE_EQ(*estimated.diagnostics.width_upper_bound,2.8);
    const auto free=std::find_if(estimated.halos.begin(),estimated.halos.end(),[](const auto & halo) {
        return halo.classification.atom_id=="21";
    });
    ASSERT_NE(free,estimated.halos.end());
    ASSERT_TRUE(free->reconstructed_parameters);
    EXPECT_NEAR((*free->reconstructed_parameters)[0],fixture.free_truth[0],2e-3);
    EXPECT_NEAR((*free->reconstructed_parameters)[1],fixture.free_truth[1],2e-3);
    EXPECT_NEAR((*free->reconstructed_parameters)[2],fixture.free_truth[2],2e-3);
    ASSERT_TRUE(estimated.diagnostics.objective);
    EXPECT_GE(*estimated.diagnostics.objective,fixture.zero_row_residual_square-1e-14);
    EXPECT_EQ(fixture.result.components.front().state->ac,target_before->ac);
    EXPECT_EQ(fixture.result.components.front().state->b,target_before->b);
    const auto zero=std::find_if(estimated.halos.begin(),estimated.halos.end(),[](const auto & halo) {
        return halo.classification.atom_id=="23";
    });
    ASSERT_NE(zero,estimated.halos.end());
    EXPECT_EQ(zero->classification.distinct_original_observations,0u);
    ASSERT_TRUE(zero->reconstructed_parameters);
    EXPECT_EQ(*zero->reconstructed_parameters,*estimated.global_target_median);
}

TEST(HaloFullABCJointEstimationTest, AllFixedHalosSkipTheFullABCOptimizer)
{
    auto fixture=MakeFullABCFixture(2);
    const auto estimated=EstimateHaloFullABCJoint(
        fixture.model,fixture.problem,fixture.result,"all-fixed-run");
    ASSERT_TRUE(estimated.available) << estimated.unavailable_reason;
    EXPECT_EQ(estimated.fixed_halo_count,3u);
    EXPECT_EQ(estimated.free_halo_count,0u);
    EXPECT_FALSE(estimated.diagnostics.solver_invoked);
    EXPECT_EQ(estimated.diagnostics.solver_status,JointCheckStatus::NotRun);
    for (const auto & halo:estimated.halos)
    {
        ASSERT_TRUE(halo.reconstructed_parameters);
        EXPECT_EQ(*halo.reconstructed_parameters,*estimated.global_target_median);
    }
}

TEST(HaloFullABCJointEstimationTest, ThreeObservationsAreFreeAndNoDonorsAreUnavailable)
{
    auto fixture=MakeFullABCFixture(3);
    auto estimated=EstimateHaloFullABCJoint(fixture.model,fixture.problem,fixture.result,"three-row-run");
    ASSERT_EQ(estimated.free_halo_count,1u);
    const auto free=std::find_if(estimated.halos.begin(),estimated.halos.end(),[](const auto & halo) {
        return halo.classification.atom_id=="21";
    });
    ASSERT_NE(free,estimated.halos.end());
    EXPECT_EQ(free->classification.distinct_original_observations,3u);
    EXPECT_EQ(free->classification.role,HaloFullABCRole::FreeFullABC);
    EXPECT_TRUE(estimated.diagnostics.solver_invoked);

    for (auto & point:fixture.result.target_evidence->point_atoms) point.local_conditioning="Warning";
    const auto unavailable=EstimateHaloFullABCJoint(
        fixture.model,fixture.problem,fixture.result,"no-donor-run");
    EXPECT_FALSE(unavailable.available);
    EXPECT_EQ(unavailable.unavailable_reason,"no-qualified-joint-target-donors");
    for (const auto & halo:unavailable.halos) EXPECT_FALSE(halo.reconstructed_parameters);
}

TEST(HaloFullABCJointEstimationTest, NoHaloCandidatesRequireNoMedianOrSolver)
{
    auto fixture=MakeFullABCFixture();
    auto input=fixture.problem.OriginalInput();
    input.atom_ids.resize(20);
    input.support.resize(20);
    input.selection_domain->pure_halo_support_exclusion=
        JointPureHaloSupportExclusionProvenance{32,32,0,0,0,0};
    JointProblem target_problem(input);
    auto result=fixture.result;
    result.atom_ids=input.atom_ids;
    result.selection_domain=input.selection_domain;
    std::vector<std::unique_ptr<AtomObject>> targets;
    for (int id=1;id<=20;++id)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(id);
        atom->SetElement(Element::CARBON);
        atom->SetPosition(4.0*static_cast<double>(id),0,0);
        targets.push_back(std::move(atom));
    }
    ModelObject target_model(std::move(targets));
    const auto estimated=EstimateHaloFullABCJoint(target_model,target_problem,result,"no-halo-run");
    EXPECT_TRUE(estimated.available);
    EXPECT_TRUE(estimated.halos.empty());
    EXPECT_FALSE(estimated.global_target_median);
    EXPECT_FALSE(estimated.diagnostics.solver_invoked);
}

TEST(HaloFullABCJointEstimationTest, ChemicalLabelsDoNotAffectTheHaloJointEstimate)
{
    auto baseline=MakeFullABCFixture();
    auto relabeled=MakeFullABCFixture();
    for (const int id:{3,7,21,22})
    {
        auto * atom=relabeled.model.FindAtomPtr(id);
        atom->SetElement(Element::OXYGEN);
        atom->SetComponentKey(static_cast<ComponentKey>(id%3+1));
        atom->SetAtomKey(static_cast<AtomKey>(id+100));
        atom->SetResidue(Residue::GLY);
    }
    const auto left=EstimateHaloFullABCJoint(
        baseline.model,baseline.problem,baseline.result,"label-invariant-run");
    const auto right=EstimateHaloFullABCJoint(
        relabeled.model,relabeled.problem,relabeled.result,"label-invariant-run");
    ASSERT_TRUE(left.available) << left.unavailable_reason;
    ASSERT_TRUE(right.available) << right.unavailable_reason;
    ASSERT_EQ(left.halos.size(),right.halos.size());
    for (std::size_t index=0;index<left.halos.size();++index)
    {
        EXPECT_EQ(left.halos[index].classification.atom_id,right.halos[index].classification.atom_id);
        EXPECT_EQ(left.halos[index].classification.role,right.halos[index].classification.role);
        EXPECT_EQ(left.halos[index].classification.distinct_original_observations,
            right.halos[index].classification.distinct_original_observations);
        EXPECT_EQ(left.halos[index].reconstructed_parameters,right.halos[index].reconstructed_parameters);
    }
}

TEST(HaloMedianReconstructionTest, ClassifiesDistinctOriginalObservationsAtThreeOnly)
{
    auto fixture=MakeObservationCountFixture();
    const auto classified=BuildHaloFullABCPreclassification(
        fixture.model,fixture.problem,fixture.result);
    ASSERT_TRUE(classified.available);
    ASSERT_TRUE(classified.global_target_median);
    EXPECT_EQ(*classified.global_target_median,(std::array<double,3>{3.0,2.0,5.0}));
    ASSERT_EQ(classified.halos.size(),5u);
    EXPECT_EQ(classified.fixed_halo_count,3u);
    EXPECT_EQ(classified.free_halo_count,2u);
    EXPECT_EQ(classified.halos[0].atom_id,"4");
    EXPECT_EQ(classified.halos[0].distinct_original_observations,1u);
    EXPECT_EQ(classified.halos[0].original_support_rows,(std::vector<std::size_t>{3}));
    EXPECT_EQ(classified.halos[0].role,HaloFullABCRole::FixedGlobalTargetMedian);
    EXPECT_EQ(classified.halos[1].distinct_original_observations,2u);
    EXPECT_EQ(classified.halos[1].role,HaloFullABCRole::FixedGlobalTargetMedian);
    EXPECT_EQ(classified.halos[2].distinct_original_observations,3u);
    EXPECT_EQ(classified.halos[2].role,HaloFullABCRole::FreeFullABC);
    EXPECT_EQ(classified.halos[3].distinct_original_observations,4u);
    EXPECT_EQ(classified.halos[3].role,HaloFullABCRole::FreeFullABC);
    EXPECT_EQ(classified.halos[4].distinct_original_observations,0u);
    EXPECT_FALSE(classified.halos[4].original_contributor_index);
    EXPECT_EQ(classified.halos[4].role,HaloFullABCRole::FixedGlobalTargetMedian);
    for(const auto & halo:classified.halos)
        EXPECT_EQ(halo.parameter_profile,*classified.global_target_median);
}

TEST(HaloMedianReconstructionTest, ThreeEqualRadiusRowsRemainFreeForFullABC)
{
    auto fixture=MakeObservationCountFixture(false,true,false);
    auto input=fixture.problem.OriginalInput();
    for(auto & support:input.support[5]) support.squared_distance=.5;
    fixture.problem=JointProblem(std::move(input));
    const auto classified=BuildHaloFullABCPreclassification(
        fixture.model,fixture.problem,fixture.result);
    ASSERT_TRUE(classified.available);
    const auto halo=std::find_if(classified.halos.begin(),classified.halos.end(),[](const auto & value) {
        return value.atom_id=="6";
    });
    ASSERT_NE(halo,classified.halos.end());
    EXPECT_EQ(halo->distinct_original_observations,3u);
    EXPECT_EQ(halo->role,HaloFullABCRole::FreeFullABC);
}

TEST(HaloMedianReconstructionTest, DuplicateSupportEntriesDoNotIncreaseDistinctObservationCount)
{
    EXPECT_EQ(DistinctOriginalSupportRows({{7,.25},{7,1.0},{9,2.0}}),
        (std::vector<std::size_t>{7,9}));
}

TEST(HaloMedianReconstructionTest, AllFixedAndAllFreeClassificationsDoNotUseRankGates)
{
    auto all_fixed=MakeObservationCountFixture(true,false,true);
    const auto fixed=BuildHaloFullABCPreclassification(all_fixed.model,all_fixed.problem,all_fixed.result);
    ASSERT_TRUE(fixed.available);
    EXPECT_EQ(fixed.free_halo_count,0u);
    EXPECT_EQ(fixed.fixed_halo_count,fixed.halos.size());

    auto all_free=MakeObservationCountFixture(false,true,false);
    const auto free=BuildHaloFullABCPreclassification(all_free.model,all_free.problem,all_free.result);
    ASSERT_TRUE(free.available);
    EXPECT_EQ(free.fixed_halo_count,0u);
    EXPECT_EQ(free.free_halo_count,free.halos.size());
}

TEST(HaloMedianReconstructionTest, ClassificationIsIndependentOfTraversalAndChemicalLabels)
{
    auto baseline=MakeObservationCountFixture();
    auto changed=MakeObservationCountFixture(false,false,true,true);
    changed.model.FindAtomPtr(5)->SetComponentKey(77);
    changed.model.FindAtomPtr(5)->SetAtomKey(1001);
    changed.model.FindAtomPtr(5)->SetElement(Element::OXYGEN);
    const auto first=BuildHaloFullABCPreclassification(baseline.model,baseline.problem,baseline.result);
    const auto second=BuildHaloFullABCPreclassification(changed.model,changed.problem,changed.result);
    EXPECT_EQ(first,second);
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
