#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <utility>

namespace {
namespace n=rhbm_gem::core::joint_component;

std::vector<n::Support> MakeSupport()
{
    std::vector<n::Support> support;
    for(Eigen::Index row=0;row<9;++row) support.push_back({row,.36*static_cast<double>(row)});
    return support;
}

struct StabilityFixture
{
    n::Domain domain;
    n::Vector y,eta,beta;
    n::EvaluationContext context;
    StabilityFixture():domain(9,{MakeSupport()}),y(9),eta(1),beta(2)
    {
        eta(0)=std::log(.8); beta<<1.35,.24;
        const double width=std::exp(eta(0));
        for(Eigen::Index row=0;row<y.size();++row)
        {
            const double square=domain.atoms[0][static_cast<std::size_t>(row)].square;
            const double radius=std::sqrt(square);
            const double gaussian=std::exp(-square/(2*width*width)) /
                std::pow(2*std::numbers::pi*width*width,1.5);
            const double charge=radius<1e-5 ? std::sqrt(2/std::numbers::pi)/width :
                std::erf(radius/(std::sqrt(2.0)*width))/radius;
            y(row)=beta(0)*gaussian+beta(1)*charge;
        }
        context=n::CreateContext(y,1);
    }
};

TEST(JointParameterStabilityTest, ExactLocalFitIsStableAndDeterministic)
{
    StabilityFixture fixture;
    const auto accepted_beta=fixture.beta;
    rhbm_gem::JointBlockParameterStability no_history;
    EXPECT_EQ(no_history.status,"NotAssessed");
    EXPECT_EQ(no_history.prediction_status,"NotAssessed");
    rhbm_gem::JointParameterStability empty_trajectory;
    EXPECT_EQ(empty_trajectory.trajectory_stability_status,"insufficient-history");
    const auto first=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        fixture.context,"test/block-1");
    const auto second=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        fixture.context,"test/block-1");
    EXPECT_EQ(first.evidence.status,"Stable");
    EXPECT_EQ(first.evidence.prediction_status,"Stable");
    EXPECT_EQ(first.evidence.objective_consistency_status,"Passed");
    ASSERT_TRUE(first.evidence.baseline_local_ac_kkt);
    EXPECT_LE(*first.evidence.baseline_local_ac_kkt,1e-10);
    ASSERT_TRUE(first.evidence.baseline_projected_width_kkt);
    EXPECT_LE(*first.evidence.baseline_projected_width_kkt,1e-12);
    EXPECT_EQ(first.evidence.baseline_independent_replay_status,"Passed");
    EXPECT_EQ(first.evidence.conditional_reprofile_status,"StrictlyQualified");
    EXPECT_EQ(first.evidence.conditional_reprofile_independent_replay_status,"Passed");
    ASSERT_TRUE(first.evidence.conditional_reprofile_objective_difference);
    EXPECT_LE(*first.evidence.conditional_reprofile_objective_difference,1e-15);
    ASSERT_TRUE(first.evidence.conditional_reprofile_scaled_ac_parameter_difference);
    EXPECT_LE(*first.evidence.conditional_reprofile_scaled_ac_parameter_difference,1e-10);
    EXPECT_TRUE(fixture.beta.isApprox(accepted_beta,0.));
    EXPECT_EQ(first.evidence.perturbation_solves,2u);
    EXPECT_EQ(first.evidence.perturbation_directions,1u);
    EXPECT_EQ(first.evidence.perturbation_seeds,2u);
    ASSERT_EQ(first.evidence.probes.size(),1u);
    EXPECT_EQ(first.evidence.probes.front().direction,"weak-singular");
    EXPECT_NE(first.evidence.probes.front().source.find("normalized-all-width-direction"),std::string::npos);
    EXPECT_LE(first.cpu_seconds,n::MaximumLocalPerturbationCpuSeconds);
    EXPECT_LE(first.evidence.perturbation_evaluations,
        first.evidence.perturbation_seeds*n::MaximumLocalPerturbationProfileEvaluations);
    ASSERT_TRUE(first.evidence.local_prediction_change_relative);
    ASSERT_TRUE(second.evidence.local_prediction_change_relative);
    EXPECT_DOUBLE_EQ(*first.evidence.local_prediction_change_relative,
        *second.evidence.local_prediction_change_relative);
    ASSERT_TRUE(first.evidence.log_width_parameter_change_inf);
    ASSERT_TRUE(second.evidence.log_width_parameter_change_inf);
    EXPECT_DOUBLE_EQ(*first.evidence.log_width_parameter_change_inf,
        *second.evidence.log_width_parameter_change_inf);
    ASSERT_TRUE(first.evidence.scaled_ac_parameter_change_inf);
    EXPECT_LE(*first.evidence.log_width_parameter_change_inf,n::LocalParameterChangeWarning);
    EXPECT_LE(*first.evidence.scaled_ac_parameter_change_inf,n::LocalParameterChangeWarning);
    EXPECT_LE(*first.evidence.local_prediction_change_relative,n::LocalPredictionChangeWarning);
}

TEST(JointParameterStabilityTest, UsesWeakLocalWidthDirectionsWithinTheBoundedBudget)
{
    std::vector<n::Support> first_support,second_support;
    for(Eigen::Index row=0;row<12;++row)
    {
        first_support.push_back({row,.22*static_cast<double>(row)});
        const double radius=.12*static_cast<double>(row)+.31;
        second_support.push_back({row,radius*radius});
    }
    n::Domain domain(12,{first_support,second_support});
    n::Vector eta(2),beta(4),y(12);
    eta<<std::log(.72),std::log(1.08);
    beta<<1.25,.19,.83,.11;
    for(Eigen::Index row=0;row<y.size();++row)
    {
        double value{};
        for(Eigen::Index atom=0;atom<2;++atom)
        {
            const auto & support=domain.atoms[static_cast<std::size_t>(atom)];
            const double radius=std::sqrt(support[static_cast<std::size_t>(row)].square);
            const double width=std::exp(eta(atom));
            const double gaussian=std::exp(-radius*radius/(2*width*width)) /
                std::pow(2*std::numbers::pi*width*width,1.5);
            const double charge=radius<1e-5 ? std::sqrt(2/std::numbers::pi)/width :
                std::erf(radius/(std::sqrt(2.0)*width))/radius;
            value+=beta(2*atom)*gaussian+beta(2*atom+1)*charge;
        }
        y(row)=value;
    }
    const auto context=n::CreateContext(y,1);
    const auto first=n::AssessLocalParameterStability(domain,y,eta,beta,context,"two-atom-block");
    ASSERT_FALSE(first.evidence.probes.empty());
    EXPECT_EQ(first.evidence.probes.front().direction,"weak-singular");
    EXPECT_GE(first.evidence.perturbation_directions,2u);
    EXPECT_LE(first.evidence.perturbation_directions,n::MaximumLocalPerturbationDirectionsPerBlock);
    EXPECT_LE(first.evidence.perturbation_seeds,
        n::LocalPerturbationSeedsPerDirection*n::MaximumLocalPerturbationDirectionsPerBlock);
    EXPECT_LE(first.evidence.perturbation_evaluations,
        first.evidence.perturbation_seeds*n::MaximumLocalPerturbationProfileEvaluations);
    ASSERT_TRUE(first.evidence.log_width_parameter_change_inf);
    ASSERT_TRUE(first.evidence.scaled_ac_parameter_change_inf);
    ASSERT_TRUE(first.evidence.local_prediction_change_relative);
    if(first.evidence.status!="Unavailable")
    {
        if(*first.evidence.local_prediction_change_relative>n::LocalPredictionChangeWarning)
            EXPECT_EQ(first.evidence.status,"Unstable");
        else if(*first.evidence.log_width_parameter_change_inf>n::LocalParameterChangeWarning ||
            *first.evidence.scaled_ac_parameter_change_inf>n::LocalParameterChangeWarning)
            EXPECT_EQ(first.evidence.status,"Sensitive");
        else EXPECT_EQ(first.evidence.status,"Stable");
    }
    for(const auto & probe:first.evidence.probes)
    {
        if(probe.unit_direction.empty()) continue;
        ASSERT_EQ(probe.unit_direction.size(),2u);
        EXPECT_NEAR(std::hypot(probe.unit_direction[0],probe.unit_direction[1]),1.,1e-12);
        EXPECT_LE(probe.evaluations,2*n::MaximumLocalPerturbationProfileEvaluations);
        if(probe.status=="Assessed")
        {
            ASSERT_TRUE(probe.parameter_change);
            ASSERT_TRUE(probe.prediction_change);
            ASSERT_TRUE(probe.objective_change);
        }
    }
    const auto second=n::AssessLocalParameterStability(domain,y,eta,beta,context,"two-atom-block");
    ASSERT_EQ(first.evidence.probes.size(),second.evidence.probes.size());
    for(std::size_t index=0;index<first.evidence.probes.size();++index)
    {
        EXPECT_EQ(first.evidence.probes[index].direction,second.evidence.probes[index].direction);
        EXPECT_EQ(first.evidence.probes[index].unit_direction,second.evidence.probes[index].unit_direction);
    }
}

TEST(JointParameterStabilityTest, AuditsPracticalBaselineThatFailsStrictLocalAcKkt)
{
    StabilityFixture fixture;
    n::Vector accepted_beta;
    bool found{};
    for(int angle=0;angle<720 && !found;++angle)
    {
        const double radians=std::numbers::pi*static_cast<double>(angle)/360.;
        n::Vector direction(2); direction<<std::cos(radians),std::sin(radians);
        for(const double magnitude:{1e-9,1e-8,1e-7,1e-6})
        {
            const n::Vector beta=fixture.beta+magnitude*direction;
            const auto state=n::EvaluateState(fixture.domain,fixture.y,fixture.eta,beta,fixture.context);
            const auto width_kkt=n::ProfileGradientInfinityNorm(state,fixture.context);
            if(!state.valid || !state.certificate.available || !state.certificate.feasible ||
                !(state.certificate.projected_kkt>1e-10 && state.certificate.projected_kkt<=1e-7) ||
                !width_kkt || *width_kkt>1e-9) continue;
            const auto replay=n::CheckIndependentReplay(fixture.domain,fixture.y,state,fixture.context,
                nullptr,nullptr,n::PracticalJointConvergenceCriteria);
            if(!replay.passed) continue;
            accepted_beta=beta; found=true; break;
        }
    }
    ASSERT_TRUE(found);
    const auto assessment=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,accepted_beta,
        fixture.context,"practical-baseline/block-1");
    EXPECT_EQ(assessment.evidence.status,"Unavailable");
    EXPECT_EQ(assessment.evidence.reason,"baseline-local-kkt-failed");
    ASSERT_TRUE(assessment.evidence.baseline_local_ac_kkt);
    EXPECT_GT(*assessment.evidence.baseline_local_ac_kkt,1e-10);
    EXPECT_LE(*assessment.evidence.baseline_local_ac_kkt,1e-7);
    ASSERT_TRUE(assessment.evidence.baseline_projected_width_kkt);
    EXPECT_LE(*assessment.evidence.baseline_projected_width_kkt,1e-9);
    EXPECT_EQ(assessment.evidence.baseline_independent_replay_status,"Passed");
    EXPECT_EQ(assessment.evidence.conditional_reprofile_status,"StrictlyQualified");
    EXPECT_EQ(assessment.evidence.conditional_reprofile_independent_replay_status,"Passed");
    EXPECT_EQ(assessment.evidence.perturbation_solves,0u);
}

TEST(JointParameterStabilityTest, ProfileSearchHonorsCpuBudget)
{
    StabilityFixture fixture;
    auto context=fixture.context;
    context.profile_budget=static_cast<int>(n::MaximumLocalPerturbationProfileEvaluations);
    context.update_budget=6;
    context.profile_cpu_budget_seconds=1e-7;
    const n::Vector perturbed_widths=(fixture.eta.array()+n::LocalWidthPerturbation).exp();
    const auto result=n::SearchProfile(fixture.domain,fixture.y,perturbed_widths,context);
    EXPECT_EQ(result.stop_reason,"cpu-time-budget");
    EXPECT_TRUE(result.stopped);
    EXPECT_LE(result.evaluations,static_cast<int>(n::MaximumLocalPerturbationProfileEvaluations));
}

TEST(JointParameterStabilityTest, ResourceBudgetExhaustionIsUnavailableWithoutSolves)
{
    StabilityFixture fixture;
    const auto result=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        fixture.context,"test/block-1",n::JointMemoryPlanningBudgetBytes);
    EXPECT_EQ(result.evidence.status,"Unavailable");
    EXPECT_EQ(result.evidence.reason.rfind("resource-limit:parameter-stability:",0),0u);
    EXPECT_EQ(result.evidence.perturbation_solves,0u);
    EXPECT_EQ(result.evidence.perturbation_evaluations,0u);
}

TEST(JointParameterStabilityTest, ProbePriorityIsSensitivityBasedAndIndependentOfInputOrder)
{
    const auto block=[](std::string id,std::string classification,double condition) {
        rhbm_gem::JointBlockConditioning value; value.block_id=std::move(id);
        value.conditioning_classification=std::move(classification);
        value.free_ac_condition=condition; value.free_ac_rank=2; value.free_ac_columns=2;
        return value;
    };
    std::vector<rhbm_gem::JointBlockConditioning> forward{
        block("block-a-well","WellConditioned",1e2),
        block("block-b-moderate","IllConditioned",1e7),
        block("block-z-highly-sensitive","IllConditioned",1e10),
        block("block-c-unavailable","Unavailable",0.)};
    const auto identities=[](const auto & blocks,const auto & order) {
        std::vector<std::string> out;
        for(const auto index:order) out.push_back(blocks.at(index).block_id);
        return out;
    };
    const auto forward_order=n::RankFixedNeighborStabilityProbesForTesting(forward);
    ASSERT_FALSE(forward_order.empty());
    EXPECT_EQ(forward.at(forward_order.front()).block_id,"block-z-highly-sensitive");
    EXPECT_EQ(std::find_if(forward_order.begin(),forward_order.end(),[&](std::size_t index) {
        return forward[index].block_id=="block-c-unavailable";
    }),forward_order.end());
    const auto expected_order=identities(forward,forward_order);
    std::reverse(forward.begin(),forward.end());
    const auto reverse_order=n::RankFixedNeighborStabilityProbesForTesting(forward);
    EXPECT_EQ(identities(forward,reverse_order),expected_order);
}

TEST(JointParameterStabilityTest, ComponentFitHonorsPerturbationClaimBudgetAndPreservesPoint)
{
    const rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    n::FixedNeighborSearchPolicy policy; policy.core_atoms=1;
    std::size_t claims{};
    policy.claim_parameter_stability_block=[&] {return ++claims<=1;};
    policy.configure_search_for_testing=[](n::FixedNeighborPolicy & search) {
        search.assess_parameter_stability=true;
    };
    const auto fit=n::FitFixedNeighborComponents(problem,{.55,.55},policy);
    n::FixedNeighborSearchPolicy baseline_policy; baseline_policy.core_atoms=1;
    baseline_policy.configure_search_for_testing=[](n::FixedNeighborPolicy & search) {
        search.assess_parameter_stability=false;
    };
    const auto baseline=n::FitFixedNeighborComponents(problem,{.55,.55},baseline_policy);
    ASSERT_EQ(fit.components.size(),1u);
    ASSERT_EQ(baseline.components.size(),1u);
    const auto & component=fit.components.front();
    const auto & baseline_component=baseline.components.front();
    ASSERT_TRUE(component.state) << component.stop_reason;
    ASSERT_TRUE(baseline_component.state) << baseline_component.stop_reason;
    ASSERT_TRUE(component.parameter_stability);
    EXPECT_EQ(component.state->ac,baseline_component.state->ac);
    EXPECT_EQ(component.state->b,baseline_component.state->b);
    EXPECT_EQ(component.state->log_b,baseline_component.state->log_b);
    EXPECT_DOUBLE_EQ(component.state->objective,baseline_component.state->objective);
    EXPECT_EQ(component.RuntimeConvergence(),baseline_component.RuntimeConvergence());
    const auto & stability=*component.parameter_stability;
    EXPECT_EQ(stability.blocks.size(),2u);
    EXPECT_LE(stability.perturbation_solves,2u);
    EXPECT_LE(claims,2u);
    EXPECT_LE(stability.recent_sweep_count,3u);
    if(stability.recent_sweep_count<3)
        EXPECT_NE(stability.trajectory_stability_status,"Stable");
}
}
