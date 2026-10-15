#include <gtest/gtest.h>
#include "core/detail/joint_component/ParameterStability.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <cmath>
#include <numbers>

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
    EXPECT_EQ(first.evidence.perturbation_solves,2u);
    EXPECT_LE(first.cpu_seconds,n::MaximumLocalPerturbationCpuSeconds);
    EXPECT_LE(first.evidence.perturbation_evaluations,
        2*n::MaximumLocalPerturbationProfileEvaluations+2);
    ASSERT_TRUE(first.evidence.local_prediction_change_relative);
    ASSERT_TRUE(second.evidence.local_prediction_change_relative);
    EXPECT_DOUBLE_EQ(*first.evidence.local_prediction_change_relative,
        *second.evidence.local_prediction_change_relative);
    ASSERT_TRUE(first.evidence.log_width_parameter_change_inf);
    ASSERT_TRUE(second.evidence.log_width_parameter_change_inf);
    EXPECT_DOUBLE_EQ(*first.evidence.log_width_parameter_change_inf,
        *second.evidence.log_width_parameter_change_inf);
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
