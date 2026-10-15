#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/LocalConditioning.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <cmath>
#include <limits>

namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;

n::LocalConditioningAssessment Assess(const n::Domain & domain,n::VectorRef y,const n::Vector & eta,
    const n::Vector & beta,const n::EvaluationContext & context)
{
    return n::AssessFixedNeighborLocalConditioning(domain,y,eta,beta,context,"test/block-1");
}

TEST(JointLocalConditioningTest, FullRankBlockUsesCompleteProfileDerivative)
{
    n::Domain domain(6,{{{0,0.0},{1,.25},{2,1.0},{3,2.25},{4,4.0},{5,6.25}}});
    n::Vector y(6); y<<.8,-.2,.5,1.1,-.6,.3;
    auto context=n::CreateContext(y,1);
    n::Vector eta(1),beta(2); eta<<std::log(.8); beta<<1.2,.3;
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_EQ(result.evidence.status,"assessed");
    EXPECT_EQ(result.evidence.free_ac_columns,2u);
    ASSERT_TRUE(result.evidence.free_ac_rank);
    EXPECT_EQ(*result.evidence.free_ac_rank,2u);
    ASSERT_TRUE(result.evidence.free_ac_condition);
    EXPECT_TRUE(std::isfinite(*result.evidence.free_ac_condition));
    EXPECT_EQ(result.evidence.free_ac_status,"full-rank");
    ASSERT_TRUE(result.evidence.projected_width_rank);
    EXPECT_EQ(*result.evidence.projected_width_rank,1u);
    ASSERT_TRUE(result.evidence.projected_width_condition);
    ASSERT_TRUE(result.evidence.profile_condition);
    EXPECT_EQ(result.evidence.profile_rank,1u);

    const auto state=n::EvaluateState(domain,y,eta,beta,context);
    ASSERT_TRUE(state.valid);
    const auto complete=n::PrepareDerivative(state,context.scale,&context);
    const auto projected=n::PrepareProjectedWidthDerivative(state,context.scale,&context);
    ASSERT_TRUE(complete.valid); ASSERT_TRUE(projected.valid);
    EXPECT_GT(complete.correction.norm(),1e-12);
    n::Matrix width_component,profile_jacobian;
    complete.Rows(0,domain.rows,width_component,profile_jacobian);
    EXPECT_GT((width_component-profile_jacobian).norm(),1e-12);
    EXPECT_GT(result.evidence.workspace_peak_memory_estimate_bytes.value_or(0u),0u);
}

TEST(JointLocalConditioningTest, ExactSingularityIsReportedWithoutInventingAConditionNumber)
{
    const std::vector<std::vector<n::Support>> supports{
        {{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}};
    n::Domain domain(5,supports);
    n::Vector y=n::Vector::Zero(5),eta=n::Vector::Constant(2,std::log(.7)),beta(4);
    beta<<1.0,.2,1.4,-.1;
    auto context=n::CreateContext(y,2);
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_EQ(result.evidence.status,"partial");
    ASSERT_TRUE(result.evidence.free_ac_rank);
    EXPECT_LT(*result.evidence.free_ac_rank,result.evidence.free_ac_columns);
    EXPECT_FALSE(result.evidence.free_ac_condition);
    EXPECT_EQ(result.evidence.free_ac_status,"rank-deficient");
    ASSERT_TRUE(result.evidence.projected_width_rank);
    EXPECT_EQ(*result.evidence.projected_width_rank,2u);
    EXPECT_TRUE(result.evidence.projected_width_condition);
    EXPECT_NE(result.evidence.profile_status,"full-rank");
}

TEST(JointLocalConditioningTest, NearSingularityRetainsTheMeasuredSpectrum)
{
    constexpr double offset=1e-3;
    n::Domain domain(5,{{{0,.1},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,.1+offset},{1,.5-offset},{2,1.0+offset},{3,1.5-offset},{4,2.0+offset}}});
    n::Vector y=n::Vector::Zero(5),eta=n::Vector::Constant(2,std::log(.7)),beta(4);
    beta<<1.0,.2,1.4,-.1;
    auto context=n::CreateContext(y,2);
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_EQ(result.evidence.status,"assessed");
    ASSERT_TRUE(result.evidence.free_ac_rank);
    ASSERT_TRUE(result.evidence.free_ac_condition);
    EXPECT_EQ(*result.evidence.free_ac_rank,result.evidence.free_ac_columns);
    EXPECT_GT(*result.evidence.free_ac_condition,1e3);
}

TEST(JointLocalConditioningTest, ActiveAmplitudeBoundaryIsKeptAsConditionalEvidence)
{
    n::Domain domain(5,{{{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}});
    n::Vector y=n::Vector::Zero(5),eta(1),beta(2); eta<<std::log(.7); beta<<0.,.2;
    auto context=n::CreateContext(y,1);
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_EQ(result.evidence.status,"assessed");
    EXPECT_EQ(result.evidence.active_amplitude_count,1u);
    EXPECT_EQ(result.evidence.free_ac_columns,1u);
    EXPECT_EQ(result.evidence.profile_status,"active-set-boundary");
    EXPECT_EQ(result.evidence.profile_reason,"profile-derivative-may-switch-active-face");
    EXPECT_EQ(result.evidence.projected_width_status,"active-face-full-rank");
}

TEST(JointLocalConditioningTest, ZeroSupportAndNonfiniteInputAreUnavailable)
{
    n::Vector empty_y(0); auto empty_context=n::CreateContext(empty_y,1);
    n::Domain empty_domain(0,{{}});
    n::Vector eta(1),beta(2); eta<<0.; beta<<1.,.2;
    const auto empty=Assess(empty_domain,empty_y,eta,beta,empty_context);
    EXPECT_EQ(empty.evidence.status,"unavailable");
    EXPECT_EQ(empty.evidence.reason,"zero-support");

    n::Domain domain(5,{{{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}});
    n::Vector y=n::Vector::Zero(5); auto context=n::CreateContext(y,1);
    eta(0)=std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite=Assess(domain,y,eta,beta,context);
    EXPECT_EQ(nonfinite.evidence.status,"unavailable");
    EXPECT_EQ(nonfinite.evidence.reason,"nonfinite-local-state");

    eta(0)=0.; context.rank.rows=-1;
    const auto decomposition=Assess(domain,y,eta,beta,context);
    EXPECT_EQ(decomposition.evidence.status,"partial");
    EXPECT_EQ(decomposition.evidence.free_ac_status,"full-rank");
    EXPECT_EQ(decomposition.evidence.projected_width_status,"unavailable");
    EXPECT_EQ(decomposition.evidence.projected_width_reason,"nonfinite-derivative");
}

TEST(JointLocalConditioningTest, MultiBlockDiagnosticsLeaveTheAcceptedFitUnchanged)
{
    const c::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto run=[&](bool diagnostics) {
        n::FixedNeighborSearchPolicy policy; policy.core_atoms=1;
        policy.configure_search_for_testing=[diagnostics](n::FixedNeighborPolicy & search) {
            search.assess_local_conditioning=diagnostics;
            search.assess_parameter_stability=false;
        };
        return n::FitFixedNeighborComponents(problem,{.55,.55},policy);
    };
    const auto with_diagnostics=run(true),without_diagnostics=run(false);
    ASSERT_EQ(with_diagnostics.components.size(),1u);
    ASSERT_EQ(without_diagnostics.components.size(),1u);
    const auto & measured=with_diagnostics.components.front();
    const auto & baseline=without_diagnostics.components.front();
    ASSERT_TRUE(measured.state); ASSERT_TRUE(baseline.state);
    ASSERT_EQ(measured.local_conditioning.size(),2u);
    EXPECT_TRUE(baseline.local_conditioning.empty());
    EXPECT_EQ(measured.state->ac,baseline.state->ac);
    EXPECT_EQ(measured.state->b,baseline.state->b);
    EXPECT_EQ(measured.state->log_b,baseline.state->log_b);
    EXPECT_DOUBLE_EQ(measured.state->objective,baseline.state->objective);
    EXPECT_EQ(measured.RuntimeConvergence(),baseline.RuntimeConvergence());
}
}
