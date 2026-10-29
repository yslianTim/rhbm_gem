#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointTestContext.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>
#include <utility>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace tc=second_stage_test::matched::joint_abc;

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
        context=tc::MakeContext(y,1);
        context.eta_upper_bound=std::log(2.0);
    }
};

TEST(JointParameterStabilityTest, ExactLocalFitIsStableAndDeterministic)
{
    StabilityFixture fixture;
    const auto accepted_beta=fixture.beta;
    rhbm_gem::JointBlockParameterStability no_history;
    EXPECT_EQ(no_history.status,"NotAssessed");
    EXPECT_EQ(no_history.prediction_status,"NotAssessed");
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
    EXPECT_EQ(first.evidence.baseline_numerical_validity_status,"Valid");
    EXPECT_EQ(first.evidence.baseline_conditional_optimality_status,"StrictlyQualified");
    EXPECT_EQ(first.evidence.baseline_eligibility_status,"Eligible");
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
    EXPECT_EQ(first.evidence.probes.front().status,"Assessed");
    EXPECT_EQ(first.evidence.probes.front().minus_seed_status,"Assessed");
    EXPECT_EQ(first.evidence.probes.front().plus_seed_status,"Assessed");
    ASSERT_TRUE(first.evidence.probes.front().minus_effective_displacement_inf);
    EXPECT_NEAR(*first.evidence.probes.front().minus_effective_displacement_inf,n::LocalWidthPerturbation,1e-15);
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

TEST(JointParameterStabilityTest, PureClassifierRespectsIndependentThresholdBoundaries)
{
    const double parameter_threshold=n::LocalParameterChangeWarning;
    const auto stable=n::ClassifyLocalParameterStability(parameter_threshold,parameter_threshold,
        n::LocalPredictionChangeWarning,true,true);
    EXPECT_EQ(stable.status,"Stable");
    EXPECT_EQ(stable.prediction_status,"Stable");

    const auto sensitive=n::ClassifyLocalParameterStability(std::nextafter(parameter_threshold,INFINITY),
        parameter_threshold,parameter_threshold,true,true);
    EXPECT_EQ(sensitive.status,"Sensitive");
    EXPECT_EQ(sensitive.prediction_status,"Stable");

    const auto unstable=n::ClassifyLocalParameterStability(parameter_threshold,parameter_threshold,
        std::nextafter(n::LocalPredictionChangeWarning,INFINITY),true,true);
    EXPECT_EQ(unstable.status,"Unstable");
    EXPECT_EQ(unstable.prediction_status,"Changed");

    const auto incomplete=n::ClassifyLocalParameterStability(parameter_threshold,parameter_threshold,
        parameter_threshold,true,false,1e-6,1e-6);
    EXPECT_EQ(incomplete.status,"Unavailable");
    EXPECT_EQ(incomplete.prediction_status,"Stable");

    const auto unavailable=n::ClassifyLocalParameterStability(0.,0.,0.,false,false,1e-6,1e-6);
    EXPECT_EQ(unavailable.status,"Unavailable");
    EXPECT_EQ(unavailable.prediction_status,"Unavailable");

    const auto nonfinite=n::ClassifyLocalParameterStability(
        std::numeric_limits<double>::quiet_NaN(),0.,0.,true,true,1e-6,1e-6);
    EXPECT_EQ(nonfinite.status,"Unavailable");
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
    auto context=tc::MakeContext(y,1);
    context.eta_upper_bound=eta(1);
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
        if(probe.status=="Assessed" || probe.status=="OneSidedAssessed")
        {
            ASSERT_TRUE(probe.parameter_change);
            ASSERT_TRUE(probe.prediction_change);
            ASSERT_TRUE(probe.objective_change);
        }
        EXPECT_TRUE(probe.minus_seed_status=="Assessed" || probe.minus_seed_status=="BoundInfeasibleDirection" ||
            probe.minus_seed_status=="ZeroDisplacement" || probe.minus_seed_status=="LocalAcKktFailed" ||
            probe.minus_seed_status=="IndependentReplayFailed" || probe.minus_seed_status=="WidthStationarityFailed" ||
            probe.minus_seed_status=="SolverFailed" || probe.minus_seed_status=="BudgetExhausted");
        EXPECT_TRUE(probe.plus_seed_status=="Assessed" || probe.plus_seed_status=="BoundInfeasibleDirection" ||
            probe.plus_seed_status=="ZeroDisplacement" || probe.plus_seed_status=="LocalAcKktFailed" ||
            probe.plus_seed_status=="IndependentReplayFailed" || probe.plus_seed_status=="WidthStationarityFailed" ||
            probe.plus_seed_status=="SolverFailed" || probe.plus_seed_status=="BudgetExhausted");
    }
    EXPECT_TRUE(std::any_of(first.evidence.probes.begin(),first.evidence.probes.end(),[](const auto & probe) {
        return probe.minus_seed_status=="BoundInfeasibleDirection" ||
            probe.plus_seed_status=="BoundInfeasibleDirection";
    }));
    const auto second=n::AssessLocalParameterStability(domain,y,eta,beta,context,"two-atom-block");
    ASSERT_EQ(first.evidence.probes.size(),second.evidence.probes.size());
    for(std::size_t index=0;index<first.evidence.probes.size();++index)
    {
        EXPECT_EQ(first.evidence.probes[index].direction,second.evidence.probes[index].direction);
        EXPECT_EQ(first.evidence.probes[index].unit_direction,second.evidence.probes[index].unit_direction);
    }
}

TEST(JointParameterStabilityTest, WidthBoundUsesOneSidedInwardProbeAndRejectsZeroDisplacement)
{
    StabilityFixture fixture;
    auto context=fixture.context;
    context.eta_upper_bound=fixture.eta(0);
    const auto assessment=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        context,"bound-active/block-1");
    ASSERT_EQ(assessment.evidence.baseline_eligibility_status,"Eligible");
    ASSERT_EQ(assessment.evidence.probes.size(),1u);
    const auto & probe=assessment.evidence.probes.front();
    EXPECT_EQ(probe.status,"OneSidedAssessed");
    EXPECT_EQ(probe.minus_seed_status,"Assessed");
    EXPECT_EQ(probe.plus_seed_status,"ZeroDisplacement");
    EXPECT_EQ(probe.plus_seed_reason,"upper-bound-projection-produced-zero-displacement");
    ASSERT_TRUE(probe.plus_effective_displacement_inf);
    EXPECT_DOUBLE_EQ(*probe.plus_effective_displacement_inf,0.);
    EXPECT_EQ(assessment.evidence.perturbation_solves,1u);
}

TEST(JointParameterStabilityTest, IndependentReplayFailureCannotBeClassifiedStable)
{
    StabilityFixture fixture;
    auto context=fixture.context;
    context.eta_upper_bound=fixture.eta(0)-.1;
    const auto assessment=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        context,"infeasible-baseline/block-1");
    EXPECT_EQ(assessment.evidence.baseline_independent_replay_status,"Failed");
    EXPECT_EQ(assessment.evidence.baseline_eligibility_status,"Unavailable");
    EXPECT_TRUE(assessment.evidence.baseline_eligibility_reason=="baseline-independent-replay-failed" ||
        assessment.evidence.baseline_eligibility_reason=="baseline-width-practical-stationarity-failed" ||
        assessment.evidence.baseline_eligibility_reason=="baseline-width-bound-infeasible");
    EXPECT_EQ(assessment.evidence.status,"Unavailable");
    EXPECT_EQ(assessment.evidence.perturbation_solves,0u);
}

TEST(JointParameterStabilityTest, LocalAcSolveFailureCannotProduceFalseStableResult)
{
    StabilityFixture fixture;
    auto context=fixture.context;
    context.linear.rank_relative=1.;
    const auto assessment=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        context,"rank-rejected/block-1");
    EXPECT_EQ(assessment.evidence.baseline_eligibility_status,"Eligible");
    EXPECT_NE(assessment.evidence.status,"Stable");
    EXPECT_TRUE(std::any_of(assessment.evidence.probes.begin(),assessment.evidence.probes.end(),[](const auto & probe) {
        return probe.minus_seed_status=="LocalAcKktFailed" || probe.plus_seed_status=="LocalAcKktFailed" ||
            probe.minus_seed_status=="SolverFailed" || probe.plus_seed_status=="SolverFailed";
    }));
}

TEST(JointParameterStabilityTest, ComponentFitHonorsPerturbationClaimBudgetAndPreservesPoint)
{
    const rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    std::size_t claims{};
    const auto run=[&](bool assess_parameter_stability) {
        n::FixedNeighborPolicy policy; policy.core_atoms=1;
        policy.maximum_sweeps=n::ProductionFixedNeighborMaximumSweeps;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2);
        policy.assess_parameter_stability=assess_parameter_stability;
        if(assess_parameter_stability)
            policy.claim_parameter_stability_block=[&](int) {return ++claims<=1;};
        const n::Vector eta=n::Vector::Constant(2,std::log(.55));
        return n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    };
    const auto measured=run(true),baseline=run(false);
    ASSERT_TRUE(measured.state.eta.allFinite())<<measured.reason;
    ASSERT_TRUE(baseline.state.eta.allFinite())<<baseline.reason;
    ASSERT_TRUE(measured.parameter_stability);
    EXPECT_EQ(measured.state.eta,baseline.state.eta);
    EXPECT_EQ(measured.state.beta,baseline.state.beta);
    EXPECT_DOUBLE_EQ(measured.state.objective,baseline.state.objective);
    EXPECT_EQ(measured.search_converged,baseline.search_converged);
    EXPECT_EQ(measured.reason,baseline.reason);
    const auto & stability=*measured.parameter_stability;
    EXPECT_EQ(stability.blocks.size(),2u);
    EXPECT_LE(stability.perturbation_solves,2u);
    EXPECT_LE(claims,2u);
}
}
