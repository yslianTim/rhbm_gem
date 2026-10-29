#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/LocalConditioning.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "core/detail/joint_component/SplitDiagnostics.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointTestContext.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>
#include <utility>

namespace {
namespace c=rhbm_gem::core;
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

struct SplitMetricFixture
{
    n::Vector eta0{2},beta0{4},prediction0{3},eta_null{2},beta_null{4},prediction_null{3},column_norms{4};
    double scale{2.};
    SplitMetricFixture()
    {
        eta0<<0.,0.; beta0<<0.,0.,0.,0.; prediction0<<3.,4.,0.;
        eta_null=eta0; beta_null=beta0; prediction_null=prediction0;
        column_norms<<2.,1.,4.,3.;
    }
    n::SplitDiagnosticsAccumulator Accumulator(bool qualified=true,std::string reason={},std::size_t directions=1)
    {
        return {eta0,beta0,prediction0,eta_null,beta_null,prediction_null,column_norms,
            scale,qualified,std::move(reason),directions};
    }
};

c::JointState TargetState(const n::BlockCoordinateState & accepted)
{
    c::JointState state;
    state.ac.assign(accepted.beta.data(),accepted.beta.data()+accepted.beta.size());
    state.log_b.assign(accepted.eta.data(),accepted.eta.data()+accepted.eta.size());
    state.b.resize(static_cast<std::size_t>(accepted.eta.size()));
    for(Eigen::Index atom=0;atom<accepted.eta.size();++atom)
        state.b[static_cast<std::size_t>(atom)]=std::exp(accepted.eta(atom));
    state.objective=accepted.objective;
    return state;
}

TEST(JointSplitDiagnosticsTest, ExactFitNullIdentityHasZeroDriftAndSeedDifference)
{
    SplitMetricFixture fixture;
    auto diagnostics=fixture.Accumulator();
    diagnostics.AddEndpoint("weak-singular","weak-singular:-",fixture.eta0,fixture.beta0,
        fixture.prediction0,true);
    diagnostics.FinishDirection("weak-singular",true);
    const auto result=diagnostics.Finish();
    EXPECT_EQ(result.status,"Qualified");
    EXPECT_DOUBLE_EQ(*result.reoptimization_prediction_drift,0.);
    EXPECT_DOUBLE_EQ(*result.reoptimization_log_width_difference_inf,0.);
    EXPECT_DOUBLE_EQ(*result.reoptimization_scaled_ac_difference,0.);
    EXPECT_DOUBLE_EQ(*result.seed_prediction_robustness,0.);
    EXPECT_DOUBLE_EQ(*result.seed_log_width_difference_inf,0.);
    EXPECT_DOUBLE_EQ(*result.seed_scaled_ac_difference,0.);
    EXPECT_DOUBLE_EQ(*result.vector_decomposition_max_error,0.);
}

TEST(JointSplitDiagnosticsTest, ReoptimizationDriftAndSeedEndpointDifferenceUseDifferentReferences)
{
    SplitMetricFixture fixture;
    fixture.eta_null<<1e-3,-2e-3;
    fixture.beta_null<<1.,-2.,3.,-4.;
    fixture.prediction_null<<3.00000001,4.,0.;
    auto diagnostics=fixture.Accumulator();
    auto seed_eta=fixture.eta_null; seed_eta(0)+=2e-12;
    auto seed_beta=fixture.beta_null; seed_beta(2)+=1e-12;
    auto seed_prediction=fixture.prediction_null; seed_prediction(1)+=1e-12;
    diagnostics.AddEndpoint("common-mode","common-mode:+",seed_eta,seed_beta,seed_prediction,true);
    diagnostics.FinishDirection("common-mode",true);
    const auto result=diagnostics.Finish();
    ASSERT_EQ(result.status,"Qualified");
    EXPECT_NEAR(*result.reoptimization_prediction_drift,2e-9,1e-16);
    EXPECT_DOUBLE_EQ(*result.reoptimization_log_width_difference_inf,2e-3);
    EXPECT_DOUBLE_EQ(*result.reoptimization_scaled_ac_difference,6.);
    EXPECT_NEAR(*result.seed_prediction_robustness,2e-13,1e-16);
    EXPECT_NEAR(*result.seed_log_width_difference_inf,2e-12,1e-20);
    EXPECT_NEAR(*result.seed_scaled_ac_difference,2e-12,3e-16);
    EXPECT_EQ(result.seed_prediction_endpoint,"common-mode:+");
}

TEST(JointSplitDiagnosticsTest, T0SmallDriftRegressionKeepsSeedDifferenceSeparate)
{
    SplitMetricFixture fixture;
    fixture.prediction0.setZero();
    fixture.prediction_null.setZero();
    fixture.eta_null(0)=8.83e-9; fixture.beta_null(0)=3.01e-9;
    fixture.prediction_null(0)=3.2813e-9;
    auto diagnostics=fixture.Accumulator();
    auto seed_prediction=fixture.prediction_null; seed_prediction(0)+=5.47e-14;
    diagnostics.AddEndpoint("weak-singular","weak-singular:+",fixture.eta_null,fixture.beta_null,
        seed_prediction,true);
    diagnostics.FinishDirection("weak-singular",true);
    const auto result=diagnostics.Finish();
    EXPECT_EQ(result.status,"Qualified");
    EXPECT_NEAR(*result.reoptimization_prediction_drift,3.2813e-9,1e-20);
    EXPECT_NEAR(*result.seed_prediction_robustness,5.47e-14,1e-20);
    EXPECT_NEAR(*result.reoptimization_scaled_ac_difference,3.01e-9,1e-20);
}

TEST(JointSplitDiagnosticsTest, T4LargeDriftAndSmallSeedRegressionDoesNotChangeLegacyClassification)
{
    SplitMetricFixture fixture;
    fixture.prediction0.setZero();
    fixture.prediction_null.setZero();
    fixture.eta_null(0)=5.2135492445426479e-7;
    fixture.beta_null(0)=1.7762220265378486e-7;
    fixture.prediction_null(0)=1.9376222241792277e-7;
    auto diagnostics=fixture.Accumulator();
    auto seed_prediction=fixture.prediction_null; seed_prediction(0)+=1.20e-13;
    diagnostics.AddEndpoint("weak-singular","weak-singular:-",fixture.eta_null,fixture.beta_null,
        seed_prediction,true);
    diagnostics.FinishDirection("weak-singular",true);
    const auto result=diagnostics.Finish();
    EXPECT_NEAR(*result.reoptimization_prediction_drift,1.9376222241792277e-7,1e-20);
    EXPECT_NEAR(*result.seed_prediction_robustness,1.20e-13,1e-20);
    EXPECT_EQ(n::ClassifyLocalParameterStability(5.2135492445426479e-7,1.7762220265378486e-7,
        1.9376222241792277e-7,true,true).status,"Unstable");
}

TEST(JointSplitDiagnosticsTest, LargeComponentReferenceBlocksUseTheSameNormalizedMetrics)
{
    const auto check=[](std::size_t rows,double drift,double seed,double width,double ac) {
        n::Vector eta0= n::Vector::Zero(1),eta_null=n::Vector::Constant(1,width);
        n::Vector beta0=n::Vector::Zero(2),beta_null=n::Vector::Zero(2);
        n::Vector prediction0=n::Vector::Zero(static_cast<Eigen::Index>(rows));
        n::Vector prediction_null=prediction0; prediction_null(0)=drift;
        n::Vector column_norms=n::Vector::Constant(2,1.);
        beta_null(0)=ac;
        n::SplitDiagnosticsAccumulator diagnostics(eta0,beta0,prediction0,eta_null,beta_null,
            prediction_null,column_norms,1.,true,{},1);
        auto seed_prediction=prediction_null; seed_prediction(0)+=seed;
        diagnostics.AddEndpoint("weak-singular","weak-singular:+",eta_null,beta_null,seed_prediction,true);
        diagnostics.FinishDirection("weak-singular",true);
        const auto result=diagnostics.Finish();
        EXPECT_EQ(result.status,"Qualified")<<rows;
        EXPECT_NEAR(*result.reoptimization_prediction_drift,drift,1e-21)<<rows;
        EXPECT_NEAR(*result.seed_prediction_robustness,seed,1e-21)<<rows;
        EXPECT_DOUBLE_EQ(*result.reoptimization_log_width_difference_inf,width)<<rows;
        EXPECT_DOUBLE_EQ(*result.reoptimization_scaled_ac_difference,ac)<<rows;
    };
    check(300,1.9387e-7,3.84e-13,5.22e-7,5.62e-8);
    check(1000,1.9380e-7,5.40e-14,5.22e-7,3.07e-8);
}

TEST(JointSplitDiagnosticsTest, UnqualifiedNullReferencesNeverProduceSplitMetrics)
{
    SplitMetricFixture fixture;
    for(const auto & failure:std::array<std::string,2>{"null-reference-kkt-failed","null-reference-budget-exhausted"})
    {
        auto diagnostics=fixture.Accumulator(false,failure);
        diagnostics.AddEndpoint("weak-singular","weak-singular:+",fixture.eta0,fixture.beta0,
            fixture.prediction0,true);
        const auto result=diagnostics.Finish();
        EXPECT_EQ(result.status,"Unavailable");
        EXPECT_EQ(result.reason,failure);
        EXPECT_FALSE(result.reoptimization_prediction_drift);
        EXPECT_FALSE(result.seed_prediction_robustness);
    }
}

TEST(JointSplitDiagnosticsTest, ReplayFailureAndIncompleteDirectionCannotClaimSeedRobustness)
{
    SplitMetricFixture fixture;
    auto diagnostics=fixture.Accumulator();
    diagnostics.AddEndpoint("weak-singular","weak-singular:+",fixture.eta0,fixture.beta0,
        fixture.prediction0,false,"independent-replay-failed");
    diagnostics.FinishDirection("weak-singular",false,"independent-replay-failed");
    const auto result=diagnostics.Finish();
    EXPECT_EQ(result.status,"Incomplete");
    EXPECT_EQ(result.incomplete_direction_count,1u);
    EXPECT_FALSE(result.seed_prediction_robustness);
    EXPECT_EQ(result.directions.front().reason,"independent-replay-failed");
}

TEST(JointSplitDiagnosticsTest, ActiveWidthOneSidedProbeCanQualifyWithItsEndpointIdentity)
{
    SplitMetricFixture fixture;
    auto diagnostics=fixture.Accumulator(true,{},1);
    auto eta=fixture.eta0; eta(0)-=1e-3;
    auto prediction=fixture.prediction0; prediction(0)+=2e-8;
    diagnostics.AddEndpoint("weak-singular","weak-singular:-",eta,fixture.beta0,prediction,true);
    diagnostics.FinishDirection("weak-singular",true);
    const auto result=diagnostics.Finish();
    EXPECT_EQ(result.status,"Qualified");
    EXPECT_EQ(result.seed_prediction_endpoint,"weak-singular:-");
}

TEST(JointSplitDiagnosticsTest, ZeroEffectiveDisplacementLeavesDirectionIncomplete)
{
    SplitMetricFixture fixture;
    auto diagnostics=fixture.Accumulator();
    diagnostics.AddEndpoint("coordinate","coordinate:+",fixture.eta0,fixture.beta0,
        fixture.prediction0,false,"zero-effective-displacement");
    diagnostics.FinishDirection("coordinate",false,"zero-effective-displacement");
    const auto result=diagnostics.Finish();
    EXPECT_EQ(result.status,"Incomplete");
    EXPECT_EQ(result.reason,"incomplete-perturbation-directions");
    EXPECT_FALSE(result.seed_prediction_robustness);
}

TEST(JointSplitDiagnosticsTest, NonfiniteReferenceInputIsUnavailable)
{
    SplitMetricFixture fixture;
    fixture.prediction_null(0)=std::numeric_limits<double>::infinity();
    auto diagnostics=fixture.Accumulator();
    const auto result=diagnostics.Finish();
    EXPECT_EQ(result.status,"Unavailable");
    EXPECT_EQ(result.reason,"invalid-or-nonfinite-split-input");
    EXPECT_FALSE(result.reoptimization_prediction_drift);
}

TEST(JointSplitDiagnosticsTest, ScaledAcDifferenceUsesTheNullReferenceColumnNorms)
{
    SplitMetricFixture fixture;
    fixture.column_norms<<3.,4.,5.,6.; fixture.scale=2.;
    fixture.beta_null<<2.,-3.,1.,-1.;
    auto diagnostics=fixture.Accumulator();
    diagnostics.FinishDirection("coordinate",true);
    const auto result=diagnostics.Finish();
    EXPECT_DOUBLE_EQ(*result.reoptimization_scaled_ac_difference,6.);
}

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

TEST(JointParameterStabilityTest, NullIdentityAndProfiledReferencePreserveAnExactLocalFit)
{
    StabilityFixture fixture;
    const auto input_eta=fixture.eta,input_beta=fixture.beta;
    const auto n0=n::EvaluateState(fixture.domain,fixture.y,input_eta,input_beta,fixture.context,
        n::ProfileEvaluationRole::AcceptedEndpoint);
    const auto identity=n::EvaluateState(fixture.domain,fixture.y,input_eta,input_beta,fixture.context,
        n::ProfileEvaluationRole::AcceptedEndpoint);
    const auto n1=n::EvaluateProfile(fixture.domain,fixture.y,input_eta,&fixture.context,
        n::ProfileEvaluationRole::AcceptedEndpoint);
    auto strict_context=fixture.context;
    strict_context.profile_budget=200; strict_context.update_budget=100;
    const auto search=n::SearchProfile(fixture.domain,fixture.y,input_eta.array().exp(),strict_context);
    const auto reference=n::EvaluateProfile(fixture.domain,fixture.y,search.eta,&strict_context,
        n::ProfileEvaluationRole::AcceptedEndpoint);

    ASSERT_TRUE(n0.valid); ASSERT_TRUE(identity.valid); ASSERT_TRUE(n1.valid); ASSERT_TRUE(reference.valid);
    EXPECT_EQ(identity.eta,n0.eta);
    EXPECT_EQ(identity.beta,n0.beta);
    const auto prediction=[](const n::Evaluation & state) {return (state.x*state.beta).eval();};
    EXPECT_DOUBLE_EQ((prediction(identity)-prediction(n0)).norm(),0.);
    EXPECT_DOUBLE_EQ(identity.certificate.objective,n0.certificate.objective);
    EXPECT_LE((prediction(n1)-prediction(n0)).norm(),1e-14);
    EXPECT_LE((prediction(reference)-prediction(n1)).norm(),1e-14);
    EXPECT_LE((reference.eta-input_eta).lpNorm<Eigen::Infinity>(),1e-12);
    EXPECT_TRUE(search.initial.beta.isApprox(n1.beta,1e-14));
    EXPECT_TRUE(n::CheckIndependentReplay(fixture.domain,fixture.y,n0,fixture.context).passed);
    EXPECT_TRUE(n::CheckIndependentReplay(fixture.domain,fixture.y,n1,fixture.context).passed);
    EXPECT_TRUE(n::CheckIndependentReplay(fixture.domain,fixture.y,reference,strict_context).passed);
    EXPECT_EQ(fixture.eta,input_eta);
    EXPECT_EQ(fixture.beta,input_beta);
}

TEST(JointParameterStabilityTest, MinimalAndWorkspaceReuseCandidatesQualifyTheSameNullAndSeedEndpoints)
{
    StabilityFixture fixture;
    const auto legacy=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        fixture.context,"split-candidate/block-1");
    const auto minimal=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        fixture.context,"split-candidate/block-1",0,n::SplitDiagnosticsMode::Minimal);
    const auto reused=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        fixture.context,"split-candidate/block-1",0,n::SplitDiagnosticsMode::ReuseWorkspace);
    ASSERT_TRUE(minimal.split_diagnostics);
    ASSERT_TRUE(reused.split_diagnostics);
    EXPECT_EQ(minimal.split_diagnostics->null_reference_status,"Qualified")
        << minimal.split_diagnostics->null_reference_failure_reason;
    EXPECT_EQ(reused.split_diagnostics->null_reference_status,"Qualified")
        << reused.split_diagnostics->null_reference_failure_reason;
    EXPECT_EQ(minimal.split_diagnostics->status,"Qualified");
    EXPECT_EQ(reused.split_diagnostics->status,"Qualified");
    EXPECT_EQ(minimal.split_diagnostics->qualified_seed_count,2u);
    EXPECT_EQ(minimal.split_diagnostics->null_replay_checks,2u);
    EXPECT_LE(*minimal.split_diagnostics->high_budget_width_difference_inf,1e-12);
    EXPECT_LE(*minimal.split_diagnostics->high_budget_scaled_ac_difference,1e-10);
    EXPECT_LE(*minimal.split_diagnostics->high_budget_prediction_difference_relative,1e-10);
    EXPECT_EQ(minimal.evidence.status,legacy.evidence.status);
    EXPECT_EQ(reused.evidence.status,legacy.evidence.status);
    EXPECT_EQ(minimal.evidence.prediction_status,legacy.evidence.prediction_status);
    EXPECT_EQ(reused.evidence.prediction_status,legacy.evidence.prediction_status);
    EXPECT_EQ(minimal.evidence.probes.size(),legacy.evidence.probes.size());
    EXPECT_EQ(reused.evidence.probes.size(),legacy.evidence.probes.size());
    EXPECT_NEAR(*minimal.split_diagnostics->reoptimization_prediction_drift,
        *reused.split_diagnostics->reoptimization_prediction_drift,1e-15);
    EXPECT_NEAR(*minimal.split_diagnostics->seed_prediction_robustness,
        *reused.split_diagnostics->seed_prediction_robustness,1e-15);
    EXPECT_LE(*minimal.split_diagnostics->vector_decomposition_max_error,1e-14);
}

TEST(JointParameterStabilityTest, SplitCandidateBaselineFailureIsUnavailableAndDoesNotChangeLegacyStatus)
{
    StabilityFixture fixture;
    auto context=fixture.context;
    context.eta_upper_bound=fixture.eta(0)-.1;
    const auto legacy=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        context,"split-baseline/block-1");
    const auto candidate=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        context,"split-baseline/block-1",0,n::SplitDiagnosticsMode::Minimal);
    ASSERT_TRUE(candidate.split_diagnostics);
    EXPECT_EQ(candidate.split_diagnostics->null_reference_status,"Unavailable");
    EXPECT_EQ(candidate.split_diagnostics->status,"Unavailable");
    EXPECT_EQ(candidate.split_diagnostics->null_reference_failure_reason,"baseline-width-bound-infeasible");
    EXPECT_FALSE(candidate.split_diagnostics->seed_prediction_robustness);
    EXPECT_EQ(candidate.evidence.status,legacy.evidence.status);
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

TEST(JointParameterStabilityTest, ReclassifiesSameQualifiedMetricsAcrossThresholdGrid)
{
    constexpr double width_change=5.2135492445426479e-7;
    constexpr double ac_change=1.7762220265378486e-7;
    constexpr double prediction_change=1.9376222241792277e-7;
    const std::array<double,3> thresholds{1e-8,1e-7,1e-6};
    const std::array<std::array<std::string,3>,3> expected{{
        {{"Unstable","Unstable","Sensitive"}},
        {{"Unstable","Unstable","Sensitive"}},
        {{"Unstable","Unstable","Stable"}}
    }};
    for(std::size_t p=0;p<thresholds.size();++p)
        for(std::size_t y=0;y<thresholds.size();++y)
        {
            const auto result=n::ClassifyLocalParameterStability(width_change,ac_change,prediction_change,
                true,true,thresholds[p],thresholds[y]);
            EXPECT_EQ(result.status,expected[p][y]);
        }
    for(const auto parameter_threshold:thresholds)
        for(const auto prediction_threshold:thresholds)
        {
            const auto result=n::ClassifyLocalParameterStability(width_change,ac_change,prediction_change,
                false,false,parameter_threshold,prediction_threshold);
            EXPECT_EQ(result.status,"Unavailable");
        }
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
    const auto candidate=n::AssessLocalParameterStability(domain,y,eta,beta,context,"two-atom-block",0,
        n::SplitDiagnosticsMode::Minimal);
    ASSERT_TRUE(candidate.split_diagnostics);
    EXPECT_EQ(candidate.evidence.status,first.evidence.status);
    EXPECT_EQ(candidate.evidence.prediction_status,first.evidence.prediction_status);
    EXPECT_NE(candidate.split_diagnostics->status,"Unavailable")
        << candidate.split_diagnostics->reason;
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

    const auto split=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        context,"bound-active/block-1",0,n::SplitDiagnosticsMode::Minimal);
    ASSERT_TRUE(split.split_diagnostics);
    EXPECT_EQ(split.evidence.status,assessment.evidence.status);
    EXPECT_EQ(split.evidence.prediction_status,assessment.evidence.prediction_status);
    EXPECT_EQ(split.split_diagnostics->null_reference_status,"Qualified")
        << split.split_diagnostics->null_reference_failure_reason;
    EXPECT_EQ(split.split_diagnostics->status,"Qualified");
    ASSERT_EQ(split.split_diagnostics->directions.size(),1u);
    EXPECT_EQ(split.split_diagnostics->qualified_seed_count,1u);
    EXPECT_TRUE(split.split_diagnostics->directions.front().prediction_endpoint.ends_with(":-"));
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
    const auto candidate=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        context,"rank-rejected/block-1",0,n::SplitDiagnosticsMode::Minimal);
    ASSERT_TRUE(candidate.split_diagnostics);
    EXPECT_EQ(candidate.evidence.status,assessment.evidence.status);
    EXPECT_EQ(candidate.evidence.prediction_status,assessment.evidence.prediction_status);
    EXPECT_NE(candidate.split_diagnostics->status,"Qualified");
    EXPECT_FALSE(candidate.split_diagnostics->seed_prediction_robustness);
}

TEST(JointParameterStabilityTest, ConditionalReprofileStatusDoesNotOverrideWidthReplayFailure)
{
    const rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",30));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    n::FixedNeighborPolicy policy;
    policy.core_atoms=12;
    policy.maximum_sweeps=n::ProductionFixedNeighborMaximumSweeps;
    policy.convergence_criteria={1e-7,1e-5};
    policy.eta_upper_bound=std::log(2.2);
    policy.assess_local_conditioning=true;
    policy.assess_parameter_stability=true;
    policy.claim_parameter_stability_block=[](int) {return true;};
    const n::Vector eta=n::Vector::Constant(30,std::log(.55));

    const auto result=n::SearchFixedNeighborComponent(
        data.domain,data.layout,data.y,eta,data.context,policy);
    ASSERT_TRUE(result.search_converged)<<result.reason;
    ASSERT_TRUE(result.parameter_stability);
    const auto found=std::find_if(result.parameter_stability->blocks.begin(),
        result.parameter_stability->blocks.end(),[](const auto & block) {
            return block.block_id=="component/block-1";
        });
    ASSERT_NE(found,result.parameter_stability->blocks.end());
    EXPECT_EQ(found->baseline_eligibility_status,"Unavailable");
    EXPECT_EQ(found->baseline_eligibility_reason,"baseline-width-practical-stationarity-failed");
    EXPECT_EQ(found->conditional_reprofile_status,"StrictlyQualified");
    ASSERT_TRUE(found->conditional_reprofile_local_ac_kkt);
    ASSERT_TRUE(found->conditional_reprofile_projected_width_kkt);
    EXPECT_LE(*found->conditional_reprofile_local_ac_kkt,1e-10);
    EXPECT_GT(*found->conditional_reprofile_projected_width_kkt,1e-9);
    EXPECT_EQ(found->conditional_reprofile_independent_replay_status,"Failed");
    EXPECT_EQ(found->conditional_reprofile_independent_replay_reason,
        "primary-bounded-width-kkt-failed");
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

TEST(JointParameterStabilityTest, SplitCandidatesPreserveTheCompleteAcceptedFitAndLegacyEvidence)
{
    const rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto run=[&](n::SplitDiagnosticsMode mode) {
        n::FixedNeighborPolicy policy;
        policy.core_atoms=1;
        policy.maximum_sweeps=n::ProductionFixedNeighborMaximumSweeps;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2);
        policy.assess_local_conditioning=true;
        policy.assess_parameter_stability=true;
        policy.claim_parameter_stability_block=[](int) {return true;};
        policy.split_diagnostics_mode=mode;
        const auto eta=n::Vector::Constant(2,std::log(.55));
        return n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    };
    const auto legacy=run(n::SplitDiagnosticsMode::Disabled);
    const auto minimal=run(n::SplitDiagnosticsMode::Minimal);
    const auto reused=run(n::SplitDiagnosticsMode::ReuseWorkspace);
    ASSERT_TRUE(legacy.certification && minimal.certification && reused.certification);
    EXPECT_EQ(legacy.state.eta,minimal.state.eta);
    EXPECT_EQ(legacy.state.eta,reused.state.eta);
    EXPECT_EQ(legacy.state.beta,minimal.state.beta);
    EXPECT_EQ(legacy.state.beta,reused.state.beta);
    EXPECT_EQ(legacy.state.prediction,minimal.state.prediction);
    EXPECT_EQ(legacy.state.prediction,reused.state.prediction);
    EXPECT_EQ(legacy.state.residual,minimal.state.residual);
    EXPECT_EQ(legacy.state.residual,reused.state.residual);
    EXPECT_DOUBLE_EQ(legacy.state.objective,minimal.state.objective);
    EXPECT_DOUBLE_EQ(legacy.state.objective,reused.state.objective);
    EXPECT_EQ(legacy.search_converged,minimal.search_converged);
    EXPECT_EQ(legacy.search_converged,reused.search_converged);
    EXPECT_EQ(legacy.sweep_count,minimal.sweep_count);
    EXPECT_EQ(legacy.sweep_count,reused.sweep_count);
    ASSERT_TRUE(legacy.parameter_stability && minimal.parameter_stability && reused.parameter_stability);
    ASSERT_EQ(legacy.parameter_stability->blocks.size(),minimal.parameter_stability->blocks.size());
    ASSERT_EQ(legacy.parameter_stability->blocks.size(),reused.parameter_stability->blocks.size());
    ASSERT_EQ(minimal.split_diagnostics.size(),reused.split_diagnostics.size());
    for(std::size_t index=0;index<legacy.parameter_stability->blocks.size();++index)
    {
        const auto & baseline=legacy.parameter_stability->blocks[index];
        const auto & b1=minimal.parameter_stability->blocks[index];
        const auto & b2=reused.parameter_stability->blocks[index];
        EXPECT_EQ(baseline.status,b1.status);
        EXPECT_EQ(baseline.status,b2.status);
        EXPECT_EQ(baseline.prediction_status,b1.prediction_status);
        EXPECT_EQ(baseline.prediction_status,b2.prediction_status);
        EXPECT_EQ(baseline.reason,b1.reason);
        EXPECT_EQ(baseline.reason,b2.reason);
    }
    for(std::size_t index=0;index<minimal.split_diagnostics.size();++index)
    {
        const auto & b1=minimal.split_diagnostics[index];
        const auto & b2=reused.split_diagnostics[index];
        EXPECT_EQ(b1.block_id,b2.block_id);
        EXPECT_EQ(b1.diagnostics.status,b2.diagnostics.status);
        EXPECT_EQ(b1.diagnostics.null_reference_status,b2.diagnostics.null_reference_status);
        EXPECT_EQ(b1.diagnostics.null_reference_failure_reason,b2.diagnostics.null_reference_failure_reason);
        EXPECT_EQ(b1.diagnostics.qualified_seed_count,b2.diagnostics.qualified_seed_count);
        EXPECT_EQ(b1.diagnostics.seed_prediction_endpoint,b2.diagnostics.seed_prediction_endpoint);
        EXPECT_EQ(b1.diagnostics.seed_log_width_endpoint,b2.diagnostics.seed_log_width_endpoint);
        EXPECT_EQ(b1.diagnostics.seed_scaled_ac_endpoint,b2.diagnostics.seed_scaled_ac_endpoint);
    }
}

TEST(JointParameterStabilityTest, NoiseOnlyBudgetExhaustionNeverProducesSeedRobustness)
{
    const c::JointProblem problem(second_stage_test::SyntheticJointWorkload("noise-only",3));
    const auto & data=c::JointProblemAccess::Get(problem);
    const auto run=[&](n::SplitDiagnosticsMode mode) {
        n::FixedNeighborPolicy policy;
        policy.core_atoms=12;
        policy.maximum_sweeps=n::ProductionFixedNeighborMaximumSweeps;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2);
        policy.assess_local_conditioning=true;
        policy.assess_parameter_stability=true;
        policy.claim_parameter_stability_block=[](int) {return true;};
        policy.split_diagnostics_mode=mode;
        const n::Vector eta=n::Vector::Constant(3,std::log(.55));
        return n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    };
    const auto legacy=run(n::SplitDiagnosticsMode::Disabled);
    const auto minimal=run(n::SplitDiagnosticsMode::Minimal);
    const auto reused=run(n::SplitDiagnosticsMode::ReuseWorkspace);
    ASSERT_TRUE(legacy.certification && minimal.certification && reused.certification);
    ASSERT_TRUE(legacy.certification->qualified() && minimal.certification->qualified() && reused.certification->qualified());
    EXPECT_EQ(legacy.state.eta,minimal.state.eta);
    EXPECT_EQ(legacy.state.eta,reused.state.eta);
    EXPECT_EQ(legacy.state.beta,minimal.state.beta);
    EXPECT_EQ(legacy.state.beta,reused.state.beta);
    EXPECT_EQ(legacy.state.prediction,minimal.state.prediction);
    EXPECT_EQ(legacy.state.prediction,reused.state.prediction);
    EXPECT_EQ(legacy.state.residual,minimal.state.residual);
    EXPECT_EQ(legacy.state.residual,reused.state.residual);
    EXPECT_DOUBLE_EQ(legacy.state.objective,minimal.state.objective);
    EXPECT_DOUBLE_EQ(legacy.state.objective,reused.state.objective);
    ASSERT_TRUE(legacy.parameter_stability && minimal.parameter_stability && reused.parameter_stability);
    ASSERT_EQ(legacy.parameter_stability->blocks.size(),minimal.parameter_stability->blocks.size());
    ASSERT_EQ(legacy.parameter_stability->blocks.size(),reused.parameter_stability->blocks.size());
    for(std::size_t index=0;index<legacy.parameter_stability->blocks.size();++index)
    {
        EXPECT_EQ(legacy.parameter_stability->blocks[index].status,minimal.parameter_stability->blocks[index].status);
        EXPECT_EQ(legacy.parameter_stability->blocks[index].status,reused.parameter_stability->blocks[index].status);
        EXPECT_EQ(legacy.parameter_stability->blocks[index].prediction_status,
            minimal.parameter_stability->blocks[index].prediction_status);
        EXPECT_EQ(legacy.parameter_stability->blocks[index].prediction_status,
            reused.parameter_stability->blocks[index].prediction_status);
    }
    ASSERT_EQ(minimal.split_diagnostics.size(),1u);
    ASSERT_EQ(reused.split_diagnostics.size(),1u);
    const auto & b1=minimal.split_diagnostics.front().diagnostics;
    const auto & b2=reused.split_diagnostics.front().diagnostics;
    EXPECT_EQ(b1.status,"Unavailable");
    EXPECT_EQ(b1.null_reference_status,"Unavailable");
    EXPECT_EQ(b1.null_reference_failure_reason,"null-reference-budget-exhausted:accepted-update-budget");
    EXPECT_EQ(b1.reason,b2.reason);
    EXPECT_EQ(b1.null_reference_status,b2.null_reference_status);
    EXPECT_EQ(b1.qualified_seed_count,0u);
    EXPECT_FALSE(b1.reoptimization_prediction_drift);
    EXPECT_FALSE(b1.seed_prediction_robustness);
    EXPECT_FALSE(b2.seed_prediction_robustness);
    EXPECT_GT(b1.null_profile_evaluations,0u);
    const auto evidence=[&](const auto & result) {
        return n::AssessTargetPointState(problem,data.layout,TargetState(result.state),
            rhbm_gem::JointCheckStatus::Passed,result.local_conditioning,
            result.parameter_stability->blocks,rhbm_gem::JointEvidenceScope::AssembledGlobal,
            data.input->observations.size());
    };
    const auto legacy_evidence=evidence(legacy),minimal_evidence=evidence(minimal),reused_evidence=evidence(reused);
    ASSERT_EQ(legacy_evidence.point_atoms.size(),minimal_evidence.point_atoms.size());
    ASSERT_EQ(legacy_evidence.point_atoms.size(),reused_evidence.point_atoms.size());
    for(std::size_t atom=0;atom<legacy_evidence.point_atoms.size();++atom)
    {
        EXPECT_EQ(legacy_evidence.point_atoms[atom].parameter_stability,
            minimal_evidence.point_atoms[atom].parameter_stability);
        EXPECT_EQ(legacy_evidence.point_atoms[atom].point_quality,minimal_evidence.point_atoms[atom].point_quality);
        EXPECT_EQ(legacy_evidence.point_atoms[atom].point_quality,reused_evidence.point_atoms[atom].point_quality);
    }
}

TEST(JointParameterStabilityTest, NearSingularDesignKeepsConditioningWarningWithSplitEvidence)
{
    c::JointProblemInput input;
    input.atom_ids={"a","b"};
    input.row_ids={"0","1","2","3","4"};
    input.observations.assign(5,0.0);
    input.selection_domain.emplace();
    input.selection_domain->target_indices={0,1};
    input.support={{{0,.1},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,.1+1e-8},{1,.5-1e-8},{2,1.0+1e-8},{3,1.5-1e-8},{4,2.0+1e-8}}};
    n::Vector eta=n::Vector::Constant(2,std::log(.7)),beta(4);
    beta<<1.0,.2,1.4,-.1;
    for(std::size_t row=0;row<input.observations.size();++row)
        for(std::size_t atom=0;atom<input.support.size();++atom)
        {
            const auto kernel=n::EvaluateKernel(input.support[atom][row].squared_distance,.7,2.5);
            input.observations[row]+=beta(2*static_cast<Eigen::Index>(atom))*kernel.gaussian+
                beta(2*static_cast<Eigen::Index>(atom)+1)*kernel.charge;
        }
    const auto original_rows=input.row_ids.size();
    const c::JointProblem problem(std::move(input));
    const auto & data=c::JointProblemAccess::Get(problem);
    auto context=data.context;
    context.eta_upper_bound=std::log(2.0);
    auto conditioning=n::AssessFixedNeighborLocalConditioning(data.domain,data.y,eta,beta,context,
        "near-singular/block-1");
    ASSERT_EQ(conditioning.evidence.conditioning_classification,"IllConditioned");
    ASSERT_TRUE(conditioning.evidence.free_ac_condition);
    EXPECT_GT(*conditioning.evidence.free_ac_condition,1e7);
    conditioning.evidence.core_atoms={0,1};
    const auto legacy=n::AssessLocalParameterStability(data.domain,data.y,eta,beta,context,
        "near-singular/block-1");
    const auto candidate=n::AssessLocalParameterStability(data.domain,data.y,eta,beta,context,
        "near-singular/block-1",0,n::SplitDiagnosticsMode::Minimal);
    ASSERT_TRUE(candidate.split_diagnostics);
    EXPECT_EQ(candidate.evidence.baseline_eligibility_status,"Eligible")
        << candidate.evidence.baseline_eligibility_reason;
    EXPECT_EQ(candidate.evidence.status,legacy.evidence.status);
    EXPECT_EQ(candidate.evidence.prediction_status,legacy.evidence.prediction_status);
    EXPECT_EQ(candidate.split_diagnostics->null_reference_status,"Qualified")
        << candidate.split_diagnostics->null_reference_failure_reason;
    EXPECT_TRUE(candidate.split_diagnostics->status=="Qualified" ||
        candidate.split_diagnostics->status=="Incomplete");
    if(candidate.split_diagnostics->status=="Qualified")
    {
        ASSERT_TRUE(candidate.split_diagnostics->reoptimization_prediction_drift);
        ASSERT_TRUE(candidate.split_diagnostics->seed_prediction_robustness);
        EXPECT_LE(*candidate.split_diagnostics->reoptimization_prediction_drift,1e-8);
        EXPECT_LE(*candidate.split_diagnostics->seed_prediction_robustness,1e-8);
    }
    c::JointState state;
    state.ac.assign(beta.data(),beta.data()+beta.size());
    state.log_b.assign(eta.data(),eta.data()+eta.size());
    state.b.resize(static_cast<std::size_t>(eta.size()));
    for(Eigen::Index atom=0;atom<eta.size();++atom)
        state.b[static_cast<std::size_t>(atom)]=std::exp(eta(atom));
    state.objective=0.;
    auto legacy_stability=legacy.evidence;
    auto candidate_stability=candidate.evidence;
    legacy_stability.core_atoms={0,1};
    candidate_stability.core_atoms={0,1};
    const auto legacy_point=n::AssessTargetPointState(problem,data.layout,state,
        rhbm_gem::JointCheckStatus::Passed,{conditioning.evidence},{legacy_stability},
        rhbm_gem::JointEvidenceScope::AssembledGlobal,original_rows);
    const auto candidate_point=n::AssessTargetPointState(problem,data.layout,state,
        rhbm_gem::JointCheckStatus::Passed,{conditioning.evidence},{candidate_stability},
        rhbm_gem::JointEvidenceScope::AssembledGlobal,original_rows);
    ASSERT_EQ(legacy_point.point_atoms.size(),candidate_point.point_atoms.size());
    for(std::size_t atom=0;atom<legacy_point.point_atoms.size();++atom)
    {
        EXPECT_EQ(legacy_point.point_atoms[atom].parameter_stability,
            candidate_point.point_atoms[atom].parameter_stability);
        EXPECT_EQ(legacy_point.point_atoms[atom].point_quality,candidate_point.point_atoms[atom].point_quality);
        EXPECT_EQ(candidate_point.point_atoms[atom].point_quality,"Warning");
    }
}

TEST(JointParameterStabilityTest, QualifiedMultiStartEndpointsAgreeOnExactSyntheticFixture)
{
    StabilityFixture fixture;
    const n::Vector low_width=n::Vector::Constant(1,.4),high_width=n::Vector::Constant(1,1.1);
    const auto low=n::SearchProfile(fixture.domain,fixture.y,low_width,fixture.context);
    const auto high=n::SearchProfile(fixture.domain,fixture.y,high_width,fixture.context);
    ASSERT_EQ(low.eta.size(),1)<<low.stop_reason;
    ASSERT_EQ(high.eta.size(),1)<<high.stop_reason;
    const auto low_endpoint=n::EvaluateProfile(fixture.domain,fixture.y,low.eta,&fixture.context);
    const auto high_endpoint=n::EvaluateProfile(fixture.domain,fixture.y,high.eta,&fixture.context);
    ASSERT_TRUE(low_endpoint.valid && high_endpoint.valid)<<low.stop_reason<<'/'<<high.stop_reason;
    const auto low_replay=n::CheckIndependentReplay(fixture.domain,fixture.y,low_endpoint,fixture.context,
        nullptr,nullptr,n::StrictJointConvergenceCriteria);
    const auto high_replay=n::CheckIndependentReplay(fixture.domain,fixture.y,high_endpoint,fixture.context,
        nullptr,nullptr,n::StrictJointConvergenceCriteria);
    EXPECT_TRUE(low_replay.passed)<<low_replay.reason;
    EXPECT_TRUE(high_replay.passed)<<high_replay.reason;
    EXPECT_LE((low.eta-high.eta).lpNorm<Eigen::Infinity>(),1e-8);
    EXPECT_LE((low_endpoint.x*low_endpoint.beta-high_endpoint.x*high_endpoint.beta).norm(),1e-12);
    const auto split=n::AssessLocalParameterStability(fixture.domain,fixture.y,fixture.eta,fixture.beta,
        fixture.context,"multi-start/block-1",0,n::SplitDiagnosticsMode::Minimal);
    ASSERT_TRUE(split.split_diagnostics);
    EXPECT_EQ(split.split_diagnostics->null_reference_status,"Qualified")
        << split.split_diagnostics->null_reference_failure_reason;
    EXPECT_EQ(split.split_diagnostics->status,"Qualified");
}

TEST(JointParameterStabilityTest, SplitMetricsQualifyAcrossTruthWidthsAndAcAmplitudes)
{
    for(const double truth_width:std::array<double,3>{.4,.8,1.0})
        for(const double amplitude:std::array<double,3>{.01,1.0,10.0})
        {
            n::Domain domain(9,{MakeSupport()});
            n::Vector eta(1),beta(2),y(9);
            eta(0)=std::log(truth_width);
            beta<<1.35*amplitude,.24*amplitude;
            for(Eigen::Index row=0;row<y.size();++row)
            {
                const double square=domain.atoms[0][static_cast<std::size_t>(row)].square;
                const auto kernel=n::EvaluateKernel(square,truth_width,2.5);
                y(row)=beta(0)*kernel.gaussian+beta(1)*kernel.charge;
            }
            auto context=tc::MakeContext(y,1);
            context.eta_upper_bound=std::log(2.0);
            const auto legacy=n::AssessLocalParameterStability(domain,y,eta,beta,context,"scale-block");
            const auto candidate=n::AssessLocalParameterStability(domain,y,eta,beta,context,"scale-block",0,
                n::SplitDiagnosticsMode::Minimal);
            ASSERT_EQ(candidate.evidence.baseline_eligibility_status,"Eligible")
                <<"width="<<truth_width<<" amplitude="<<amplitude<<' '<<candidate.evidence.baseline_eligibility_reason;
            ASSERT_TRUE(candidate.split_diagnostics);
            EXPECT_EQ(candidate.evidence.status,legacy.evidence.status);
            EXPECT_EQ(candidate.evidence.prediction_status,legacy.evidence.prediction_status);
            EXPECT_EQ(candidate.split_diagnostics->status,"Qualified")
                <<"width="<<truth_width<<" amplitude="<<amplitude<<' '<<candidate.split_diagnostics->reason;
            ASSERT_TRUE(candidate.split_diagnostics->reoptimization_prediction_drift);
            ASSERT_TRUE(candidate.split_diagnostics->seed_prediction_robustness);
            ASSERT_TRUE(candidate.split_diagnostics->reoptimization_scaled_ac_difference);
            ASSERT_TRUE(candidate.split_diagnostics->seed_scaled_ac_difference);
            EXPECT_TRUE(std::isfinite(*candidate.split_diagnostics->reoptimization_prediction_drift));
            EXPECT_TRUE(std::isfinite(*candidate.split_diagnostics->seed_prediction_robustness));
            EXPECT_TRUE(std::isfinite(*candidate.split_diagnostics->reoptimization_scaled_ac_difference));
            EXPECT_TRUE(std::isfinite(*candidate.split_diagnostics->seed_scaled_ac_difference));
        }
}

TEST(JointParameterStabilityTest, RankDeficientDesignCannotProduceSplitMetricsOrGoodPointQuality)
{
    const std::vector<n::Support> support{{{0,.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}};
    n::Domain domain(5,{support,support});
    n::Vector y=n::Vector::Zero(5),eta=n::Vector::Constant(2,std::log(.7)),beta=n::Vector::Zero(4);
    auto context=tc::MakeContext(y,2);
    context.eta_upper_bound=std::log(2.0);
    const auto conditioning=n::AssessFixedNeighborLocalConditioning(domain,y,eta,beta,context,
        "rank-deficient/block-1");
    EXPECT_EQ(conditioning.evidence.conditioning_classification,"Singular");
    ASSERT_TRUE(conditioning.evidence.free_ac_rank);
    EXPECT_LT(*conditioning.evidence.free_ac_rank,conditioning.evidence.free_ac_columns);
    const auto legacy=n::AssessLocalParameterStability(domain,y,eta,beta,context,"rank-deficient/block-1");
    const auto candidate=n::AssessLocalParameterStability(domain,y,eta,beta,context,
        "rank-deficient/block-1",0,n::SplitDiagnosticsMode::Minimal);
    ASSERT_TRUE(candidate.split_diagnostics);
    EXPECT_EQ(candidate.evidence.status,legacy.evidence.status);
    EXPECT_EQ(candidate.evidence.prediction_status,legacy.evidence.prediction_status);
    EXPECT_EQ(candidate.split_diagnostics->status,"Unavailable");
    EXPECT_EQ(candidate.split_diagnostics->null_reference_failure_reason,
        "null-reference-numerical-validity-failed");
    EXPECT_FALSE(candidate.split_diagnostics->seed_prediction_robustness);

    rhbm_gem::JointBlockConditioning warning=conditioning.evidence;
    warning.core_atoms={0,1};
    auto stability=candidate.evidence;
    stability.core_atoms={0,1};
    c::JointState state=TargetState({eta,beta,n::Vector::Zero(5),y,0.});
    const c::JointProblemInput input=[&] {
        c::JointProblemInput problem_input;
        problem_input.atom_ids={"a","b"};
        for(std::size_t row=0;row<5;++row) problem_input.row_ids.push_back(std::to_string(row));
        problem_input.observations.assign(5,0.0);
        problem_input.support.resize(2);
        for(std::size_t atom=0;atom<2;++atom)
            for(const auto & entry:support)
                problem_input.support[atom].push_back({static_cast<std::size_t>(entry.row),entry.square});
        problem_input.selection_domain.emplace();
        problem_input.selection_domain->target_indices={0,1};
        return problem_input;
    }();
    const c::JointProblem problem(input);
    const auto & data=c::JointProblemAccess::Get(problem);
    const auto point=n::AssessTargetPointState(problem,data.layout,state,rhbm_gem::JointCheckStatus::Passed,
        {warning},{stability},rhbm_gem::JointEvidenceScope::AssembledGlobal,5);
    ASSERT_EQ(point.point_atoms.size(),2u);
    for(const auto & atom:point.point_atoms)
    {
        EXPECT_EQ(atom.local_conditioning,"Warning");
        EXPECT_EQ(atom.point_quality,"Warning");
    }
}

TEST(JointParameterStabilityTest, SplitCandidatePreservesSyntheticTruthAccuracyAcrossNoiseScales)
{
    for(const double noise_multiplier:std::array<double,2>{1.0,10.0})
    {
        auto input=second_stage_test::SyntheticJointWorkload("chain",3);
        std::vector<double> truth(input.observations.size());
        for(std::size_t atom=0;atom<input.support.size();++atom)
            for(const auto & entry:input.support[atom])
            {
                const auto kernel=n::EvaluateKernel(entry.squared_distance,.5,2.5);
                truth[entry.row]+=2.*kernel.gaussian+.2*kernel.charge;
            }
        for(std::size_t row=0;row<input.observations.size();++row)
            input.observations[row]=truth[row]+noise_multiplier*(input.observations[row]-truth[row]);
        const c::JointProblem problem(std::move(input));
        const auto & data=c::JointProblemAccess::Get(problem);
        const auto run=[&](n::SplitDiagnosticsMode mode) {
            n::FixedNeighborPolicy policy;
            policy.core_atoms=12;
            policy.maximum_sweeps=n::ProductionFixedNeighborMaximumSweeps;
            policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
            policy.eta_upper_bound=std::log(2.2);
            policy.assess_local_conditioning=true;
            policy.assess_parameter_stability=true;
            policy.claim_parameter_stability_block=[](int) {return true;};
            policy.split_diagnostics_mode=mode;
            const n::Vector eta=n::Vector::Constant(3,std::log(.55));
            return n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
        };
        const auto legacy=run(n::SplitDiagnosticsMode::Disabled);
        const auto candidate=run(n::SplitDiagnosticsMode::Minimal);
        ASSERT_TRUE(legacy.certification && candidate.certification);
        ASSERT_TRUE(legacy.certification->qualified() && candidate.certification->qualified());
        EXPECT_EQ(legacy.state.eta,candidate.state.eta);
        EXPECT_EQ(legacy.state.beta,candidate.state.beta);
        EXPECT_EQ(legacy.state.prediction,candidate.state.prediction);
        EXPECT_EQ(legacy.state.residual,candidate.state.residual);
        EXPECT_DOUBLE_EQ(legacy.state.objective,candidate.state.objective);
        ASSERT_TRUE(legacy.parameter_stability && candidate.parameter_stability);
        ASSERT_EQ(legacy.parameter_stability->blocks.size(),candidate.parameter_stability->blocks.size());
        for(std::size_t block=0;block<legacy.parameter_stability->blocks.size();++block)
            EXPECT_EQ(legacy.parameter_stability->blocks[block].status,
                candidate.parameter_stability->blocks[block].status);
        double max_a_error{},max_c_error{},max_b_error{};
        for(Eigen::Index atom=0;atom<candidate.state.eta.size();++atom)
        {
            max_a_error=std::max(max_a_error,std::abs(candidate.state.beta(2*atom)-2.));
            max_c_error=std::max(max_c_error,std::abs(candidate.state.beta(2*atom+1)-.2));
            max_b_error=std::max(max_b_error,std::abs(std::exp(candidate.state.eta(atom))-.5));
        }
        const Eigen::Map<const n::Vector> truth_vector(truth.data(),static_cast<Eigen::Index>(truth.size()));
        const double prediction_rmse=(candidate.state.prediction-truth_vector).norm()/
            std::sqrt(static_cast<double>(truth.size()));
        EXPECT_TRUE(std::isfinite(prediction_rmse));
        EXPECT_LE(max_a_error,1e-2)<<"noise_multiplier="<<noise_multiplier;
        EXPECT_LE(max_c_error,1e-2)<<"noise_multiplier="<<noise_multiplier;
        EXPECT_LE(max_b_error,1e-2)<<"noise_multiplier="<<noise_multiplier;
        EXPECT_TRUE(candidate.final_sweep.has_value());
        EXPECT_LE(candidate.final_sweep->global_ac_kkt,
            n::PracticalJointConvergenceCriteria.global_ac_kkt);
        EXPECT_LE(candidate.final_sweep->stationarity_width_gradient_inf_norm,
            n::PracticalJointConvergenceCriteria.global_width_gradient);
    }
}

TEST(JointParameterStabilityTest, SplitMetricValuesMatchAcrossMinimalAndWorkspaceReuseOnThirtyTargets)
{
    const rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",30));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    const auto run=[&](n::SplitDiagnosticsMode mode) {
        n::FixedNeighborPolicy policy;
        policy.core_atoms=12;
        policy.maximum_sweeps=n::ProductionFixedNeighborMaximumSweeps;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(2.2);
        policy.assess_local_conditioning=true;
        policy.assess_parameter_stability=true;
        policy.claim_parameter_stability_block=[](int) {return true;};
        policy.split_diagnostics_mode=mode;
        const auto eta=n::Vector::Constant(30,std::log(.55));
        return n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    };
    const auto minimal=run(n::SplitDiagnosticsMode::Minimal);
    const auto reused=run(n::SplitDiagnosticsMode::ReuseWorkspace);
    ASSERT_TRUE(minimal.certification && reused.certification);
    EXPECT_EQ(minimal.state.eta,reused.state.eta);
    EXPECT_EQ(minimal.state.beta,reused.state.beta);
    EXPECT_EQ(minimal.state.prediction,reused.state.prediction);
    EXPECT_EQ(minimal.state.residual,reused.state.residual);
    EXPECT_DOUBLE_EQ(minimal.state.objective,reused.state.objective);
    ASSERT_TRUE(minimal.parameter_stability && reused.parameter_stability);
    ASSERT_EQ(minimal.parameter_stability->blocks.size(),reused.parameter_stability->blocks.size());
    ASSERT_EQ(minimal.split_diagnostics.size(),reused.split_diagnostics.size());
    const auto expect_same=[](const std::optional<double> & left,const std::optional<double> & right) {
        ASSERT_EQ(left.has_value(),right.has_value());
        if(left) EXPECT_DOUBLE_EQ(*left,*right);
    };
    for(std::size_t index=0;index<minimal.parameter_stability->blocks.size();++index)
    {
        EXPECT_EQ(minimal.parameter_stability->blocks[index].status,reused.parameter_stability->blocks[index].status);
        EXPECT_EQ(minimal.parameter_stability->blocks[index].prediction_status,
            reused.parameter_stability->blocks[index].prediction_status);
    }
    for(std::size_t index=0;index<minimal.split_diagnostics.size();++index)
    {
        const auto & left=minimal.split_diagnostics[index].diagnostics;
        const auto & right=reused.split_diagnostics[index].diagnostics;
        EXPECT_EQ(minimal.split_diagnostics[index].block_id,reused.split_diagnostics[index].block_id);
        EXPECT_EQ(left.status,right.status);
        EXPECT_EQ(left.reason,right.reason);
        EXPECT_EQ(left.null_reference_status,right.null_reference_status);
        EXPECT_EQ(left.null_reference_failure_reason,right.null_reference_failure_reason);
        EXPECT_EQ(left.high_budget_reference_status,right.high_budget_reference_status);
        EXPECT_EQ(left.qualified_seed_count,right.qualified_seed_count);
        EXPECT_EQ(left.assessed_direction_count,right.assessed_direction_count);
        EXPECT_EQ(left.incomplete_direction_count,right.incomplete_direction_count);
        EXPECT_EQ(left.null_profile_evaluations,right.null_profile_evaluations);
        EXPECT_EQ(left.null_accepted_updates,right.null_accepted_updates);
        expect_same(left.reoptimization_prediction_drift,right.reoptimization_prediction_drift);
        expect_same(left.reoptimization_log_width_difference_inf,right.reoptimization_log_width_difference_inf);
        expect_same(left.reoptimization_scaled_ac_difference,right.reoptimization_scaled_ac_difference);
        expect_same(left.seed_prediction_robustness,right.seed_prediction_robustness);
        expect_same(left.seed_log_width_difference_inf,right.seed_log_width_difference_inf);
        expect_same(left.seed_scaled_ac_difference,right.seed_scaled_ac_difference);
        expect_same(left.high_budget_width_difference_inf,right.high_budget_width_difference_inf);
        expect_same(left.high_budget_scaled_ac_difference,right.high_budget_scaled_ac_difference);
        expect_same(left.high_budget_prediction_difference_relative,right.high_budget_prediction_difference_relative);
        ASSERT_EQ(left.directions.size(),right.directions.size());
        for(std::size_t direction=0;direction<left.directions.size();++direction)
        {
            EXPECT_EQ(left.directions[direction].direction,right.directions[direction].direction);
            EXPECT_EQ(left.directions[direction].status,right.directions[direction].status);
            EXPECT_EQ(left.directions[direction].reason,right.directions[direction].reason);
            EXPECT_EQ(left.directions[direction].prediction_endpoint,right.directions[direction].prediction_endpoint);
            EXPECT_EQ(left.directions[direction].log_width_endpoint,right.directions[direction].log_width_endpoint);
            EXPECT_EQ(left.directions[direction].scaled_ac_endpoint,right.directions[direction].scaled_ac_endpoint);
        }
    }
}
}
