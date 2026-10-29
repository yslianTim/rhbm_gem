#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "core/detail/joint_component/SplitDiagnostics.hpp"
#include "core/detail/joint_component/Problem.hpp"
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
}
