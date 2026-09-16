#include <gtest/gtest.h>
#include "support/MDPDEExperiment.hpp"
#include "support/ForwardModelExperiment.hpp"
#include "support/EndpointRefinementExperiment.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <cmath>
#include <filesystem>

namespace {
second_stage_test::ShapeFixture Fixture(double alpha = 0.1)
{
    second_stage_test::ShapeFixture f; f.alpha = alpha;
    f.dataset.X.resize(40,2); f.dataset.y.resize(40);
    for (int i = 0; i < 40; ++i)
    {
        const double r{static_cast<double>(i)/40.0};
        f.dataset.X.row(i) << 1.0, -0.5*r*r;
        f.dataset.y(i) = 1.2 - 2.0*r*r + 0.01*std::sin(3.0*i+0.2);
    }
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(alpha, f.dataset, f.options);
    return f;
}
}

TEST(MDPDEExperimentTest, AlphaZeroMatchesClosedFormVarianceAndNormalEquations)
{
    const auto f{Fixture(0.0)};
    const auto beta{second_stage_test::MDPDETestBeta(f.dataset, Eigen::VectorXd::Ones(40), "svd")};
    const double v{(f.dataset.y-f.dataset.X*beta).squaredNorm()/40.0};
    const auto e{second_stage_test::EvaluateMDPDEEquations(f.dataset,0.0,beta,v,1e-8)};
    ASSERT_TRUE(e.valid); EXPECT_LT(e.scaled.lpNorm<Eigen::Infinity>(),1e-10);
    EXPECT_DOUBLE_EQ(e.denominator,40.0);
}

TEST(MDPDEExperimentTest, RootAndFixedPointAgreeOnFiniteVarianceFixture)
{
    const auto f{Fixture()}; const auto result{second_stage_test::CompareMDPDE(f)};
    ASSERT_TRUE(result.at("exact_replay").as_bool());
    double reference_v{};
    for (const auto & value : result.at("methods").as_array())
    {
        const auto & m{value.as_object()};
        if (m.at("method") == "production") continue;
        ASSERT_TRUE(m.at("reference_pass").as_bool()) << m.at("method");
        const double v{m.at("variance").as_double()};
        if (reference_v == 0.0) reference_v = v;
        EXPECT_NEAR(v,reference_v,1e-11);
    }
}

TEST(MDPDEExperimentTest, FreshRobustEquationsMatchIndependentScalarCalculation)
{
    const auto f{Fixture()}; Eigen::Vector2d beta; beta << 1.19,3.95;
    constexpr double v{0.0001}, floor{0.01}, alpha{0.5};
    const auto e{second_stage_test::EvaluateMDPDEEquations(f.dataset,alpha,beta,v,floor)};
    Eigen::Vector3d expected{Eigen::Vector3d::Zero()};
    double denominator{};
    for (Eigen::Index i=0;i<f.dataset.y.size();++i)
    {
        const double r{f.dataset.y(i)-f.dataset.X.row(i).dot(beta)};
        const double w{std::max(floor,std::exp(-alpha*r*r/(2*v)))};
        expected.head<2>() += f.dataset.X.row(i).transpose()*w*r;
        expected(2) += w*(r*r/v-1)/static_cast<double>(f.dataset.y.size());
        denominator += w;
        EXPECT_DOUBLE_EQ(e.weights(i),w);
    }
    const double correction{alpha/std::pow(1+alpha,1.5)};
    expected(2) += correction; denominator -= static_cast<double>(f.dataset.y.size())*correction;
    EXPECT_TRUE(e.raw.isApprox(expected,1e-12)); EXPECT_NEAR(e.denominator,denominator,1e-12);
}

TEST(MDPDEExperimentTest, BoundaryRankFloorAndDenominatorAreExplicit)
{
    auto f{Fixture()};
    const auto beta{f.expected.beta_ols};
    EXPECT_EQ(second_stage_test::EvaluateMDPDEEquations(f.dataset,.1,beta,0.0,1e-8).reason,"variance-boundary");
    f.dataset.X.col(1) = f.dataset.X.col(0);
    EXPECT_EQ(second_stage_test::EvaluateMDPDEEquations(f.dataset,.1,beta,1.0,1e-8).reason,"rank-deficient");
    f = Fixture(); f.dataset.y.array() += 100.0;
    const auto e{second_stage_test::EvaluateMDPDEEquations(f.dataset,.1,beta,1e-6,1e-8)};
    EXPECT_EQ(e.floor_count,40); EXPECT_EQ(e.reason,"invalid-denominator");
}

TEST(MDPDEExperimentTest, ExactFitIsNotReportedAsPositiveVarianceReference)
{
    auto f{Fixture()}; Eigen::Vector2d beta; beta << 1.0,4.0;
    f.dataset.y = f.dataset.X*beta;
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto result{second_stage_test::CompareMDPDE(f)};
    EXPECT_EQ(result.at("classification"),"roundoff-exact-fit-boundary");
    EXPECT_EQ(result.at("methods").as_array().size(),1u);
}

TEST(MDPDEExperimentTest, NearZeroFiniteNoiseRemainsAnEquationTest)
{
    auto f{Fixture()}; Eigen::Vector2d beta; beta << 1.0,4.0;
    for (int i=0;i<40;++i) f.dataset.y(i)=f.dataset.X.row(i).dot(beta)+1e-7*std::sin(3.0*i);
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto result{second_stage_test::CompareMDPDE(f,false)};
    EXPECT_FALSE(result.contains("classification"));
    EXPECT_EQ(result.at("methods").as_array().size(),8u);
    for (const auto & method : result.at("methods").as_array())
        if (method.at("method").as_string().starts_with("root-"))
            EXPECT_LE(method.at("equation_evaluations").as_int64() +
                method.at("verification_equation_evaluations").as_int64(),2000);
}

TEST(MDPDEExperimentTest, InvalidShapeStartingPointIsRecordedWithoutRootIterations)
{
    auto f{Fixture()};
    f.dataset.y = -f.dataset.y;
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto result{second_stage_test::CompareMDPDE(f,false)};
    for (const auto & method : result.at("methods").as_array())
    {
        if (!method.at("method").as_string().starts_with("root-")) continue;
        EXPECT_EQ(method.at("stop"),"invalid-start");
        EXPECT_EQ(method.at("iterations"),0);
        EXPECT_EQ(method.at("equation_evaluations"),0);
        EXPECT_FALSE(method.at("equation_pass").as_bool());
    }
}

TEST(MDPDEExperimentTest, SamplingWrapperUsesGridNodesAndExistingBoundaryClamping)
{
    rhbm_gem::MapObject map({4,4,4},{0.1,0.1,0.1},{-0.1,-0.1,-0.1});
    auto values{std::make_unique<double[]>(64)};
    for (std::size_t i=0;i<64;++i) values[i]=static_cast<double>(i);
    map.SetMapValueArray(std::move(values));
    const auto samples{second_stage_test::SampleExperimentPoints(map,{{0.0,{-0.1,-0.1,-0.1},true},{0.0,{0.0,0.0,0.0},true}})};
    EXPECT_DOUBLE_EQ(samples[0].response,map.GetMapValue(0,0,0));
    EXPECT_DOUBLE_EQ(samples[1].response,map.GetMapValue(1,1,1));
}

TEST(EndpointRefinementTest, AcceptedResultHasFreshWeightsAndExistingCovarianceFormula)
{
    for (const double alpha : {0.0,0.1,0.2})
    {
        const auto f{Fixture(alpha)};
        const auto refined{second_stage_test::RefineMDPDEEndpoint(f,2000)};
        ASSERT_TRUE(refined.accepted) << refined.evidence;
        const auto & result{refined.result};
        const auto fresh{second_stage_test::EvaluateMDPDEEquations(f.dataset,alpha,result.beta_mdpde,
            result.sigma_square,f.options.data_weight_min)};
        ASSERT_TRUE(fresh.valid); EXPECT_LE(fresh.scaled.lpNorm<Eigen::Infinity>(),1e-8);
        EXPECT_TRUE((result.data_weight.diagonal().array() == fresh.weights.array()).all());
        const double trace{fresh.weights.cwiseInverse().sum()};
        for (Eigen::Index i=0;i<fresh.weights.size();++i)
            EXPECT_DOUBLE_EQ(result.data_covariance.diagonal()(i),static_cast<double>(fresh.weights.size())*result.sigma_square/fresh.weights(i)/trace);
        EXPECT_TRUE((result.beta_ols.array() == f.expected.beta_ols.array()).all());
        EXPECT_EQ(result.diagnostics.iterations,f.expected.diagnostics.iterations);
        EXPECT_LE(refined.evidence.at("candidate_equation_evaluations").as_int64(),2000);
    }
}

TEST(EndpointRefinementTest, BudgetRejectionPreservesNativeHistoryButIsUnqualified)
{
    const auto f{Fixture()}; ASSERT_EQ(f.expected.status,rhbm_gem::RHBMEstimationStatus::SUCCESS);
    const auto refined{second_stage_test::RefineMDPDEEndpoint(f,11)};
    EXPECT_FALSE(refined.accepted); EXPECT_EQ(refined.evidence.at("reason"),"budget-exhausted");
    EXPECT_EQ(refined.result.status,f.expected.status);
    EXPECT_EQ(refined.result.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
    EXPECT_TRUE((refined.result.beta_mdpde.array() == f.expected.beta_mdpde.array()).all());
    EXPECT_DOUBLE_EQ(refined.result.sigma_square,f.expected.sigma_square);
    EXPECT_TRUE((refined.result.data_weight.diagonal().array() == f.expected.data_weight.diagonal().array()).all());
    EXPECT_TRUE((refined.result.data_covariance.diagonal().array() == f.expected.data_covariance.diagonal().array()).all());
    EXPECT_LE(refined.evidence.at("candidate_equation_evaluations").as_int64(),11);
}

TEST(EndpointRefinementTest, KnownOtherRootFailsBranchComparison)
{
    const auto path{std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() /
        "fixtures/mdpde/shape-maximum-iterations-final-32-33-99.txt"};
    auto f{second_stage_test::ReadShapeFixture(path.string())};
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto refined{second_stage_test::RefineMDPDEEndpoint(f,2000)};
    ASSERT_TRUE(refined.accepted) << refined.evidence;
    const double other_v{2.18540209421e-5}, other_a{6.99969719603}, other_b{0.499945192644};
    Eigen::Vector2d other; other << std::log(other_a/std::pow(2*std::acos(-1.0)*other_b*other_b,1.5)),1/(other_b*other_b);
    const auto e{second_stage_test::EvaluateMDPDEEquations(f.dataset,f.alpha,other,other_v,f.options.data_weight_min)};
    ASSERT_TRUE(e.valid); EXPECT_LT(e.scaled.lpNorm<Eigen::Infinity>(),1e-8);
    const auto reference{second_stage_test::EvaluateMDPDEEquations(f.dataset,f.alpha,refined.result.beta_mdpde,
        refined.result.sigma_square,f.options.data_weight_min)};
    const auto branch{second_stage_test::CompareMDPDEBranches(other,other_v,e,refined.result.beta_mdpde,
        refined.result.sigma_square,reference,f.options.data_weight_min)};
    EXPECT_FALSE(branch.at("pass").as_bool()); EXPECT_FALSE(branch.at("floor_masks_equal").as_bool());
}

TEST(EndpointRefinementTest, InvalidEndpointsAreNeverPromoted)
{
    auto f{Fixture()}; f.expected.sigma_square=0;
    auto result{second_stage_test::RefineMDPDEEndpoint(f,64)};
    EXPECT_FALSE(result.accepted); EXPECT_EQ(result.evidence.at("reason"),"variance-boundary");
    f=Fixture(); f.dataset.X.col(1)=f.dataset.X.col(0);
    result=second_stage_test::RefineMDPDEEndpoint(f,64);
    EXPECT_FALSE(result.accepted); EXPECT_EQ(result.evidence.at("reason"),"rank-deficient");
    f=Fixture(); f.dataset.y.array()+=100;
    result=second_stage_test::RefineMDPDEEndpoint(f,64);
    EXPECT_FALSE(result.accepted); EXPECT_EQ(result.evidence.at("reason"),"invalid-denominator");
    f=Fixture(); f.dataset.X.conservativeResize(1,2); f.dataset.y.conservativeResize(1);
    f.expected=rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto insufficient{rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(f.dataset,f.alpha,f.options,f.expected)};
    EXPECT_EQ(insufficient.status,rhbm_gem::RHBMEstimationStatus::INSUFFICIENT_DATA);
    EXPECT_EQ(insufficient.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
    EXPECT_EQ(insufficient.refinement->reference_updates,0);
    EXPECT_EQ(insufficient.refinement->candidate_equation_evaluations,1);
}

TEST(EndpointRefinementTest, ExactFitAndNearZeroNoiseRemainDistinct)
{
    auto f{Fixture()}; Eigen::Vector2d beta; beta << 1.0,4.0; f.dataset.y=f.dataset.X*beta;
    f.expected=rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto exact{second_stage_test::RefineMDPDEEndpoint(f,128)};
    EXPECT_FALSE(exact.accepted); EXPECT_EQ(exact.evidence.at("reason"),"roundoff-exact-fit-boundary");
    for (int i=0;i<40;++i) f.dataset.y(i)+=1e-7*std::sin(3.0*i);
    f.expected=rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto noisy{second_stage_test::RefineMDPDEEndpoint(f,128)};
    EXPECT_NE(noisy.evidence.at("reason"),"roundoff-exact-fit-boundary");
    EXPECT_LE(noisy.evidence.at("candidate_equation_evaluations").as_int64(),128);
    if (noisy.accepted) EXPECT_TRUE(noisy.evidence.at("branch").at("pass").as_bool());
    else EXPECT_EQ(noisy.result.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
}

TEST(EndpointRefinementTest, PoliciesDoNotConfuseNativeSuccessWithFreshQualification)
{
    using second_stage_test::EndpointPolicy; using second_stage_test::ShouldRefineEndpoint;
    using rhbm_gem::RHBMEstimationStatus;
    for (const bool fresh : {false,true})
    {
        EXPECT_FALSE(ShouldRefineEndpoint(EndpointPolicy::Legacy,RHBMEstimationStatus::MAX_ITERATIONS_REACHED,fresh));
        EXPECT_FALSE(ShouldRefineEndpoint(EndpointPolicy::FailedOnly,RHBMEstimationStatus::SUCCESS,fresh));
        EXPECT_TRUE(ShouldRefineEndpoint(EndpointPolicy::FailedOnly,RHBMEstimationStatus::MAX_ITERATIONS_REACHED,fresh));
        EXPECT_TRUE(ShouldRefineEndpoint(EndpointPolicy::FreshResidual,RHBMEstimationStatus::MAX_ITERATIONS_REACHED,fresh));
        EXPECT_EQ(ShouldRefineEndpoint(EndpointPolicy::FreshResidual,RHBMEstimationStatus::SUCCESS,fresh),!fresh);
    }
}

TEST(EndpointRefinementTest, FailedOnlyLeavesNativeSuccessUntouchedDespiteFreshResidual)
{
    auto f{Fixture()};
    f.options.tolerance = 0.1;
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    ASSERT_EQ(f.expected.status,rhbm_gem::RHBMEstimationStatus::SUCCESS);
    const auto fresh{second_stage_test::EvaluateMDPDEEquations(f.dataset,f.alpha,f.expected.beta_mdpde,
        f.expected.sigma_square,f.options.data_weight_min)};
    ASSERT_GT(fresh.scaled.lpNorm<Eigen::Infinity>(),1e-8);
    const auto result{rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(f.dataset,f.alpha,f.options,f.expected)};
    EXPECT_FALSE(result.refinement);
    EXPECT_EQ(result.Qualification(),rhbm_gem::RHBMSolveQualification::NativeSuccess);
    EXPECT_TRUE((result.beta_mdpde.array() == f.expected.beta_mdpde.array()).all());
    EXPECT_DOUBLE_EQ(result.sigma_square,f.expected.sigma_square);
    EXPECT_TRUE((result.data_weight.diagonal().array() == f.expected.data_weight.diagonal().array()).all());
    EXPECT_TRUE((result.data_covariance.diagonal().array() == f.expected.data_covariance.diagonal().array()).all());
}

TEST(EndpointRefinementTest, ProductionRefinesCapturedFailuresWithoutRewritingNativeStatus)
{
    const auto directory{std::filesystem::path(__FILE__).parent_path().parent_path().parent_path()/"fixtures/mdpde"};
    int count{};
    for (const auto & entry : std::filesystem::directory_iterator(directory))
    {
        if (!entry.path().filename().string().starts_with("shape-") || entry.path().extension() != ".txt") continue;
        auto f{second_stage_test::ReadShapeFixture(entry.path().string())};
        f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
        ASSERT_EQ(f.expected.status,rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED);
        const auto result{rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(f.dataset,f.alpha,f.options,f.expected)};
        ASSERT_TRUE(result.refinement);
        ASSERT_TRUE(result.refinement->accepted) << result.refinement->reason;
        EXPECT_EQ(result.status,f.expected.status);
        EXPECT_EQ(result.Qualification(),rhbm_gem::RHBMSolveQualification::RefinedSuccess);
        EXPECT_LE(result.refinement->candidate_equation_evaluations,128);
        EXPECT_LE(*result.refinement->candidate_residual,1e-8);
        EXPECT_LE(*result.refinement->reference_residual,1e-10);
        EXPECT_TRUE(*result.refinement->floor_masks_equal);
        EXPECT_EQ(result.diagnostics.iterations,f.expected.diagnostics.iterations);
        const auto rejected{rhbm_gem::mdpde_detail::RefineMDPDEEndpoint(f.dataset,f.alpha,f.options,f.expected,11)};
        EXPECT_FALSE(rejected.refinement->accepted);
        EXPECT_EQ(rejected.status,f.expected.status);
        EXPECT_EQ(rejected.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
        EXPECT_TRUE((rejected.beta_mdpde.array() == f.expected.beta_mdpde.array()).all());
        EXPECT_DOUBLE_EQ(rejected.sigma_square,f.expected.sigma_square);
        EXPECT_TRUE((rejected.data_weight.diagonal().array() == f.expected.data_weight.diagonal().array()).all());
        EXPECT_TRUE((rejected.data_covariance.diagonal().array() == f.expected.data_covariance.diagonal().array()).all());
        ++count;
    }
    EXPECT_EQ(count,3);
}

TEST(EndpointRefinementTest, ReferenceBudgetFailureDoesNotUseTheCandidateOrReferenceAsFallback)
{
    auto f{Fixture()};
    f.expected.status = rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    f.expected.diagnostics.iterations = 10000;
    const auto result{rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(f.dataset,f.alpha,f.options,f.expected)};
    ASSERT_TRUE(result.refinement);
    EXPECT_EQ(result.refinement->reason,"reference-unqualified");
    EXPECT_EQ(result.refinement->reference_stop,"budget-exhausted");
    EXPECT_EQ(result.refinement->reference_updates,0);
    EXPECT_EQ(result.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
    EXPECT_TRUE((result.beta_mdpde.array() == f.expected.beta_mdpde.array()).all());
    EXPECT_DOUBLE_EQ(result.sigma_square,f.expected.sigma_square);
}

TEST(EndpointRefinementTest, ReferenceStagnationIsRejectedEvenWithQualifiedRootEquations)
{
    auto f{Fixture(0.0)};
    for (int i=0;i<40;++i)
    {
        const double r{static_cast<double>(i)/40.0};
        f.dataset.y(i) = 1.0 - 2.0*r*r + 1e-6*std::sin(3.0*i+0.2);
    }
    f.expected = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options);
    const auto result{rhbm_gem::mdpde_detail::RefineMDPDEEndpoint(f.dataset,f.alpha,f.options,f.expected)};
    ASSERT_TRUE(result.refinement);
    EXPECT_LE(result.refinement->candidate_residual.value_or(INFINITY),1e-8);
    EXPECT_EQ(result.refinement->reference_stop,"stalled");
    EXPECT_EQ(result.refinement->reason,"reference-unqualified");
    EXPECT_EQ(result.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
    EXPECT_EQ(result.status,f.expected.status);
    EXPECT_TRUE((result.beta_mdpde.array() == f.expected.beta_mdpde.array()).all());
    EXPECT_TRUE((result.data_covariance.diagonal().array() == f.expected.data_covariance.diagonal().array()).all());
}

TEST(EndpointRefinementTest, FloorMaskMismatchRejectsOtherwiseCloseBranches)
{
    const auto f{Fixture()};
    const auto reference{second_stage_test::EvaluateMDPDEEquations(f.dataset,f.alpha,
        f.expected.beta_mdpde,f.expected.sigma_square,f.options.data_weight_min)};
    auto candidate{reference};
    constexpr double floor{1e-8};
    candidate.weights(0) = floor;
    auto other{reference}; other.weights(0) = floor + 1e-10;
    const auto branch{rhbm_gem::mdpde_detail::CompareMDPDEBranches(f.expected.beta_mdpde,f.expected.sigma_square,
        candidate,f.expected.beta_mdpde,f.expected.sigma_square,other,floor)};
    ASSERT_TRUE(branch.weight_max_difference);
    EXPECT_LT(*branch.weight_max_difference,1e-6);
    EXPECT_FALSE(*branch.floor_masks_equal);
    EXPECT_FALSE(branch.pass);
}
