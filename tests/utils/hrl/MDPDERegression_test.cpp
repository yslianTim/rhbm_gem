#include <gtest/gtest.h>
#include "support/MDPDETestSupport.hpp"
#include "support/EndpointRefinementTestSupport.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>

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

TEST(MDPDENumericsTest, AlphaZeroMatchesClosedFormVarianceAndNormalEquations)
{
    const auto f{Fixture(0.0)};
    const auto beta{second_stage_test::MDPDETestBeta(f.dataset, Eigen::VectorXd::Ones(40), "svd")};
    const double v{(f.dataset.y-f.dataset.X*beta).squaredNorm()/40.0};
    const auto e{second_stage_test::EvaluateMDPDEEquations(f.dataset,0.0,beta,v,1e-8)};
    ASSERT_TRUE(e.valid); EXPECT_LT(e.scaled.lpNorm<Eigen::Infinity>(),1e-10);
    EXPECT_DOUBLE_EQ(e.denominator,40.0);
    const auto production{rhbm_gem::rhbm_helper::EstimateBetaMDPDE(0.0,f.dataset,f.options)};
    EXPECT_EQ(production.status,rhbm_gem::RHBMEstimationStatus::SUCCESS);
    EXPECT_TRUE(production.beta_mdpde.isApprox(beta,1e-12));
    EXPECT_NEAR(production.sigma_square,v,1e-12);
}

TEST(MDPDENumericsTest, ProductionFixedPointMatchesIndependentQRReference)
{
    const auto f{Fixture()};
    const auto production{rhbm_gem::rhbm_helper::EstimateBetaMDPDE(f.alpha,f.dataset,f.options)};
    Eigen::VectorXd beta{f.expected.beta_ols};
    double variance{(f.dataset.y-f.dataset.X*beta).squaredNorm()/(f.dataset.y.size()-1)};
    bool converged{};
    for (int i=0;i<f.options.max_iterations;++i)
    {
        const auto equations{second_stage_test::EvaluateMDPDEEquations(
            f.dataset,f.alpha,beta,variance,f.options.data_weight_min)};
        ASSERT_TRUE(equations.valid) << equations.reason;
        const auto previous_beta{beta}; const double previous_variance{variance};
        beta=second_stage_test::MDPDETestBeta(f.dataset,equations.weights,"qr");
        variance=second_stage_test::MDPDETestVariance(f.dataset,f.alpha,equations.weights,beta);
        const double beta_change{(beta-previous_beta).squaredNorm()};
        const double variance_change{std::abs(variance-previous_variance)/
            std::max({std::abs(variance),std::abs(previous_variance),f.options.data_weight_min})};
        if(beta_change<f.options.tolerance && variance_change<f.options.tolerance)
        {converged=true;break;}
    }
    ASSERT_TRUE(converged);
    ASSERT_EQ(production.status,rhbm_gem::RHBMEstimationStatus::SUCCESS);
    EXPECT_LT((production.beta_mdpde-beta).norm(),1e-8);
    EXPECT_NEAR(production.sigma_square,variance,1e-11);
    const auto fresh{second_stage_test::EvaluateMDPDEEquations(f.dataset,f.alpha,
        production.beta_mdpde,production.sigma_square,f.options.data_weight_min)};
    ASSERT_TRUE(fresh.valid);
}

TEST(MDPDENumericsTest, FreshRobustEquationsMatchIndependentScalarCalculation)
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

TEST(MDPDENumericsTest, BoundaryRankFloorAndDenominatorAreExplicit)
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
    const auto branch{rhbm_gem::mdpde_detail::CompareMDPDEBranches(other,other_v,e,refined.result.beta_mdpde,
        refined.result.sigma_square,reference,f.options.data_weight_min)};
    EXPECT_FALSE(branch.pass); ASSERT_TRUE(branch.floor_masks_equal); EXPECT_FALSE(*branch.floor_masks_equal);
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

TEST(EndpointRefinementTest, FailedOnlyDoesNotPromoteAnIneligibleRankDeficientEndpoint)
{
    auto f{Fixture()};
    f.dataset.X.col(1) = f.dataset.X.col(0);
    auto endpoint{f.expected};
    endpoint.status = rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    const auto result{rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(
        f.dataset, f.alpha, f.options, endpoint)};
    ASSERT_TRUE(result.refinement);
    EXPECT_FALSE(result.refinement->accepted);
    EXPECT_EQ(result.refinement->reason,"rank-deficient");
    EXPECT_EQ(result.refinement->reference_updates,0);
    EXPECT_EQ(result.status,endpoint.status);
    EXPECT_EQ(result.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
    EXPECT_TRUE((result.beta_mdpde.array() == endpoint.beta_mdpde.array()).all());
    EXPECT_DOUBLE_EQ(result.sigma_square,endpoint.sigma_square);
}

TEST(EndpointRefinementTest, NegativeAndNonfiniteVarianceAreNeverPromoted)
{
    for (const double variance : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::infinity()})
    {
        auto f{Fixture()};
        f.expected.sigma_square = variance;
        const auto result{second_stage_test::RefineMDPDEEndpoint(f,128)};
        EXPECT_FALSE(result.accepted);
        EXPECT_EQ(result.result.Qualification(),rhbm_gem::RHBMSolveQualification::Unqualified);
        if (std::isnan(variance)) EXPECT_TRUE(std::isnan(result.result.sigma_square));
        else EXPECT_DOUBLE_EQ(result.result.sigma_square,variance);
    }
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
    EXPECT_EQ(result.status,f.expected.status);
    EXPECT_EQ(result.Qualification(),rhbm_gem::RHBMSolveQualification::NativeSuccess);
    EXPECT_TRUE((result.beta_mdpde.array() == f.expected.beta_mdpde.array()).all());
    EXPECT_DOUBLE_EQ(result.sigma_square,f.expected.sigma_square);
    EXPECT_EQ(result.diagnostics.iterations,f.expected.diagnostics.iterations);
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

TEST(EndpointRefinementTest, Captured6Z6UFirstStageFailuresHaveReproducibleRecoveryEvidence)
{
    const auto path{std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() /
        "fixtures/mdpde/6z6u-fixed-halo-first-stage.json"};
    std::ifstream input(path);
    ASSERT_TRUE(input) << path;
    boost::json::parse_options parse_options;
    parse_options.numbers = boost::json::number_precision::precise;
    const auto fixture = boost::json::parse(
        std::string(std::istreambuf_iterator<char>(input), {}), {}, parse_options);
    ASSERT_EQ(fixture.at("fixture_contract").as_string(), "6z6u-first-stage-mdpde-diagnostic-v1");
    const auto & halos = fixture.at("failed_halos").as_array();
    ASSERT_EQ(halos.size(), 4u);

    for (const auto & halo : halos)
    {
        const auto & input_data = halo.at("dataset");
        const auto & rows = input_data.at("design_matrix").as_array();
        const auto & response = input_data.at("log_response").as_array();
        ASSERT_GT(rows.size(), 2u);
        ASSERT_EQ(rows.size(), response.size());
        rhbm_gem::RHBMMemberDataset dataset;
        dataset.X.resize(static_cast<Eigen::Index>(rows.size()), 2);
        dataset.y.resize(static_cast<Eigen::Index>(response.size()));
        for (std::size_t row = 0; row < rows.size(); ++row)
        {
            const auto & values = rows[row].as_array();
            ASSERT_EQ(values.size(), 2u);
            dataset.X(static_cast<Eigen::Index>(row), 0) = boost::json::value_to<double>(values[0]);
            dataset.X(static_cast<Eigen::Index>(row), 1) = boost::json::value_to<double>(values[1]);
            dataset.y(static_cast<Eigen::Index>(row)) = boost::json::value_to<double>(response[row]);
        }
        const double alpha = boost::json::value_to<double>(halo.at("local_alpha"));
        rhbm_gem::RHBMExecutionOptions options;
        options.thread_size = 1;
        options.tolerance = boost::json::value_to<double>(
            halo.at("native_solver_options").at("tolerance"));
        options.data_weight_min = boost::json::value_to<double>(
            halo.at("native_solver_options").at("data_weight_min"));

        options.max_iterations = 100;
        const auto baseline = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(alpha, dataset, options);
        EXPECT_EQ(baseline.status, rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED)
            << halo.at("atom_id").as_string();
        EXPECT_EQ(baseline.Qualification(), rhbm_gem::RHBMSolveQualification::Unqualified);
        const auto & expected_baseline = halo.at("runs").at("native_100");
        const auto & expected_beta = expected_baseline.at("beta_mdpde").as_array();
        ASSERT_EQ(expected_beta.size(), 2u);
        for (std::size_t k = 0; k < 2; ++k)
            EXPECT_NEAR(baseline.beta_mdpde(static_cast<Eigen::Index>(k)),
                boost::json::value_to<double>(expected_beta[k]), 1e-12);

        for (const int budget : {300, 1000})
        {
            options.max_iterations = budget;
            const auto extended = rhbm_gem::rhbm_helper::EstimateBetaMDPDE(alpha, dataset, options);
            const auto key = "native_" + std::to_string(budget);
            const auto & expected = halo.at("runs").at(key);
            EXPECT_EQ(extended.status, rhbm_gem::RHBMEstimationStatus::SUCCESS)
                << halo.at("atom_id").as_string() << " at " << budget;
            EXPECT_EQ(extended.Qualification(), rhbm_gem::RHBMSolveQualification::NativeSuccess);
            EXPECT_EQ(extended.diagnostics.iterations,
                boost::json::value_to<int>(expected.at("iterations")));
            EXPECT_LE(*extended.diagnostics.relative_variance_change, options.tolerance);
            const auto & expected_extended_beta = expected.at("beta_mdpde").as_array();
            for (std::size_t k = 0; k < 2; ++k)
                EXPECT_NEAR(extended.beta_mdpde(static_cast<Eigen::Index>(k)),
                    boost::json::value_to<double>(expected_extended_beta[k]), 1e-12);
            if (budget == 300)
            {
                const auto native_success_unchanged = rhbm_gem::core::detail::RecoverFailedFirstStageEstimate(
                    dataset, alpha, options, extended);
                EXPECT_EQ(native_success_unchanged.Qualification(), rhbm_gem::RHBMSolveQualification::NativeSuccess);
                EXPECT_FALSE(native_success_unchanged.failed_only_recovery);
                EXPECT_EQ(native_success_unchanged.diagnostics.iterations, extended.diagnostics.iterations);
                EXPECT_TRUE(native_success_unchanged.beta_mdpde.isApprox(extended.beta_mdpde, 0.0));
            }
            if (budget == 1000)
            {
                const auto & trajectory = halo.at("native_fixed_point_trajectory_to_1000");
                EXPECT_TRUE(boost::json::value_to<bool>(trajectory.at("converged")));
                EXPECT_EQ(boost::json::value_to<int>(trajectory.at("iterations")),
                    extended.diagnostics.iterations);
                EXPECT_EQ(boost::json::value_to<int>(trajectory.at("nonfinite_intermediates")), 0);
                EXPECT_FALSE(boost::json::value_to<bool>(halo.at("initial_ols").at("roundoff_exact_fit")));
            }
        }

        options.max_iterations = 100;
        const auto recovered = rhbm_gem::mdpde_detail::ApplyFailedOnlyRefinement(
            dataset, alpha, options, baseline);
        const auto & expected_recovery = halo.at("runs").at("existing_failed_only_refinement");
        const auto & recovery = expected_recovery.at("recovery");
        ASSERT_TRUE(recovered.refinement);
        EXPECT_EQ(recovered.status, baseline.status);
        EXPECT_EQ(recovered.refinement->accepted, boost::json::value_to<bool>(recovery.at("accepted")));
        EXPECT_EQ(recovered.refinement->reason, recovery.at("reason").as_string());
        if (recovered.refinement->accepted)
        {
            EXPECT_EQ(recovered.Qualification(), rhbm_gem::RHBMSolveQualification::RefinedSuccess);
            EXPECT_LE(*recovered.refinement->candidate_residual, 1e-8);
            EXPECT_LE(*recovered.refinement->reference_residual, 1e-10);
            EXPECT_TRUE(*recovered.refinement->floor_masks_equal);
            EXPECT_LE(recovered.refinement->candidate_equation_evaluations,
                rhbm_gem::mdpde_detail::kEndpointEquationBudget);
        }
        else
        {
            EXPECT_EQ(recovered.Qualification(), rhbm_gem::RHBMSolveQualification::Unqualified);
            EXPECT_TRUE((recovered.beta_mdpde.array() == baseline.beta_mdpde.array()).all());
            EXPECT_DOUBLE_EQ(recovered.sigma_square, baseline.sigma_square);
        }

        const auto production_recovery = rhbm_gem::core::detail::RecoverFailedFirstStageEstimate(
            dataset, alpha, options, baseline);
        EXPECT_EQ(production_recovery.status, baseline.status)
            << halo.at("atom_id").as_string();
        ASSERT_TRUE(production_recovery.failed_only_recovery);
        ASSERT_TRUE(production_recovery.failed_only_recovery->accepted)
            << halo.at("atom_id").as_string() << ": "
            << production_recovery.failed_only_recovery->reason;
        EXPECT_EQ(production_recovery.Qualification(), rhbm_gem::RHBMSolveQualification::RefinedSuccess);
        if (halo.at("atom_id").as_string() == "16664")
        {
            EXPECT_FALSE(production_recovery.failed_only_recovery->endpoint_refinement->accepted);
            ASSERT_TRUE(production_recovery.failed_only_recovery->bounded_native_retry);
            EXPECT_EQ(production_recovery.failed_only_recovery->bounded_native_retry->max_iterations, 300);
            EXPECT_EQ(production_recovery.failed_only_recovery->bounded_native_retry->native_status,
                static_cast<int>(rhbm_gem::RHBMEstimationStatus::SUCCESS));
            EXPECT_EQ(production_recovery.failed_only_recovery->bounded_native_retry->iterations, 268);
            ASSERT_TRUE(production_recovery.failed_only_recovery->bounded_native_retry->squared_beta_change);
            EXPECT_LT(*production_recovery.failed_only_recovery->bounded_native_retry->squared_beta_change,
                options.tolerance);
            ASSERT_TRUE(production_recovery.failed_only_recovery->bounded_native_retry->relative_variance_change);
            EXPECT_LT(*production_recovery.failed_only_recovery->bounded_native_retry->relative_variance_change,
                options.tolerance);
            ASSERT_TRUE(production_recovery.failed_only_recovery->bounded_native_retry->equation_residual);
            EXPECT_GT(*production_recovery.failed_only_recovery->bounded_native_retry->equation_residual, 1e-8);
        }
        else
        {
            ASSERT_TRUE(production_recovery.failed_only_recovery->endpoint_refinement);
            EXPECT_TRUE(production_recovery.failed_only_recovery->endpoint_refinement->accepted);
            EXPECT_FALSE(production_recovery.failed_only_recovery->bounded_native_retry);
        }

        rhbm_gem::JointInitializationAtom initialization_atom;
        initialization_atom.id = halo.at("atom_id").as_string().c_str();
        initialization_atom.target = false;
        initialization_atom.mdpde = {
            production_recovery.beta_mdpde(0), production_recovery.beta_mdpde(1), 0.0};
        initialization_atom.native_status = static_cast<int>(production_recovery.status);
        initialization_atom.first_stage_fit_qualified = true;
        initialization_atom.first_stage_qualification = production_recovery.Qualification();
        initialization_atom.first_stage_recovery = production_recovery.failed_only_recovery;
        std::string failure_reason;
        const auto fixed_halo = rhbm_gem::core::detail::MakeFixedHaloParameter(
            initialization_atom, 0, initialization_atom.id, failure_reason);
        ASSERT_TRUE(fixed_halo) << failure_reason;
        EXPECT_EQ(fixed_halo->first_stage_status,
            static_cast<int>(rhbm_gem::RHBMEstimationStatus::MAX_ITERATIONS_REACHED));
        EXPECT_DOUBLE_EQ(fixed_halo->charge, 0.0);
        EXPECT_TRUE(rhbm_gem::IsQualifiedFirstStageFixedHalo(*fixed_halo));
    }
}
