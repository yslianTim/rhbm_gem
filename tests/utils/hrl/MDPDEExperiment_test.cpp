#include <gtest/gtest.h>
#include "support/MDPDEExperiment.hpp"
#include "support/ForwardModelExperiment.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <cmath>

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
