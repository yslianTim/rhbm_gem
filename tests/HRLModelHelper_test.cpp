#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <cmath>

#include "HRLModelHelper.hpp"
#include "Logger.hpp"

namespace {

using Eigen::VectorXd;
using DataTuple = std::tuple<std::vector<Eigen::VectorXd>, std::string>;

DataTuple CreateMember(const std::vector<std::pair<double, double>>& xy, const std::string& name)
{
    std::vector<VectorXd> data;
    data.reserve(xy.size());
    for (const auto& [x, y] : xy)
    {
        VectorXd v(3);
        v << 1.0, x, y;
        data.emplace_back(v);
    }
    return {data, name};
}

}

TEST(HRLModelHelperTest, DefaultInitializationValues)
{
    const int basis_size{ 2 };
    const int member_size{ 2 };

    HRLModelHelper helper{ basis_size, member_size };
    Eigen::MatrixXd expected_lambda{ Eigen::MatrixXd::Identity(basis_size, basis_size) };
    EXPECT_TRUE(helper.GetCapitalLambdaMatrix().isApprox(expected_lambda));

    Eigen::VectorXd expected_mu{ Eigen::VectorXd::Ones(basis_size) };
    EXPECT_TRUE(helper.GetMuVectorPrior().isApprox(expected_mu));
    EXPECT_TRUE(helper.GetMuVectorMDPDE().isApprox(expected_mu));
    EXPECT_TRUE(helper.GetMuVectorMean().isApprox(expected_mu));

    Eigen::VectorXd expected_beta{ Eigen::VectorXd::Ones(basis_size) };
    for (int member = 0; member < member_size; member++)
    {
        EXPECT_TRUE(helper.GetBetaMatrixOLS(member).isApprox(expected_beta));
        EXPECT_TRUE(helper.GetBetaMatrixMDPDE(member).isApprox(expected_beta));
        EXPECT_TRUE(helper.GetBetaMatrixPosterior(member).isApprox(expected_beta));
        EXPECT_DOUBLE_EQ(0.0, helper.GetStatisticalDistance(member));
        EXPECT_FALSE(helper.GetOutlierFlag(member));
    }
}

TEST(HRLModelHelperTest, AcceptsProperlySizedData)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<Eigen::VectorXd> member_data;
    Eigen::VectorXd sample1(2);
    sample1 << 1.0, 2.0;
    Eigen::VectorXd sample2(2);
    sample2 << 2.0, 5.0;
    member_data.emplace_back(sample1);
    member_data.emplace_back(sample2);
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member_data, "member1");
    ASSERT_NO_THROW(helper.SetDataArray(data_array));
    helper.RunEstimation(0.0, 0.0);

    Eigen::VectorXd beta{ helper.GetBetaMatrixMDPDE(0) };
    ASSERT_EQ(beta.size(), 1);
    EXPECT_NEAR(beta(0), 2.4, 1e-9);

    Eigen::VectorXd mu{ helper.GetMuVectorMDPDE() };
    ASSERT_EQ(mu.size(), 1);
    EXPECT_NEAR(mu(0), 2.4, 1e-9);

    Eigen::VectorXd mu_prior{ helper.GetMuVectorPrior() };
    ASSERT_EQ(mu_prior.size(), 1);
    EXPECT_NEAR(mu_prior(0), 2.4, 1e-9);
}

TEST(HRLModelHelperTest, ThrowsWhenMemberSizeMismatch)
{
    HRLModelHelper helper{ 1, 1 };
    Eigen::VectorXd sample(2);
    sample << 1.0, 1.0;
    std::vector<Eigen::VectorXd> member1{ sample };
    std::vector<Eigen::VectorXd> member2{ sample };
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member1, "m1");
    data_array.emplace_back(member2, "m2");

    EXPECT_THROW(helper.SetDataArray(data_array), std::runtime_error);
}

TEST(HRLModelHelperTest, ThrowsWhenSampleVectorSizeMismatch)
{
    HRLModelHelper helper{ 1, 1 };
    // Create a sample with size not equal to basis_size + 1 (expected 2)
    Eigen::VectorXd bad_sample(3);
    bad_sample << 1.0, 2.0, 3.0;
    std::vector<Eigen::VectorXd> member{ bad_sample };
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member, "member1");
    EXPECT_THROW(helper.SetDataArray(data_array), std::runtime_error);
}

TEST(HRLModelHelperTest, EstimationOnSmallSyntheticData)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_1") };
    std::vector<DataTuple> data_array{member0, member1};

    HRLModelHelper helper(2, 2);
    helper.SetDataArray(data_array);
    helper.RunEstimation(0.0, 0.0);

    Eigen::Vector2d expected;
    expected << 5.0 / 6.0, 2.5; // Analytic OLS solution

    for (int i = 0; i < 2; ++i)
    {
        auto beta_post{ helper.GetBetaMatrixPosterior(i) };
        auto beta_mdpde{ helper.GetBetaMatrixMDPDE(i) };
        auto beta_ols{ helper.GetBetaMatrixOLS(i) };
        ASSERT_NEAR(expected(0), beta_post(0), 1.0e-6);
        ASSERT_NEAR(expected(1), beta_post(1), 1.0e-6);
        ASSERT_NEAR(expected(0), beta_mdpde(0), 1.0e-6);
        ASSERT_NEAR(expected(1), beta_mdpde(1), 1.0e-6);
        ASSERT_NEAR(expected(0), beta_ols(0), 1.0e-6);
        ASSERT_NEAR(expected(1), beta_ols(1), 1.0e-6);
        EXPECT_TRUE(std::isfinite(helper.GetStatisticalDistance(i)));
        EXPECT_FALSE(helper.GetOutlierFlag(i));
    }

    auto mu_prior{ helper.GetMuVectorPrior() };
    auto mu_mdpde{ helper.GetMuVectorMDPDE() };
    auto mu_mean{ helper.GetMuVectorMean() };
    ASSERT_NEAR(expected(0), mu_prior(0), 1.0e-6);
    ASSERT_NEAR(expected(1), mu_prior(1), 1.0e-6);
    ASSERT_NEAR(expected(0), mu_mdpde(0), 1.0e-6);
    ASSERT_NEAR(expected(1), mu_mdpde(1), 1.0e-6);
    ASSERT_NEAR(expected(0), mu_mean(0), 1.0e-6);
    ASSERT_NEAR(expected(1), mu_mean(1), 1.0e-6);
}

TEST(HRLModelHelperTest, OutlierFlagsWithDivergentData)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 10.0}, {1.0, 12.0}, {2.0, 15.0}}, "member_1") };
    std::vector<DataTuple> data_array{member0, member1};

    HRLModelHelper helper(2, 2);
    helper.SetDataArray(data_array);
    helper.RunEstimation(0.0, 0.0);

    Eigen::Vector2d expected0, expected1;
    expected0 << 5.0 / 6.0, 2.5;
    expected1 << 59.0 / 6.0, 2.5;

    auto beta_ols0{ helper.GetBetaMatrixOLS(0) };
    auto beta_ols1{ helper.GetBetaMatrixOLS(1) };
    ASSERT_NEAR(expected0(0), beta_ols0(0), 1.0e-6);
    ASSERT_NEAR(expected0(1), beta_ols0(1), 1.0e-6);
    ASSERT_NEAR(expected1(0), beta_ols1(0), 1.0e-6);
    ASSERT_NEAR(expected1(1), beta_ols1(1), 1.0e-6);

    auto beta_mdpde0{ helper.GetBetaMatrixMDPDE(0) };
    auto beta_mdpde1{ helper.GetBetaMatrixMDPDE(1) };
    ASSERT_NEAR(expected0(0), beta_mdpde0(0), 1.0e-6);
    ASSERT_NEAR(expected0(1), beta_mdpde0(1), 1.0e-6);
    ASSERT_NEAR(expected1(0), beta_mdpde1(0), 1.0e-6);
    ASSERT_NEAR(expected1(1), beta_mdpde1(1), 1.0e-6);

    auto mu_prior{ helper.GetMuVectorPrior() };
    auto mu_mdpde{ helper.GetMuVectorMDPDE() };
    auto mu_mean{ helper.GetMuVectorMean() };
    Eigen::Vector2d expected_mu;
    expected_mu << (5.0 / 6.0 + 59.0 / 6.0) / 2.0, (2.5 + 2.5) / 2.0;
    ASSERT_NEAR(expected_mu(0), mu_prior(0), 1.0e-6);
    ASSERT_NEAR(expected_mu(1), mu_prior(1), 1.0e-6);
    ASSERT_NEAR(expected_mu(0), mu_mdpde(0), 1.0e-6);
    ASSERT_NEAR(expected_mu(1), mu_mdpde(1), 1.0e-6);
    ASSERT_NEAR(expected_mu(0), mu_mean(0), 1.0e-6);
    ASSERT_NEAR(expected_mu(1), mu_mean(1), 1.0e-6);

    for (int i = 0; i < 2; ++i)
    {
        EXPECT_TRUE(std::isfinite(helper.GetStatisticalDistance(i)));
        EXPECT_FALSE(helper.GetOutlierFlag(i));
    }
}

class HRLModelHelperFixture : public ::testing::Test
{
protected:
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> m_data;
    LogLevel m_prev_level;

    void SetUp(void) override
    {
        m_prev_level = Logger::GetLogLevel();
        Logger::SetLogLevel(LogLevel::Info);
        std::vector<Eigen::VectorXd> member;
        Eigen::VectorXd d1(2);
        d1 << 1.0, 1.0;
        member.emplace_back(d1);
        Eigen::VectorXd d2(2);
        d2 << 2.0, 2.0;
        member.emplace_back(d2);
        Eigen::VectorXd d3(2);
        d3 << 3.0, 3.0;
        member.emplace_back(d3);
        Eigen::VectorXd d4(2);
        d4 << 4.0, 30.0;
        member.emplace_back(d4);
        m_data.emplace_back(member, "m1");
    }

    void TearDown(void) override
    {
        Logger::SetLogLevel(m_prev_level);
    }
};

TEST_F(HRLModelHelperFixture, HonorsMaximumIteration)
{
    // Baseline with large iteration count for reference
    HRLModelHelper baseline(1, 1);
    baseline.SetDataArray(m_data);
    baseline.SetMaximumIteration(50);
    baseline.SetTolerance(0.0);
    baseline.RunEstimation(0.5, 0.0);
    double beta_baseline{baseline.GetBetaMatrixMDPDE(0)(0)};

    // Run with single iteration
    HRLModelHelper limited1(1, 1);
    limited1.SetDataArray(m_data);
    limited1.SetMaximumIteration(1);
    limited1.SetTolerance(0.0);
    testing::internal::CaptureStderr();
    limited1.RunEstimation(0.5, 0.0);
    const std::string log1{ testing::internal::GetCapturedStderr() };
    EXPECT_NE(std::string::npos, log1.find("Reach maximum iterations"));
    double beta1{limited1.GetBetaMatrixMDPDE(0)(0)};

    // Run with five iterations
    HRLModelHelper limited5(1, 1);
    limited5.SetDataArray(m_data);
    limited5.SetMaximumIteration(5);
    limited5.SetTolerance(0.0);
    limited5.RunEstimation(0.5, 0.0);
    double beta5{ limited5.GetBetaMatrixMDPDE(0)(0) };

    // Verify more iterations move result toward baseline
    EXPECT_GT(std::abs(beta1 - beta_baseline), std::abs(beta5 - beta_baseline));
}

TEST_F(HRLModelHelperFixture, RespectsToleranceSetting)
{
    // Reference solution with strict tolerance and many iterations
    HRLModelHelper baseline(1, 1);
    baseline.SetDataArray(m_data);
    baseline.SetMaximumIteration(100);
    baseline.SetTolerance(1.0e-12);
    baseline.RunEstimation(0.5, 0.0);
    double beta_full{ baseline.GetBetaMatrixMDPDE(0)(0) };

    // Strict tolerance with limited iterations should hit maximum
    HRLModelHelper strict(1, 1);
    strict.SetDataArray(m_data);
    strict.SetMaximumIteration(5);
    strict.SetTolerance(1.0e-12);
    testing::internal::CaptureStderr();
    strict.RunEstimation(0.5, 0.0);
    const std::string log_strict{ testing::internal::GetCapturedStderr() };
    EXPECT_NE(std::string::npos, log_strict.find("Reach maximum iterations"));
    double beta_strict{ strict.GetBetaMatrixMDPDE(0)(0) };

    // Relaxed tolerance should converge early and avoid warning
    HRLModelHelper relaxed(1, 1);
    relaxed.SetDataArray(m_data);
    relaxed.SetMaximumIteration(5);
    relaxed.SetTolerance(1.0e-1);
    testing::internal::CaptureStderr();
    relaxed.RunEstimation(0.5, 0.0);
    const std::string log_relaxed{ testing::internal::GetCapturedStderr() };
    EXPECT_EQ(std::string::npos, log_relaxed.find("Reach maximum iterations"));
    double beta_relaxed{ relaxed.GetBetaMatrixMDPDE(0)(0) };

    // Relaxed tolerance stops early, yielding a result further from baseline
    EXPECT_LT(std::abs(beta_strict - beta_full), std::abs(beta_relaxed - beta_full));
}
