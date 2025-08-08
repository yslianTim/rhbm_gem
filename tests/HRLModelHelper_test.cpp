#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>

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

TEST(HRLModelHelperTest, ThrowsOnNonPositiveBasisSize)
{
    EXPECT_THROW(HRLModelHelper(0, 1), std::invalid_argument);
    EXPECT_THROW(HRLModelHelper(-1, 1), std::invalid_argument);
}

TEST(HRLModelHelperTest, ThrowsOnNonPositiveMemberSize)
{
    EXPECT_THROW(HRLModelHelper(1, 0), std::invalid_argument);
    EXPECT_THROW(HRLModelHelper(1, -1), std::invalid_argument);
}

TEST(HRLModelHelperTest, SetMaximumIterationRejectsZero)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.SetMaximumIteration(0), std::invalid_argument);
}

TEST(HRLModelHelperTest, SetMaximumIterationRejectsTooLarge)
{
    HRLModelHelper helper{ 1, 1 };
    unsigned int too_large{ static_cast<unsigned int>(std::numeric_limits<int>::max()) + 1u };
    EXPECT_THROW(helper.SetMaximumIteration(too_large), std::out_of_range);
}

TEST(HRLModelHelperTest, ThrowsOnNegativeTolerance)
{
    HRLModelHelper helper{1, 1};
    EXPECT_THROW(helper.SetTolerance(-1.0), std::invalid_argument);
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
    EXPECT_TRUE(helper.GetCapitalLambdaMatrix().isApprox(Eigen::MatrixXd::Zero(1, 1)));

    Eigen::VectorXd beta{ helper.GetBetaMatrixMDPDE(0) };
    ASSERT_EQ(beta.size(), 1);
    EXPECT_NEAR(beta(0), 2.4, 1e-9);

    Eigen::VectorXd mu{ helper.GetMuVectorMDPDE() };
    ASSERT_EQ(mu.size(), 1);
    EXPECT_NEAR(mu(0), 2.4, 1e-9);
    EXPECT_TRUE(mu.isApprox(beta));

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

    EXPECT_THROW(helper.SetDataArray(data_array), std::invalid_argument);
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
    EXPECT_THROW(helper.SetDataArray(data_array), std::invalid_argument);
}

TEST(HRLModelHelperTest, DISABLED_ThrowsWhenSampleContainsNaN)
{
    HRLModelHelper helper{ 1, 1 };
    Eigen::VectorXd bad_sample(2);
    bad_sample << 1.0, std::numeric_limits<double>::quiet_NaN();
    std::vector<Eigen::VectorXd> member{ bad_sample };
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member, "member1");
    EXPECT_THROW(helper.SetDataArray(data_array), std::invalid_argument);
}

TEST(HRLModelHelperTest, ThrowsWhenMemberDataEmpty)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<Eigen::VectorXd> empty_member;
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(empty_member, "member1");
    EXPECT_THROW(helper.SetDataArray(data_array), std::invalid_argument);
}

TEST(HRLModelHelperTest, ThrowsWhenDataArrayNotSet)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(0.0, 0.0), std::runtime_error);
}

TEST(HRLModelHelperTest, ThrowsWhenAlphaRIsNegative)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(-0.1, 0.0), std::invalid_argument);
}

TEST(HRLModelHelperTest, ThrowsWhenAlphaGIsNegative)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(0.0, -0.1), std::invalid_argument);
}

TEST(HRLModelHelperTest, ThrowsWhenAlphaRIsNotFinite)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(std::numeric_limits<double>::infinity(), 0.0),
                 std::invalid_argument);
    EXPECT_THROW(helper.RunEstimation(std::numeric_limits<double>::quiet_NaN(), 0.0),
                 std::invalid_argument);
}

TEST(HRLModelHelperTest, DISABLED_ThrowsWhenDataVarianceDenominatorNonPositive)
{
    // Craft data with a huge outlier so that all weights collapse to zero
    auto member{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 1.0e9}}, "degenerate") };
    std::vector<DataTuple> data_array{ member };

    HRLModelHelper helper{ 2, 1 };
    helper.SetDataArray(data_array);
    // Large alpha_r forces denominator to be non-positive
    EXPECT_THROW(helper.RunEstimation(1.0e9, 0.0), std::runtime_error);
}

TEST(HRLModelHelperTest, GettersThrowOnInvalidId)
{
    const int basis_size{ 1 };
    const int member_size{ 1 };
    HRLModelHelper helper{ basis_size, member_size };

    const int invalid_id{ member_size };
    EXPECT_THROW(helper.GetOutlierFlag(invalid_id), std::out_of_range);
    EXPECT_THROW(helper.GetStatisticalDistance(invalid_id), std::out_of_range);
    EXPECT_THROW(helper.GetCapitalSigmaMatrixPosterior(invalid_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixPosterior(invalid_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixMDPDE(invalid_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixOLS(invalid_id), std::out_of_range);

    const int negative_id{ -1 };
    EXPECT_THROW(helper.GetOutlierFlag(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetStatisticalDistance(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetCapitalSigmaMatrixPosterior(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixPosterior(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixMDPDE(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixOLS(negative_id), std::out_of_range);
}

TEST(HRLModelHelperTest, EstimationOnSmallSyntheticData)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_1") };
    std::vector<DataTuple> data_array{member0, member1};

    HRLModelHelper helper(2, 2);
    helper.SetDataArray(data_array);
    helper.SetMaximumIteration(500);
    helper.SetTolerance(0.0);
    helper.SetMaximumIteration(500);
    helper.SetTolerance(0.0);
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

TEST(HRLModelHelperTest, PosteriorSigmaMatrixIsSymmetricPositive)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_1") };
    std::vector<DataTuple> data_array{ member0, member1 };
    const int basis_size{ 2 };
    const int member_size{ 2 };
    HRLModelHelper helper(basis_size, member_size);
    helper.SetDataArray(data_array);
    helper.RunEstimation(0.0, 0.0);
    for (int i = 0; i < member_size; ++i)
    {
        const auto & sigma{ helper.GetCapitalSigmaMatrixPosterior(i) };
        ASSERT_EQ(basis_size, sigma.rows());
        ASSERT_EQ(basis_size, sigma.cols());
        EXPECT_TRUE(sigma.isApprox(sigma.transpose()));
        for (int j = 0; j < basis_size; ++j)
        {
            EXPECT_GT(sigma(j, j), 0.0);
        }
    }
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

TEST(HRLModelHelperTest, AlphaGReweightsOutliers)
{
    auto inlier0{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}}, "inlier_0") };
    auto inlier1{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}}, "inlier_1") };
    auto outlier{ CreateMember({{0.0, 10.0}, {1.0, 20.0}, {2.0, 30.0}}, "outlier") };
    std::vector<DataTuple> data_array{inlier0, inlier1, outlier};

    HRLModelHelper helper(2, 3);
    helper.SetDataArray(data_array);
    helper.RunEstimation(0.0, 0.0);
    Eigen::VectorXd mu_prior_initial{ helper.GetMuVectorPrior() };
    helper.RunEstimation(0.0, 0.5);
    Eigen::VectorXd mu_prior_reweighted{ helper.GetMuVectorPrior() };
    Eigen::VectorXd inlier_mean(2);
    inlier_mean << 1.0, 1.0;

    double dist_initial{ (mu_prior_initial - inlier_mean).norm() };
    double dist_reweighted{ (mu_prior_reweighted - inlier_mean).norm() };
    EXPECT_LT(dist_reweighted, dist_initial);
}

TEST(HRLModelHelperTest, RobustEstimationDownweightsOutlier)
{
    // Majority member follows a clean linear trend y = 2x + 1
    auto majority{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 7.0}}, "majority") };
    // Second member contains an outlier at the last point
    auto outlier{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 20.0}}, "outlier") };
    std::vector<DataTuple> data_array{ majority, outlier };

    HRLModelHelper helper(2, 2);
    helper.SetDataArray(data_array);

    // First run with no robustness (alpha = 0)
    helper.RunEstimation(0.0, 0.0);
    Eigen::Vector2d majority_beta{ helper.GetBetaMatrixMDPDE(0) };
    Eigen::Vector2d outlier_beta_ols{ helper.GetBetaMatrixMDPDE(1) };

    // Second run with robustness enabled (alpha = 0.5)
    helper.RunEstimation(0.5, 0.5);
    Eigen::Vector2d outlier_beta_robust{ helper.GetBetaMatrixMDPDE(1) };

    // The robust run should yield coefficients closer to the majority trend
    double diff_ols{ (outlier_beta_ols - majority_beta).norm() };
    double diff_robust{ (outlier_beta_robust - majority_beta).norm() };
    EXPECT_LT(diff_robust, diff_ols);
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

TEST(HRLModelHelperTest, ExtremeResidualsYieldFinitePosteriorCovariance)
{
    Eigen::VectorXd sample1(2);
    sample1 << 1.0, 0.0;
    Eigen::VectorXd sample2(2);
    sample2 << 1.0, 1.0e8;
    std::vector<Eigen::VectorXd> member{ sample1, sample2 };
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member, "extreme");

    HRLModelHelper helper{ 1, 1 };
    helper.SetDataArray(data_array);
    helper.RunEstimation(1.0e20, 0.0);

    const auto & sigma{ helper.GetCapitalSigmaMatrixPosterior(0) };
    ASSERT_EQ(1, sigma.rows());
    ASSERT_EQ(1, sigma.cols());
    double val{ sigma(0, 0) };
    EXPECT_TRUE(std::isfinite(val));
    EXPECT_GT(val, 0.0);
}
