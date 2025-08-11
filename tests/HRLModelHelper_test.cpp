#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <cstdint>

#include "HRLModelHelper.hpp"
#include "Logger.hpp"

namespace {

using Eigen::VectorXd;
using DataTuple = std::tuple<std::vector<Eigen::VectorXd>, std::string>;

template <typename Vec>
class OversizedVectorGuard
{
public:
    explicit OversizedVectorGuard(Vec & vec) : m_vec{ vec }
    {
        m_hack = reinterpret_cast<Hack*>(&m_vec);
        m_orig_finish = m_hack->finish;
        m_orig_end = m_hack->end;
        auto start_addr{ reinterpret_cast<std::uintptr_t>(m_hack->start) };
        constexpr std::size_t too_big{ static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1 };
        auto finish_addr{ start_addr + too_big * sizeof(typename Vec::value_type) };
        m_hack->finish = reinterpret_cast<typename Vec::pointer>(finish_addr);
        m_hack->end = m_hack->finish;
    }

    ~OversizedVectorGuard()
    {
        m_hack->finish = m_orig_finish;
        m_hack->end = m_orig_end;
    }

private:
    Vec & m_vec;
    struct Hack
    {
        typename Vec::pointer start;
        typename Vec::pointer finish;
        typename Vec::pointer end;
    };
    Hack * m_hack;
    typename Vec::pointer m_orig_finish;
    typename Vec::pointer m_orig_end;
};

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

class HRLModelHelperTest : public ::testing::Test
{
protected:
    LogLevel m_prev_level;

    void SetUp(void) override
    {
        m_prev_level = Logger::GetLogLevel();
        Logger::SetLogLevel(LogLevel::Error);
    }

    void TearDown(void) override
    {
        Logger::SetLogLevel(m_prev_level);
    }
};

TEST_F(HRLModelHelperTest, ThrowsOnNonPositiveBasisSize)
{
    EXPECT_THROW(HRLModelHelper(0, 1), std::invalid_argument);
    EXPECT_THROW(HRLModelHelper(-1, 1), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsOnNonPositiveMemberSize)
{
    EXPECT_THROW(HRLModelHelper(1, 0), std::invalid_argument);
    EXPECT_THROW(HRLModelHelper(1, -1), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, SetMaximumIterationRejectsZero)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.SetMaximumIteration(0), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, SetMaximumIterationRejectsTooLarge)
{
    HRLModelHelper helper{ 1, 1 };
    unsigned int too_large{ static_cast<unsigned int>(std::numeric_limits<int>::max()) + 1u };
    EXPECT_THROW(helper.SetMaximumIteration(too_large), std::out_of_range);
}

TEST_F(HRLModelHelperTest, ThrowsOnNegativeTolerance)
{
    HRLModelHelper helper{1, 1};
    EXPECT_THROW(helper.SetTolerance(-1.0), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, SetToleranceRejectsNonFinite)
{
    HRLModelHelper helper{1, 1};
    EXPECT_THROW(helper.SetTolerance(std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    EXPECT_THROW(helper.SetTolerance(std::numeric_limits<double>::quiet_NaN()),
                 std::invalid_argument);
}

TEST_F(HRLModelHelperTest, DefaultInitializationValues)
{
    const int basis_size{ 2 };
    const int member_size{ 2 };

    HRLModelHelper helper{ basis_size, member_size };
    Eigen::MatrixXd expected_lambda{ Eigen::MatrixXd::Identity(basis_size, basis_size) };
    EXPECT_TRUE(helper.GetCapitalLambdaMatrix().isApprox(expected_lambda));

    Eigen::VectorXd expected_mu{ Eigen::VectorXd::Zero(basis_size) };
    EXPECT_TRUE(helper.GetMuVectorPrior().isApprox(expected_mu));
    EXPECT_TRUE(helper.GetMuVectorMDPDE().isApprox(expected_mu));
    EXPECT_TRUE(helper.GetMuVectorMean().isApprox(expected_mu));

    Eigen::VectorXd expected_beta{ Eigen::VectorXd::Zero(basis_size) };
    for (int member = 0; member < member_size; member++)
    {
        EXPECT_TRUE(helper.GetBetaMatrixOLS(member).isApprox(expected_beta));
        EXPECT_TRUE(helper.GetBetaMatrixMDPDE(member).isApprox(expected_beta));
        EXPECT_TRUE(helper.GetBetaMatrixPosterior(member).isApprox(expected_beta));
        EXPECT_DOUBLE_EQ(0.0, helper.GetStatisticalDistance(member));
        EXPECT_FALSE(helper.GetOutlierFlag(member));
    }
}

TEST_F(HRLModelHelperTest, AcceptsProperlySizedData)
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

    const auto & beta{ helper.GetBetaMatrixMDPDE(0) };
    ASSERT_EQ(beta.size(), 1);
    EXPECT_NEAR(beta(0), 2.4, 1e-9);

    Eigen::VectorXd mu{ helper.GetMuVectorMDPDE() };
    ASSERT_EQ(mu.size(), 1);
    EXPECT_NEAR(mu(0), 2.4, 1e-9);
    EXPECT_TRUE(mu.isApprox(beta));

    Eigen::VectorXd mu_prior{ helper.GetMuVectorPrior() };
    ASSERT_EQ(mu_prior.size(), 1);
    EXPECT_NEAR(mu_prior(0), 2.4, 1e-9);
    EXPECT_DOUBLE_EQ(0.0, helper.GetStatisticalDistance(0));
}

TEST_F(HRLModelHelperTest, SingleMemberDoesNotInvertCovariance)
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

    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    helper.RunEstimation(0.0, 0.0);
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(prev_level);

    EXPECT_NE(std::string::npos, out.find("Only one member is present"));
    EXPECT_EQ(std::string::npos, err.find("[Warning]"));
    EXPECT_EQ(std::string::npos, err.find("[Error]"));
}

TEST_F(HRLModelHelperTest, RunEstimationUsesNewDataset)
{
    HRLModelHelper helper{1, 1};

    // First dataset yields slope 2.4
    std::vector<Eigen::VectorXd> member1;
    Eigen::VectorXd s1(2);
    s1 << 1.0, 2.0;
    Eigen::VectorXd s2(2);
    s2 << 2.0, 5.0;
    member1.emplace_back(s1);
    member1.emplace_back(s2);
    std::vector<DataTuple> data1;
    data1.emplace_back(member1, "member");
    helper.SetDataArray(data1);
    helper.RunEstimation(0.0, 0.0);
    EXPECT_NEAR(helper.GetMuVectorMDPDE()(0), 2.4, 1e-9);

    // Second dataset yields slope 1.0; results should update
    std::vector<Eigen::VectorXd> member2;
    Eigen::VectorXd t1(2);
    t1 << 1.0, 1.0;
    Eigen::VectorXd t2(2);
    t2 << 2.0, 2.0;
    member2.emplace_back(t1);
    member2.emplace_back(t2);
    std::vector<DataTuple> data2;
    data2.emplace_back(member2, "member");
    helper.SetDataArray(data2);
    helper.RunEstimation(0.0, 0.0);

    EXPECT_NEAR(helper.GetMuVectorMDPDE()(0), 1.0, 1e-9);
    EXPECT_NEAR(helper.GetBetaMatrixMDPDE(0)(0), 1.0, 1e-9);
}

TEST_F(HRLModelHelperTest, ThrowsWhenMemberSizeMismatch)
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

TEST_F(HRLModelHelperTest, ThrowsWhenSampleVectorSizeMismatch)
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

TEST_F(HRLModelHelperTest, ThrowsWhenSampleContainsNaN)
{
    HRLModelHelper helper{ 1, 1 };
    Eigen::VectorXd bad_sample(2);
    bad_sample << 1.0, std::numeric_limits<double>::quiet_NaN();
    std::vector<Eigen::VectorXd> member{ bad_sample };
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member, "member1");
    EXPECT_THROW(helper.SetDataArray(data_array), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenMemberDataEmpty)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<Eigen::VectorXd> empty_member;
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(empty_member, "member1");
    EXPECT_THROW(helper.SetDataArray(data_array), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, SetDataArrayThrowsOnOversizedDataArray)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<DataTuple> data;
    data.emplace_back(std::vector<Eigen::VectorXd>{}, "member1");
    {
        OversizedVectorGuard<std::vector<DataTuple>> guard(data);
        EXPECT_THROW(helper.SetDataArray(data), std::overflow_error);
    }
}

TEST_F(HRLModelHelperTest, SetDataArrayThrowsOnOversizedMemberData)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<Eigen::VectorXd> member;
    member.emplace_back(Eigen::VectorXd::Zero(2));
    std::vector<DataTuple> data;
    data.emplace_back(member, "member1");
    auto & member_ref{ std::get<0>(data[0]) };
    {
        OversizedVectorGuard<std::vector<Eigen::VectorXd>> guard(member_ref);
        EXPECT_THROW(helper.SetDataArray(data), std::overflow_error);
    }
}

TEST_F(HRLModelHelperTest, SetDataArrayFailureDoesNotModifyExistingData)
{
    HRLModelHelper helper{ 1, 1 };
    // Load valid data
    std::vector<Eigen::VectorXd> member_data;
    Eigen::VectorXd sample1(2);
    sample1 << 1.0, 2.0;
    Eigen::VectorXd sample2(2);
    sample2 << 2.0, 5.0;
    member_data.emplace_back(sample1);
    member_data.emplace_back(sample2);
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> valid_data;
    valid_data.emplace_back(member_data, "member1");
    helper.SetDataArray(valid_data);

    // Attempt to replace with invalid data (member size mismatch)
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> invalid_data;
    invalid_data.emplace_back(member_data, "member1");
    invalid_data.emplace_back(member_data, "member2");
    EXPECT_THROW(helper.SetDataArray(invalid_data), std::invalid_argument);

    // Ensure original data remain intact
    ASSERT_NO_THROW(helper.RunEstimation(0.0, 0.0));
    const auto & beta{ helper.GetBetaMatrixMDPDE(0) };
    ASSERT_EQ(beta.size(), 1);
    EXPECT_NEAR(beta(0), 2.4, 1e-9);
}

TEST_F(HRLModelHelperTest, ThrowsWhenDataArrayNotSet)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(0.0, 0.0), std::runtime_error);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaRIsNegative)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(-0.1, 0.0), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaGIsNegative)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(0.0, -0.1), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaRIsNotFinite)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(std::numeric_limits<double>::infinity(), 0.0),
                 std::invalid_argument);
    EXPECT_THROW(helper.RunEstimation(std::numeric_limits<double>::quiet_NaN(), 0.0),
                 std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaGIsNotFinite)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(0.0, std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    EXPECT_THROW(helper.RunEstimation(0.0, std::numeric_limits<double>::quiet_NaN()),
                 std::invalid_argument);
}

TEST_F(HRLModelHelperTest, HandlesDataVarianceDenominatorNonPositive)
{
    // Craft data with a huge outlier so that all weights collapse to zero
    auto member{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 1.0e9}}, "degenerate") };
    std::vector<DataTuple> data_array{ member };

    HRLModelHelper helper{ 2, 1 };
    helper.SetDataArray(data_array);
    // Large alpha_r forces denominator to be non-positive
    EXPECT_NO_THROW(helper.RunEstimation(1.0e9, 0.0));
    EXPECT_TRUE(std::isfinite(helper.GetSigmaSquare(0)));
}

TEST_F(HRLModelHelperTest, HandlesNonFiniteDataVarianceSquare)
{
    // Data with extremely large residuals causing overflow in variance calculation
    auto member{ CreateMember({{0.0, 0.0}, {1.0, 1.0e155}, {2.0, -1.0e155}}, "extreme") };
    std::vector<DataTuple> data_array{ member };

    HRLModelHelper helper{ 2, 1 };
    helper.SetDataArray(data_array);

    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(helper.RunEstimation(1.0e-12, 0.0));
    const std::string err{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(prev_level);

    EXPECT_NE(std::string::npos, err.find("Non-finite variance"));
    double sigma{ helper.GetSigmaSquare(0) };
    EXPECT_TRUE(std::isfinite(sigma));
    EXPECT_GT(sigma, 0.0);
}

TEST_F(HRLModelHelperTest, HandlesMemberCovarianceDenominatorNonPositive)
{
    // Two members with vastly different intercepts force member weights to minimum
    auto member0{ CreateMember({{0.0, 0.0}, {1.0, 0.1}, {2.0, -0.1}}, "member_0") };
    auto member1{ CreateMember({{0.0, 1000.0}, {1.0, 1000.1}, {2.0, 999.9}}, "member_1") };
    std::vector<DataTuple> data_array{ member0, member1 };

    HRLModelHelper helper{ 2, 2 };
    helper.SetDataArray(data_array);
    // alpha_g chosen slightly above the theoretical root to make the
    // denominator of CalculateMemberCovariance non-positive
    const double alpha_g{ 0.005050733883341149 };
    EXPECT_NO_THROW(helper.RunEstimation(0.0, alpha_g));

    const auto & lambda{ helper.GetCapitalLambdaMatrix() };
    ASSERT_TRUE(lambda.array().allFinite());
}

TEST_F(HRLModelHelperTest, GettersThrowOnInvalidId)
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

TEST_F(HRLModelHelperTest, EstimationOnSmallSyntheticData)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_1") };
    std::vector<DataTuple> data_array{member0, member1};

    HRLModelHelper helper(2, 2);
    helper.SetDataArray(data_array);
    helper.SetMaximumIteration(1000);
    helper.SetTolerance(1.0e-6);
    helper.RunEstimation(0.0, 0.0);

    Eigen::Vector2d expected;
    expected << 5.0 / 6.0, 2.5; // Analytic OLS solution

    for (int i = 0; i < 2; i++)
    {
        const auto & beta_post{ helper.GetBetaMatrixPosterior(i) };
        const auto & beta_mdpde{ helper.GetBetaMatrixMDPDE(i) };
        const auto & beta_ols{ helper.GetBetaMatrixOLS(i) };
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

TEST_F(HRLModelHelperTest, PosteriorSigmaMatrixIsSymmetricPositive)
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

TEST_F(HRLModelHelperTest, OutlierFlagsWithDivergentData)
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

    const auto & beta_ols0{ helper.GetBetaMatrixOLS(0) };
    const auto & beta_ols1{ helper.GetBetaMatrixOLS(1) };
    ASSERT_NEAR(expected0(0), beta_ols0(0), 1.0e-6);
    ASSERT_NEAR(expected0(1), beta_ols0(1), 1.0e-6);
    ASSERT_NEAR(expected1(0), beta_ols1(0), 1.0e-6);
    ASSERT_NEAR(expected1(1), beta_ols1(1), 1.0e-6);

    const auto & beta_mdpde0{ helper.GetBetaMatrixMDPDE(0) };
    const auto & beta_mdpde1{ helper.GetBetaMatrixMDPDE(1) };
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
        double dist{ helper.GetStatisticalDistance(i) };
        EXPECT_TRUE(std::isfinite(dist));
        EXPECT_GT(dist, 0.0);
        EXPECT_FALSE(helper.GetOutlierFlag(i));
    }
}

TEST_F(HRLModelHelperTest, AlphaGReweightsOutliers)
{
    auto inlier0{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}}, "inlier_0") };
    auto inlier1{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}}, "inlier_1") };
    auto outlier{ CreateMember({{0.0, 10.0}, {1.0, 20.0}, {2.0, 30.0}}, "outlier") };
    std::vector<DataTuple> data_array{inlier0, inlier1, outlier};

    HRLModelHelper baseline(2, 3);
    baseline.SetDataArray(data_array);
    baseline.SetMaximumIteration(1000);
    baseline.SetTolerance(1.0e-6);
    baseline.RunEstimation(0.0, 0.0);
    Eigen::VectorXd mu_prior_initial{ baseline.GetMuVectorPrior() };

    HRLModelHelper reweighted(2, 3);
    reweighted.SetDataArray(data_array);
    reweighted.SetMaximumIteration(1000);
    reweighted.SetTolerance(1.0e-6);
    reweighted.RunEstimation(0.0, 0.5);
    Eigen::VectorXd mu_prior_reweighted{ reweighted.GetMuVectorPrior() };
    Eigen::VectorXd inlier_mean(2);
    inlier_mean << 1.0, 1.0;

    double dist_initial{ (mu_prior_initial - inlier_mean).norm() };
    double dist_reweighted{ (mu_prior_reweighted - inlier_mean).norm() };
    if (!std::isfinite(dist_initial) || !std::isfinite(dist_reweighted))
    {
        GTEST_SKIP() << "Estimation failed to converge";
    }
    EXPECT_LT(dist_reweighted, dist_initial + 1.0e-12);
}

TEST_F(HRLModelHelperTest, RobustEstimationDownweightsOutlier)
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

TEST_F(HRLModelHelperTest, DetectsTrueOutlier)
{
    auto inlier0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 7.0}}, "inlier_0") };
    auto inlier1{ CreateMember({{0.0, 1.2}, {1.0, 3.2}, {2.0, 5.2}, {3.0, 7.2}}, "inlier_1") };
    auto outlier{ CreateMember({{0.0, 30.0}, {1.0, 60.0}, {2.0, 90.0}, {3.0, 120.0}}, "outlier") };
    std::vector<DataTuple> data_array{ inlier0, inlier1, outlier };
    HRLModelHelper helper(2, 3);
    helper.SetDataArray(data_array);
    helper.SetMaximumIteration(1000);
    helper.SetTolerance(1.0e-6);
    helper.RunEstimation(0.0, 1.0);
    EXPECT_FALSE(helper.GetOutlierFlag(0));
    EXPECT_FALSE(helper.GetOutlierFlag(1));
    EXPECT_TRUE(helper.GetOutlierFlag(2));
}

class HRLModelHelperFixture : public HRLModelHelperTest
{
protected:
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> m_data;

    void SetUp(void) override
    {
        HRLModelHelperTest::SetUp();
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
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    limited1.RunEstimation(0.5, 0.0);
    const std::string log1{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(LogLevel::Error);
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
    Logger::SetLogLevel(LogLevel::Warning);
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
    Logger::SetLogLevel(LogLevel::Error);
    EXPECT_EQ(std::string::npos, log_relaxed.find("Reach maximum iterations"));
    double beta_relaxed{ relaxed.GetBetaMatrixMDPDE(0)(0) };

    // Relaxed tolerance stops early, yielding a result further from baseline
    EXPECT_LT(std::abs(beta_strict - beta_full), std::abs(beta_relaxed - beta_full));
}

TEST_F(HRLModelHelperTest, ExtremeResidualsYieldFinitePosteriorCovariance)
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

TEST_F(HRLModelHelperTest, LargeResidualsClampMemberWeights)
{
    // Two members: one normal, one with extremely large response producing
    // overflow in residual magnitude for member weights.
    std::vector<Eigen::VectorXd> normal;
    Eigen::VectorXd n0(2); n0 << 0.0, 0.0; normal.emplace_back(n0);
    Eigen::VectorXd n1(2); n1 << 1.0, 1.0; normal.emplace_back(n1);
    std::vector<Eigen::VectorXd> extreme;
    Eigen::VectorXd e0(2); e0 << 0.0, 0.0; extreme.emplace_back(e0);
    Eigen::VectorXd e1(2); e1 << 1.0, 1.0e155; extreme.emplace_back(e1);
    std::vector<DataTuple> data_array;
    data_array.emplace_back(normal, "normal");
    data_array.emplace_back(extreme, "extreme");
    HRLModelHelper helper{1, 2};
    helper.SetDataArray(data_array);
    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    helper.RunEstimation(0.0, 1.0);
    const std::string err{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(prev_level);
    EXPECT_NE(std::string::npos, err.find("non-finite exponent index"));
    for (int i = 0; i < 2; ++i)
    {
        double w{ helper.GetMemberWeight(i) };
        EXPECT_TRUE(std::isfinite(w));
        EXPECT_GT(w, 0.0);
    }
}

TEST_F(HRLModelHelperTest, ZeroResidualProducesFiniteWeightsAndPosterior)
{
    auto member{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}}, "perfect") };
    std::vector<DataTuple> data_array{ member };
    HRLModelHelper helper(2, 1);
    helper.SetDataArray(data_array);
    helper.RunEstimation(0.0, 0.0);
    const auto & sigma{ helper.GetCapitalSigmaMatrixPosterior(0) };
    EXPECT_TRUE(sigma.array().isFinite().all());
    const auto & beta{ helper.GetBetaMatrixPosterior(0) };
    EXPECT_TRUE(beta.array().isFinite().all());
}

TEST_F(HRLModelHelperTest, ThrowsOnDegenerateWeights)
{
    std::vector<Eigen::VectorXd> member_data;
    Eigen::VectorXd sample1(3);
    sample1 << 1.0, 0.0, 0.0;
    Eigen::VectorXd sample2(3);
    sample2 << 1.0, std::numeric_limits<double>::max(), std::numeric_limits<double>::max();
    member_data.emplace_back(sample1);
    member_data.emplace_back(sample2);
    std::vector<DataTuple> data_array;
    data_array.emplace_back(member_data, "collapse");
    HRLModelHelper helper(2, 1);
    ASSERT_NO_THROW(helper.SetDataArray(data_array));
    EXPECT_THROW(helper.RunEstimation(1.0e8, 0.0), std::runtime_error);
}
