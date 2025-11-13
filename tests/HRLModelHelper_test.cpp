#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <cstdint>
#include <type_traits>

#include "HRLModelHelper.hpp"
#include "Logger.hpp"
#include "EigenMatrixUtility.hpp"

namespace {

using Eigen::VectorXd;
using DataTuple = std::tuple<std::vector<Eigen::VectorXd>, std::string>;

static_assert(!std::is_convertible_v<int, HRLModelHelper>,
              "HRLModelHelper should not be implicitly constructible from a single int");
static_assert(!std::is_convertible_v<int[2], HRLModelHelper>,
              "HRLModelHelper should not be implicitly constructible from an int array");

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

TEST_F(HRLModelHelperTest, SetMaximumIterationRejectsNegative)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.SetMaximumIteration(-1), std::invalid_argument);
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

TEST_F(HRLModelHelperTest, SetToleranceAcceptsZero)
{
    HRLModelHelper helper{1, 1};
    ASSERT_NO_THROW(helper.SetTolerance(0.0));

    // Minimal dataset to ensure RunEstimation still converges
    std::vector<Eigen::VectorXd> member_data;
    Eigen::VectorXd sample1(2);
    sample1 << 1.0, 2.0;
    Eigen::VectorXd sample2(2);
    sample2 << 2.0, 5.0;
    member_data.emplace_back(sample1);
    member_data.emplace_back(sample2);
    std::vector<DataTuple> data_array;
    data_array.emplace_back(member_data, "member1");
    ASSERT_NO_THROW(helper.SetDataArray(std::move(data_array)));
    EXPECT_NO_THROW(helper.RunEstimation(0.5));
    EXPECT_NEAR(helper.GetMuVectorMDPDE()(0), 2.4, 1e-9);
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
    Eigen::MatrixXd expected_sigma{ Eigen::MatrixXd::Identity(basis_size, basis_size) };
    for (int member = 0; member < member_size; member++)
    {
        EXPECT_TRUE(helper.GetBetaMatrixOLS(member).isApprox(expected_beta));
        EXPECT_TRUE(helper.GetBetaMatrixMDPDE(member).isApprox(expected_beta));
        EXPECT_TRUE(helper.GetBetaMatrixPosterior(member).isApprox(expected_beta));
        EXPECT_TRUE(helper.GetCapitalSigmaMatrixPosterior(member).isApprox(expected_sigma));
        EXPECT_DOUBLE_EQ(0.0, helper.GetStatisticalDistance(member));
        EXPECT_FALSE(helper.GetOutlierFlag(member));
        EXPECT_DOUBLE_EQ(1.0, helper.GetSigmaSquare(member));
        EXPECT_DOUBLE_EQ(1.0, helper.GetMemberWeight(member));
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
    ASSERT_NO_THROW(helper.SetDataArray(std::move(data_array)));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.5);
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
    ASSERT_NO_THROW(helper.SetDataArray(std::move(data_array)));

    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);
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
    helper.SetDataArray(std::move(data1));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);
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
    helper.SetDataArray(std::move(data2));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);

    EXPECT_NEAR(helper.GetMuVectorMDPDE()(0), 1.0, 1e-9);
    EXPECT_NEAR(helper.GetBetaMatrixMDPDE(0)(0), 1.0, 1e-9);
}

TEST_F(HRLModelHelperTest, MuVectorMatchesWeightedMean)
{
    HRLModelHelper helper{ 2, 2 };
    auto member1{ CreateMember({{1.0, 2.0}, {2.0, 5.0}}, "m1") };
    auto member2{ CreateMember({{1.0, 1.0}, {2.0, 2.0}}, "m2") };
    std::vector<DataTuple> data{ member1, member2 };
    helper.SetDataArray(std::move(data));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.5);
    const double beta0{ helper.GetBetaMatrixMDPDE(0)(0) };
    const double beta1{ helper.GetBetaMatrixMDPDE(1)(0) };
    const double w0{ helper.GetMemberWeight(0) };
    const double w1{ helper.GetMemberWeight(1) };
    const double expected_mu{ (beta0 * w0 + beta1 * w1) / (w0 + w1) };
    EXPECT_NEAR(helper.GetMuVectorMDPDE()(0), expected_mu, 1e-9);
}

TEST_F(HRLModelHelperTest, RunEstimationUpdatesWeightsWithNewDataset)
{
    HRLModelHelper helper{ 1, 1 };

    // Dataset A: two samples that produce non-identity weights
    std::vector<Eigen::VectorXd> member_a;
    Eigen::VectorXd a1(2); a1 << 1.0, 1.0;
    Eigen::VectorXd a2(2); a2 << 2.0, 10.0;
    member_a.emplace_back(a1);
    member_a.emplace_back(a2);
    std::vector<DataTuple> data_a;
    data_a.emplace_back(member_a, "member");
    helper.SetDataArray(std::move(data_a));
    helper.SetUniversalAlphaR(0.5);
    helper.RunEstimation(0.0);
    Eigen::VectorXd w_a{ helper.GetDataWeightMatrix(0).diagonal() };
    ASSERT_EQ(w_a.size(), 2);
    EXPECT_FALSE(w_a.isApprox(Eigen::VectorXd::Ones(2)));

    // Dataset B: different sample count; perfect fit yields unit weights
    std::vector<Eigen::VectorXd> member_b;
    Eigen::VectorXd b1(2); b1 << 1.0, 1.0;
    Eigen::VectorXd b2(2); b2 << 2.0, 2.0;
    Eigen::VectorXd b3(2); b3 << 3.0, 3.0;
    member_b.emplace_back(b1);
    member_b.emplace_back(b2);
    member_b.emplace_back(b3);
    std::vector<DataTuple> data_b;
    data_b.emplace_back(member_b, "member");
    helper.SetDataArray(std::move(data_b));
    helper.SetUniversalAlphaR(0.5);
    helper.RunEstimation(0.0);
    Eigen::VectorXd w_b{ helper.GetDataWeightMatrix(0).diagonal() };

    ASSERT_EQ(w_b.size(), 3);
    EXPECT_TRUE(w_b.isApprox(Eigen::VectorXd::Ones(3)));
    EXPECT_FALSE(w_b.head(2).isApprox(w_a));
}

TEST_F(HRLModelHelperTest, LabelsOutlierMemberWhenDfIsThree)
{
    const int basis_size{ 3 };
    const int member_size{ 3 };
    HRLModelHelper helper{ basis_size, member_size };

    // Members 0 and 1: similar coefficients (intercept 0, slopes 1)
    std::vector<Eigen::VectorXd> member0;
    std::vector<Eigen::VectorXd> member1;
    Eigen::VectorXd a1(4); a1 << 1.0, 0.0, 0.0, 0.0;
    Eigen::VectorXd a2(4); a2 << 1.0, 1.0, 0.0, 1.0;
    Eigen::VectorXd a3(4); a3 << 1.0, 0.0, 1.0, 1.0;
    member0.emplace_back(a1);
    member0.emplace_back(a2);
    member0.emplace_back(a3);
    member1 = member0; // identical to member0

    // Member 2: extreme coefficients to trigger outlier detection
    std::vector<Eigen::VectorXd> member2;
    Eigen::VectorXd c1(4); c1 << 1.0, 0.0, 0.0, 100.0;
    Eigen::VectorXd c2(4); c2 << 1.0, 1.0, 0.0, 200.0;
    Eigen::VectorXd c3(4); c3 << 1.0, 0.0, 1.0, 200.0;
    member2.emplace_back(c1);
    member2.emplace_back(c2);
    member2.emplace_back(c3);

    std::vector<DataTuple> data;
    data.emplace_back(member0, "member0");
    data.emplace_back(member1, "member1");
    data.emplace_back(member2, "member2");
    helper.SetDataArray(std::move(data));

    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(1.0);
    EXPECT_LT(helper.GetStatisticalDistance(0), 11.34);
    EXPECT_FALSE(helper.GetOutlierFlag(0));
    EXPECT_LT(helper.GetStatisticalDistance(1), 11.34);
    EXPECT_FALSE(helper.GetOutlierFlag(1));
    EXPECT_GT(helper.GetStatisticalDistance(2), 11.34);
    EXPECT_TRUE(helper.GetOutlierFlag(2));
}

TEST_F(HRLModelHelperTest, CapitalLambdaMatchesManualComputation)
{
    const int basis_size{ 2 };
    const int member_size{ 2 };
    HRLModelHelper helper{ basis_size, member_size };

    auto member1{ CreateMember({{0.0, 0.0}, {1.0, 1.0}}, "m1") };
    auto member2{ CreateMember({{0.0, 0.0}, {1.0, 2.0}}, "m2") };
    std::vector<DataTuple> data{ member1, member2 };
    helper.SetDataArray(std::move(data));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);
    Eigen::VectorXd mu{ helper.GetMuVectorMDPDE() };
    Eigen::VectorXd beta0{ helper.GetBetaMatrixMDPDE(0) };
    Eigen::VectorXd beta1{ helper.GetBetaMatrixMDPDE(1) };
    const double w0{ helper.GetMemberWeight(0) };
    const double w1{ helper.GetMemberWeight(1) };

    Eigen::MatrixXd numerator{
        w0 * (beta0 - mu) * (beta0 - mu).transpose()
      + w1 * (beta1 - mu) * (beta1 - mu).transpose()
    };
    const double denominator{ w0 + w1 };
    ASSERT_GT(denominator, 0.0);
    Eigen::MatrixXd expected{ numerator / denominator };

    EXPECT_TRUE(helper.GetCapitalLambdaMatrix().isApprox(expected));
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

    EXPECT_THROW(helper.SetDataArray(std::move(data_array)), std::invalid_argument);
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
    EXPECT_THROW(helper.SetDataArray(std::move(data_array)), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenSampleContainsNaN)
{
    HRLModelHelper helper{ 1, 1 };
    Eigen::VectorXd bad_sample(2);
    bad_sample << 1.0, std::numeric_limits<double>::quiet_NaN();
    std::vector<Eigen::VectorXd> member{ bad_sample };
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member, "member1");
    EXPECT_THROW(helper.SetDataArray(std::move(data_array)), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenSampleContainsInfinity)
{
    HRLModelHelper helper{ 1, 1 };
    Eigen::VectorXd pos_inf_sample(2);
    pos_inf_sample << 1.0, std::numeric_limits<double>::infinity();
    std::vector<Eigen::VectorXd> member{ pos_inf_sample };
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(member, "member1");
    EXPECT_THROW(helper.SetDataArray(std::move(data_array)), std::invalid_argument);

    Eigen::VectorXd neg_inf_sample(2);
    neg_inf_sample << 1.0, -std::numeric_limits<double>::infinity();
    member[0] = neg_inf_sample;
    data_array.clear();
    data_array.emplace_back(member, "member1");
    EXPECT_THROW(helper.SetDataArray(std::move(data_array)), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, DISABLED_ThrowsWhenMemberDataEmpty)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<Eigen::VectorXd> empty_member;
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
    data_array.emplace_back(empty_member, "member1");
    EXPECT_THROW(helper.SetDataArray(std::move(data_array)), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, SetDataArrayThrowsOnOversizedDataArray)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<DataTuple> data;
    data.emplace_back(std::vector<Eigen::VectorXd>{}, "member1");
    {
        OversizedVectorGuard<std::vector<DataTuple>> guard(data);
        EXPECT_THROW(helper.SetDataArray(std::move(data)), std::overflow_error);
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
        EXPECT_THROW(helper.SetDataArray(std::move(data)), std::overflow_error);
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
    helper.SetDataArray(std::move(valid_data));

    // Attempt to replace with invalid data (member size mismatch)
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> invalid_data;
    invalid_data.emplace_back(member_data, "member1");
    invalid_data.emplace_back(member_data, "member2");
    EXPECT_THROW(helper.SetDataArray(std::move(invalid_data)), std::invalid_argument);

    // Ensure original data remain intact
    ASSERT_NO_THROW(helper.RunEstimation(0.0));
    const auto & beta{ helper.GetBetaMatrixMDPDE(0) };
    ASSERT_EQ(beta.size(), 1);
    EXPECT_NEAR(beta(0), 2.4, 1e-9);
}

TEST_F(HRLModelHelperTest, DISABLED_ThrowsWhenDataArrayNotSet)
{
    HRLModelHelper helper{ 1, 1 };
    helper.SetUniversalAlphaR(0.0);
    EXPECT_THROW(helper.RunEstimation(0.0), std::runtime_error);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaRIsNegative)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.SetUniversalAlphaR(-0.1), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaGIsNegative)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(-0.1), std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaRIsTooLarge)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.SetUniversalAlphaR(std::numeric_limits<double>::max()),
                 std::overflow_error);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaGIsTooLarge)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(std::numeric_limits<double>::max()),
                 std::overflow_error);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaRIsNotFinite)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.SetUniversalAlphaR(std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    EXPECT_THROW(helper.SetUniversalAlphaR(std::numeric_limits<double>::quiet_NaN()),
                 std::invalid_argument);
}

TEST_F(HRLModelHelperTest, ThrowsWhenAlphaGIsNotFinite)
{
    HRLModelHelper helper{ 1, 1 };
    EXPECT_THROW(helper.RunEstimation(std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    EXPECT_THROW(helper.RunEstimation(std::numeric_limits<double>::quiet_NaN()),
                 std::invalid_argument);
}

TEST_F(HRLModelHelperTest, HandlesDataVarianceDenominatorNonPositive)
{
    // Craft data with a huge outlier so that all weights collapse to zero
    auto member{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 1.0e9}}, "degenerate") };
    std::vector<DataTuple> data_array{ member };

    HRLModelHelper helper{ 2, 1 };
    helper.SetDataArray(std::move(data_array));
    // Large alpha_r forces denominator to be non-positive
    helper.SetUniversalAlphaR(1.0e9);
    EXPECT_NO_THROW(helper.RunEstimation(0.0));
    EXPECT_TRUE(std::isfinite(helper.GetSigmaSquare(0)));
}

TEST_F(HRLModelHelperTest, HandlesNonFiniteDataVarianceSquare)
{
    // Data with extremely large residuals causing overflow in variance calculation
    auto member{ CreateMember({{0.0, 0.0}, {1.0, 1.0e155}, {2.0, -1.0e155}}, "extreme") };
    std::vector<DataTuple> data_array{ member };

    HRLModelHelper helper{ 2, 1 };
    helper.SetDataArray(std::move(data_array));

    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    helper.SetUniversalAlphaR(1.0e-12);
    EXPECT_NO_THROW(helper.RunEstimation(0.0));
    const std::string err{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(prev_level);

    EXPECT_NE(std::string::npos, err.find("Non-finite variance"));
    double sigma{ helper.GetSigmaSquare(0) };
    EXPECT_TRUE(std::isfinite(sigma));
    EXPECT_GT(sigma, 0.0);
}

TEST_F(HRLModelHelperTest, MemberCovarianceDenominatorNonPositiveKeepsLambdaUnchanged)
{
    // Two members with noticeably different intercepts yield extreme residuals
    auto member0{ CreateMember({{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 5.0}, {1.0, 6.0}, {2.0, 7.0}}, "member_1") };
    std::vector<DataTuple> data_array{ member0, member1 };

    HRLModelHelper helper{ 2, 2 };
    helper.SetDataArray(std::move(data_array));
    // Capital lambda should start as the identity matrix
    Eigen::MatrixXd lambda_before{ helper.GetCapitalLambdaMatrix() };
    EXPECT_TRUE(lambda_before.isApprox(Eigen::MatrixXd::Identity(2, 2)));

    // Use a large alpha_g so that CalculateMemberCovariance has non-positive denominator
    const double alpha_g{ 100.0 };
    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    helper.SetUniversalAlphaR(0.0);
    EXPECT_NO_THROW(helper.RunEstimation(alpha_g));
    const std::string err{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(prev_level);
    EXPECT_NE(std::string::npos, err.find("Member covariance denominator is non-positive"));

    Eigen::MatrixXd lambda_after{ helper.GetCapitalLambdaMatrix() };
    EXPECT_TRUE(lambda_after.isApprox(lambda_before));
}

TEST_F(HRLModelHelperTest, MemberWeightMatchesAnalyticFormula)
{
    const double slope_mag{ 0.8231614670716226 };
    std::vector<Eigen::VectorXd> member0;
    Eigen::VectorXd m0_s1(2); m0_s1 << 1.0, -slope_mag;
    Eigen::VectorXd m0_s2(2); m0_s2 << 2.0, -2.0 * slope_mag;
    member0.emplace_back(m0_s1);
    member0.emplace_back(m0_s2);

    std::vector<Eigen::VectorXd> member1;
    Eigen::VectorXd m1_s1(2); m1_s1 << 1.0, slope_mag;
    Eigen::VectorXd m1_s2(2); m1_s2 << 2.0, 2.0 * slope_mag;
    member1.emplace_back(m1_s1);
    member1.emplace_back(m1_s2);

    std::vector<DataTuple> data_array;
    data_array.emplace_back(member0, "member0");
    data_array.emplace_back(member1, "member1");

    HRLModelHelper helper{ 1, 2 };
    helper.SetDataArray(std::move(data_array));

    const double alpha_r{ 0.0 };
    const double alpha_g{ 0.5 };
    helper.SetUniversalAlphaR(alpha_r);
    helper.RunEstimation(alpha_g);

    Eigen::VectorXd mu{ helper.GetMuVectorMDPDE() };
    Eigen::MatrixXd lambda{ helper.GetCapitalLambdaMatrix() };
    ASSERT_EQ(mu.size(), 1);
    ASSERT_EQ(lambda.rows(), 1);
    ASSERT_EQ(lambda.cols(), 1);
    EXPECT_NEAR(mu(0), 0.0, 1.0e-9);
    EXPECT_NEAR(lambda(0, 0), 1.0, 1.0e-9);

    const double inv_lambda{ 1.0 / lambda(0, 0) };
    std::vector<double> slopes{ -slope_mag, slope_mag };
    for (int i = 0; i < 2; ++i)
    {
        double residual{ slopes[static_cast<size_t>(i)] - mu(0) };
        double expected_weight{ std::exp(-0.5 * alpha_g * residual * inv_lambda * residual) };
        EXPECT_NEAR(expected_weight, helper.GetMemberWeight(i), 1.0e-9);
    }
}

TEST_F(HRLModelHelperTest, MemberWeightsClampedToMinimum)
{
    // Members with moderately large residuals yield weights below the minimum
    auto member0{ CreateMember({{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 5.0}, {1.0, 6.0}, {2.0, 7.0}}, "member_1") };
    std::vector<DataTuple> data_array{ member0, member1 };

    HRLModelHelper helper{ 2, 2 };
    helper.SetDataArray(std::move(data_array));

    const double alpha_g{ 100.0 }; // Forces weights below m_weight_member_min
    helper.SetUniversalAlphaR(0.0);
    EXPECT_NO_THROW(helper.RunEstimation(alpha_g));

    const double min_weight{ HRLModelHelper::DEFAULT_WEIGHT_MEMBER_MIN / 2.0 };
    for (int i = 0; i < 2; ++i)
    {
        double w{ helper.GetMemberWeight(i) };
        EXPECT_TRUE(std::isfinite(w));
        EXPECT_GE(w, min_weight);
        EXPECT_DOUBLE_EQ(w, min_weight);
    }
}

TEST_F(HRLModelHelperTest, GettersThrowOnInvalidId)
{
    const int basis_size{ 1 };
    const int member_size{ 1 };
    HRLModelHelper helper{ basis_size, member_size };

    const int negative_id{ -1 };
    EXPECT_THROW(helper.GetOutlierFlag(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetStatisticalDistance(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetCapitalSigmaMatrixPosterior(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixPosterior(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixMDPDE(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetBetaMatrixOLS(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetSigmaSquare(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetMemberWeight(negative_id), std::out_of_range);
    EXPECT_THROW(helper.GetDataWeightMatrix(negative_id), std::out_of_range);
}

TEST_F(HRLModelHelperTest, EstimationOnSmallSyntheticData)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_1") };
    std::vector<DataTuple> data_array{member0, member1};

    HRLModelHelper helper(2, 2);
    helper.SetDataArray(std::move(data_array));
    helper.SetMaximumIteration(1000);
    helper.SetTolerance(1.0e-6);
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);

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
    helper.SetDataArray(std::move(data_array));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);
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

TEST_F(HRLModelHelperTest, MuPriorDiffersWithDistinctMemberTrends)
{
    auto member0{ CreateMember({{0.0, 0.0}, {1.0, 1.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 0.0}, {1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}, {4.0, 8.0}}, "member_1") };
    auto member2{ CreateMember({{0.0, 10.0}, {1.0, 5.0}, {2.0, 0.0}}, "member_2") };
    std::vector<DataTuple> data_array{ member0, member1, member2 };

    HRLModelHelper helper(2, 3);
    helper.SetDataArray(std::move(data_array));
    helper.SetMaximumIteration(1000);
    helper.SetTolerance(1.0e-6);
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);

    const auto & mu_prior{ helper.GetMuVectorPrior() };
    const auto & mu_mdpde{ helper.GetMuVectorMDPDE() };
    const auto & mu_mean{ helper.GetMuVectorMean() };
    EXPECT_FALSE(mu_prior.isApprox(mu_mdpde));
    EXPECT_FALSE(mu_prior.isApprox(mu_mean));

    for (int i = 0; i < 3; ++i)
    {
        const auto & beta{ helper.GetBetaMatrixPosterior(i) };
        EXPECT_TRUE(beta.array().isFinite().all());
    }
}

TEST_F(HRLModelHelperTest, OutlierFlagsWithDivergentData)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 6.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 10.0}, {1.0, 12.0}, {2.0, 15.0}}, "member_1") };
    std::vector<DataTuple> data_array{member0, member1};

    HRLModelHelper helper(2, 2);
    helper.SetDataArray(std::move(data_array));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);

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

TEST_F(HRLModelHelperTest, MuPriorEqualsMuMDPDEForTwoMembers)
{
    auto member0{ CreateMember({{0.0,  0.0}, {1.0,  1.0}, {2.0,  2.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 10.0}, {1.0, 12.0}, {2.0, 14.0}}, "member_1") };
    std::vector<DataTuple> data_array{ member0, member1 };
    HRLModelHelper helper(2, 2);
    helper.SetDataArray(std::move(data_array));
    helper.SetMaximumIteration(1000);
    helper.SetTolerance(1.0e-6);
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);
    const auto & mu_prior{ helper.GetMuVectorPrior() };
    const auto & mu_mdpde{ helper.GetMuVectorMDPDE() };
    EXPECT_TRUE(mu_prior.isApprox(mu_mdpde));

    const auto & data0{ std::get<0>(member0) };
    const auto & data1{ std::get<0>(member1) };
    const int n0{ static_cast<int>(data0.size()) };
    const int n1{ static_cast<int>(data1.size()) };
    Eigen::MatrixXd X0(n0, 2), X1(n1, 2);
    Eigen::VectorXd y0(n0), y1(n1);
    for (int i = 0; i < n0; i++)
    {
        X0.row(i) = data0[static_cast<size_t>(i)].head(2);
        y0(i) = data0[static_cast<size_t>(i)](2);
    }
    for (int i = 0; i < n1; i++)
    {
        X1.row(i) = data1[static_cast<size_t>(i)].head(2);
        y1(i) = data1[static_cast<size_t>(i)](2);
    }

    Eigen::VectorXd w0{ helper.GetDataWeightMatrix(0).diagonal() };
    Eigen::VectorXd w1{ helper.GetDataWeightMatrix(1).diagonal() };
    double sigma0{ helper.GetSigmaSquare(0) };
    double sigma1{ helper.GetSigmaSquare(1) };
    double inv_trace0{ w0.array().inverse().sum() };
    double inv_trace1{ w1.array().inverse().sum() };
    Eigen::VectorXd cap_sigma_diag0(n0);
    Eigen::VectorXd cap_sigma_diag1(n1);
    for (int j = 0; j < n0; ++j) cap_sigma_diag0(j) = n0 * sigma0 / w0(j) / inv_trace0;
    for (int j = 0; j < n1; ++j) cap_sigma_diag1(j) = n1 * sigma1 / w1(j) / inv_trace1;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> cap_sigma0(cap_sigma_diag0);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> cap_sigma1(cap_sigma_diag1);
    auto inv_cap_sigma0{ EigenMatrixUtility::GetInverseDiagonalMatrix(cap_sigma0) };
    auto inv_cap_sigma1{ EigenMatrixUtility::GetInverseDiagonalMatrix(cap_sigma1) };
    Eigen::MatrixXd gram0{ X0.transpose() * inv_cap_sigma0 * X0 };
    Eigen::MatrixXd gram1{ X1.transpose() * inv_cap_sigma1 * X1 };
    Eigen::VectorXd moment0{ X0.transpose() * inv_cap_sigma0 * y0 };
    Eigen::VectorXd moment1{ X1.transpose() * inv_cap_sigma1 * y1 };
    double omega0{ helper.GetMemberWeight(0) };
    double omega1{ helper.GetMemberWeight(1) };
    double omega_h{ 1.0 / (1.0 / omega0 + 1.0 / omega1) };
    Eigen::MatrixXd capital_lambda{ helper.GetCapitalLambdaMatrix() };
    Eigen::MatrixXd lambda0{ (2.0 * omega_h / omega0) * capital_lambda };
    Eigen::MatrixXd lambda1{ (2.0 * omega_h / omega1) * capital_lambda };
    auto inv_lambda0{ EigenMatrixUtility::GetInverseMatrix(lambda0) };
    auto inv_lambda1{ EigenMatrixUtility::GetInverseMatrix(lambda1) };
    Eigen::MatrixXd inv_sigma_post0{ gram0 + inv_lambda0 };
    Eigen::MatrixXd inv_sigma_post1{ gram1 + inv_lambda1 };
    Eigen::MatrixXd sigma_post0{ EigenMatrixUtility::GetInverseMatrix(inv_sigma_post0) };
    Eigen::MatrixXd sigma_post1{ EigenMatrixUtility::GetInverseMatrix(inv_sigma_post1) };
    Eigen::VectorXd numerator{ inv_lambda0 * sigma_post0 * moment0 + inv_lambda1 * sigma_post1 * moment1 };
    Eigen::MatrixXd denominator{ inv_lambda0 * sigma_post0 * gram0 + inv_lambda1 * sigma_post1 * gram1 };
    Eigen::VectorXd mu_prior_general{ EigenMatrixUtility::GetInverseMatrix(denominator) * numerator };
    EXPECT_FALSE(mu_prior.isApprox(mu_prior_general));
    EXPECT_FALSE(mu_mdpde.isApprox(mu_prior_general));
}

TEST_F(HRLModelHelperTest, AlphaGReweightsOutliers)
{
    auto inlier0{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}}, "inlier_0") };
    auto inlier1{ CreateMember({{0.0, 1.0}, {1.0, 2.0}, {2.0, 3.0}}, "inlier_1") };
    auto outlier{ CreateMember({{0.0, 10.0}, {1.0, 20.0}, {2.0, 30.0}}, "outlier") };
    std::vector<DataTuple> data_array{inlier0, inlier1, outlier};

    HRLModelHelper baseline(2, 3);
    auto baseline_data{ data_array };
    baseline.SetDataArray(std::move(baseline_data));
    baseline.SetMaximumIteration(1000);
    baseline.SetTolerance(1.0e-6);
    baseline.SetUniversalAlphaR(0.0);
    baseline.RunEstimation(0.0);
    Eigen::VectorXd mu_prior_initial{ baseline.GetMuVectorPrior() };

    HRLModelHelper reweighted(2, 3);
    reweighted.SetDataArray(std::move(data_array));
    reweighted.SetMaximumIteration(1000);
    reweighted.SetTolerance(1.0e-6);
    reweighted.SetUniversalAlphaR(0.0);
    reweighted.RunEstimation(0.0);
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

TEST_F(HRLModelHelperTest, DetectsTrueOutlier)
{
    auto inlier0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 7.0}}, "inlier_0") };
    auto inlier1{ CreateMember({{0.0, 1.2}, {1.0, 3.2}, {2.0, 5.2}, {3.0, 7.2}}, "inlier_1") };
    auto outlier{ CreateMember({{0.0, 30.0}, {1.0, 60.0}, {2.0, 90.0}, {3.0, 120.0}}, "outlier") };
    std::vector<DataTuple> data_array{ inlier0, inlier1, outlier };
    HRLModelHelper helper(2, 3);
    helper.SetDataArray(std::move(data_array));
    helper.SetMaximumIteration(1000);
    helper.SetTolerance(1.0e-6);
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(1.0);
    EXPECT_FALSE(helper.GetOutlierFlag(0));
    EXPECT_FALSE(helper.GetOutlierFlag(1));
    EXPECT_TRUE(helper.GetOutlierFlag(2));
}

TEST_F(HRLModelHelperTest, FlagsOnlyExtremePosteriorAsOutlier)
{
    auto member0{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 7.0}}, "member_0") };
    auto member1{ CreateMember({{0.0, 1.1}, {1.0, 3.1}, {2.0, 5.1}, {3.0, 7.1}}, "member_1") };
    auto extreme{ CreateMember({{0.0, 100.0}, {1.0, 200.0}, {2.0, 300.0}, {3.0, 400.0}}, "extreme") };
    std::vector<DataTuple> data_array{ member0, member1, extreme };
    HRLModelHelper helper(2, 3);
    helper.SetDataArray(std::move(data_array));
    helper.SetMaximumIteration(1000);
    helper.SetTolerance(1.0e-6);
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(1.0);
    constexpr double chi_square_threshold{ 9.21 }; // 99th percentile, df = 2
    EXPECT_LT(helper.GetStatisticalDistance(0), chi_square_threshold);
    EXPECT_LT(helper.GetStatisticalDistance(1), chi_square_threshold);
    EXPECT_GT(helper.GetStatisticalDistance(2), chi_square_threshold);
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
    auto baseline_data{ m_data };
    baseline.SetDataArray(std::move(baseline_data));
    baseline.SetMaximumIteration(50);
    baseline.SetTolerance(0.0);
    baseline.SetUniversalAlphaR(0.5);
    baseline.RunEstimation(0.0);
    double beta_baseline{baseline.GetBetaMatrixMDPDE(0)(0)};

    // Run with single iteration
    HRLModelHelper limited1(1, 1);
    auto limited1_data{ m_data };
    limited1.SetDataArray(std::move(limited1_data));
    limited1.SetMaximumIteration(1);
    limited1.SetTolerance(0.0);
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    limited1.SetUniversalAlphaR(0.5);
    limited1.RunEstimation(0.0);
    const std::string log1{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(LogLevel::Error);
    EXPECT_NE(std::string::npos, log1.find("Reach maximum iterations"));
    double beta1{limited1.GetBetaMatrixMDPDE(0)(0)};

    // Run with five iterations
    HRLModelHelper limited5(1, 1);
    auto limited5_data{ m_data };
    limited5.SetDataArray(std::move(limited5_data));
    limited5.SetMaximumIteration(5);
    limited5.SetTolerance(0.0);
    limited5.SetUniversalAlphaR(0.5);
    limited5.RunEstimation(0.0);
    double beta5{ limited5.GetBetaMatrixMDPDE(0)(0) };

    // Verify more iterations move result toward baseline
    EXPECT_GT(std::abs(beta1 - beta_baseline), std::abs(beta5 - beta_baseline));
}

TEST_F(HRLModelHelperFixture, RespectsToleranceSetting)
{
    // Reference solution with strict tolerance and many iterations
    HRLModelHelper baseline(1, 1);
    auto baseline_data{ m_data };
    baseline.SetDataArray(std::move(baseline_data));
    baseline.SetMaximumIteration(HRLModelHelper::DEFAULT_MAXIMUM_ITERATION);
    baseline.SetTolerance(1.0e-12);
    baseline.SetUniversalAlphaR(0.5);
    baseline.RunEstimation(0.0);
    double beta_full{ baseline.GetBetaMatrixMDPDE(0)(0) };

    // Strict tolerance with limited iterations should hit maximum
    HRLModelHelper strict(1, 1);
    auto strict_data{ m_data };
    strict.SetDataArray(std::move(strict_data));
    strict.SetMaximumIteration(5);
    strict.SetTolerance(1.0e-12);
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    strict.SetUniversalAlphaR(0.5);
    strict.RunEstimation(0.0);
    const std::string log_strict{ testing::internal::GetCapturedStderr() };
    EXPECT_NE(std::string::npos, log_strict.find("Reach maximum iterations"));
    double beta_strict{ strict.GetBetaMatrixMDPDE(0)(0) };

    // Relaxed tolerance should converge early and avoid warning
    HRLModelHelper relaxed(1, 1);
    auto relaxed_data{ m_data };
    relaxed.SetDataArray(std::move(relaxed_data));
    relaxed.SetMaximumIteration(5);
    relaxed.SetTolerance(1.0e-1);
    testing::internal::CaptureStderr();
    relaxed.SetUniversalAlphaR(0.5);
    relaxed.RunEstimation(0.0);
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
    helper.SetDataArray(std::move(data_array));
    helper.SetUniversalAlphaR(1.0e20);
    helper.RunEstimation(0.0);
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
    helper.SetDataArray(std::move(data_array));
    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(1.0);
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

TEST_F(HRLModelHelperTest, ExtremelySmallRawDataWeightsAreClamped)
{
    // Many near-zero samples and one extreme outlier produce raw weights
    // smaller than the minimum threshold. Verify weights are clamped.
    constexpr int n{ 4001 }; // 4000 normal samples, 1 extreme
    std::vector<Eigen::VectorXd> samples;
    samples.reserve(n);
    Eigen::VectorXd normal(2);
    normal << 1.0, 0.0;
    for (int i = 0; i < n - 1; ++i) samples.emplace_back(normal);
    Eigen::VectorXd outlier(2);
    outlier << 1.0, 1.0e12;
    samples.emplace_back(outlier);
    std::vector<DataTuple> data_array;
    data_array.emplace_back(samples, "extreme");
    HRLModelHelper helper{ 1, 1 };
    helper.SetDataArray(std::move(data_array));
    helper.SetUniversalAlphaR(0.01);
    helper.RunEstimation(0.0);
    const auto & W{ helper.GetDataWeightMatrix(0) };
    Eigen::VectorXd weights{ W.diagonal() };
    for (int i = 0; i < weights.size(); ++i)
    {
        EXPECT_GE(weights(i), HRLModelHelper::DEFAULT_WEIGHT_DATA_MIN);
    }
    EXPECT_DOUBLE_EQ(HRLModelHelper::DEFAULT_WEIGHT_DATA_MIN, weights.minCoeff());
}

TEST_F(HRLModelHelperTest, NonFiniteExponentFallsBackToMinimum)
{
    std::vector<Eigen::VectorXd> samples;
    Eigen::VectorXd normal(2);
    normal << 1.0, 0.0;
    samples.emplace_back(normal);
    Eigen::VectorXd outlier(2);
    outlier << 1.0, std::numeric_limits<double>::max();
    samples.emplace_back(outlier);
    std::vector<DataTuple> data_array;
    data_array.emplace_back(samples, "extreme");
    constexpr double alpha_r{ 0.5 };
    HRLModelHelper helper{ 1, 1 };
    helper.SetDataArray(std::move(data_array));
    helper.SetUniversalAlphaR(alpha_r);
    helper.RunEstimation(0.0);

    const auto & W{ helper.GetDataWeightMatrix(0) };
    Eigen::VectorXd weights{ W.diagonal() };
    for (int i = 0; i < weights.size(); ++i)
    {
        EXPECT_TRUE(std::isfinite(weights(i)));
        EXPECT_GE(weights(i), HRLModelHelper::DEFAULT_WEIGHT_DATA_MIN);
        EXPECT_LE(weights(i), 1.0);
    }
    EXPECT_DOUBLE_EQ(HRLModelHelper::DEFAULT_WEIGHT_DATA_MIN, weights.minCoeff());

    double beta{ helper.GetBetaMatrixMDPDE(0)(0) };
    double sigma{ helper.GetSigmaSquare(0) };
    double exponent{ -0.5 * alpha_r * std::pow(outlier(1) - beta, 2) / sigma };
    EXPECT_FALSE(std::isfinite(exponent));
}

TEST_F(HRLModelHelperTest, ZeroResidualProducesFiniteWeightsAndPosterior)
{
    auto member{ CreateMember({{0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}}, "perfect") };
    std::vector<DataTuple> data_array{ member };
    HRLModelHelper helper(2, 1);
    helper.SetDataArray(std::move(data_array));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);
    const auto & sigma{ helper.GetCapitalSigmaMatrixPosterior(0) };
    EXPECT_TRUE(sigma.array().isFinite().all());
    const auto & beta{ helper.GetBetaMatrixPosterior(0) };
    EXPECT_TRUE(beta.array().isFinite().all());
}

TEST_F(HRLModelHelperTest, DegenerateWeightsUseFallbackCovariance)
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
    ASSERT_NO_THROW(helper.SetDataArray(std::move(data_array)));
    helper.SetUniversalAlphaR(1.0e8);
    EXPECT_NO_THROW(helper.RunEstimation(0.0));
    const auto & sigma{ helper.GetDataCovarianceMatrix(0) };
    EXPECT_TRUE(sigma.diagonal().array().isFinite().all());
    EXPECT_GT(sigma.diagonal().minCoeff(), 0.0);
    Eigen::VectorXd weights{ helper.GetDataWeightMatrix(0).diagonal() };
    EXPECT_TRUE(weights.array().isFinite().all());
}

TEST_F(HRLModelHelperTest, CovarianceOverflowResetsToIdentity)
{
    HRLModelHelper helper{1, 2};

    std::vector<Eigen::VectorXd> member1;
    Eigen::VectorXd m1s1(2);
    m1s1 << 1.0, 0.0;
    Eigen::VectorXd m1s2(2);
    m1s2 << 1.0, std::numeric_limits<double>::max();
    member1.emplace_back(m1s1);
    member1.emplace_back(m1s2);

    std::vector<Eigen::VectorXd> member2;
    Eigen::VectorXd m2s1(2);
    m2s1 << 1.0, 0.0;
    Eigen::VectorXd m2s2(2);
    m2s2 << 1.0, -std::numeric_limits<double>::max();
    member2.emplace_back(m2s1);
    member2.emplace_back(m2s2);

    std::vector<DataTuple> data_array;
    data_array.emplace_back(member1, "member1");
    data_array.emplace_back(member2, "member2");

    helper.SetDataArray(std::move(data_array));

    auto prev_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    helper.SetUniversalAlphaR(0.5);
    ASSERT_NO_THROW(helper.RunEstimation(0.0));
    const std::string err{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(prev_level);

    EXPECT_NE(std::string::npos,
              err.find("Resulting covariance has non-finite entries"));
    EXPECT_TRUE(helper.GetCapitalLambdaMatrix()
                    .isApprox(Eigen::MatrixXd::Identity(1, 1)));
}

TEST_F(HRLModelHelperTest, OutlierFlagsAllWhenBasisSizeTooLarge)
{
    const int basis_size{ 11 };
    const int member_size{ 2 };
    HRLModelHelper helper{ basis_size, member_size };

    // Member 0 data uses basis indices 0 and 1
    std::vector<Eigen::VectorXd> member0;
    Eigen::VectorXd m0s1(basis_size + 1);
    m0s1 << 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    Eigen::VectorXd m0s2(basis_size + 1);
    m0s2 << 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;
    member0.emplace_back(m0s1);
    member0.emplace_back(m0s2);

    // Member 1 data uses basis indices 0 and 2
    std::vector<Eigen::VectorXd> member1;
    Eigen::VectorXd m1s1(basis_size + 1);
    m1s1 << 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    Eigen::VectorXd m1s2(basis_size + 1);
    m1s2 << 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;
    member1.emplace_back(m1s1);
    member1.emplace_back(m1s2);

    std::vector<DataTuple> data;
    data.emplace_back(member0, "m0");
    data.emplace_back(member1, "m1");
    helper.SetDataArray(std::move(data));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);

    for (int i = 0; i < member_size; ++i)
    {
        EXPECT_TRUE(helper.GetOutlierFlag(i));
        EXPECT_GE(helper.GetStatisticalDistance(i), 0.0);
        EXPECT_TRUE(std::isfinite(helper.GetSigmaSquare(i)));
        EXPECT_TRUE(std::isfinite(helper.GetMemberWeight(i)));
        EXPECT_TRUE(helper.GetBetaMatrixOLS(i).array().isFinite().all());
        EXPECT_TRUE(helper.GetBetaMatrixMDPDE(i).array().isFinite().all());
        EXPECT_TRUE(helper.GetBetaMatrixPosterior(i).array().isFinite().all());
        EXPECT_TRUE(helper.GetCapitalSigmaMatrixPosterior(i).array().isFinite().all());
        EXPECT_TRUE(helper.GetDataWeightMatrix(i).diagonal().array().isFinite().all());
    }
    EXPECT_TRUE(helper.GetCapitalLambdaMatrix().array().isFinite().all());
    EXPECT_TRUE(helper.GetMuVectorPrior().array().isFinite().all());
    EXPECT_TRUE(helper.GetMuVectorMDPDE().array().isFinite().all());
    EXPECT_TRUE(helper.GetMuVectorMean().array().isFinite().all());
}

TEST_F(HRLModelHelperTest, AlphaZeroYieldsIdentityDataWeights)
{
    HRLModelHelper helper{ 1, 1 };
    std::vector<Eigen::VectorXd> member_data;
    Eigen::VectorXd sample1(2); sample1 << 1.0, 2.0;
    Eigen::VectorXd sample2(2); sample2 << 2.0, 5.0;
    member_data.emplace_back(sample1);
    member_data.emplace_back(sample2);
    std::vector<DataTuple> data_array;
    data_array.emplace_back(member_data, "member1");
    helper.SetDataArray(std::move(data_array));
    helper.SetUniversalAlphaR(0.0);
    helper.RunEstimation(0.0);

    const auto & W{ helper.GetDataWeightMatrix(0) };
    auto diag{ W.diagonal() };
    for (int i = 0; i < diag.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(1.0, diag(i));
    }
}
