#include <gtest/gtest.h>

#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/HRLGroupEstimator.hpp>

namespace
{
Eigen::VectorXd MakeVector(std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index{ 0 };
    for (double value : values)
    {
        result(index++) = value;
    }
    return result;
}

HRLDiagonalMatrix MakeDiagonal(std::initializer_list<double> values)
{
    HRLDiagonalMatrix result(static_cast<Eigen::Index>(values.size()));
    result.diagonal() = MakeVector(values);
    return result;
}

SeriesPoint MakeSeriesPoint(std::initializer_list<double> values)
{
    std::vector<double> row(values);
    const auto response{ row.back() };
    row.pop_back();
    return SeriesPoint{ std::move(row), response };
}

HRLMemberDataset MakeDataset(std::initializer_list<SeriesPoint> rows)
{
    return rhbm_gem::rhbm_helper::BuildMemberDataset(SeriesPointList(rows));
}

HRLBetaEstimateResult MakeFitResult(
    std::initializer_list<double> beta_values,
    double sigma_square,
    std::initializer_list<double> weight_values,
    std::initializer_list<double> covariance_values,
    HRLEstimationStatus status = HRLEstimationStatus::SUCCESS)
{
    HRLBetaEstimateResult fit_result;
    fit_result.status = status;
    fit_result.beta_mdpde = MakeVector(beta_values);
    fit_result.sigma_square = sigma_square;
    fit_result.data_weight = MakeDiagonal(weight_values);
    fit_result.data_covariance = MakeDiagonal(covariance_values);
    return fit_result;
}
} // namespace

TEST(HRLGroupEstimatorTest, EstimateSingleMemberUsesFallbackResult)
{
    const auto input{
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 1.0 }), MakeSeriesPoint({ 1.0, 1.0, 3.0 }) })
            },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 })
            }
        )
    };

    const auto result{ HRLGroupEstimator().Estimate(input, 0.0) };

    EXPECT_EQ(HRLEstimationStatus::SINGLE_MEMBER, result.status);
    EXPECT_TRUE(result.mu_mean.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(result.mu_prior.isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    EXPECT_TRUE(result.beta_posterior_matrix.col(0).isApprox(MakeVector({ 1.0, 2.0 }), 1e-12));
    ASSERT_EQ(1u, result.capital_sigma_posterior_list.size());
    EXPECT_TRUE(result.capital_sigma_posterior_list[0].isApprox(Eigen::MatrixXd::Zero(2, 2), 1e-12));
}

TEST(HRLGroupEstimatorTest, EstimateRejectsMissingMemberFitResults)
{
    HRLGroupEstimationInput input;
    input.basis_size = 2;
    input.member_datasets.push_back(
        { Eigen::MatrixXd::Identity(2, 2), MakeVector({ 1.0, 2.0 }), MakeVector({ 1.0, 1.0 }) });

    EXPECT_THROW(HRLGroupEstimator().Estimate(input, 0.0), std::invalid_argument);
}

TEST(HRLGroupEstimatorTest, EstimateRejectsMismatchedWeightSize)
{
    HRLGroupEstimationInput input;
    input.basis_size = 2;
    input.member_datasets.push_back({
        Eigen::MatrixXd::Identity(2, 2),
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 1.0, 1.0 })
    });

    auto fit_result{ MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0 }, { 0.25, 0.25 }) };
    input.member_fit_results.push_back(fit_result);

    EXPECT_THROW(HRLGroupEstimator().Estimate(input, 0.0), std::invalid_argument);
}

TEST(HRLGroupEstimatorTest, EstimateRejectsInconsistentScoreSize)
{
    HRLGroupEstimationInput input;
    input.basis_size = 2;
    input.member_datasets.push_back({
        Eigen::MatrixXd::Identity(2, 2),
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 1.0 })
    });
    input.member_fit_results.push_back(
        MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 }));

    EXPECT_THROW(HRLGroupEstimator().Estimate(input, 0.0), std::invalid_argument);
}

TEST(HRLGroupEstimatorTest, EstimateRejectsBasisSizeMismatchWithMemberDatasets)
{
    HRLGroupEstimationInput input;
    input.basis_size = 3;
    input.member_datasets.push_back({
        Eigen::MatrixXd::Identity(2, 2),
        MakeVector({ 1.0, 2.0 }),
        MakeVector({ 1.0, 1.0 })
    });
    input.member_fit_results.push_back(
        MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 }));

    EXPECT_THROW(HRLGroupEstimator().Estimate(input, 0.0), std::invalid_argument);
}

TEST(HRLGroupEstimatorTest, EstimateTwoMembersPinsMuPriorToMuMDPDE)
{
    const auto input{
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 1.0 }), MakeSeriesPoint({ 1.0, 1.0, 3.0 }) }),
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 2.0 }), MakeSeriesPoint({ 1.0, 1.0, 4.0 }) })
            },
            {
                MakeFitResult({ 1.0, 2.0 }, 0.25, { 1.0, 1.0 }, { 0.25, 0.25 }),
                MakeFitResult({ 2.0, 2.0 }, 0.5, { 1.0, 1.0 }, { 0.5, 0.5 })
            }
        )
    };

    const auto result{ HRLGroupEstimator().Estimate(input, 0.0) };

    EXPECT_EQ(HRLEstimationStatus::SUCCESS, result.status);
    EXPECT_TRUE(result.mu_prior.isApprox(result.mu_mdpde, 1e-12));
    EXPECT_EQ(2, result.beta_posterior_matrix.cols());
    ASSERT_EQ(2u, result.capital_sigma_posterior_list.size());
    EXPECT_EQ(2, result.statistical_distance_array.size());
    EXPECT_EQ(2, result.outlier_flag_array.size());
}

TEST(HRLGroupEstimatorTest, EstimateAllowsConvergenceFallbackStatuses)
{
    const auto input{
        rhbm_gem::rhbm_helper::BuildGroupInput(
            {
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 1.0 }), MakeSeriesPoint({ 1.0, 1.0, 3.0 }) }),
                MakeDataset({ MakeSeriesPoint({ 1.0, 0.0, 2.0 }), MakeSeriesPoint({ 1.0, 1.0, 4.0 }) })
            },
            {
                MakeFitResult(
                    { 1.0, 2.0 },
                    0.25,
                    { 1.0, 1.0 },
                    { 0.25, 0.25 },
                    HRLEstimationStatus::MAX_ITERATIONS_REACHED),
                MakeFitResult(
                    { 2.0, 2.0 },
                    0.5,
                    { 1.0, 1.0 },
                    { 0.5, 0.5 },
                    HRLEstimationStatus::NUMERICAL_FALLBACK)
            }
        )
    };

    const auto result{ HRLGroupEstimator().Estimate(input, 0.0) };

    EXPECT_EQ(HRLEstimationStatus::SUCCESS, result.status);
    EXPECT_EQ(2, result.beta_posterior_matrix.cols());
}
