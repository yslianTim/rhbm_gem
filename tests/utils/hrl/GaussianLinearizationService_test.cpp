#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>

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

rhbm_gem::GaussianLinearizationService MakeDatasetService()
{
    return rhbm_gem::GaussianLinearizationService{
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset()
    };
}

} // namespace

TEST(GaussianLinearizationServiceTest, BuildLinearModelSeriesTransformsPositiveWeightedSamplesOnly)
{
    const LocalPotentialSampleList sampling_entries{
        {0.1f, 4.0f, 0.5f},
        {0.2f, -2.0f, 7.0f},
        {0.3f, 8.0f, 2.5f},
        {0.4f, 9.0f, 0.0f},
    };

    const auto series{ MakeDatasetService().BuildLinearModelSeries(sampling_entries) };

    ASSERT_EQ(series.size(), 2U);
    EXPECT_NEAR(-0.5 * 0.1 * 0.1, series.at(0).GetBasisValue(0), 1.0e-7);
    EXPECT_NEAR(std::log(4.0), series.at(0).response, 1.0e-7);
    EXPECT_NEAR(-0.5 * 0.3 * 0.3, series.at(1).GetBasisValue(0), 1.0e-7);
    EXPECT_NEAR(std::log(8.0), series.at(1).response, 1.0e-7);
}

TEST(GaussianLinearizationServiceTest, EncodeGaussianToBetaMatchesClosedForm)
{
    const auto beta{
        MakeDatasetService().EncodeGaussianToBeta(MakeVector({ 2.0, 0.5, 0.0 }))
    };
    const auto expected_beta0{
        std::log(2.0) - 1.5 * std::log(Constants::two_pi * 0.5 * 0.5)
    };
    const auto expected_beta1{ 1.0 / (0.5 * 0.5) };

    ASSERT_EQ(beta.size(), 2);
    EXPECT_NEAR(expected_beta0, beta(0), 1.0e-12);
    EXPECT_NEAR(expected_beta1, beta(1), 1.0e-12);
}

TEST(GaussianLinearizationServiceTest, DecodeGroupEstimateMatchesClosedForm)
{
    const auto estimate{
        MakeDatasetService().DecodeGroupEstimate(MakeVector({ std::log(2.0), 4.0 }))
    };
    const auto expected_amplitude{ 2.0 * std::pow(Constants::two_pi / 4.0, 1.5) };
    const auto expected_width{ 1.0 / std::sqrt(4.0) };

    EXPECT_NEAR(expected_amplitude, estimate.amplitude, 1.0e-12);
    EXPECT_NEAR(expected_width, estimate.width, 1.0e-12);
}

TEST(GaussianLinearizationServiceTest, DecodePosteriorEstimateMatchesCurrentVarianceFormula)
{
    const auto linear_model{ MakeVector({ 0.5, 2.0 }) };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.05,
                         0.05, 0.2;

    const auto posterior{
        MakeDatasetService().DecodePosteriorEstimate(linear_model, covariance_matrix)
    };
    const auto expected_amplitude{
        std::exp(0.5) * std::pow(Constants::two_pi / 2.0, 1.5)
    };
    const auto expected_width{ 1.0 / std::sqrt(2.0) };
    const auto expected_var_amplitude{
        std::pow(Constants::two_pi, 3) * std::exp(1.0) * std::pow(2.0, -5) *
        (std::pow(2.0, 2) * 0.1 + 2.25 * 0.2 - 3.0 * 2.0 * 0.05)
    };
    const auto expected_var_width{ 0.25 * std::pow(2.0, -3) * 0.2 };

    EXPECT_NEAR(expected_amplitude, posterior.estimate.amplitude, 1.0e-12);
    EXPECT_NEAR(expected_width, posterior.estimate.width, 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_amplitude), posterior.variance.amplitude, 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_width), posterior.variance.width, 1.0e-12);
}
