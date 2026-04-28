#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

namespace
{
namespace ls = rhbm_gem::linearization_service;

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

const ls::LinearizationSpec & AtomLinearizationSpec()
{
    static const auto spec{ ls::LinearizationSpec::AtomDecode() };
    return spec;
}

const ls::LinearizationSpec & BondLinearizationSpec()
{
    static const auto spec{ ls::LinearizationSpec::BondDecode() };
    return spec;
}

} // namespace

TEST(LinearizationServiceTest, BuildDatasetSeriesTransformsEffectiveSamplesWithinRange)
{
    const LocalPotentialSampleList sampling_entries{
        {0.1f, 4.0f, 0.5f},
        {0.2f, -2.0f, 7.0f},
        {0.3f, 8.0f, 2.5f},
        {0.4f, 9.0f, 0.0f},
        {0.8f, 16.0f, 3.0f},
    };

    const auto series{
        ls::BuildDatasetSeries(sampling_entries, ls::LinearizationRange{ 0.0, 0.5 })
    };

    ASSERT_EQ(series.size(), 2U);
    EXPECT_NEAR(1.0, series.at(0).GetBasisValue(0), 1.0e-7);
    EXPECT_NEAR(-0.5 * 0.1 * 0.1, series.at(0).GetBasisValue(1), 1.0e-7);
    EXPECT_NEAR(std::log(4.0), series.at(0).response, 1.0e-7);
    EXPECT_FLOAT_EQ(0.5f, series.at(0).score);
    EXPECT_NEAR(1.0, series.at(1).GetBasisValue(0), 1.0e-7);
    EXPECT_NEAR(-0.5 * 0.3 * 0.3, series.at(1).GetBasisValue(1), 1.0e-7);
    EXPECT_NEAR(std::log(8.0), series.at(1).response, 1.0e-7);
    EXPECT_FLOAT_EQ(2.5f, series.at(1).score);
}

TEST(LinearizationServiceTest, BuildDatasetSeriesRejectsInvalidRange)
{
    const LocalPotentialSampleList sampling_entries{
        {0.1f, 4.0f, 0.5f},
    };

    EXPECT_THROW(
        ls::BuildDatasetSeries(sampling_entries, ls::LinearizationRange{ 1.0, 0.0 }),
        std::invalid_argument);
    EXPECT_THROW(
        ls::BuildDatasetSeries(
            sampling_entries,
            ls::LinearizationRange{ std::numeric_limits<double>::quiet_NaN(), 1.0 }),
        std::invalid_argument);
}

TEST(LinearizationServiceTest, EncodeGaussianToParameterVectorMatchesClosedForm)
{
    const rhbm_gem::GaussianModel3D model{ 2.0, 0.5, 0.25 };
    const auto beta{ ls::EncodeGaussianToParameterVector(AtomLinearizationSpec(), model) };
    const auto expected_beta0{
        std::log(2.0) - 1.5 * std::log(Constants::two_pi * 0.5 * 0.5)
    };
    const auto expected_beta1{ 1.0 / (0.5 * 0.5) };

    ASSERT_EQ(beta.size(), 2);
    EXPECT_NEAR(expected_beta0, beta(0), 1.0e-12);
    EXPECT_NEAR(expected_beta1, beta(1), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorMatchesClosedForm)
{
    const auto estimate{
        ls::DecodeParameterVector(AtomLinearizationSpec(), MakeVector({ std::log(2.0), 4.0 }))
    };
    const auto expected_amplitude{ 2.0 * std::pow(Constants::two_pi / 4.0, 1.5) };
    const auto expected_width{ 1.0 / std::sqrt(4.0) };

    EXPECT_NEAR(expected_amplitude, estimate.GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(expected_width, estimate.GetWidth(), 1.0e-12);
    EXPECT_NEAR(0.0, estimate.GetIntercept(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorReturnsStandardDeviation)
{
    const auto linear_model{ MakeVector({ 0.5, 2.0 }) };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.05,
                         0.05, 0.2;

    const auto gaussian{
        ls::DecodeParameterVector(AtomLinearizationSpec(), linear_model, covariance_matrix)
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

    EXPECT_NEAR(expected_amplitude, gaussian.GetModel().GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(expected_width, gaussian.GetModel().GetWidth(), 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_amplitude), gaussian.GetStandardDeviationModel().GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_width), gaussian.GetStandardDeviationModel().GetWidth(), 1.0e-12);
    EXPECT_NEAR(0.0, gaussian.GetStandardDeviationModel().GetIntercept(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorForBondUsesTwoDimensionalUncertainty)
{
    const auto linear_model{ MakeVector({ 0.5, 2.0 }) };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.05,
                         0.05, 0.2;

    const auto gaussian{
        ls::DecodeParameterVector(BondLinearizationSpec(), linear_model, covariance_matrix)
    };
    const auto expected_amplitude{ std::exp(0.5) * Constants::two_pi / 2.0 };
    const auto expected_width{ 1.0 / std::sqrt(2.0) };
    const auto expected_var_amplitude{
        expected_amplitude * expected_amplitude *
        (0.1 + 0.2 / (2.0 * 2.0) - 2.0 * 0.05 / 2.0)
    };
    const auto expected_var_width{ 0.25 * std::pow(2.0, -3) * 0.2 };

    EXPECT_NEAR(expected_amplitude, gaussian.GetModel().GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(expected_width, gaussian.GetModel().GetWidth(), 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_amplitude), gaussian.GetStandardDeviationModel().GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_width), gaussian.GetStandardDeviationModel().GetWidth(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorRejectsMismatchedBasisSize)
{
    EXPECT_THROW(
        ls::DecodeParameterVector(AtomLinearizationSpec(), MakeVector({ 0.5, 2.0, 0.25 })),
        std::invalid_argument);

    Eigen::MatrixXd covariance_matrix{ Eigen::MatrixXd::Identity(2, 2) };
    EXPECT_THROW(
        ls::DecodeParameterVector(
            AtomLinearizationSpec(),
            MakeVector({ 0.5, 2.0, 0.25 }),
            covariance_matrix),
        std::invalid_argument);
}
