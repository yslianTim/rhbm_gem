#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
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

const ls::LinearizationSpec & DatasetLinearizationSpec()
{
    static const auto spec{ ls::LinearizationSpec::DefaultDataset() };
    return spec;
}

ls::LinearizationSpec DirectGaussianModelSpec()
{
    auto spec{ ls::LinearizationSpec::DefaultDataset() };
    spec.basis_size = rhbm_gem::GaussianModel3D::ParameterSize();
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
        ls::BuildDatasetSeries(DatasetLinearizationSpec(), sampling_entries, 0.0, 0.5)
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

TEST(LinearizationServiceTest, BuildDatasetSeriesRejectsNonLogQuadraticBasisSize)
{
    auto spec{ DatasetLinearizationSpec() };
    spec.basis_size = rhbm_gem::GaussianModel3D::ParameterSize();
    const LocalPotentialSampleList sampling_entries{
        {0.1f, 4.0f, 0.5f},
    };

    EXPECT_THROW(
        ls::BuildDatasetSeries(spec, sampling_entries, 0.0, 0.5),
        std::invalid_argument);
}

TEST(LinearizationServiceTest, EncodeGaussianToParameterVectorMatchesClosedForm)
{
    const rhbm_gem::GaussianModel3D model{ 2.0, 0.5 };
    const auto beta{ ls::EncodeGaussianToParameterVector(DatasetLinearizationSpec(), model) };
    const auto expected_beta0{
        std::log(2.0) - 1.5 * std::log(Constants::two_pi * 0.5 * 0.5)
    };
    const auto expected_beta1{ 1.0 / (0.5 * 0.5) };

    ASSERT_EQ(beta.size(), 2);
    EXPECT_NEAR(expected_beta0, beta(0), 1.0e-12);
    EXPECT_NEAR(expected_beta1, beta(1), 1.0e-12);
}

TEST(LinearizationServiceTest, EncodeGaussianEstimateToParameterVectorMatchesClosedForm)
{
    const rhbm_gem::GaussianModel3D estimate{ 2.0, 0.5 };
    const auto beta{ ls::EncodeGaussianToParameterVector(DatasetLinearizationSpec(), estimate) };
    const auto expected_beta0{
        std::log(estimate.GetAmplitude())
            - 1.5 * std::log(Constants::two_pi * estimate.GetWidth() * estimate.GetWidth())
    };
    const auto expected_beta1{ 1.0 / (estimate.GetWidth() * estimate.GetWidth()) };

    ASSERT_EQ(beta.size(), 2);
    EXPECT_NEAR(expected_beta0, beta(0), 1.0e-12);
    EXPECT_NEAR(expected_beta1, beta(1), 1.0e-12);
}

TEST(LinearizationServiceTest, EncodeGaussianModelToParameterVectorMatchesRoundTripModel)
{
    const rhbm_gem::GaussianModel3D model{ 2.0, 0.5, 0.25 };
    const auto spec{ DirectGaussianModelSpec() };

    const auto typed_beta{ ls::EncodeGaussianToParameterVector(spec, model) };
    const auto round_trip_beta{
        ls::EncodeGaussianToParameterVector(spec, rhbm_gem::GaussianModel3D::FromVector(model.ToVector()))
    };

    ASSERT_EQ(typed_beta.size(), rhbm_gem::GaussianModel3D::ParameterSize());
    EXPECT_TRUE(typed_beta.isApprox(round_trip_beta, 1.0e-12));
    EXPECT_NEAR(0.25, typed_beta(rhbm_gem::GaussianModel3D::InterceptIndex()), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorMatchesClosedForm)
{
    const auto estimate{
        ls::DecodeParameterVector(DatasetLinearizationSpec(), MakeVector({ std::log(2.0), 4.0 }))
    };
    const auto expected_amplitude{ 2.0 * std::pow(Constants::two_pi / 4.0, 1.5) };
    const auto expected_width{ 1.0 / std::sqrt(4.0) };

    EXPECT_NEAR(expected_amplitude, estimate.GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(expected_width, estimate.GetWidth(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorPreservesDirect3DModelIntercept)
{
    const auto model{
        ls::DecodeParameterVector(DirectGaussianModelSpec(), MakeVector({ 2.0, 0.5, 0.25 }))
    };

    EXPECT_NEAR(2.0, model.GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(0.5, model.GetWidth(), 1.0e-12);
    EXPECT_NEAR(0.25, model.GetIntercept(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorReturnsStandardDeviation)
{
    const auto linear_model{ MakeVector({ 0.5, 2.0 }) };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.05,
                         0.05, 0.2;

    const auto gaussian{
        ls::DecodeParameterVector(DatasetLinearizationSpec(), linear_model, covariance_matrix)
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
}

TEST(LinearizationServiceTest, DecodeParameterVectorFor2BetaModelAddsZeroIntercept)
{
    const auto linear_model{ MakeVector({ 0.5, 2.0 }) };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.05,
                         0.05, 0.2;

    const auto gaussian{
        ls::DecodeParameterVector(
            DatasetLinearizationSpec(),
            linear_model,
            covariance_matrix)
    };

    EXPECT_NEAR(0.0, gaussian.GetModel().GetIntercept(), 1.0e-12);
    EXPECT_NEAR(0.0, gaussian.GetStandardDeviationModel().GetIntercept(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeParameterVectorForDirect3DModelUsesDiagonalStandardDeviation)
{
    const auto linear_model{ MakeVector({ 2.0, 0.5, 0.25 }) };
    Eigen::MatrixXd covariance_matrix{ Eigen::MatrixXd::Zero(3, 3) };
    covariance_matrix(0, 0) = 0.04;
    covariance_matrix(1, 1) = 0.09;
    covariance_matrix(2, 2) = 0.16;

    const auto gaussian{
        ls::DecodeParameterVector(
            DirectGaussianModelSpec(),
            linear_model,
            covariance_matrix)
    };

    EXPECT_NEAR(2.0, gaussian.GetModel().GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(0.5, gaussian.GetModel().GetWidth(), 1.0e-12);
    EXPECT_NEAR(0.25, gaussian.GetModel().GetIntercept(), 1.0e-12);
    EXPECT_NEAR(0.2, gaussian.GetStandardDeviationModel().GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(0.3, gaussian.GetStandardDeviationModel().GetWidth(), 1.0e-12);
    EXPECT_NEAR(0.4, gaussian.GetStandardDeviationModel().GetIntercept(), 1.0e-12);
}
