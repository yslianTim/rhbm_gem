#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>

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

TEST(LinearizationServiceTest, BuildLinearModelSeriesTransformsPositiveWeightedSamplesOnly)
{
    const LocalPotentialSampleList sampling_entries{
        {0.1f, 4.0f, 0.5f},
        {0.2f, -2.0f, 7.0f},
        {0.3f, 8.0f, 2.5f},
        {0.4f, 9.0f, 0.0f},
    };

    const auto series{ ls::BuildLinearModelSeries(DatasetLinearizationSpec(), sampling_entries) };

    ASSERT_EQ(series.size(), 2U);
    EXPECT_NEAR(-0.5 * 0.1 * 0.1, series.at(0).GetBasisValue(0), 1.0e-7);
    EXPECT_NEAR(std::log(4.0), series.at(0).response, 1.0e-7);
    EXPECT_NEAR(-0.5 * 0.3 * 0.3, series.at(1).GetBasisValue(0), 1.0e-7);
    EXPECT_NEAR(std::log(8.0), series.at(1).response, 1.0e-7);
}

TEST(LinearizationServiceTest, EncodeGaussianToBetaMatchesClosedForm)
{
    const auto beta{
        ls::EncodeGaussianToBeta(DatasetLinearizationSpec(), MakeVector({ 2.0, 0.5, 0.0 }))
    };
    const auto expected_beta0{
        std::log(2.0) - 1.5 * std::log(Constants::two_pi * 0.5 * 0.5)
    };
    const auto expected_beta1{ 1.0 / (0.5 * 0.5) };

    ASSERT_EQ(beta.size(), 2);
    EXPECT_NEAR(expected_beta0, beta(0), 1.0e-12);
    EXPECT_NEAR(expected_beta1, beta(1), 1.0e-12);
}

TEST(LinearizationServiceTest, EncodeGaussianEstimateToBetaMatchesClosedForm)
{
    const rhbm_gem::GaussianEstimate estimate{ 2.0, 0.5 };
    const auto beta{ ls::EncodeGaussianToBeta(DatasetLinearizationSpec(), estimate) };
    const auto expected_beta0{
        std::log(estimate.amplitude)
            - 1.5 * std::log(Constants::two_pi * estimate.width * estimate.width)
    };
    const auto expected_beta1{ 1.0 / (estimate.width * estimate.width) };

    ASSERT_EQ(beta.size(), 2);
    EXPECT_NEAR(expected_beta0, beta(0), 1.0e-12);
    EXPECT_NEAR(expected_beta1, beta(1), 1.0e-12);
}

TEST(LinearizationServiceTest, EncodeGaussianModelToBetaMatchesVectorOverload)
{
    const rhbm_gem::GaussianModel3D model{ 2.0, 0.5, 0.25 };
    const auto spec{ DirectGaussianModelSpec() };

    const auto typed_beta{ ls::EncodeGaussianToBeta(spec, model) };
    const auto vector_beta{ ls::EncodeGaussianToBeta(spec, model.ToVector()) };

    ASSERT_EQ(typed_beta.size(), rhbm_gem::GaussianModel3D::ParameterSize());
    EXPECT_TRUE(typed_beta.isApprox(vector_beta, 1.0e-12));
    EXPECT_NEAR(0.25, typed_beta(rhbm_gem::GaussianModel3D::InterceptIndex()), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeGroupEstimateMatchesClosedForm)
{
    const auto estimate{
        ls::DecodeGroupEstimate(DatasetLinearizationSpec(), MakeVector({ std::log(2.0), 4.0 }))
    };
    const auto expected_amplitude{ 2.0 * std::pow(Constants::two_pi / 4.0, 1.5) };
    const auto expected_width{ 1.0 / std::sqrt(4.0) };

    EXPECT_NEAR(expected_amplitude, estimate.amplitude, 1.0e-12);
    EXPECT_NEAR(expected_width, estimate.width, 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeGroupBetaPreservesDirect3DModelIntercept)
{
    const auto gaussian{
        ls::DecodeGroupBeta(DirectGaussianModelSpec(), MakeVector({ 2.0, 0.5, 0.25 }))
    };

    ASSERT_EQ(gaussian.size(), rhbm_gem::GaussianModel3D::ParameterSize());
    EXPECT_NEAR(2.0, gaussian(rhbm_gem::GaussianModel3D::AmplitudeIndex()), 1.0e-12);
    EXPECT_NEAR(0.5, gaussian(rhbm_gem::GaussianModel3D::WidthIndex()), 1.0e-12);
    EXPECT_NEAR(0.25, gaussian(rhbm_gem::GaussianModel3D::InterceptIndex()), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeGroupModel3DReturnsTypedModel)
{
    const auto model{
        ls::DecodeGroupModel3D(DirectGaussianModelSpec(), MakeVector({ 2.0, 0.5, 0.25 }))
    };

    EXPECT_NEAR(2.0, model.GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(0.5, model.GetWidth(), 1.0e-12);
    EXPECT_NEAR(0.25, model.GetIntercept(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeLocalModel3DUsesTypedContext)
{
    auto spec{ DirectGaussianModelSpec() };
    spec.requires_local_context = true;
    const rhbm_gem::GaussianModel3D local_model{ 2.0, 0.5, 0.25 };
    const auto context{ ls::LinearizationContext::FromModel(local_model) };

    const auto decoded{
        ls::DecodeLocalModel3D(spec, MakeVector({ std::log(3.0), 0.25, 0.5 }), context)
    };

    EXPECT_NEAR(6.0, decoded.GetAmplitude(), 1.0e-12);
    EXPECT_NEAR(std::exp(0.25) + 0.5, decoded.GetWidth(), 1.0e-12);
    EXPECT_NEAR(0.75, decoded.GetIntercept(), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodeGaussianEstimateWithUncertaintyReturnsStandardDeviation)
{
    const auto linear_model{ MakeVector({ 0.5, 2.0 }) };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.05,
                         0.05, 0.2;

    const auto gaussian{
        ls::DecodeGaussianEstimateWithUncertainty(DatasetLinearizationSpec(), linear_model, covariance_matrix)
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

    EXPECT_NEAR(expected_amplitude, gaussian.estimate.amplitude, 1.0e-12);
    EXPECT_NEAR(expected_width, gaussian.estimate.width, 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_amplitude), gaussian.standard_deviation.amplitude, 1.0e-12);
    EXPECT_NEAR(std::sqrt(expected_var_width), gaussian.standard_deviation.width, 1.0e-12);
}

TEST(LinearizationServiceTest, DecodePosteriorFor2Beta3DModelAddsZeroIntercept)
{
    const auto linear_model{ MakeVector({ 0.5, 2.0 }) };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.05,
                         0.05, 0.2;

    const auto decoded{
        ls::DecodePosterior(DatasetLinearizationSpec(), linear_model, covariance_matrix)
    };
    const auto gaussian{ std::get<0>(decoded) };
    const auto standard_deviation{ std::get<1>(decoded) };

    ASSERT_EQ(gaussian.size(), rhbm_gem::GaussianModel3D::ParameterSize());
    ASSERT_EQ(standard_deviation.size(), rhbm_gem::GaussianModel3D::ParameterSize());
    EXPECT_NEAR(0.0, gaussian(rhbm_gem::GaussianModel3D::InterceptIndex()), 1.0e-12);
    EXPECT_NEAR(0.0, standard_deviation(rhbm_gem::GaussianModel3D::InterceptIndex()), 1.0e-12);
}

TEST(LinearizationServiceTest, DecodePosteriorForDirect3DModelUsesDiagonalStandardDeviation)
{
    const auto linear_model{ MakeVector({ 2.0, 0.5, 0.25 }) };
    Eigen::MatrixXd covariance_matrix{ Eigen::MatrixXd::Zero(3, 3) };
    covariance_matrix(0, 0) = 0.04;
    covariance_matrix(1, 1) = 0.09;
    covariance_matrix(2, 2) = 0.16;

    const auto decoded{
        ls::DecodePosterior(DirectGaussianModelSpec(), linear_model, covariance_matrix)
    };
    const auto gaussian{ std::get<0>(decoded) };
    const auto standard_deviation{ std::get<1>(decoded) };

    EXPECT_TRUE(gaussian.isApprox(linear_model, 1.0e-12));
    EXPECT_NEAR(0.2, standard_deviation(rhbm_gem::GaussianModel3D::AmplitudeIndex()), 1.0e-12);
    EXPECT_NEAR(0.3, standard_deviation(rhbm_gem::GaussianModel3D::WidthIndex()), 1.0e-12);
    EXPECT_NEAR(0.4, standard_deviation(rhbm_gem::GaussianModel3D::InterceptIndex()), 1.0e-12);
}
