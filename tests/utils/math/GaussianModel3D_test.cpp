#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rg = rhbm_gem;

namespace
{

Eigen::VectorXd MakeVector(std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index{ 0 };
    for (const auto value : values)
    {
        result(index++) = value;
    }
    return result;
}

} // namespace

TEST(GaussianModel3DTest, RoundTripsThroughCanonicalVector)
{
    const rg::GaussianModel3D model{ 9.0, 1.5, 0.25 };

    const auto parameters{ model.ToVector() };
    const auto round_trip{ rg::GaussianModel3D::FromVector(parameters) };

    ASSERT_EQ(parameters.size(), rg::GaussianModel3D::ParameterSize());
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::AmplitudeIndex()), model.GetAmplitude());
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::WidthIndex()), model.GetWidth());
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::InterceptIndex()), model.GetIntercept());
    EXPECT_DOUBLE_EQ(round_trip.GetAmplitude(), model.GetAmplitude());
    EXPECT_DOUBLE_EQ(round_trip.GetWidth(), model.GetWidth());
    EXPECT_DOUBLE_EQ(round_trip.GetIntercept(), model.GetIntercept());
}

TEST(GaussianModel3DTest, DefaultsToZeroIntercept)
{
    const rg::GaussianModel3D model{};

    EXPECT_DOUBLE_EQ(model.GetAmplitude(), 0.0);
    EXPECT_DOUBLE_EQ(model.GetWidth(), 1.0);
    EXPECT_DOUBLE_EQ(model.GetIntercept(), 0.0);
    EXPECT_DOUBLE_EQ(model.GetModelParameter(rg::GaussianModel3D::InterceptIndex()), 0.0);
}

TEST(GaussianModel3DTest, FromVectorRequiresExactlyThreeFiniteEntries)
{
    EXPECT_NO_THROW(rg::GaussianModel3D::FromVector(MakeVector({ 1.0, 0.5, 0.0 })));
    EXPECT_THROW(rg::GaussianModel3D::FromVector(MakeVector({ 1.0, 0.5 })),
                 std::invalid_argument);
    EXPECT_THROW(rg::GaussianModel3D::FromVector(MakeVector({ 1.0, 0.5, 0.0, 2.0 })),
                 std::invalid_argument);
    EXPECT_THROW(rg::GaussianModel3D::FromVector(
                     MakeVector({ 1.0, std::numeric_limits<double>::quiet_NaN(), 0.0 })),
                 std::invalid_argument);
}

TEST(GaussianModel3DTest, FromVectorPrefixAllowsLongerLegacyVectors)
{
    const auto model{
        rg::GaussianModel3D::FromVectorPrefix(MakeVector({ 1.0, 0.5, 0.25, 8.0 })) };

    EXPECT_DOUBLE_EQ(model.GetAmplitude(), 1.0);
    EXPECT_DOUBLE_EQ(model.GetWidth(), 0.5);
    EXPECT_DOUBLE_EQ(model.GetIntercept(), 0.25);
}

TEST(GaussianModel3DTest, RejectsInvalidSamplingModels)
{
    EXPECT_NO_THROW(rg::GaussianModel3D::RequireFinitePositiveWidthModel(
        rg::GaussianModel3D{ 1.0, 0.5, 0.0 }));
    EXPECT_THROW(
        rg::GaussianModel3D::RequireFinitePositiveWidthModel(
            rg::GaussianModel3D{ 1.0, 0.0, 0.0 }),
        std::invalid_argument);
    EXPECT_THROW(
        rg::GaussianModel3D::RequireFinitePositiveWidthModel(
            rg::GaussianModel3D{ 1.0, -0.5, 0.0 }),
        std::invalid_argument);
    EXPECT_THROW(
        rg::GaussianModel3D::RequireFinitePositiveWidthModel(
            rg::GaussianModel3D{ 1.0, std::numeric_limits<double>::infinity(), 0.0 }),
        std::invalid_argument);
}

TEST(GaussianModel3DTest, WithMethodsReturnAdjustedCopies)
{
    const rg::GaussianModel3D model{ 9.0, 1.5, 0.25 };

    EXPECT_DOUBLE_EQ(model.WithAmplitude(8.0).GetAmplitude(), 8.0);
    EXPECT_DOUBLE_EQ(model.WithWidth(0.75).GetWidth(), 0.75);
    EXPECT_DOUBLE_EQ(model.WithIntercept(0.5).GetIntercept(), 0.5);
    EXPECT_DOUBLE_EQ(model.GetAmplitude(), 9.0);
    EXPECT_DOUBLE_EQ(model.GetWidth(), 1.5);
    EXPECT_DOUBLE_EQ(model.GetIntercept(), 0.25);
}

TEST(GaussianModel3DTest, IntensityAndResponseMatchClosedForm)
{
    const rg::GaussianModel3D model{ 9.0, 1.5, 0.25 };
    constexpr double distance{ 0.75 };
    const auto expected_intensity{
        model.GetAmplitude()
        * std::pow(Constants::two_pi * model.GetWidth() * model.GetWidth(), -1.5) };
    const auto expected_response{
        expected_intensity
            * std::exp(-0.5 * distance * distance / (model.GetWidth() * model.GetWidth())) +
        model.GetIntercept() };

    EXPECT_DOUBLE_EQ(model.Intensity(), expected_intensity);
    EXPECT_DOUBLE_EQ(model.ResponseAtDistance(distance), expected_response);
}

TEST(GaussianModel3DTest, ZeroWidthKeepsExistingFallbackBehavior)
{
    const rg::GaussianModel3D model{ 9.0, 0.0, 0.25 };

    EXPECT_DOUBLE_EQ(model.Intensity(), 0.0);
    EXPECT_DOUBLE_EQ(model.ResponseAtDistance(0.75), model.GetIntercept());
}

TEST(GaussianModel3DTest, DisplayParameterKeepsIntensitySeparateFromModelIntercept)
{
    const rg::GaussianModel3D model{ 9.0, 1.5, 0.25 };

    EXPECT_DOUBLE_EQ(
        model.GetModelParameter(rg::GaussianModel3D::InterceptIndex()),
        model.GetIntercept());
    EXPECT_DOUBLE_EQ(
        model.GetDisplayParameter(rg::GaussianModel3D::InterceptIndex()),
        model.Intensity());
}

TEST(GaussianModel3DTest, UncertaintyRoundTripsThroughCanonicalVector)
{
    const rg::GaussianModel3DUncertainty uncertainty{ 0.5, 0.2, 0.1 };

    const auto parameters{ uncertainty.ToVector() };

    ASSERT_EQ(parameters.size(), rg::GaussianModel3D::ParameterSize());
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::AmplitudeIndex()), 0.5);
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::WidthIndex()), 0.2);
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::InterceptIndex()), 0.1);
}

TEST(GaussianModel3DTest, UncertaintyRequiresFiniteNonNegativeValues)
{
    EXPECT_NO_THROW(rg::GaussianModel3DUncertainty::RequireFiniteNonNegativeUncertainty(
        rg::GaussianModel3DUncertainty{ 0.5, 0.2, 0.0 }));
    EXPECT_THROW(
        rg::GaussianModel3DUncertainty::RequireFiniteNonNegativeUncertainty(
            rg::GaussianModel3DUncertainty{ -0.5, 0.2, 0.0 }),
        std::invalid_argument);
    EXPECT_THROW(
        rg::GaussianModel3DUncertainty::RequireFiniteNonNegativeUncertainty(
            rg::GaussianModel3DUncertainty{
                0.5,
                std::numeric_limits<double>::infinity(),
                0.0 }),
        std::invalid_argument);
}

TEST(GaussianModel3DTest, WithUncertaintyComputesIntensityStandardDeviation)
{
    const rg::GaussianModel3D model{ 9.0, 1.5 };
    const rg::GaussianModel3DUncertainty standard_deviation{ 0.5, 0.2 };
    const rg::GaussianModel3DWithUncertainty gaussian{ model, standard_deviation };

    const auto expected_intensity{
        model.GetAmplitude()
        * std::pow(Constants::two_pi * model.GetWidth() * model.GetWidth(), -1.5) };
    const auto expected_intensity_standard_deviation{
        std::sqrt(
            std::pow(
                expected_intensity * standard_deviation.GetAmplitude() / model.GetAmplitude(),
                2)
            + std::pow(
                -3.0 * model.GetAmplitude() * std::pow(Constants::two_pi, -1.5)
                    * std::pow(model.GetWidth(), -4) * standard_deviation.GetWidth(),
                2))
    };

    EXPECT_DOUBLE_EQ(gaussian.GetModel().GetAmplitude(), model.GetAmplitude());
    EXPECT_DOUBLE_EQ(
        gaussian.GetStandardDeviationModel().GetWidth(),
        standard_deviation.GetWidth());
    EXPECT_DOUBLE_EQ(
        gaussian.GetDisplayStandardDeviation(rg::GaussianModel3D::InterceptIndex()),
        expected_intensity_standard_deviation);
}
