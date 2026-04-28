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

    ASSERT_EQ(parameters.size(), rg::GaussianModel3D::kParameterSize);
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::kAmplitudeIndex), model.amplitude);
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::kWidthIndex), model.width);
    EXPECT_DOUBLE_EQ(parameters(rg::GaussianModel3D::kInterceptIndex), model.intercept);
    EXPECT_DOUBLE_EQ(round_trip.amplitude, model.amplitude);
    EXPECT_DOUBLE_EQ(round_trip.width, model.width);
    EXPECT_DOUBLE_EQ(round_trip.intercept, model.intercept);
}

TEST(GaussianModel3DTest, DefaultsToZeroIntercept)
{
    const rg::GaussianModel3D model{};

    EXPECT_DOUBLE_EQ(model.amplitude, 0.0);
    EXPECT_DOUBLE_EQ(model.width, 1.0);
    EXPECT_DOUBLE_EQ(model.intercept, 0.0);
    EXPECT_DOUBLE_EQ(model.GetModelParameter(rg::GaussianModel3D::kInterceptIndex), 0.0);
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

    EXPECT_DOUBLE_EQ(model.amplitude, 1.0);
    EXPECT_DOUBLE_EQ(model.width, 0.5);
    EXPECT_DOUBLE_EQ(model.intercept, 0.25);
}

TEST(GaussianModel3DTest, IntensityAndResponseMatchClosedForm)
{
    const rg::GaussianModel3D model{ 9.0, 1.5, 0.25 };
    constexpr double distance{ 0.75 };
    const auto expected_intensity{
        model.amplitude * std::pow(Constants::two_pi * model.width * model.width, -1.5) };
    const auto expected_response{
        expected_intensity * std::exp(-0.5 * distance * distance / (model.width * model.width)) +
        model.intercept };

    EXPECT_DOUBLE_EQ(model.Intensity(), expected_intensity);
    EXPECT_DOUBLE_EQ(model.ResponseAtDistance(distance), expected_response);
}

TEST(GaussianModel3DTest, ZeroWidthKeepsExistingFallbackBehavior)
{
    const rg::GaussianModel3D model{ 9.0, 0.0, 0.25 };

    EXPECT_DOUBLE_EQ(model.Intensity(), 0.0);
    EXPECT_DOUBLE_EQ(model.ResponseAtDistance(0.75), model.intercept);
}
