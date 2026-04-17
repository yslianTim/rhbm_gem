#include <gtest/gtest.h>

#include <vector>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/GaussianResponseMath.hpp>

namespace
{

Eigen::VectorXd MakePoint3D(double x, double y, double z)
{
    Eigen::VectorXd point{ Eigen::VectorXd::Zero(3) };
    point << x, y, z;
    return point;
}

} // namespace

TEST(GaussianResponseMathTest, GetGaussianResponseAtDistanceMatchesClosedForm)
{
    const auto response{
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtDistance(0.5, 0.25)
    };
    const auto width_square{ 0.25 * 0.25 };
    const auto expected{
        1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * 0.5 * 0.5 / width_square)
    };
    EXPECT_NEAR(expected, response, 1.0e-12);
}

TEST(GaussianResponseMathTest, GetGaussianResponseAtPointMatchesDistanceVersionAtCenterAxis)
{
    const auto point{ MakePoint3D(0.5, 0.0, 0.0) };
    const auto center{ MakePoint3D(0.0, 0.0, 0.0) };
    const auto point_response{
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtPoint(point, center, 0.25)
    };
    const auto distance_response{
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtDistance(0.5, 0.25)
    };
    EXPECT_NEAR(distance_response, point_response, 1.0e-12);
}

TEST(GaussianResponseMathTest, GetGaussianResponseAtPointWithNeighborhoodAddsNeighborContribution)
{
    const auto point{ MakePoint3D(0.0, 0.0, 0.0) };
    const auto center{ MakePoint3D(0.0, 0.0, 0.0) };
    const std::vector<Eigen::VectorXd> neighbors{ MakePoint3D(1.0, 0.0, 0.0) };

    const auto base_response{
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtPoint(point, center, 0.5)
    };
    const auto neighborhood_response{
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtPointWithNeighborhood(
            point,
            center,
            neighbors,
            0.5)
    };
    EXPECT_GT(neighborhood_response, base_response);
}

TEST(GaussianResponseMathTest, GaussianResponseRejectsNonPositiveWidth)
{
    const auto point{ MakePoint3D(0.0, 0.0, 0.0) };
    const auto center{ MakePoint3D(0.0, 0.0, 0.0) };
    const std::vector<Eigen::VectorXd> neighbors{ MakePoint3D(1.0, 0.0, 0.0) };

    EXPECT_THROW(
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtDistance(0.0, 0.0),
        std::runtime_error);
    EXPECT_THROW(
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtPoint(point, center, -1.0),
        std::runtime_error);
    EXPECT_THROW(
        rhbm_gem::GaussianResponseMath::GetGaussianResponseAtPointWithNeighborhood(
            point,
            center,
            neighbors,
            0.0),
        std::runtime_error);
}
