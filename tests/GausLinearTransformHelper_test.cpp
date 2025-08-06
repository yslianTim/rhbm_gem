#include <gtest/gtest.h>
#include <stdexcept>
#include <cmath>

#include "GausLinearTransformHelper.hpp"
#include "Constants.hpp"

TEST(GausLinearTransformHelperTest, BuildLinearModelDataVectorNormalCase)
{
    const double gaus_x{ 2.0 };
    const double gaus_y{ 4.0 };
    auto result{ GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y) };
    ASSERT_EQ(3, result.size());
    EXPECT_DOUBLE_EQ(1.0, result(0));
    EXPECT_DOUBLE_EQ(-0.5 * gaus_x * gaus_x, result(1));
    EXPECT_DOUBLE_EQ(std::log(gaus_y), result(2));
}

TEST(GausLinearTransformHelperTest, BuildLinearModelDataVectorThrowsForNonPositiveY)
{
    const double gaus_x{ 1.0 };
    EXPECT_THROW(
        GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, 0.0),
        std::runtime_error
    );
    EXPECT_THROW(
        GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, -1.0),
        std::runtime_error
    );
}

TEST(GausLinearTransformHelperTest, BuildGausModelPositiveBeta)
{
    Eigen::VectorXd linear_model(2);
    const double beta0{ 1.2 };
    const double beta1{ 3.4 };
    linear_model << beta0, beta1;
    auto gaus_model{ GausLinearTransformHelper::BuildGausModel(linear_model) };
    const double expected_amplitude{
        std::exp(beta0) * std::pow(Constants::two_pi / beta1, 1.5)
    };
    const double expected_width{ 1.0 / std::sqrt(beta1) };
    EXPECT_NEAR(expected_amplitude, gaus_model(0), 1e-9);
    EXPECT_NEAR(expected_width, gaus_model(1), 1e-9);
}

TEST(GausLinearTransformHelperTest, BuildGausModelNonPositiveBeta)
{
    Eigen::VectorXd linear_model(2);
    linear_model << 0.0, 0.0;
    auto gaus_model{ GausLinearTransformHelper::BuildGausModel(linear_model) };
    EXPECT_NEAR(0.0, gaus_model(0), 1e-12);
    EXPECT_NEAR(0.0, gaus_model(1), 1e-12);
}

TEST(GausLinearTransformHelperTest, BuildGausModelNonPositiveBeta1)
{
    // beta1 equal to zero should yield a zero vector
    Eigen::VectorXd linear_model(2);
    linear_model << 0.5, 0.0;
    auto gaus_model{ GausLinearTransformHelper::BuildGausModel(linear_model) };
    EXPECT_DOUBLE_EQ(0.0, gaus_model(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model(1));

    // beta1 negative should also yield a zero vector
    linear_model(1) = -1.0;
    gaus_model = GausLinearTransformHelper::BuildGausModel(linear_model);
    EXPECT_DOUBLE_EQ(0.0, gaus_model(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model(1));
}

TEST(GausLinearTransformHelperTest, BuildGausModelWithVariance)
{
    // Construct example linear model [beta0, beta1] with beta1 > 0
    const double beta0{ 0.5 };
    const double beta1{ 2.0 };
    Eigen::VectorXd linear_model(2);
    linear_model << beta0, beta1;

    // Covariance matrix
    const double var_beta0{ 0.1 };
    const double var_beta1{ 0.2 };
    const double cov_beta0_beta1{ 0.05 };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << var_beta0, cov_beta0_beta1,
                          cov_beta0_beta1, var_beta1;

    // Expected Gaussian model values based on implementation formulas
    const double expected_amplitude{ std::exp(beta0) * std::pow(Constants::two_pi / beta1, 1.5) };
    const double expected_width{ 1.0 / std::sqrt(beta1) };
    const double var_amplitude{
        std::pow(Constants::two_pi, 3) * std::exp(2.0 * beta0) * std::pow(beta1, -5) *
        (std::pow(beta1, 2) * var_beta0 + 2.25 * var_beta1 - 3.0 * beta1 * cov_beta0_beta1)
    };
    const double var_width{ 0.25 * std::pow(beta1, -3) * var_beta1 };
    auto [gaus_model, gaus_model_variance]{
        GausLinearTransformHelper::BuildGausModelWithVariance(linear_model, covariance_matrix)
    };

    EXPECT_NEAR(expected_amplitude, gaus_model(0), 1e-9);
    EXPECT_NEAR(expected_width, gaus_model(1), 1e-9);
    EXPECT_NEAR(std::sqrt(var_amplitude), gaus_model_variance(0), 1e-9);
    EXPECT_NEAR(std::sqrt(var_width), gaus_model_variance(1), 1e-9);
}
