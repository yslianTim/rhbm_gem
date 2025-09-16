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

TEST(GausLinearTransformHelperTest, BuildLinearModelDataVectorNegativeX)
{
    const double gaus_x{ 2.0 };
    const double gaus_y{ 4.0 };
    auto positive_result{ GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y) };
    auto negative_result{ GausLinearTransformHelper::BuildLinearModelDataVector(-gaus_x, gaus_y) };
    ASSERT_EQ(positive_result.size(), negative_result.size());
    EXPECT_DOUBLE_EQ(positive_result(0), negative_result(0));
    EXPECT_DOUBLE_EQ(positive_result(1), negative_result(1));
    EXPECT_DOUBLE_EQ(positive_result(2), negative_result(2));
}

TEST(GausLinearTransformHelperTest, BuildLinearModelDataVectorZeroX)
{
    const double gaus_x{ 0.0 };
    const double gaus_y{ 3.0 };
    auto result{ GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y) };
    ASSERT_EQ(3, result.size());
    EXPECT_DOUBLE_EQ(1.0, result(0));
    EXPECT_DOUBLE_EQ(0.0, result(1));
    EXPECT_DOUBLE_EQ(std::log(gaus_y), result(2));
}

TEST(GausLinearTransformHelperTest, BuildLinearModelDataVectorSmallPositiveY)
{
    const double gaus_x{ 1.0 };
    const double gaus_y{ 1e-10 };
    auto result{ GausLinearTransformHelper::BuildLinearModelDataVector(gaus_x, gaus_y) };
    ASSERT_EQ(3, result.size());
    EXPECT_DOUBLE_EQ(1.0, result(0));
    EXPECT_DOUBLE_EQ(-0.5 * gaus_x * gaus_x, result(1));
    EXPECT_DOUBLE_EQ(std::log(1e-10), result(2));
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

TEST(GausLinearTransformHelperTest, BuildGaus3DModelNegativeBeta0)
{
    Eigen::VectorXd linear_model(2);
    const double beta0{ -0.5 };
    const double beta1{ 2.3 };
    linear_model << beta0, beta1;
    auto gaus_model{ GausLinearTransformHelper::BuildGaus3DModel(linear_model) };
    const double expected_amplitude{ std::exp(beta0) * std::pow(Constants::two_pi / beta1, 1.5) };
    const double expected_width{ 1.0 / std::sqrt(beta1) };
    EXPECT_NEAR(expected_amplitude, gaus_model(0), 1e-9);
    EXPECT_NEAR(expected_width, gaus_model(1), 1e-9);
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelPositiveBeta)
{
    Eigen::VectorXd linear_model(2);
    const double beta0{ 1.2 };
    const double beta1{ 3.4 };
    linear_model << beta0, beta1;
    auto gaus_model{ GausLinearTransformHelper::BuildGaus3DModel(linear_model) };
    const double expected_amplitude{
        std::exp(beta0) * std::pow(Constants::two_pi / beta1, 1.5)
    };
    const double expected_width{ 1.0 / std::sqrt(beta1) };
    EXPECT_NEAR(expected_amplitude, gaus_model(0), 1e-9);
    EXPECT_NEAR(expected_width, gaus_model(1), 1e-9);
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelNonPositiveBeta)
{
    Eigen::VectorXd linear_model(2);
    linear_model << 0.0, 0.0;
    auto gaus_model{ GausLinearTransformHelper::BuildGaus3DModel(linear_model) };
    EXPECT_NEAR(0.0, gaus_model(0), 1e-12);
    EXPECT_NEAR(0.0, gaus_model(1), 1e-12);
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelNonPositiveBeta1)
{
    // beta1 equal to zero should yield a zero vector
    Eigen::VectorXd linear_model(2);
    linear_model << 0.5, 0.0;
    auto gaus_model{ GausLinearTransformHelper::BuildGaus3DModel(linear_model) };
    EXPECT_DOUBLE_EQ(0.0, gaus_model(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model(1));

    // beta1 negative should also yield a zero vector
    linear_model(1) = -1.0;
    gaus_model = GausLinearTransformHelper::BuildGaus3DModel(linear_model);
    EXPECT_DOUBLE_EQ(0.0, gaus_model(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model(1));
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelWithVariance)
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
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, covariance_matrix)
    };

    EXPECT_NEAR(expected_amplitude, gaus_model(0), 1e-9);
    EXPECT_NEAR(expected_width, gaus_model(1), 1e-9);
    EXPECT_NEAR(std::sqrt(var_amplitude), gaus_model_variance(0), 1e-9);
    EXPECT_NEAR(std::sqrt(var_width), gaus_model_variance(1), 1e-9);
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelWithVarianceNonPositiveBeta1)
{
    Eigen::VectorXd linear_model(2);
    Eigen::MatrixXd covariance_matrix{ Eigen::MatrixXd::Zero(2, 2) };

    // beta1 equal to zero should yield zero model and variance
    linear_model << 0.5, 0.0;
    auto [gaus_model, gaus_model_variance]{
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, covariance_matrix)
    };
    EXPECT_DOUBLE_EQ(0.0, gaus_model(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model(1));
    EXPECT_DOUBLE_EQ(0.0, gaus_model_variance(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model_variance(1));

    // beta1 negative should also yield zero model and variance
    linear_model(1) = -1.0;
    std::tie(gaus_model, gaus_model_variance) =
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, covariance_matrix);
    EXPECT_DOUBLE_EQ(0.0, gaus_model(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model(1));
    EXPECT_DOUBLE_EQ(0.0, gaus_model_variance(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model_variance(1));
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelWithVarianceThrowsForNon2x2Covariance)
{
    Eigen::VectorXd linear_model(2);
    linear_model << 0.5, 2.0;
    Eigen::MatrixXd bad_matrix(3, 3);
    bad_matrix.setIdentity();
    EXPECT_THROW(
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, bad_matrix),
        std::invalid_argument
    );
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelWithVarianceNonSymmetric)
{
    Eigen::VectorXd linear_model(2);
    linear_model << 0.5, 2.0;
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << 0.1, 0.2,
                         0.3, 0.4;
    EXPECT_THROW(
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, covariance_matrix),
        std::invalid_argument
    );
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelWithVarianceZeroCovariance)
{
    // Example linear model with positive beta1 so BuildGaus3DModel returns non-zero values
    const double beta0{ 0.3 };
    const double beta1{ 1.7 };
    Eigen::VectorXd linear_model(2);
    linear_model << beta0, beta1;

    // Zero covariance matrix
    Eigen::MatrixXd covariance_matrix{ Eigen::MatrixXd::Zero(2, 2) };
    auto [gaus_model, gaus_model_variance]{
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, covariance_matrix)
    };
    auto expected_gaus_model{ GausLinearTransformHelper::BuildGaus3DModel(linear_model) };
    EXPECT_DOUBLE_EQ(expected_gaus_model(0), gaus_model(0));
    EXPECT_DOUBLE_EQ(expected_gaus_model(1), gaus_model(1));
    EXPECT_DOUBLE_EQ(0.0, gaus_model_variance(0));
    EXPECT_DOUBLE_EQ(0.0, gaus_model_variance(1));
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelWithVarianceNegativeVariance)
{
    // Linear model with beta1 > 0 ensures non-zero Gaussian model
    const double beta0{ 0.5 };
    const double beta1{ 2.0 };
    Eigen::VectorXd linear_model(2);
    linear_model << beta0, beta1;

    // Covariance matrix with negative variances on the diagonal
    const double var_beta0{ -0.1 };
    const double var_beta1{ -0.2 };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << var_beta0, 0.0, 0.0, var_beta1;
    auto [gaus_model, gaus_model_variance]{
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, covariance_matrix)
    };
    auto expected_gaus_model{ GausLinearTransformHelper::BuildGaus3DModel(linear_model) };
    EXPECT_DOUBLE_EQ(expected_gaus_model(0), gaus_model(0));
    EXPECT_DOUBLE_EQ(expected_gaus_model(1), gaus_model(1));
    EXPECT_TRUE(std::isnan(gaus_model_variance(0)));
    EXPECT_TRUE(std::isnan(gaus_model_variance(1)));
}

TEST(GausLinearTransformHelperTest, BuildGaus3DModelWithVarianceLargeCovariance)
{
    // Linear model with positive beta1
    const double beta0{ 0.5 };
    const double beta1{ 2.0 };
    Eigen::VectorXd linear_model(2);
    linear_model << beta0, beta1;

    // Covariance matrix with large off-diagonal covariance leading to negative amplitude variance
    const double var_beta0{ 0.1 };
    const double var_beta1{ 0.2 };
    const double cov_beta0_beta1{ 0.8 };
    Eigen::MatrixXd covariance_matrix(2, 2);
    covariance_matrix << var_beta0, cov_beta0_beta1,
                         cov_beta0_beta1, var_beta1;

    auto [gaus_model, gaus_model_variance]{
        GausLinearTransformHelper::BuildGaus3DModelWithVariance(linear_model, covariance_matrix)
    };
    auto expected_gaus_model{ GausLinearTransformHelper::BuildGaus3DModel(linear_model) };
    EXPECT_DOUBLE_EQ(expected_gaus_model(0), gaus_model(0));
    EXPECT_DOUBLE_EQ(expected_gaus_model(1), gaus_model(1));
    EXPECT_TRUE(std::isnan(gaus_model_variance(0)));
    const double expected_sd_width{ std::sqrt(0.25 * std::pow(beta1, -3) * var_beta1) };
    EXPECT_DOUBLE_EQ(expected_sd_width, gaus_model_variance(1));
}
