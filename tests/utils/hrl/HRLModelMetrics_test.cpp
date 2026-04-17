#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <stdexcept>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/HRLModelMetrics.hpp>

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
} // namespace

TEST(HRLModelMetricsTest, CalculateNormalizedResidualReturnsExpectedRatio)
{
    const Eigen::VectorXd linear_estimate{ MakeVector({ std::log(2.0), 4.0 }) };
    const Eigen::VectorXd gaus_truth{ MakeVector({ 1.0, 0.4, 0.0 }) };

    const auto residual{
        HRLModelMetrics::CalculateNormalizedResidual(linear_estimate, gaus_truth)
    };

    const double expected_amplitude{
        std::exp(std::log(2.0)) * std::pow(Constants::two_pi / 4.0, 1.5)
    };
    const double expected_width{ 1.0 / std::sqrt(4.0) };

    ASSERT_EQ(residual.size(), 3);
    EXPECT_NEAR((expected_amplitude - 1.0) / 1.0, residual(0), 1e-12);
    EXPECT_NEAR((expected_width - 0.4) / 0.4, residual(1), 1e-12);
    EXPECT_TRUE(std::isnan(residual(2)));
}

TEST(HRLModelMetricsTest, CalculateNormalizedResidualRejectsDimensionMismatch)
{
    const Eigen::VectorXd linear_estimate{ MakeVector({ std::log(2.0), 4.0 }) };
    const Eigen::VectorXd wrong_truth{ MakeVector({ 1.0, 0.4 }) };

    EXPECT_THROW(
        HRLModelMetrics::CalculateNormalizedResidual(linear_estimate, wrong_truth),
        std::invalid_argument
    );
}

TEST(HRLModelMetricsTest, CalculateAbsoluteGausDifferenceReturnsAbsoluteDelta)
{
    const Eigen::VectorXd linear_a{ MakeVector({ std::log(2.0), 4.0 }) };
    const Eigen::VectorXd linear_b{ MakeVector({ std::log(1.0), 1.0 }) };

    const auto difference{
        HRLModelMetrics::CalculateAbsoluteGausDifference(linear_a, linear_b)
    };

    const double amplitude_a{
        2.0 * std::pow(Constants::two_pi / 4.0, 1.5)
    };
    const double width_a{ 1.0 / std::sqrt(4.0) };
    const double amplitude_b{
        1.0 * std::pow(Constants::two_pi / 1.0, 1.5)
    };
    const double width_b{ 1.0 };

    ASSERT_EQ(difference.size(), 3);
    EXPECT_NEAR(std::abs(amplitude_a - amplitude_b), difference(0), 1e-12);
    EXPECT_NEAR(std::abs(width_a - width_b), difference(1), 1e-12);
    EXPECT_DOUBLE_EQ(0.0, difference(2));
}
