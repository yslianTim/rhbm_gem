#include <gtest/gtest.h>

#include <vector>

#include <rhbm_gem/utils/algorithm/LinearRegressionSample.hpp>
#include <rhbm_gem/utils/algorithm/RobustSlopeEstimator.hpp>
#include <rhbm_gem/utils/algorithm/RobustSlopeOptions.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

alg::RobustSlopeOptions MakeOptions()
{
    alg::RobustSlopeOptions options;
    options.maximum_iterations = 50;
    options.tolerance = 1.0e-10;
    options.scale_multiplier = 1.4826;
    options.scale_min = 1.0e-12;
    options.cutoff_multiplier = 1.345;
    return options;
}

} // namespace

TEST(RobustSlopeEstimatorTest, OrdinarySlopeThroughOriginMatchesClosedForm)
{
    const std::vector<alg::LinearRegressionSample> samples{
        { 1.0, 2.0 },
        { 2.0, 4.0 },
        { 3.0, 6.0 }
    };

    double slope{ 0.0 };

    ASSERT_TRUE(alg::RobustSlopeEstimator::EstimateOrdinarySlopeThroughOrigin(samples, slope));
    EXPECT_DOUBLE_EQ(2.0, slope);
}

TEST(RobustSlopeEstimatorTest, HuberSlopeDownweightsLargeOutlier)
{
    const std::vector<alg::LinearRegressionSample> samples{
        { 1.0, 2.0 },
        { 2.0, 4.0 },
        { 3.0, 6.0 },
        { 4.0, 100.0 }
    };
    double ordinary_slope{ 0.0 };
    double robust_slope{ 0.0 };

    ASSERT_TRUE(alg::RobustSlopeEstimator::EstimateOrdinarySlopeThroughOrigin(
        samples,
        ordinary_slope));
    ASSERT_TRUE(alg::RobustSlopeEstimator::EstimateHuberSlopeThroughOrigin(
        samples,
        MakeOptions(),
        robust_slope));

    EXPECT_LT(robust_slope, ordinary_slope);
    EXPECT_NEAR(2.0, robust_slope, 5.0);
}

TEST(RobustSlopeEstimatorTest, RejectsEmptyAndDegenerateSamples)
{
    double slope{ 0.0 };

    EXPECT_FALSE(alg::RobustSlopeEstimator::EstimateHuberSlopeThroughOrigin(
        {},
        MakeOptions(),
        slope));
    EXPECT_FALSE(alg::RobustSlopeEstimator::EstimateOrdinarySlopeThroughOrigin(
        std::vector<alg::LinearRegressionSample>{
            { 0.0, 1.0 },
            { 0.0, 2.0 }
        },
        slope));
}
