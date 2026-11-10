#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "core/estimator/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/core/EstimatorTester.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>

namespace {
constexpr double kExpectedLocalFittingRadius{1.0};

namespace rt = rhbm_gem::core;
namespace rt_detail = rhbm_gem::core::detail;
namespace tdf = rhbm_gem::core;
namespace rg = rhbm_gem;

tdf::GaussianParameterDistribution MakeDistribution(
    const rg::GaussianModel3D & mean,
    const rg::GaussianModel3DUncertainty & sigma = rg::GaussianModel3DUncertainty{ 0.05, 0.025, 0.01 })
{
    return tdf::GaussianParameterDistribution{ mean, sigma };
}

rt::LocalTestOptions MakeLocalOptions(
    double alpha_r,
    bool alpha_training)
{
    rt::LocalTestOptions options;
    options.requested_alpha_r = alpha_r;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

rt::GroupTestOptions MakeGroupOptions(
    double alpha_g,
    bool alpha_training)
{
    rt::GroupTestOptions options;
    options.requested_alpha_g = alpha_g;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

void ExpectBiasStatisticSize(const rt::BiasStatistics & bias)
{
    EXPECT_EQ(bias.mean.size(), rg::GaussianModel3D::ParameterSize());
    EXPECT_EQ(bias.sigma.size(), rg::GaussianModel3D::ParameterSize());
}

} // namespace

TEST(EstimatorTesterTest, PreparedLocalGaussianDatasetMatchesLegacyBuilder)
{
    LocalPotentialSampleList samples{
        { 0.6, SamplingPoint{ 0.0 } },
        { 0.8, SamplingPoint{ 0.25 } },
        { -1.0, SamplingPoint{ 0.5 } },
        { 0.7, SamplingPoint{ 1.0 } },
        { 0.9, SamplingPoint{ 1.25 } }
    };
    constexpr double range_min{ 0.25 };
    constexpr double range_max{ 1.0 };
    const rg::GaussianModel3D offset_model{ 1.0, 0.5, 0.1 };
    auto adjusted_samples{ samples };
    std::vector<double> response_list;
    response_list.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); i++)
    {
        const auto evaluation{
            offset_model.EvaluateAtDistance(samples.at(i).point.distance)
        };
        adjusted_samples.at(i).response = samples.at(i).response -
            (evaluation.response - evaluation.signal);
        response_list.emplace_back(samples.at(i).response);
    }

    const auto legacy_dataset{
        rhbm_gem::rhbm_helper::BuildMemberDataset(
            adjusted_samples,
            range_min,
            range_max)
    };
    const rt_detail::PreparedLocalGaussianDesign design{
        samples,
        range_min,
        range_max
    };
    const auto prepared_dataset{
        design.BuildDataset(response_list, offset_model)
    };

    EXPECT_TRUE(prepared_dataset.X.isApprox(legacy_dataset.X, 0.0));
    EXPECT_TRUE(prepared_dataset.y.isApprox(legacy_dataset.y, 0.0));

    const rt_detail::PreparedLocalGaussianDesign signal_design{
        samples, 0.0, kExpectedLocalFittingRadius
    };
    const auto signal_dataset{ signal_design.BuildDataset(response_list, offset_model) };
    ASSERT_EQ(signal_dataset.X.rows(), 3);
    EXPECT_DOUBLE_EQ(signal_dataset.X(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(signal_dataset.X(2, 1), -0.5);
    const rt::FitOptions options;
    const auto fitted{ rt::EstimateLocalGaussian(samples, 0.0, options, offset_model) };
    const auto expected{ signal_design.Estimate(response_list, 0.0, options.thread_size, offset_model) };
    EXPECT_DOUBLE_EQ(fitted.mdpde.GetModel().GetAmplitude(), expected.mdpde.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(fitted.mdpde.GetModel().GetWidth(), expected.mdpde.GetModel().GetWidth());


    std::fill(response_list.begin(), response_list.end(), -1.0);
    for (auto & sample : adjusted_samples) sample.response = -1.0;
    const auto legacy_fallback{
        rhbm_gem::rhbm_helper::BuildMemberDataset(
            adjusted_samples,
            range_min,
            range_max)
    };
    const auto prepared_fallback{
        design.BuildDataset(response_list, offset_model)
    };
    EXPECT_TRUE(prepared_fallback.X.isApprox(legacy_fallback.X, 0.0));
    EXPECT_TRUE(prepared_fallback.y.isApprox(legacy_fallback.y, 0.0));
}

TEST(EstimatorTesterTest, RunLocalEstimationTestPopulatesBiasOutputs)
{
    constexpr double alpha_r{ 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(alpha_r, true))
    };

    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_median.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_median.value(), 0.0);
}

TEST(EstimatorTesterTest, RunGroupEstimationTestPopulatesBiasOutputs)
{
    constexpr double alpha_g{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunGroupEstimationTest(test_input, MakeGroupOptions(alpha_g, true))
    };

    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_median.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_median.value(), 0.0);
}

TEST(EstimatorTesterTest, RunLocalEstimationTestSkipsTrainedAlphaWhenDisabled)
{
    constexpr double alpha_r{ 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(alpha_r, false))
    };

    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_median.has_value());
}

TEST(EstimatorTesterTest, RunGroupEstimationTestSkipsTrainedAlphaWhenDisabled)
{
    constexpr double alpha_g{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunGroupEstimationTest(test_input, MakeGroupOptions(alpha_g, false))
    };

    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_median.has_value());
}

TEST(EstimatorTesterTest, RunLocalEstimationTestRejectsNonFiniteTruth)
{
    auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };
    test_input.gaus_true = rg::GaussianModel3D{
        std::numeric_limits<double>::quiet_NaN(),
        0.5,
        0.0
    };

    EXPECT_THROW(
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(0.5, true)),
        std::invalid_argument
    );
}
