#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace {
namespace rg = rhbm_gem;
namespace rgc = rhbm_gem::core;

rgc::FitOptions MakeFitOptions()
{
    rgc::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 0.9;
    options.quiet_mode = true;
    return options;
}

LocalPotentialSampleList BuildSamples(
    const rg::GaussianModel3D & model,
    const std::vector<double> & distance_list)
{
    LocalPotentialSampleList sample_entries;
    sample_entries.reserve(distance_list.size());
    for (const auto distance : distance_list)
    {
        sample_entries.emplace_back(LocalPotentialSample{
            static_cast<float>(model.ResponseAtDistance(distance)),
            SamplingPoint{ static_cast<float>(distance) }
        });
    }
    return sample_entries;
}

LocalPotentialSampleList BuildFullRangeSamples(const rg::GaussianModel3D & model)
{
    return BuildSamples(
        model,
        {
            0.0, 0.1, 0.2, 0.3, 0.4,
            0.5, 0.6, 0.7, 0.8, 0.9,
            1.1, 1.2, 1.3, 1.4, 1.5,
            1.6, 1.7, 1.8, 1.9, 2.0
        });
}

} // namespace

TEST(GaussianEstimatorTest, EstimatesNonZeroInterceptFromNoiselessSamples)
{
    const rg::GaussianModel3D truth{ 1.2, 0.45, 0.2 };
    const auto result{
        rgc::EstimateLocalGaussianWithIntercept(
            BuildFullRangeSamples(truth),
            0.5,
            MakeFitOptions())
    };
    const auto & estimate{ result.mdpde.GetModel() };

    EXPECT_NEAR(truth.GetAmplitude(), estimate.GetAmplitude(), 1.0e-3);
    EXPECT_NEAR(truth.GetWidth(), estimate.GetWidth(), 1.0e-3);
    EXPECT_NEAR(truth.GetIntercept(), estimate.GetIntercept(), 1.0e-3);
}

TEST(GaussianEstimatorTest, UsesAlphaRForRobustInterceptEstimation)
{
    const rg::GaussianModel3D truth{ 1.2, 0.45, 0.2 };
    auto sample_entries{ BuildFullRangeSamples(truth) };
    for (auto & sample : sample_entries)
    {
        if (std::abs(static_cast<double>(sample.point.distance) - 1.5) < 1.0e-6)
        {
            sample.response += 2.0f;
        }
    }

    const auto ols_result{
        rgc::EstimateLocalGaussianWithIntercept(
            sample_entries,
            0.0,
            MakeFitOptions())
    };
    const auto robust_result{
        rgc::EstimateLocalGaussianWithIntercept(
            sample_entries,
            1.0,
            MakeFitOptions())
    };
    const auto ols_error{
        std::abs(ols_result.mdpde.GetModel().GetIntercept() - truth.GetIntercept())
    };
    const auto robust_error{
        std::abs(robust_result.mdpde.GetModel().GetIntercept() - truth.GetIntercept())
    };

    EXPECT_LT(robust_error, ols_error);
}

TEST(GaussianEstimatorTest, KeepsClampedInitialInterceptWithoutEnoughResidualSamples)
{
    const rg::GaussianModel3D truth{ 1.2, 0.45, 0.2 };
    const auto sample_entries{
        BuildSamples(truth, { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 })
    };
    const auto result{
        rgc::EstimateLocalGaussianWithIntercept(
            sample_entries,
            0.5,
            MakeFitOptions(),
            2.0)
    };

    EXPECT_DOUBLE_EQ(1.0, result.mdpde.GetModel().GetIntercept());
}

TEST(GaussianEstimatorTest, IgnoresResidualSamplesOutsideInterceptRange)
{
    const rg::GaussianModel3D truth{ 1.2, 0.45, 0.2 };
    const auto baseline_result{
        rgc::EstimateLocalGaussianWithIntercept(
            BuildFullRangeSamples(truth),
            0.5,
            MakeFitOptions())
    };
    auto contaminated_samples{ BuildFullRangeSamples(truth) };
    contaminated_samples.emplace_back(LocalPotentialSample{
        static_cast<float>(truth.ResponseAtDistance(2.5) + 100.0),
        SamplingPoint{ 2.5f }
    });
    const auto contaminated_result{
        rgc::EstimateLocalGaussianWithIntercept(
            contaminated_samples,
            0.5,
            MakeFitOptions())
    };

    EXPECT_NEAR(
        baseline_result.mdpde.GetModel().GetIntercept(),
        contaminated_result.mdpde.GetModel().GetIntercept(),
        1.0e-12);
}
