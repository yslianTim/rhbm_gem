#include <gtest/gtest.h>

#include <cmath>

#include <rhbm_gem/utils/hrl/GaussianStatistics.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

namespace rg = rhbm_gem;

TEST(GaussianStatisticsTest, LocalAndGroupIntensityStayEquivalent)
{
    const rg::GaussianEstimate estimate{ 9.0, 1.5 };
    const rg::GaussianEstimate variance{ 0.5, 0.2 };

    rg::LocalPotentialEntry local_entry;
    local_entry.SetEstimateMDPDE(estimate);
    local_entry.SetAnnotation(
        "component",
        rg::LocalPotentialAnnotation{ rg::GaussianPosterior{ estimate, variance }, false, 0.0 });

    rg::AtomGroupPotentialEntry group_entry;
    group_entry.SetGroupStatistics(42, {}, {}, estimate, variance, 0.0);
    const auto group_posterior{ group_entry.BuildPriorPosterior(42) };

    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().GetParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group_entry.GetPrior(42).Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group_entry.GetPrior(42).GetParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(
        local_entry.FindAnnotation("component")->posterior.GetVariance(2),
        group_posterior.GetVariance(2));
}

TEST(GaussianStatisticsTest, PublicGaussianStatisticsHeaderExposesStableValueMath)
{
    const rg::GaussianEstimate estimate{ 9.0, 1.5 };
    const rg::GaussianEstimate variance{ 0.5, 0.2 };
    const rg::GaussianPosterior posterior{ estimate, variance };

    const auto expected_intensity{
        estimate.amplitude
        * std::pow(Constants::two_pi * estimate.width * estimate.width, -1.5) };
    const auto beta{ estimate.ToBeta() };

    ASSERT_EQ(beta.size(), 3);
    EXPECT_DOUBLE_EQ(estimate.GetParameter(2), expected_intensity);
    EXPECT_DOUBLE_EQ(estimate.Intensity(), expected_intensity);
    EXPECT_DOUBLE_EQ(
        beta(0),
        std::log(estimate.amplitude)
            - 1.5 * std::log(Constants::two_pi * estimate.width * estimate.width));
    EXPECT_DOUBLE_EQ(beta(1), 1.0 / (estimate.width * estimate.width));
    EXPECT_DOUBLE_EQ(beta(2), 0.0);
    EXPECT_DOUBLE_EQ(posterior.GetEstimate(2), expected_intensity);
    EXPECT_DOUBLE_EQ(
        posterior.GetVariance(2),
        std::sqrt(
            std::pow(expected_intensity * variance.amplitude / estimate.amplitude, 2)
            + std::pow(
                -3.0 * estimate.amplitude * std::pow(Constants::two_pi, -1.5)
                    * std::pow(estimate.width, -4) * variance.width,
                2)));
}
