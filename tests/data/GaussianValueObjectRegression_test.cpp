#include <gtest/gtest.h>

#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>

namespace rg = rhbm_gem;

TEST(GaussianValueObjectRegressionTest, LocalAndGroupIntensityStayEquivalent)
{
    const rg::GaussianEstimate estimate{ 9.0, 1.5 };
    const rg::GaussianEstimate variance{ 0.5, 0.2 };

    rg::LocalPotentialEntry local_entry;
    local_entry.SetEstimateMDPDE(estimate);
    local_entry.SetPosterior("component", rg::GaussianPosterior{ estimate, variance });

    rg::GroupPotentialEntry group_entry;
    auto & group_bucket{ group_entry.EnsureGroup(42) };
    group_bucket.prior = estimate;
    group_bucket.prior_variance = variance;

    const auto & group{ group_entry.GetGroup(42) };
    rg::GaussianPosterior group_posterior;
    group_posterior.estimate = group.prior;
    group_posterior.variance = group.prior_variance;

    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(local_entry.GetEstimateMDPDE().GetParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group.prior.Intensity(), estimate.Intensity());
    EXPECT_DOUBLE_EQ(group.prior.GetParameter(2), estimate.Intensity());
    EXPECT_DOUBLE_EQ(
        local_entry.GetPosterior("component").GetVariance(2),
        group_posterior.GetVariance(2));
}
