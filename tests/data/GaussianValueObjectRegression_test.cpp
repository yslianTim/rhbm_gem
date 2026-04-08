#include <gtest/gtest.h>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

namespace rg = rhbm_gem;

TEST(GaussianValueObjectRegressionTest, LocalAndGroupIntensityStayEquivalent)
{
    const rg::GaussianEstimate estimate{ 9.0, 1.5 };
    const rg::GaussianEstimate variance{ 0.5, 0.2 };

    rg::LocalPotentialEntry local_entry;
    local_entry.SetEstimateMDPDE(estimate);
    local_entry.SetAnnotation(
        "component",
        rg::LocalPotentialAnnotation{ rg::GaussianPosterior{ estimate, variance }, false, 0.0 });

    rg::GroupPotentialEntry group_entry;
    group_entry.MarkAsAtomEntry();
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

TEST(GaussianValueObjectRegressionTest, GroupPotentialEntryRejectsMixedMemberKinds)
{
    rg::GroupPotentialEntry atom_entry;
    atom_entry.MarkAsAtomEntry();
    EXPECT_THROW(atom_entry.ReserveBondMembers(1, 1U), std::runtime_error);

    rg::GroupPotentialEntry bond_entry;
    bond_entry.MarkAsBondEntry();
    EXPECT_THROW(bond_entry.ReserveAtomMembers(1, 1U), std::runtime_error);
}
