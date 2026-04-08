#include <gtest/gtest.h>

#include <type_traits>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

namespace rg = rhbm_gem;

namespace {

template <typename EntryT, typename MemberT, typename = void>
struct HasTypedAddMember : std::false_type {};

template <typename EntryT, typename MemberT>
struct HasTypedAddMember<
    EntryT,
    MemberT,
    std::void_t<decltype(std::declval<EntryT &>().AddMember(
        std::declval<GroupKey>(),
        std::declval<MemberT &>()))>> : std::true_type {};

static_assert(HasTypedAddMember<rg::AtomGroupPotentialEntry, rg::AtomObject>::value);
static_assert(!HasTypedAddMember<rg::AtomGroupPotentialEntry, rg::BondObject>::value);
static_assert(HasTypedAddMember<rg::BondGroupPotentialEntry, rg::BondObject>::value);
static_assert(!HasTypedAddMember<rg::BondGroupPotentialEntry, rg::AtomObject>::value);

} // namespace

TEST(GaussianValueObjectRegressionTest, LocalAndGroupIntensityStayEquivalent)
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
