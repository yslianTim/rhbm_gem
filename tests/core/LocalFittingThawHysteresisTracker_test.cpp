#include <gtest/gtest.h>

#include "detail/LocalFittingThawHysteresisTracker.hpp"

namespace {
namespace detail = rhbm_gem::core::detail;

} // namespace

TEST(LocalFittingThawHysteresisTrackerTest, StartsAtBaseThreshold)
{
    detail::LocalFittingThawHysteresisTracker tracker{ 2, 2.0, 8.0, 0.9 };

    EXPECT_DOUBLE_EQ(0.25, tracker.GetThreshold(1, 0.25));
    EXPECT_FALSE(tracker.ShouldThaw(1, 0.249, 0.25));
    EXPECT_TRUE(tracker.ShouldThaw(1, 0.25, 0.25));
}

TEST(LocalFittingThawHysteresisTrackerTest, DoublesThresholdAfterDependencyThaw)
{
    detail::LocalFittingThawHysteresisTracker tracker{ 1, 2.0, 8.0, 0.9 };

    tracker.RecordDependencyThaw(0);

    EXPECT_DOUBLE_EQ(0.2, tracker.GetThreshold(0, 0.1));
    EXPECT_FALSE(tracker.ShouldThaw(0, 0.199, 0.1));
    EXPECT_TRUE(tracker.ShouldThaw(0, 0.2, 0.1));
}

TEST(LocalFittingThawHysteresisTrackerTest, CapsThresholdMultiplier)
{
    detail::LocalFittingThawHysteresisTracker tracker{ 1, 2.0, 8.0, 0.9 };

    tracker.RecordDependencyThaw(0);
    tracker.RecordDependencyThaw(0);
    tracker.RecordDependencyThaw(0);
    tracker.RecordDependencyThaw(0);

    EXPECT_DOUBLE_EQ(0.8, tracker.GetThreshold(0, 0.1));
}

TEST(LocalFittingThawHysteresisTrackerTest, DecaysFrozenThresholdTowardBase)
{
    detail::LocalFittingThawHysteresisTracker tracker{ 1, 2.0, 8.0, 0.9 };

    tracker.RecordDependencyThaw(0);
    tracker.RecordDependencyThaw(0);
    tracker.RecordDependencyThaw(0);
    tracker.DecayFrozen(0);

    EXPECT_DOUBLE_EQ(0.72, tracker.GetThreshold(0, 0.1));

    for (int i = 0; i < 100; i++)
    {
        tracker.DecayFrozen(0);
    }

    EXPECT_DOUBLE_EQ(0.1, tracker.GetThreshold(0, 0.1));
}
