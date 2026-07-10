#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

#include <rhbm_gem/utils/algorithm/ScaleReferenceTracker.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;

} // namespace

TEST(ScaleReferenceTrackerTest, EmptyTrackerHasNoReference)
{
    alg::ScaleReferenceTracker tracker{ 3 };

    EXPECT_FALSE(tracker.GetCommittedReference().has_value());
    EXPECT_FALSE(tracker.IsLocked());
}

TEST(ScaleReferenceTrackerTest, InitialSampleSeedsReferenceWithoutLocking)
{
    alg::ScaleReferenceTracker tracker{ 2, 4.0 };

    const auto reference{ tracker.GetCommittedReference() };

    ASSERT_TRUE(reference.has_value());
    EXPECT_DOUBLE_EQ(4.0, *reference);
    EXPECT_FALSE(tracker.IsLocked());
}

TEST(ScaleReferenceTrackerTest, ProvisionalReferenceIncludesCandidateWithoutCommitting)
{
    alg::ScaleReferenceTracker tracker{ 3, 2.0 };

    const auto provisional_reference{ tracker.GetProvisionalReference(4.0) };
    const auto committed_reference{ tracker.GetCommittedReference() };

    ASSERT_TRUE(provisional_reference.has_value());
    EXPECT_DOUBLE_EQ(3.0, *provisional_reference);
    ASSERT_TRUE(committed_reference.has_value());
    EXPECT_DOUBLE_EQ(2.0, *committed_reference);
}

TEST(ScaleReferenceTrackerTest, CommittedReferenceUsesRecentSamples)
{
    alg::ScaleReferenceTracker tracker{ 2, 1.0 };

    tracker.CommitScaleSample(3.0);
    tracker.CommitScaleSample(5.0);
    const auto reference{ tracker.GetCommittedReference() };

    ASSERT_TRUE(reference.has_value());
    EXPECT_DOUBLE_EQ(4.0, *reference);
}

TEST(ScaleReferenceTrackerTest, LocksAfterFiniteCommittedWarmupSamples)
{
    alg::ScaleReferenceTracker tracker{ 2, 1.0 };

    tracker.CommitScaleSample(3.0);
    EXPECT_FALSE(tracker.IsLocked());

    tracker.CommitScaleSample(5.0);
    EXPECT_TRUE(tracker.IsLocked());
}

TEST(ScaleReferenceTrackerTest, LockedTrackerIgnoresProvisionalAndCommittedSamples)
{
    alg::ScaleReferenceTracker tracker{ 1, 2.0 };

    tracker.CommitScaleSample(4.0);
    ASSERT_TRUE(tracker.IsLocked());
    const auto committed_reference{ tracker.GetCommittedReference() };
    const auto provisional_reference{ tracker.GetProvisionalReference(100.0) };
    tracker.CommitScaleSample(100.0);
    const auto updated_reference{ tracker.GetCommittedReference() };

    ASSERT_TRUE(committed_reference.has_value());
    ASSERT_TRUE(provisional_reference.has_value());
    ASSERT_TRUE(updated_reference.has_value());
    EXPECT_DOUBLE_EQ(*committed_reference, *provisional_reference);
    EXPECT_DOUBLE_EQ(*committed_reference, *updated_reference);
}

TEST(ScaleReferenceTrackerTest, NonFiniteSamplesDoNotCreateReferenceOrAdvanceWarmup)
{
    const auto infinity{ std::numeric_limits<double>::infinity() };
    const auto nan{ std::numeric_limits<double>::quiet_NaN() };
    alg::ScaleReferenceTracker tracker{ 1 };

    tracker.CommitScaleSample(infinity);
    EXPECT_FALSE(tracker.IsLocked());
    EXPECT_FALSE(tracker.GetProvisionalReference(nan).has_value());

    tracker.CommitScaleSample(2.0);

    EXPECT_TRUE(tracker.IsLocked());
    ASSERT_TRUE(tracker.GetCommittedReference().has_value());
    EXPECT_DOUBLE_EQ(2.0, *tracker.GetCommittedReference());
}

TEST(ScaleReferenceTrackerTest, RejectsZeroWarmupSampleCount)
{
    EXPECT_THROW(alg::ScaleReferenceTracker(0), std::invalid_argument);
}
