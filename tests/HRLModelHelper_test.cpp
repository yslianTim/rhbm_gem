#include <gtest/gtest.h>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wkeyword-macro"
#endif
#define private public
#include "HRLModelHelper.hpp"
#undef private
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

TEST(HRLModelHelperTest, DistancesFarFromThresholdAreFlaggedConsistently)
{
    Eigen::ArrayXd distances(3);
    distances << 1.0, 12.0, 100.0;

    const auto flags{ HRLModelHelper::CalculateOutlierMemberFlag(3, distances) };

    ASSERT_EQ(3, flags.size());
    EXPECT_FALSE(flags(0));
    EXPECT_TRUE(flags(1));
    EXPECT_TRUE(flags(2));
}

TEST(HRLModelHelperTest, BasisTwoBoundaryReflectsBoostOrLegacyThreshold)
{
    Eigen::ArrayXd distances(1);
    distances << 9.2102;

    const auto flags{ HRLModelHelper::CalculateOutlierMemberFlag(2, distances) };

#ifdef HAVE_BOOST
    EXPECT_FALSE(flags(0));
#else
    EXPECT_TRUE(flags(0));
#endif
}
