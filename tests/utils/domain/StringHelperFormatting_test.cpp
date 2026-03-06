#include <gtest/gtest.h>

#include <limits>

#include "StringHelper.hpp"

TEST(StringHelperFormattingTest, ToStringWithPrecisionUsesDefaultPrecision)
{
    double value{ 123.456789 };
    EXPECT_EQ("123.456789", StringHelper::ToStringWithPrecision(value));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionUsesCustomPrecision)
{
    double value{ 123.456789 };
    EXPECT_EQ("123.46", StringHelper::ToStringWithPrecision(value, 2));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionIntegerUsesDefaultPrecision)
{
    int value{ 42 };
    EXPECT_EQ("42.000000", StringHelper::ToStringWithPrecision<double>(value));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionNegativeUsesCustomPrecision)
{
    double value{ -123.456789 };
    EXPECT_EQ("-123.46", StringHelper::ToStringWithPrecision(value, 2));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionZeroPrecision)
{
    double value{ 123.456 };
    EXPECT_EQ("123", StringHelper::ToStringWithPrecision(value, 0));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionHandlesInfinity)
{
    const double pos_inf{ std::numeric_limits<double>::infinity() };
    const double neg_inf{ -std::numeric_limits<double>::infinity() };
    EXPECT_EQ("inf", StringHelper::ToStringWithPrecision(pos_inf));
    EXPECT_EQ("-inf", StringHelper::ToStringWithPrecision(neg_inf));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionHandlesNaN)
{
    const double nan_val{ std::numeric_limits<double>::quiet_NaN() };
    EXPECT_EQ("nan", StringHelper::ToStringWithPrecision(nan_val));
}

TEST(StringHelperFormattingTest, PadWithSpacesPadsShortString)
{
    std::string input{ "abc" };
    EXPECT_EQ("abc   ", StringHelper::PadWithSpaces(input, 6));
}

TEST(StringHelperFormattingTest, PadWithSpacesReturnsOriginalWhenLongEnough)
{
    std::string input{ "abcdef" };
    EXPECT_EQ("abcdef", StringHelper::PadWithSpaces(input, 3));
}
