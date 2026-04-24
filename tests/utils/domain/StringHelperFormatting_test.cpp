#include <gtest/gtest.h>

#include <limits>

#include <rhbm_gem/utils/domain/StringHelper.hpp>

TEST(StringHelperFormattingTest, ToStringWithPrecisionUsesDefaultPrecision)
{
    double value{ 123.456789 };
    EXPECT_EQ("123.456789", rhbm_gem::string_helper::ToStringWithPrecision(value));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionUsesCustomPrecision)
{
    double value{ 123.456789 };
    EXPECT_EQ("123.46", rhbm_gem::string_helper::ToStringWithPrecision(value, 2));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionIntegerUsesDefaultPrecision)
{
    int value{ 42 };
    EXPECT_EQ("42.000000", rhbm_gem::string_helper::ToStringWithPrecision<double>(value));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionNegativeUsesCustomPrecision)
{
    double value{ -123.456789 };
    EXPECT_EQ("-123.46", rhbm_gem::string_helper::ToStringWithPrecision(value, 2));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionZeroPrecision)
{
    double value{ 123.456 };
    EXPECT_EQ("123", rhbm_gem::string_helper::ToStringWithPrecision(value, 0));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionHandlesInfinity)
{
    const double pos_inf{ std::numeric_limits<double>::infinity() };
    const double neg_inf{ -std::numeric_limits<double>::infinity() };
    EXPECT_EQ("inf", rhbm_gem::string_helper::ToStringWithPrecision(pos_inf));
    EXPECT_EQ("-inf", rhbm_gem::string_helper::ToStringWithPrecision(neg_inf));
}

TEST(StringHelperFormattingTest, ToStringWithPrecisionHandlesNaN)
{
    const double nan_val{ std::numeric_limits<double>::quiet_NaN() };
    EXPECT_EQ("nan", rhbm_gem::string_helper::ToStringWithPrecision(nan_val));
}

TEST(StringHelperFormattingTest, PadWithSpacesPadsShortString)
{
    std::string input{ "abc" };
    EXPECT_EQ("abc   ", rhbm_gem::string_helper::PadWithSpaces(input, 6));
}

TEST(StringHelperFormattingTest, PadWithSpacesReturnsOriginalWhenLongEnough)
{
    std::string input{ "abcdef" };
    EXPECT_EQ("abcdef", rhbm_gem::string_helper::PadWithSpaces(input, 3));
}
