#include <gtest/gtest.h>

#include "StringHelper.hpp"

TEST(StringHelperTest, ExtractCharWithinBounds)
{
    std::string input = "ABC";
    EXPECT_EQ("B", StringHelper::ExtractCharAsString(input, 1));
}

TEST(StringHelperTest, ExtractCharOutOfBoundsReturnsSpace)
{
    std::string input = "ABC";
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, 5));
}

TEST(StringHelperTest, ConvertCharArrayToStringTrimsSpaces)
{
    const char data[] = "AB C";
    EXPECT_EQ("ABC", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperTest, StripCarriageReturnRemovesTrailingCR)
{
    std::string line = "abc\r";
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc", line);
}
