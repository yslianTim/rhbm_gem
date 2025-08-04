#include <gtest/gtest.h>

#include "StringHelper.hpp"

TEST(StringHelperTest, ExtractCharWithinBounds)
{
    std::string input{ "ABC" };
    EXPECT_EQ("B", StringHelper::ExtractCharAsString(input, 1));
}

TEST(StringHelperTest, ExtractCharOutOfBoundsReturnsSpace)
{
    std::string input{ "ABC" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, 5));
}

TEST(StringHelperTest, ExtractCharEmptyInputReturnsSpace)
{
    std::string input{ "" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, 0));
}

TEST(StringHelperTest, ExtractCharAtInputSizeReturnsSpace)
{
    std::string input{ "ABC" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, input.size()));
}

TEST(StringHelperTest, ConvertCharArrayToStringTrimsSpaces)
{
    const char data[]{ "AB C" };
    EXPECT_EQ("ABC", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperTest, ConvertCharArrayToStringNullptrReturnsEmpty)
{
    EXPECT_EQ("", StringHelper::ConvertCharArrayToString(nullptr));
}

TEST(StringHelperTest, ConvertCharArrayToStringEmptyArrayReturnsEmpty)
{
    const char data[]{ "" };
    EXPECT_EQ("", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperTest, StripCarriageReturnRemovesTrailingCR)
{
    std::string line{ "abc\r" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc", line);
}

TEST(StringHelperTest, ToStringWithPrecisionUsesDefaultPrecision)
{
    double value{ 123.456789 };
    EXPECT_EQ("123.456789", StringHelper::ToStringWithPrecision(value));
}

TEST(StringHelperTest, ToStringWithPrecisionUsesCustomPrecision)
{
    double value{ 123.456789 };
    EXPECT_EQ("123.46", StringHelper::ToStringWithPrecision(value, 2));
}

TEST(StringHelperTest, SplitStringLineAsTokensHandlesQuotedStrings)
{
    std::string line{ "one 'two words' three" };
    std::vector<std::string> expected{"one", "two words", "three"};
    EXPECT_EQ(expected, StringHelper::SplitStringLineAsTokens(line, 3));
}

TEST(StringHelperTest, SplitStringLineAsTokensUsesCustomDelimiter)
{
    std::string line{ "one \"two words\" three" };
    std::vector<std::string> expected{"one", "two words", "three"};
    EXPECT_EQ(expected, StringHelper::SplitStringLineAsTokens(line, 3, '"'));
}

TEST(StringHelperTest, ParseListOptionParsesDoubles)
{
    const std::string input{ "1.0,2.5" };
    const std::vector<double> expected{1.0, 2.5};
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input));
}

TEST(StringHelperTest, ParseListOptionParsesStrings)
{
    const std::string input{ "a,b" };
    const std::vector<std::string> expected{"a", "b"};
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input));
}

TEST(StringHelperTest, ParseListOptionEmptyInputReturnsEmptyVector)
{
    const auto values{ StringHelper::ParseListOption<double>("") };
    EXPECT_TRUE(values.empty());
}

TEST(StringHelperTest, SplitCommaSeparatedString)
{
    std::string line{ "a,b,c" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
    EXPECT_EQ("c", tokens[2]);
}

TEST(StringHelperTest, SplitStringWithCustomDelimiter)
{
    std::string line{ "a;b;c" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line, ';') };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
    EXPECT_EQ("c", tokens[2]);
}

TEST(StringHelperTest, SplitStringConsecutiveDelimitersIgnored)
{
    std::string line{ "a,,b,,,c" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
    EXPECT_EQ("c", tokens[2]);
}

TEST(StringHelperTest, StripCarriageReturnNoChangeWithoutCR)
{
    std::string line{ "abc" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc", line);
}

TEST(StringHelperTest, StripCarriageReturnKeepsEmptyStringEmpty)
{
    std::string line{ "" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_TRUE(line.empty());
}
