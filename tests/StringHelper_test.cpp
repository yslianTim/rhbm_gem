#include <gtest/gtest.h>
#include <limits>

#include "StringHelper.hpp"

TEST(StringHelperTest, ExtractCharWithinBounds)
{
    std::string input{ "ABC" };
    EXPECT_EQ("B", StringHelper::ExtractCharAsString(input, 1));
}

TEST(StringHelperTest, ExtractCharAtZeroReturnsFirstChar)
{
    std::string input{ "ABC" };
    EXPECT_EQ("A", StringHelper::ExtractCharAsString(input, 0));
}

TEST(StringHelperTest, ExtractCharOutOfBoundsReturnsSpace)
{
    std::string input{ "ABC" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, 5));
}

TEST(StringHelperTest, ExtractCharMaxIndexReturnsSpace)
{
    std::string input{ "ABC" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, std::numeric_limits<std::size_t>::max()));
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

TEST(StringHelperTest, ConvertCharArrayToStringTrimsLeadingAndTrailingSpaces)
{
    const char data[]{ " ABC " };
    EXPECT_EQ("ABC", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperTest, ConvertCharArrayToStringOnlySpacesReturnsEmpty)
{
    const char data[]{ "   " };
    EXPECT_EQ("", StringHelper::ConvertCharArrayToString(data));
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

TEST(StringHelperTest, StripCarriageReturnNoChangeWithoutCR)
{
    std::string line{ "abc" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc", line);
}

TEST(StringHelperTest, StripCarriageReturnInMiddleNoChange)
{
    std::string line{ "ab\rcd" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("ab\rcd", line);
}

TEST(StringHelperTest, StripCarriageReturnRemovesOnlyOneTrailingCR)
{
    std::string line{ "abc\r\r" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc\r", line);
}

TEST(StringHelperTest, StripCarriageReturnKeepsEmptyStringEmpty)
{
    std::string line{ "" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_TRUE(line.empty());
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

TEST(StringHelperTest, SplitStringEmptyInputReturnsEmptyVector)
{
    std::string line{ "" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    EXPECT_TRUE(tokens.empty());
}

TEST(StringHelperTest, SplitStringTrailingDelimiterIgnored)
{
    std::string line{ "a,b," };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(2u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
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

TEST(StringHelperTest, ParseListOptionParsesDoublesWithCustomDelimiter)
{
    const std::string input{ "1.0;2.5" };
    const std::vector<double> expected{1.0, 2.5};
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input, ';'));
}

TEST(StringHelperTest, ParseListOptionParsesStringsWithCustomDelimiter)
{
    const std::string input{ "a;b" };
    const std::vector<std::string> expected{"a", "b"};
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input, ';'));
}

TEST(StringHelperTest, ParseListOptionParsesDoublesWithWhitespace)
{
    const std::string input{ "1, 2" };
    const std::vector<double> expected{1.0, 2.0};
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input));
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

TEST(StringHelperTest, SplitStringLineAsTokensHandlesMultipleSpaces)
{
    std::string line{ "  one   two   three  " };
    std::vector<std::string> expected{"one", "two", "three"};
    EXPECT_EQ(expected, StringHelper::SplitStringLineAsTokens(line, 3));
}

TEST(StringHelperTest, SplitStringLineAsTokensUnclosedQuoteReturnsRemaining)
{
    std::string line{ "one 'two three" };
    auto tokens{ StringHelper::SplitStringLineAsTokens(line, 3) };
    ASSERT_EQ(2u, tokens.size());
    EXPECT_EQ("one", tokens[0]);
    EXPECT_EQ("two three", tokens[1]);
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

TEST(StringHelperTest, ToStringWithPrecisionIntegerUsesDefaultPrecision)
{
    int value{ 42 };
    EXPECT_EQ("42.000000", StringHelper::ToStringWithPrecision<double>(value));
}

TEST(StringHelperTest, ToStringWithPrecisionNegativeUsesCustomPrecision)
{
    double value{ -123.456789 };
    EXPECT_EQ("-123.46", StringHelper::ToStringWithPrecision(value, 2));
}
