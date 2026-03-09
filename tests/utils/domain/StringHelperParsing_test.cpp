#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

#include <rhbm_gem/utils/StringHelper.hpp>

TEST(StringHelperParsingTest, ExtractCharWithinBounds)
{
    std::string input{ "ABC" };
    EXPECT_EQ("B", StringHelper::ExtractCharAsString(input, 1));
}

TEST(StringHelperParsingTest, ExtractCharAtZeroReturnsFirstChar)
{
    std::string input{ "ABC" };
    EXPECT_EQ("A", StringHelper::ExtractCharAsString(input, 0));
}

TEST(StringHelperParsingTest, ExtractCharAtLastIndex)
{
    std::string input{ "ABC" };
    EXPECT_EQ("C", StringHelper::ExtractCharAsString(input, input.size() - 1));
}

TEST(StringHelperParsingTest, ExtractCharOutOfBoundsReturnsSpace)
{
    std::string input{ "ABC" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, 5));
}

TEST(StringHelperParsingTest, ExtractCharMaxIndexReturnsSpace)
{
    std::string input{ "ABC" };
    EXPECT_EQ(
        " ",
        StringHelper::ExtractCharAsString(
            input,
            std::numeric_limits<std::size_t>::max()));
}

TEST(StringHelperParsingTest, ExtractCharEmptyInputReturnsSpace)
{
    std::string input{ "" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, 0));
}

TEST(StringHelperParsingTest, ExtractCharAtInputSizeReturnsSpace)
{
    std::string input{ "ABC" };
    EXPECT_EQ(" ", StringHelper::ExtractCharAsString(input, input.size()));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringTrimsSpaces)
{
    const char data[]{ "AB C" };
    EXPECT_EQ("ABC", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringTrimsLeadingAndTrailingSpaces)
{
    const char data[]{ " ABC " };
    EXPECT_EQ("ABC", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringOnlySpacesReturnsEmpty)
{
    const char data[]{ "   " };
    EXPECT_EQ("", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringNullptrReturnsEmpty)
{
    EXPECT_EQ("", StringHelper::ConvertCharArrayToString(nullptr));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringEmptyArrayReturnsEmpty)
{
    const char data[]{ "" };
    EXPECT_EQ("", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringPreservesTabs)
{
    const char data[]{ "A\tB\n" };
    EXPECT_EQ("A\tB\n", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringStopsAtNullChar)
{
    const char data[]{ "AB\0CD" };
    EXPECT_EQ("AB", StringHelper::ConvertCharArrayToString(data));
}

TEST(StringHelperParsingTest, StripCarriageReturnRemovesTrailingCR)
{
    std::string line{ "abc\r" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc", line);
}

TEST(StringHelperParsingTest, StripCarriageReturnNoChangeWithoutCR)
{
    std::string line{ "abc" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc", line);
}

TEST(StringHelperParsingTest, StripCarriageReturnInMiddleNoChange)
{
    std::string line{ "ab\rcd" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("ab\rcd", line);
}

TEST(StringHelperParsingTest, StripCarriageReturnDoesNotRemoveCRBeforeLF)
{
    std::string line{ "abc\r\n" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc\r\n", line);
}

TEST(StringHelperParsingTest, StripCarriageReturnRemovesOnlyOneTrailingCR)
{
    std::string line{ "abc\r\r" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ("abc\r", line);
}

TEST(StringHelperParsingTest, StripCarriageReturnOnlyCRReturnsEmpty)
{
    std::string line{ "\r" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_TRUE(line.empty());
}

TEST(StringHelperParsingTest, StripCarriageReturnKeepsEmptyStringEmpty)
{
    std::string line{ "" };
    StringHelper::StripCarriageReturn(line);
    EXPECT_TRUE(line.empty());
}

TEST(StringHelperParsingTest, SplitCommaSeparatedString)
{
    std::string line{ "a,b,c" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
    EXPECT_EQ("c", tokens[2]);
}

TEST(StringHelperParsingTest, SplitStringWithCustomDelimiter)
{
    std::string line{ "a;b;c" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line, ';') };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
    EXPECT_EQ("c", tokens[2]);
}

TEST(StringHelperParsingTest, SplitStringConsecutiveDelimitersIgnored)
{
    std::string line{ "a,,b,,,c" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
    EXPECT_EQ("c", tokens[2]);
}

TEST(StringHelperParsingTest, SplitStringEmptyInputReturnsEmptyVector)
{
    std::string line{ "" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    EXPECT_TRUE(tokens.empty());
}

TEST(StringHelperParsingTest, SplitStringTrailingDelimiterIgnored)
{
    std::string line{ "a,b," };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(2u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
}

TEST(StringHelperParsingTest, SplitStringLeadingDelimiterIgnored)
{
    std::string line{ ",a,b" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(2u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ("b", tokens[1]);
}

TEST(StringHelperParsingTest, SplitStringPreservesWhitespaceInsideTokens)
{
    std::string line{ "a, b ,c" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ(" b ", tokens[1]);
    EXPECT_EQ("c", tokens[2]);
}

TEST(StringHelperParsingTest, SplitStringWhitespaceOnlyToken)
{
    std::string line{ "a, ,b" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(3u, tokens.size());
    EXPECT_EQ("a", tokens[0]);
    EXPECT_EQ(" ", tokens[1]);
    EXPECT_EQ("b", tokens[2]);
}

TEST(StringHelperParsingTest, SplitStringWithoutDelimiterReturnsWholeString)
{
    std::string line{ "abc" };
    auto tokens{ StringHelper::SplitStringLineFromDelimiter(line) };
    ASSERT_EQ(1u, tokens.size());
    EXPECT_EQ("abc", tokens[0]);
}

TEST(StringHelperParsingTest, ParseListOptionParsesDoubles)
{
    const std::string input{ "1.0,2.5" };
    const std::vector<double> expected{ 1.0, 2.5 };
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input));
}

TEST(StringHelperParsingTest, ParseListOptionParsesStrings)
{
    const std::string input{ "a,b" };
    const std::vector<std::string> expected{ "a", "b" };
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input));
}

TEST(StringHelperParsingTest, ParseListOptionEmptyInputReturnsEmptyVector)
{
    const auto values{ StringHelper::ParseListOption<double>("") };
    EXPECT_TRUE(values.empty());
}

TEST(StringHelperParsingTest, ParseListOptionParsesDoublesWithCustomDelimiter)
{
    const std::string input{ "1.0;2.5" };
    const std::vector<double> expected{ 1.0, 2.5 };
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input, ';'));
}

TEST(StringHelperParsingTest, ParseListOptionParsesStringsWithCustomDelimiter)
{
    const std::string input{ "a;b" };
    const std::vector<std::string> expected{ "a", "b" };
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input, ';'));
}

TEST(StringHelperParsingTest, ParseListOptionParsesDoublesWithWhitespace)
{
    const std::string input{ "1, 2" };
    const std::vector<double> expected{ 1.0, 2.0 };
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input));
}

TEST(StringHelperParsingTest, ParseListOptionConsecutiveDelimitersIgnored)
{
    const std::string input{ "1,,2" };
    const std::vector<double> expected{ 1.0, 2.0 };
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input));
}

TEST(StringHelperParsingTest, ParseListOptionTrailingDelimiterIgnored)
{
    const std::string input{ "1,2," };
    const std::vector<double> expected{ 1.0, 2.0 };
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input));
}

TEST(StringHelperParsingTest, ParseListOptionLeadingDelimiterIgnored)
{
    const std::string input{ ",a,b" };
    const std::vector<std::string> expected{ "a", "b" };
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input));
}

TEST(StringHelperParsingTest, ParseListOptionTrailingDelimiterIgnoredStrings)
{
    const std::string input{ "a,b," };
    const std::vector<std::string> expected{ "a", "b" };
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input));
}

TEST(StringHelperParsingTest, ParseListOptionConsecutiveDelimitersIgnoredForStrings)
{
    const std::string input{ "a,,b" };
    const std::vector<std::string> expected{ "a", "b" };
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input));
}

TEST(StringHelperParsingTest, ParseListOptionStringsWithWhitespace)
{
    const std::string input{ "a, b" };
    const std::vector<std::string> expected{ "a", " b" };
    EXPECT_EQ(expected, StringHelper::ParseListOption<std::string>(input));
}

TEST(StringHelperParsingTest, ParseListOptionWhitespaceTokenThrows)
{
    const std::string input{ "1, ,2" };
    EXPECT_THROW(StringHelper::ParseListOption<double>(input), std::invalid_argument);
}

TEST(StringHelperParsingTest, ParseListOptionInvalidDoubleThrows)
{
    const std::string input{ "1,a" };
    EXPECT_THROW(StringHelper::ParseListOption<double>(input), std::invalid_argument);
}

TEST(StringHelperParsingTest, ParseListOptionPartialDoubleThrows)
{
    const std::string input{ "1.0a,2.0" };
    EXPECT_THROW(StringHelper::ParseListOption<double>(input), std::invalid_argument);
}

TEST(StringHelperParsingTest, ParseListOptionParsesNegativeDoubles)
{
    const std::string input{ "-1.5,-2.5" };
    const std::vector<double> expected{ -1.5, -2.5 };
    EXPECT_EQ(expected, StringHelper::ParseListOption<double>(input));
}

TEST(StringHelperParsingTest, ParseListOptionHandlesVeryLargeDoubles)
{
    const std::string input{ "1e308,2e308" };
    try
    {
        const auto values{ StringHelper::ParseListOption<double>(input) };
        ASSERT_EQ(2u, values.size());
        EXPECT_DOUBLE_EQ(1e308, values[0]);
        EXPECT_TRUE(std::isinf(values[1]));
    }
    catch (const std::out_of_range &)
    {
        SUCCEED();
    }
}
