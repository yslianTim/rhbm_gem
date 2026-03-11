#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <rhbm_gem/utils/domain/StringHelper.hpp>

namespace
{

struct ExtractCharCase
{
    const char * name;
    std::string input;
    size_t index;
    std::string expected;
};

class ExtractCharParsingTest : public ::testing::TestWithParam<ExtractCharCase>
{
};

TEST_P(ExtractCharParsingTest, ReturnsExpectedCharacterOrSpace)
{
    const auto & test_case{ GetParam() };
    EXPECT_EQ(
        test_case.expected,
        StringHelper::ExtractCharAsString(test_case.input, test_case.index));
}

INSTANTIATE_TEST_SUITE_P(
    ExtractCharCases,
    ExtractCharParsingTest,
    ::testing::Values(
        ExtractCharCase{ "WithinBounds", "ABC", 1, "B" },
        ExtractCharCase{ "AtZero", "ABC", 0, "A" },
        ExtractCharCase{ "AtLastIndex", "ABC", 2, "C" },
        ExtractCharCase{ "AtInputSize", "ABC", 3, " " },
        ExtractCharCase{ "OutOfBounds", "ABC", 5, " " },
        ExtractCharCase{ "MaxIndex", "ABC", std::numeric_limits<size_t>::max(), " " },
        ExtractCharCase{ "EmptyInput", "", 0, " " }),
    [](const ::testing::TestParamInfo<ExtractCharCase> & info)
    {
        return info.param.name;
    });

struct ConvertCharArrayCase
{
    const char * name;
    const char * input;
    std::string expected;
};

class ConvertCharArrayParsingTest : public ::testing::TestWithParam<ConvertCharArrayCase>
{
};

TEST_P(ConvertCharArrayParsingTest, ConvertsAsExpected)
{
    const auto & test_case{ GetParam() };
    EXPECT_EQ(test_case.expected, StringHelper::ConvertCharArrayToString(test_case.input));
}

INSTANTIATE_TEST_SUITE_P(
    ConvertCharArrayCases,
    ConvertCharArrayParsingTest,
    ::testing::Values(
        ConvertCharArrayCase{ "TrimsSpaces", "AB C", "ABC" },
        ConvertCharArrayCase{ "TrimsLeadingAndTrailingSpaces", " ABC ", "ABC" },
        ConvertCharArrayCase{ "OnlySpacesReturnsEmpty", "   ", "" },
        ConvertCharArrayCase{ "EmptyArrayReturnsEmpty", "", "" },
        ConvertCharArrayCase{ "PreservesTabsAndNewlines", "A\tB\n", "A\tB\n" }),
    [](const ::testing::TestParamInfo<ConvertCharArrayCase> & info)
    {
        return info.param.name;
    });

TEST(StringHelperParsingTest, ConvertCharArrayToStringNullptrReturnsEmpty)
{
    EXPECT_EQ("", StringHelper::ConvertCharArrayToString(nullptr));
}

TEST(StringHelperParsingTest, ConvertCharArrayToStringStopsAtNullChar)
{
    const char data[]{ "AB\0CD" };
    EXPECT_EQ("AB", StringHelper::ConvertCharArrayToString(data));
}

struct StripCarriageReturnCase
{
    const char * name;
    std::string input;
    std::string expected;
};

class StripCarriageReturnParsingTest
    : public ::testing::TestWithParam<StripCarriageReturnCase>
{
};

TEST_P(StripCarriageReturnParsingTest, ProducesExpectedString)
{
    auto line{ GetParam().input };
    StringHelper::StripCarriageReturn(line);
    EXPECT_EQ(GetParam().expected, line);
}

INSTANTIATE_TEST_SUITE_P(
    StripCarriageReturnCases,
    StripCarriageReturnParsingTest,
    ::testing::Values(
        StripCarriageReturnCase{ "RemovesTrailingCR", "abc\r", "abc" },
        StripCarriageReturnCase{ "NoChangeWithoutCR", "abc", "abc" },
        StripCarriageReturnCase{ "MiddleCRUnchanged", "ab\rcd", "ab\rcd" },
        StripCarriageReturnCase{ "CRBeforeLFUnchanged", "abc\r\n", "abc\r\n" },
        StripCarriageReturnCase{ "RemovesOnlyOneTrailingCR", "abc\r\r", "abc\r" },
        StripCarriageReturnCase{ "OnlyCRReturnsEmpty", "\r", "" },
        StripCarriageReturnCase{ "KeepsEmptyStringEmpty", "", "" }),
    [](const ::testing::TestParamInfo<StripCarriageReturnCase> & info)
    {
        return info.param.name;
    });

struct SplitCase
{
    const char * name;
    std::string input;
    char delimiter;
    std::vector<std::string> expected;
};

class SplitStringParsingTest : public ::testing::TestWithParam<SplitCase>
{
};

TEST_P(SplitStringParsingTest, SplitsAsExpected)
{
    const auto & test_case{ GetParam() };
    EXPECT_EQ(
        test_case.expected,
        StringHelper::SplitStringLineFromDelimiter(test_case.input, test_case.delimiter));
}

INSTANTIATE_TEST_SUITE_P(
    SplitStringCases,
    SplitStringParsingTest,
    ::testing::Values(
        SplitCase{ "CommaSeparated", "a,b,c", ',', { "a", "b", "c" } },
        SplitCase{ "CustomDelimiter", "a;b;c", ';', { "a", "b", "c" } },
        SplitCase{ "ConsecutiveDelimitersIgnored", "a,,b,,,c", ',', { "a", "b", "c" } },
        SplitCase{ "EmptyInputReturnsEmpty", "", ',', {} },
        SplitCase{ "TrailingDelimiterIgnored", "a,b,", ',', { "a", "b" } },
        SplitCase{ "LeadingDelimiterIgnored", ",a,b", ',', { "a", "b" } },
        SplitCase{ "PreservesWhitespaceInsideTokens", "a, b ,c", ',', { "a", " b ", "c" } },
        SplitCase{ "WhitespaceOnlyTokenPreserved", "a, ,b", ',', { "a", " ", "b" } },
        SplitCase{ "WithoutDelimiterReturnsWholeString", "abc", ',', { "abc" } }),
    [](const ::testing::TestParamInfo<SplitCase> & info)
    {
        return info.param.name;
    });

struct ParseDoubleCase
{
    const char * name;
    std::string input;
    char delimiter;
    std::vector<double> expected;
};

class ParseDoubleListOptionTest : public ::testing::TestWithParam<ParseDoubleCase>
{
};

TEST_P(ParseDoubleListOptionTest, ParsesExpectedValues)
{
    const auto & test_case{ GetParam() };
    EXPECT_EQ(
        test_case.expected,
        StringHelper::ParseListOption<double>(test_case.input, test_case.delimiter));
}

INSTANTIATE_TEST_SUITE_P(
    ParseDoubleCases,
    ParseDoubleListOptionTest,
    ::testing::Values(
        ParseDoubleCase{ "Basic", "1.0,2.5", ',', { 1.0, 2.5 } },
        ParseDoubleCase{ "CustomDelimiter", "1.0;2.5", ';', { 1.0, 2.5 } },
        ParseDoubleCase{ "WhitespaceIsAccepted", "1, 2", ',', { 1.0, 2.0 } },
        ParseDoubleCase{ "NegativeValues", "-1.5,-2.5", ',', { -1.5, -2.5 } }),
    [](const ::testing::TestParamInfo<ParseDoubleCase> & info)
    {
        return info.param.name;
    });

struct ParseStringCase
{
    const char * name;
    std::string input;
    char delimiter;
    std::vector<std::string> expected;
};

class ParseStringListOptionTest : public ::testing::TestWithParam<ParseStringCase>
{
};

TEST_P(ParseStringListOptionTest, ParsesExpectedValues)
{
    const auto & test_case{ GetParam() };
    EXPECT_EQ(
        test_case.expected,
        StringHelper::ParseListOption<std::string>(test_case.input, test_case.delimiter));
}

INSTANTIATE_TEST_SUITE_P(
    ParseStringCases,
    ParseStringListOptionTest,
    ::testing::Values(
        ParseStringCase{ "Basic", "a,b", ',', { "a", "b" } },
        ParseStringCase{ "CustomDelimiter", "a;b", ';', { "a", "b" } },
        ParseStringCase{ "WhitespaceIsPreserved", "a, b", ',', { "a", " b" } }),
    [](const ::testing::TestParamInfo<ParseStringCase> & info)
    {
        return info.param.name;
    });

TEST(StringHelperParsingTest, ParseListOptionEmptyInputReturnsEmptyVector)
{
    const auto values{ StringHelper::ParseListOption<double>("") };
    EXPECT_TRUE(values.empty());
}

struct InvalidDoubleListCase
{
    const char * name;
    std::string input;
};

class ParseInvalidDoubleListOptionTest
    : public ::testing::TestWithParam<InvalidDoubleListCase>
{
};

TEST_P(ParseInvalidDoubleListOptionTest, ThrowsInvalidArgument)
{
    EXPECT_THROW(
        StringHelper::ParseListOption<double>(GetParam().input),
        std::invalid_argument);
}

INSTANTIATE_TEST_SUITE_P(
    InvalidDoubleCases,
    ParseInvalidDoubleListOptionTest,
    ::testing::Values(
        InvalidDoubleListCase{ "WhitespaceToken", "1, ,2" },
        InvalidDoubleListCase{ "InvalidDoubleToken", "1,a" },
        InvalidDoubleListCase{ "PartialDoubleToken", "1.0a,2.0" }),
    [](const ::testing::TestParamInfo<InvalidDoubleListCase> & info)
    {
        return info.param.name;
    });

TEST(StringHelperParsingTest, ParseListOptionOverflowThrowsOutOfRange)
{
    EXPECT_THROW(
        StringHelper::ParseListOption<double>("1e308,2e308"),
        std::out_of_range);
}

} // namespace
