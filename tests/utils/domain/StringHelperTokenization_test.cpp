#include <gtest/gtest.h>

#include <rhbm_gem/utils/domain/StringHelper.hpp>

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensHandlesQuotedStrings)
{
    std::string line{ "one 'two words' three" };
    std::vector<std::string> expected{ "one", "two words", "three" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3));
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensHandlesQuotedCommas)
{
    std::string line{ "one 'two,three' four" };
    std::vector<std::string> expected{ "one", "two,three", "four" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3));
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensUsesCustomDelimiter)
{
    std::string line{ "one \"two words\" three" };
    std::vector<std::string> expected{ "one", "two words", "three" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3, '"'));
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensHandlesMultipleSpaces)
{
    std::string line{ "  one   two   three  " };
    std::vector<std::string> expected{ "one", "two", "three" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3));
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensUnclosedQuoteReturnsRemaining)
{
    std::string line{ "one 'two three" };
    auto tokens{ rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3) };
    ASSERT_EQ(2u, tokens.size());
    EXPECT_EQ("one", tokens[0]);
    EXPECT_EQ("two three", tokens[1]);
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensEmptyInputReturnsEmptyVector)
{
    std::string line{ "" };
    auto tokens{ rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3) };
    EXPECT_TRUE(tokens.empty());
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensHandlesEmptyQuotedToken)
{
    std::string line{ "one '' three" };
    std::vector<std::string> expected{ "one", "", "three" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3));
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensIgnoresZeroTokenCount)
{
    std::string line{ "one two" };
    std::vector<std::string> expected{ "one", "two" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 0));
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensIgnoresSmallerTokenCount)
{
    std::string line{ "one two three" };
    std::vector<std::string> expected{ "one", "two", "three" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 2));
}

TEST(StringHelperTokenizationTest, SplitStringLineAsTokensHandlesTabsAndNewlines)
{
    std::string line{ "one\ttwo\nthree" };
    std::vector<std::string> expected{ "one", "two", "three" };
    EXPECT_EQ(expected, rhbm_gem::string_helper::SplitStringLineAsTokens(line, 3));
}
