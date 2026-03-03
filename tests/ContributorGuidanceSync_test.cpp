#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

TEST(ContributorGuidanceSyncTest, DevelopmentGuidelinesMatchCatalogBasedRegistrationModel)
{
    const auto doc_path{
        std::filesystem::path(__FILE__).parent_path().parent_path()
        / "docs" / "developer" / "development-guidelines.md"
    };

    std::ifstream doc_stream(doc_path);
    ASSERT_TRUE(doc_stream.is_open());
    const std::string doc_content{
        std::istreambuf_iterator<char>(doc_stream),
        std::istreambuf_iterator<char>()
    };

    EXPECT_NE(
        doc_content.find("New built-in commands should be added through `BuiltInCommandCatalog()`"),
        std::string::npos);
    EXPECT_NE(
        doc_content.find(
            "Reserve `CommandRegistrar<...>` for extension-style or test-only registration paths"),
        std::string::npos);
    EXPECT_EQ(
        doc_content.find(
            "Continue using the existing `CommandRegistrar<...>` registration model for new commands"),
        std::string::npos);
}
