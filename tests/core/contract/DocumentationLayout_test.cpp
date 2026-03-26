#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <string>
#include <vector>

#include "support/CommandTestHelpers.hpp"

namespace fs = std::filesystem;

namespace {

std::string ReadFileContent(const fs::path & path)
{
    std::ifstream input{ path };
    if (!input.is_open())
    {
        return {};
    }
    return std::string(
        std::istreambuf_iterator<char>{ input },
        std::istreambuf_iterator<char>{});
}

std::vector<std::string> CollectProjectRelativeMarkdownLinks(const std::string & content)
{
    const std::regex link_pattern{ R"(\[[^\]]+\]\(([^)]+)\))" };

    std::vector<std::string> links;
    for (std::sregex_iterator iter{ content.begin(), content.end(), link_pattern };
        iter != std::sregex_iterator{};
        ++iter)
    {
        std::string target{ (*iter)[1].str() };
        const auto fragment_pos{ target.find('#') };
        if (fragment_pos != std::string::npos)
        {
            target.erase(fragment_pos);
        }

        if (target.empty()
            || target.rfind("http://", 0) == 0
            || target.rfind("https://", 0) == 0
            || target.rfind("mailto:", 0) == 0
            || target.front() == '#')
        {
            continue;
        }

        links.push_back(target);
    }

    return links;
}

std::vector<fs::path> DocumentationFilesToValidate()
{
    const auto project_root{ command_test::ProjectRootPath() };

    std::vector<fs::path> files{
        project_root / "README.md",
        project_root / "resources" / "README.md",
    };

    for (const auto & entry : fs::recursive_directory_iterator(project_root / "docs"))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".md")
        {
            files.push_back(entry.path());
        }
    }

    return files;
}

} // namespace

TEST(DocumentationLayoutTest, MarkdownLinksResolveWithinRepository)
{
    for (const auto & markdown_path : DocumentationFilesToValidate())
    {
        SCOPED_TRACE(markdown_path.string());

        const auto content{ ReadFileContent(markdown_path) };
        ASSERT_FALSE(content.empty()) << markdown_path;

        for (const auto & relative_target : CollectProjectRelativeMarkdownLinks(content))
        {
            const auto resolved_target{
                (markdown_path.parent_path() / relative_target).lexically_normal()
            };
            EXPECT_TRUE(fs::exists(resolved_target))
                << "Broken markdown link target: " << relative_target
                << " from " << markdown_path.string();
        }
    }
}

TEST(DocumentationLayoutTest, CommandDocumentationIndexExistsForScaffoldOutput)
{
    const auto project_root{ command_test::ProjectRootPath() };
    const auto command_docs_dir{ project_root / "docs" / "developer" / "commands" };
    const auto command_docs_index{ command_docs_dir / "README.md" };

    EXPECT_TRUE(fs::is_directory(command_docs_dir));
    EXPECT_TRUE(fs::exists(command_docs_index));
}

TEST(DocumentationLayoutTest, ResourcesTreeKeepsNarrativeDocumentationInDocsDirectory)
{
    const auto project_root{ command_test::ProjectRootPath() };
    EXPECT_FALSE(fs::exists(project_root / "resources" / "docs"));
}
