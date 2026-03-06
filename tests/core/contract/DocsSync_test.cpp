#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include "BuiltInCommandCatalogInternal.hpp"
#include "CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

std::string ReadFileContent(const std::filesystem::path & path)
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        return {};
    }
    return std::string(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
}

std::string ExtractGeneratedBlock(
    const std::string & content,
    std::string_view begin_marker,
    std::string_view end_marker)
{
    const auto begin_pos{ content.find(begin_marker) };
    if (begin_pos == std::string::npos)
    {
        return {};
    }

    auto block_begin{ content.find('\n', begin_pos) };
    if (block_begin == std::string::npos)
    {
        return {};
    }
    ++block_begin;

    const auto end_pos{ content.find(end_marker, block_begin) };
    if (end_pos == std::string::npos)
    {
        return {};
    }

    return content.substr(block_begin, end_pos - block_begin);
}

std::string NormalizeBlock(std::string text)
{
    text.erase(std::remove(text.begin(), text.end(), '\r'), text.end());
    const auto is_space{
        [](unsigned char c) { return std::isspace(c) != 0; }
    };

    while (!text.empty() && is_space(static_cast<unsigned char>(text.front())))
    {
        text.erase(text.begin());
    }
    while (!text.empty() && is_space(static_cast<unsigned char>(text.back())))
    {
        text.pop_back();
    }

    return text;
}

std::string BuildExpectedManifestBlock()
{
    std::ostringstream output;
    const auto & catalog{ rg::BuiltInCommandCatalog() };
    for (std::size_t index = 0; index < catalog.size(); ++index)
    {
        output << (index + 1) << ". `" << catalog[index].name << "`\n";
    }
    return output.str();
}

std::string BuildExpectedSurfaceMatrixBlock()
{
    std::ostringstream output;
    output << "| Command | Uses database at runtime | Uses output folder |\n";
    output << "| --- | --- | --- |\n";
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        output << "| `" << descriptor.name << "` | "
               << (rg::UsesDatabaseAtRuntime(descriptor.common_options) ? "yes" : "no")
               << " | "
               << (rg::UsesOutputFolder(descriptor.common_options) ? "yes" : "no")
               << " |\n";
    }
    return output.str();
}

std::string BuildExpectedPythonSurfaceBlock()
{
    std::ostringstream output;
    output << "### Built-in Python command classes\n";
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        output << "- `" << descriptor.python_binding_name << "`\n";
    }

    output << "\n### Shared diagnostics types\n";
    output << "- `LogLevel`\n";
    output << "- `ValidationPhase`\n";
    output << "- `ValidationIssue`\n";

    output << "\n### Shared diagnostics methods on built-in Python commands\n";
    output << "- `PrepareForExecution()`\n";
    output << "- `HasValidationErrors()`\n";
    output << "- `GetValidationIssues()`\n";

    return output.str();
}

} // namespace

TEST(DocsSyncTest, CommandArchitectureGeneratedSectionsMatchBuiltInCatalog)
{
    const auto project_root{ command_test::ProjectRootPath() };
    const auto document_path{
        project_root / "docs" / "developer" / "architecture" / "command-architecture.md"
    };
    const auto document{ ReadFileContent(document_path) };
    ASSERT_FALSE(document.empty()) << document_path.string();

    const auto actual_manifest{
        ExtractGeneratedBlock(
            document,
            "<!-- BEGIN GENERATED: built-in-command-manifest -->",
            "<!-- END GENERATED: built-in-command-manifest -->")
    };
    const auto actual_matrix{
        ExtractGeneratedBlock(
            document,
            "<!-- BEGIN GENERATED: command-surface-matrix -->",
            "<!-- END GENERATED: command-surface-matrix -->")
    };
    const auto actual_python_surface{
        ExtractGeneratedBlock(
            document,
            "<!-- BEGIN GENERATED: built-in-python-command-surface -->",
            "<!-- END GENERATED: built-in-python-command-surface -->")
    };

    ASSERT_FALSE(actual_manifest.empty());
    ASSERT_FALSE(actual_matrix.empty());
    ASSERT_FALSE(actual_python_surface.empty());

    EXPECT_EQ(
        NormalizeBlock(actual_manifest),
        NormalizeBlock(BuildExpectedManifestBlock()));
    EXPECT_EQ(
        NormalizeBlock(actual_matrix),
        NormalizeBlock(BuildExpectedSurfaceMatrixBlock()));
    EXPECT_EQ(
        NormalizeBlock(actual_python_surface),
        NormalizeBlock(BuildExpectedPythonSurfaceBlock()));
}
