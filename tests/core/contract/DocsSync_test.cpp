#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "internal/CommandCatalog.hpp"
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
    const auto & catalog{ rg::CommandCatalog() };
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
    for (const auto & descriptor : rg::CommandCatalog())
    {
        const auto common_options{ rg::CommonOptionsForCommand(descriptor) };
        output << "| `" << descriptor.name << "` | "
               << (rg::HasCommonOption(
                    common_options,
                    rg::CommonOption::Database) ? "yes" : "no")
               << " | "
               << (rg::HasCommonOption(
                    common_options,
                    rg::CommonOption::OutputFolder) ? "yes" : "no")
               << " |\n";
    }
    return output.str();
}

std::string BuildExpectedPythonSurfaceBlock()
{
    std::vector<std::string_view> command_types;
#define RHBM_GEM_COMMAND(COMMAND_ID, COMMAND_TYPE, CLI_NAME, DESCRIPTION, PROFILE) \
    command_types.emplace_back(#COMMAND_TYPE);
#include "internal/CommandList.def"
#undef RHBM_GEM_COMMAND

    std::ostringstream output;
    output << "### Python request types\n";
    output << "- `CommonCommandRequest`\n";
    for (const auto & command_type : command_types)
    {
        const std::string type_text{ command_type };
        const std::string stem{
            type_text.size() >= 7 && type_text.substr(type_text.size() - 7) == "Command"
                ? type_text.substr(0, type_text.size() - 7)
                : type_text
        };
        output << "- `" << stem << "Request`\n";
    }

    output << "\n### Python run functions\n";
    for (const auto & command_type : command_types)
    {
        const std::string type_text{ command_type };
        const std::string stem{
            type_text.size() >= 7 && type_text.substr(type_text.size() - 7) == "Command"
                ? type_text.substr(0, type_text.size() - 7)
                : type_text
        };
        output << "- `Run" << stem << "(...)`\n";
    }

    output << "\n### Shared diagnostics types\n";
    output << "- `LogLevel`\n";
    output << "- `ValidationPhase`\n";
    output << "- `ValidationIssue`\n";
    output << "- `ExecutionReport`\n";

    return output.str();
}

} // namespace

TEST(DocsSyncTest, CommandArchitectureGeneratedSectionsMatchCatalog)
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
            "<!-- BEGIN GENERATED: command-manifest -->",
            "<!-- END GENERATED: command-manifest -->")
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
            "<!-- BEGIN GENERATED: command-python-surface -->",
            "<!-- END GENERATED: command-python-surface -->")
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

    EXPECT_EQ(document.find("RegisterCLIOptionsExtend"), std::string::npos);
    EXPECT_EQ(document.find("RegisterCLIOptionsBasic"), std::string::npos);
}
