#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <string>
#include <vector>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandSystem.hpp>

namespace rg = rhbm_gem;

namespace {

std::vector<std::string> BuildExpectedCommandNames()
{
    std::vector<std::string> expected;
    for (const auto & command : rg::ListCommands())
    {
        expected.emplace_back(command.name);
    }
    return expected;
}

std::string ReadManifestContent()
{
    const auto manifest_path{
        command_test::ProjectRootPath() / "include" / "rhbm_gem" / "core" / "command"
            / "CommandManifest.def"
    };
    std::ifstream manifest_stream{ manifest_path };
    if (!manifest_stream.is_open())
    {
        return {};
    }

    return std::string{
        std::istreambuf_iterator<char>{ manifest_stream },
        std::istreambuf_iterator<char>{}
    };
}

std::vector<std::string> ExtractCommandNames(const std::string & manifest)
{
    const std::regex command_pattern{
        R"manifest(RHBM_GEM_COMMAND\(\s*[A-Za-z_][A-Za-z0-9_]*\s*,\s*"([^"]+)")manifest"
    };

    std::vector<std::string> command_names;
    for (std::sregex_iterator iter{ manifest.begin(), manifest.end(), command_pattern };
        iter != std::sregex_iterator{};
        ++iter)
    {
        command_names.push_back((*iter)[1].str());
    }
    return command_names;
}

std::vector<std::string> ParseCommandNamesFromManifest()
{
    auto manifest{ ReadManifestContent() };
    if (manifest.empty())
    {
        return {};
    }

#ifndef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    const std::regex experimental_block_pattern{
        R"(#ifdef\s+RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE[\s\S]*?#endif\s*)"
    };
    manifest = std::regex_replace(manifest, experimental_block_pattern, "");
#endif

    return ExtractCommandNames(manifest);
}

std::vector<std::string> ParseExperimentalCommandNamesFromManifest()
{
    const auto manifest{ ReadManifestContent() };
    if (manifest.empty())
    {
        return {};
    }

    const std::regex experimental_block_pattern{
        R"(#ifdef\s+RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE\s*([\s\S]*?)#endif)"
    };
    std::smatch match;
    if (!std::regex_search(manifest, match, experimental_block_pattern) || match.size() < 2)
    {
        return {};
    }

    return ExtractCommandNames(match[1].str());
}

} // namespace

TEST(CommandCatalogTest, ListCommandsMatchesManifestOrder)
{
    EXPECT_EQ(ParseCommandNamesFromManifest(), BuildExpectedCommandNames());
}

TEST(CommandCatalogTest, ManifestKeepsExperimentalEntriesBehindFeatureGuard)
{
    const std::vector<std::string> expected_experimental_names{
        "map_visualization",
        "position_estimation",
    };
    EXPECT_EQ(ParseExperimentalCommandNamesFromManifest(), expected_experimental_names);
}

TEST(CommandCatalogTest, RunCommandCLIReturnsSuccessForHelpRequest)
{
    char program[]{ "RHBM-GEM" };
    char help_flag[]{ "--help" };
    char * argv[]{ program, help_flag };
    const int argc{ static_cast<int>(std::size(argv)) };

    EXPECT_EQ(rg::RunCommandCLI(argc, argv), 0);
}

TEST(CommandCatalogTest, RunCommandCLIReturnsFailureForParseErrors)
{
    char program[]{ "RHBM-GEM" };
    char command[]{ "map_simulation" };
    char * argv[]{ program, command };
    const int argc{ static_cast<int>(std::size(argv)) };

    EXPECT_NE(rg::RunCommandCLI(argc, argv), 0);
}
