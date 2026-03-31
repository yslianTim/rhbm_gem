#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/CommandContract.hpp>
#include <CLI/CLI.hpp>

namespace rg = rhbm_gem;

namespace {

std::vector<rg::CommandDescriptor> BuildExpectedCommandMetadata()
{
    const auto catalog{ rg::GetCommandCatalog() };
    return { catalog.begin(), catalog.end() };
}

std::vector<std::pair<std::string_view, rg::CommandId>> BuildExpectedCommandIdTokens()
{
    std::vector<std::pair<std::string_view, rg::CommandId>> expected;
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    expected.emplace_back(#COMMAND_ID, rg::CommandId::COMMAND_ID);
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
    return expected;
}

bool CommandUsesDatabase(std::string_view cli_name)
{
    return cli_name == "potential_analysis"
        || cli_name == "potential_display"
        || cli_name == "result_dump";
}

std::string ReadManifestContent()
{
    const auto manifest_path{
        command_test::ProjectRootPath() / "include" / "rhbm_gem" / "core" / "command"
        / "CommandList.def"
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

std::vector<std::string> ExtractCommandIdTokens(const std::string & manifest)
{
    const std::regex command_pattern{
        R"(RHBM_GEM_COMMAND\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*,)"
    };

    std::vector<std::string> command_ids;
    for (std::sregex_iterator iter{ manifest.begin(), manifest.end(), command_pattern };
        iter != std::sregex_iterator{};
        ++iter)
    {
        command_ids.push_back((*iter)[1].str());
    }
    return command_ids;
}

std::vector<std::string> ParseCommandIdTokensFromManifest()
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

    return ExtractCommandIdTokens(manifest);
}

std::vector<std::string> ParseExperimentalCommandIdTokensFromManifest()
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

    return ExtractCommandIdTokens(match[1].str());
}

} // namespace

TEST(CommandCatalogTest, ConfigureCommandCliBuildsOneSubcommandPerManifestEntry)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    const auto subcommands{
        app.get_subcommands([](CLI::App * subcommand)
        {
            return subcommand != nullptr && !subcommand->get_name().empty();
        })
    };
    const auto expected_metadata{ BuildExpectedCommandMetadata() };
    ASSERT_EQ(subcommands.size(), expected_metadata.size());
    std::unordered_set<int> unique_ids;
    std::unordered_set<std::string> unique_names;
    for (std::size_t index = 0; index < expected_metadata.size(); ++index)
    {
        EXPECT_EQ(subcommands[index]->get_name(), expected_metadata[index].cli_name);
        EXPECT_EQ(subcommands[index]->get_description(), expected_metadata[index].description);

        unique_ids.insert(static_cast<int>(expected_metadata[index].id));
        unique_names.emplace(expected_metadata[index].cli_name);
    }

    EXPECT_EQ(unique_ids.size(), expected_metadata.size());
    EXPECT_EQ(unique_names.size(), expected_metadata.size());
}

TEST(CommandCatalogTest, CommandIdEnumMatchesManifestOrderAndIndexing)
{
    const auto expected_command_ids{ BuildExpectedCommandIdTokens() };
    const auto manifest_ids{ ParseCommandIdTokensFromManifest() };
    ASSERT_FALSE(manifest_ids.empty());
    ASSERT_EQ(manifest_ids.size(), expected_command_ids.size());

    for (std::size_t index = 0; index < expected_command_ids.size(); ++index)
    {
        const auto & expected{ expected_command_ids[index] };
        EXPECT_EQ(manifest_ids[index], expected.first);
        EXPECT_EQ(static_cast<int>(expected.second), static_cast<int>(index));
    }
}

TEST(CommandCatalogTest, ManifestKeepsExperimentalEntriesBehindFeatureGuard)
{
    const std::vector<std::string> expected_experimental_ids{
        "MapVisualization",
        "PositionEstimation",
    };
    EXPECT_EQ(ParseExperimentalCommandIdTokensFromManifest(), expected_experimental_ids);
}

TEST(CommandCatalogTest, ConfigureCommandCliSharedOptionsMatchCommandMetadata)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    for (const auto & descriptor : BuildExpectedCommandMetadata())
    {
        auto * subcommand{ app.get_subcommand(std::string(descriptor.cli_name)) };
        ASSERT_NE(subcommand, nullptr) << descriptor.cli_name;

        const std::string help_text{
            subcommand->help(subcommand->get_name(), CLI::AppFormatMode::Sub)
        };
        EXPECT_EQ(
            help_text.find("--database") != std::string::npos,
            CommandUsesDatabase(descriptor.cli_name))
            << descriptor.cli_name;
        EXPECT_NE(help_text.find("--folder"), std::string::npos) << descriptor.cli_name;
    }
}

TEST(CommandCatalogTest, ConfigureCommandCliPropagatesCommandFailureAsRuntimeError)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    EXPECT_THROW(
        app.parse("map_simulation --model missing.cif --blurring-width 1.0", false),
        CLI::RuntimeError);
}

TEST(CommandCatalogTest, ConfigureCommandCliAcceptsRepeatedReferenceGroups)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    EXPECT_THROW(
        app.parse(
            "potential_display --painter model --model-keylist model_a "
            "--ref-group A=m1,m2 --ref-group B=m3",
            false),
        CLI::RuntimeError);
}

TEST(CommandCatalogTest, ConfigureCommandCliEnforcesRequiredOptionsBeforeCallback)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    EXPECT_THROW(app.parse("map_simulation", false), CLI::RequiredError);
    EXPECT_THROW(app.parse("potential_display --model-keylist model_a", false), CLI::RequiredError);
}

TEST(CommandCatalogTest, ConfigureCommandCliAppliesEnumTransformBeforeCommandExecution)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    EXPECT_THROW(
        app.parse("potential_display --painter MODEL --model-keylist model_a", false),
        CLI::RuntimeError);
    EXPECT_THROW(
        app.parse("potential_display --painter not-a-painter --model-keylist model_a", false),
        CLI::ValidationError);
}

TEST(CommandCatalogTest, ConfigureCommandCliRejectsMalformedReferenceGroupsAtParseTime)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    EXPECT_THROW(
        app.parse(
            "potential_display --painter model --model-keylist model_a "
            "--ref-group =m1,m2",
            false),
        CLI::ValidationError);
    EXPECT_THROW(
        app.parse(
            "potential_display --painter model --model-keylist model_a "
            "--ref-group Group=",
            false),
        CLI::ValidationError);
}
