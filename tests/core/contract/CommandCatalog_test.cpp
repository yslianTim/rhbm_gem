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

#include "internal/command/CommandCli.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <CLI/CLI.hpp>

namespace rg = rhbm_gem;

namespace {

std::vector<rg::CommandInfo> BuildExpectedCommandMetadata()
{
    const auto & catalog{ rg::ListCommands() };
    return { catalog.begin(), catalog.end() };
}

std::vector<std::string> BuildExpectedCommandNames()
{
    std::vector<std::string> expected;
    for (const auto & command : rg::ListCommands())
    {
        expected.emplace_back(command.name);
    }
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
        command_test::ProjectRootPath() / "src" / "core" / "command" / "CommandManifest.def"
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
    std::unordered_set<std::string> unique_names;
    for (std::size_t index = 0; index < expected_metadata.size(); ++index)
    {
        EXPECT_EQ(subcommands[index]->get_name(), expected_metadata[index].name);
        EXPECT_EQ(subcommands[index]->get_description(), expected_metadata[index].description);

        unique_names.emplace(expected_metadata[index].name);
    }

    EXPECT_EQ(unique_names.size(), expected_metadata.size());
}

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

TEST(CommandCatalogTest, ConfigureCommandCliSharedOptionsMatchCommandMetadata)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    for (const auto & descriptor : BuildExpectedCommandMetadata())
    {
        auto * subcommand{ app.get_subcommand(std::string(descriptor.name)) };
        ASSERT_NE(subcommand, nullptr) << descriptor.name;

        const std::string help_text{
            subcommand->help(subcommand->get_name(), CLI::AppFormatMode::Sub)
        };
        EXPECT_EQ(
            help_text.find("--database") != std::string::npos,
            CommandUsesDatabase(descriptor.name))
            << descriptor.name;
        EXPECT_NE(help_text.find("--folder"), std::string::npos) << descriptor.name;
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
