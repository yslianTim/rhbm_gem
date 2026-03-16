#include <gtest/gtest.h>

#include <algorithm>
#include <array>
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
#include "internal/CommandCatalog.hpp"
#include <CLI/CLI.hpp>

namespace rg = rhbm_gem;

namespace {

struct ExpectedCommandMetadata
{
    rg::CommandId id;
    std::string_view name;
    rg::CommonOptionProfile profile;
    bool uses_database;
};

std::vector<ExpectedCommandMetadata> BuildExpectedCommandMetadata()
{
    std::vector<ExpectedCommandMetadata> expected;
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    expected.push_back(ExpectedCommandMetadata{                                                 \
        rg::CommandId::COMMAND_ID,                                                              \
        CLI_NAME,                                                                               \
        rg::CommonOptionProfile::PROFILE,                                                       \
        rg::CommonOptionProfile::PROFILE == rg::CommonOptionProfile::DatabaseWorkflow});
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
    return expected;
}

std::vector<std::pair<std::string_view, rg::CommandId>> BuildExpectedCommandIdTokens()
{
    std::vector<std::pair<std::string_view, rg::CommandId>> expected;
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    expected.emplace_back(#COMMAND_ID, rg::CommandId::COMMAND_ID);
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
    return expected;
}

std::vector<std::string> ParseCommandIdTokensFromManifest()
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

    const std::string manifest{
        std::istreambuf_iterator<char>{ manifest_stream },
        std::istreambuf_iterator<char>{}
    };
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

std::string ReadFileContent(const std::filesystem::path & path)
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

} // namespace

TEST(CommandCatalogTest, CatalogMatchesExpectedMetadataAndOrder)
{
    const auto expected_metadata{ BuildExpectedCommandMetadata() };
    const auto & catalog{ rg::CommandCatalog() };
    ASSERT_EQ(catalog.size(), expected_metadata.size());

    std::unordered_set<int> unique_ids;
    std::unordered_set<std::string> unique_names;
    for (std::size_t index = 0; index < expected_metadata.size(); ++index)
    {
        const auto & descriptor{ catalog[index] };
        const auto & expected{ expected_metadata[index] };
        const auto common_options{ rg::CommonOptionsForCommand(descriptor) };

        EXPECT_EQ(descriptor.id, expected.id);
        EXPECT_EQ(std::string_view{ descriptor.name }, expected.name);
        EXPECT_EQ(descriptor.profile, expected.profile);
        EXPECT_EQ(
            rg::HasCommonOption(common_options, rg::CommonOption::Database),
            expected.uses_database);
        EXPECT_TRUE(rg::HasCommonOption(common_options, rg::CommonOption::Threading));
        EXPECT_TRUE(rg::HasCommonOption(common_options, rg::CommonOption::Verbose));
        EXPECT_TRUE(rg::HasCommonOption(common_options, rg::CommonOption::OutputFolder));

        unique_ids.insert(static_cast<int>(descriptor.id));
        unique_names.emplace(descriptor.name);
    }

    EXPECT_EQ(unique_ids.size(), catalog.size());
    EXPECT_EQ(unique_names.size(), catalog.size());
}

TEST(CommandCatalogTest, RegisterCommandSubcommandsBuildsOneSubcommandPerManifestEntry)
{
    CLI::App app{"catalog"};
    rg::RegisterCommandSubcommands(app);

    const auto subcommands{
        app.get_subcommands([](CLI::App * subcommand)
        {
            return subcommand != nullptr && !subcommand->get_name().empty();
        })
    };
    const auto expected_metadata{ BuildExpectedCommandMetadata() };
    ASSERT_EQ(subcommands.size(), expected_metadata.size());
    for (std::size_t index = 0; index < expected_metadata.size(); ++index)
    {
        EXPECT_EQ(subcommands[index]->get_name(), expected_metadata[index].name);
    }
}

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
    for (std::size_t index = 0; index < expected_metadata.size(); ++index)
    {
        EXPECT_EQ(subcommands[index]->get_name(), expected_metadata[index].name);
    }
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

TEST(CommandCatalogTest, RunSurfaceAndPythonBindingsMatchManifestOrder)
{
    const auto project_root{ command_test::ProjectRootPath() };
    const auto command_api{
        ReadFileContent(
            project_root / "include" / "rhbm_gem" / "core" / "command" / "CommandApi.hpp")
    };
    const auto bindings{
        ReadFileContent(project_root / "bindings" / "CommandApiBindings.cpp")
    };

    ASSERT_FALSE(command_api.empty());
    ASSERT_FALSE(bindings.empty());
    EXPECT_NE(
        command_api.find("ExecutionReport Run##COMMAND_ID"),
        std::string::npos);
    EXPECT_NE(
        command_api.find("#include <rhbm_gem/core/command/CommandList.def>"),
        std::string::npos);
    EXPECT_NE(
        bindings.find("module.def(\"Run\" #COMMAND_ID"),
        std::string::npos);
    EXPECT_NE(
        bindings.find("#include <rhbm_gem/core/command/CommandList.def>"),
        std::string::npos);
}

TEST(CommandCatalogTest, PythonBindingsExposeRequestAndReportSurface)
{
    const auto bindings{
        ReadFileContent(command_test::ProjectRootPath() / "bindings" / "CommandApiBindings.cpp")
    };

    ASSERT_FALSE(bindings.empty());

    for (const std::string_view request_name : {
             "CommonCommandRequest",
             "PotentialAnalysisRequest",
             "PotentialDisplayRequest",
             "ResultDumpRequest",
             "MapSimulationRequest",
             "MapVisualizationRequest",
             "PositionEstimationRequest",
             "HRLModelTestRequest",
         })
    {
        EXPECT_NE(
            bindings.find(
                "py::class_<" + std::string(request_name) + ">(module, \""
                + std::string(request_name) + "\")"),
            std::string::npos)
            << request_name;
    }

    for (const std::string_view common_field : {
             "thread_size",
             "verbose_level",
             "database_path",
             "folder_path",
         })
    {
        EXPECT_NE(
            bindings.find(".def_readwrite(\"" + std::string(common_field) + "\""),
            std::string::npos)
            << common_field;
    }

    for (const std::string_view report_field : {
             "prepared",
             "executed",
             "validation_issues",
         })
    {
        EXPECT_NE(
            bindings.find(".def_readonly(\"" + std::string(report_field) + "\""),
            std::string::npos)
            << report_field;
    }
}

TEST(CommandCatalogTest, ModelTestCliNameMapsToHRLModelTestCommandId)
{
    const auto & catalog{ rg::CommandCatalog() };
    const auto iter = std::find_if(
        catalog.begin(),
        catalog.end(),
        [](const rg::CommandDescriptor & descriptor)
        {
            return std::string_view{ descriptor.name } == "model_test";
        });

    ASSERT_NE(iter, catalog.end());
    EXPECT_EQ(iter->id, rg::CommandId::HRLModelTest);
}

TEST(CommandCatalogTest, ConfigureCommandCliSharedOptionsMatchCommandMetadata)
{
    CLI::App app{"RHBM-GEM"};
    rg::ConfigureCommandCli(app);

    for (const auto & descriptor : rg::CommandCatalog())
    {
        auto * subcommand{ app.get_subcommand(std::string(descriptor.name)) };
        ASSERT_NE(subcommand, nullptr) << descriptor.name;
        const auto common_options{ rg::CommonOptionsForCommand(descriptor) };

        const std::string help_text{
            subcommand->help(subcommand->get_name(), CLI::AppFormatMode::Sub)
        };
        EXPECT_EQ(
            help_text.find("--database") != std::string::npos,
            rg::HasCommonOption(common_options, rg::CommonOption::Database))
            << descriptor.name;
        EXPECT_EQ(
            help_text.find("--folder") != std::string::npos,
            rg::HasCommonOption(common_options, rg::CommonOption::OutputFolder))
            << descriptor.name;
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

TEST(CommandCatalogTest, CommandTestsDoNotIncludeCommandPrivateWorkflowHeaders)
{
    const auto command_test_dir{
        command_test::ProjectRootPath() / "tests" / "core" / "command"
    };
    const std::array<std::string_view, 1> forbidden_includes{
        "internal/workflow/",
    };

    for (const auto & entry : std::filesystem::directory_iterator{ command_test_dir })
    {
        if (!entry.is_regular_file() || entry.path().extension() != ".cpp")
        {
            continue;
        }

        const auto content{ ReadFileContent(entry.path()) };
        ASSERT_FALSE(content.empty()) << entry.path().string();
        for (const auto forbidden : forbidden_includes)
        {
            EXPECT_EQ(content.find(forbidden), std::string::npos)
                << entry.path().filename().string();
        }
    }
}
