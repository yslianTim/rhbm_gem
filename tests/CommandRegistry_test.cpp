#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <unordered_set>
#include <vector>

#include <CLI/CLI.hpp>

#include "Application.hpp"
#include "CommandRegistry.hpp"
#include "RegisterBuiltInCommands.hpp"

namespace rg = rhbm_gem;

TEST(CommandRegistryTest, ContainsAllKnownCommandsWithoutDuplicates)
{
    rg::RegisterBuiltInCommands();
    const auto & commands{ rg::CommandRegistry::Instance().GetCommands() };
    std::vector<std::string> names;
    names.reserve(commands.size());
    for (const auto & info : commands)
    {
        names.push_back(info.name);
    }

    const std::unordered_set<std::string> unique_names(names.begin(), names.end());
    EXPECT_EQ(unique_names.size(), names.size());

    const std::vector<std::string> expected_names{
        "potential_analysis",
        "potential_display",
        "result_dump",
        "map_simulation",
        "map_visualization",
        "position_estimation",
        "model_test"
    };
    for (const auto & expected_name : expected_names)
    {
        EXPECT_NE(
            std::find(names.begin(), names.end(), expected_name),
            names.end());
    }
}

TEST(CommandRegistryTest, ApplicationUsesRegistryOrderForSubcommands)
{
    rg::RegisterBuiltInCommands();
    CLI::App app{"RHBM-GEM"};
    rg::Application controller(app);

    const auto & commands{ rg::CommandRegistry::Instance().GetCommands() };
    const auto subcommands{
        app.get_subcommands([](CLI::App * subcommand)
        {
            return subcommand != nullptr && !subcommand->get_name().empty();
        })
    };

    ASSERT_EQ(subcommands.size(), commands.size());
    for (std::size_t index = 0; index < commands.size(); ++index)
    {
        EXPECT_EQ(subcommands[index]->get_name(), commands[index].name);
    }
}
