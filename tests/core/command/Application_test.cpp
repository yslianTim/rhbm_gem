#include <gtest/gtest.h>

#include <string>

#include <CLI/CLI.hpp>

#include <rhbm_gem/core/command/Application.hpp>
#include "internal/CommandCatalogInternal.hpp"

namespace rg = rhbm_gem;

TEST(ApplicationTest, HelpOrderMatchesCommandCatalog)
{
    CLI::App app{"RHBM-GEM"};
    rg::Application controller(app);

    const auto subcommands{
        app.get_subcommands([](CLI::App * subcommand)
        {
            return subcommand != nullptr && !subcommand->get_name().empty();
        })
    };

    const auto & catalog{ rg::CommandCatalog() };
    ASSERT_EQ(subcommands.size(), catalog.size());
    for (std::size_t index = 0; index < catalog.size(); ++index)
    {
        EXPECT_EQ(subcommands[index]->get_name(), catalog[index].name);
    }
}

TEST(ApplicationTest, SharedOptionsMatchCommandMetadata)
{
    CLI::App app{"RHBM-GEM"};
    rg::Application controller(app);

    for (const auto & descriptor : rg::CommandCatalog())
    {
        auto * subcommand{ app.get_subcommand(std::string(descriptor.name)) };
        ASSERT_NE(subcommand, nullptr) << descriptor.name;

        const std::string help_text{
            subcommand->help(subcommand->get_name(), CLI::AppFormatMode::Sub)
        };
        EXPECT_EQ(
            help_text.find("--database") != std::string::npos,
            rg::UsesDatabaseAtRuntime(descriptor.common_options))
            << descriptor.name;
        EXPECT_EQ(
            help_text.find("--folder") != std::string::npos,
            rg::UsesOutputFolder(descriptor.common_options))
            << descriptor.name;
    }
}

TEST(ApplicationTest, CommandFailurePropagatesAsRuntimeError)
{
    CLI::App app{"RHBM-GEM"};
    rg::Application controller(app);

    EXPECT_THROW(
        app.parse("map_simulation --model missing.cif --blurring-width 1.0", false),
        CLI::RuntimeError);
}
