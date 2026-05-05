#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "command/detail/CommandCatalog.hpp"
#include <rhbm_gem/core/command/CommandSystem.hpp>

namespace rg = rhbm_gem;

namespace {

std::vector<std::string> BuildExpectedCommandNames()
{
    std::vector<std::string> expected{
        "potential_analysis",
        "potential_display",
        "result_dump",
        "map_simulation",
        "rhbm_test",
    };
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    expected.emplace_back("map_visualization");
    expected.emplace_back("position_estimation");
#endif
    return expected;
}

std::vector<std::string> BuildCommandNames()
{
    std::vector<std::string> names;
    rg::command_internal::VisitCommandCatalog([&](const auto & entry)
    {
        names.emplace_back(entry.cli_name);
    });
    return names;
}

} // namespace

TEST(CommandCatalogTest, CommandCatalogKeepsExpectedOrder)
{
    EXPECT_EQ(BuildExpectedCommandNames(), BuildCommandNames());
}

TEST(CommandCatalogTest, ExperimentalCommandVisibilityFollowsFeatureGuard)
{
    const std::vector<std::string> expected_experimental_names{
        "map_visualization",
        "position_estimation",
    };
    const auto command_names{ BuildCommandNames() };

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    ASSERT_GE(command_names.size(), expected_experimental_names.size());
    const std::vector<std::string> trailing_names{
        command_names.end() - static_cast<std::ptrdiff_t>(expected_experimental_names.size()),
        command_names.end()
    };
    EXPECT_EQ(trailing_names, expected_experimental_names);
#else
    for (const auto & experimental_name : expected_experimental_names)
    {
        EXPECT_EQ(
            std::find(command_names.begin(), command_names.end(), experimental_name),
            command_names.end());
    }
#endif
}

TEST(CommandCatalogTest, RunCommandCLIReturnsSuccessForHelpRequest)
{
    char program[]{ "RHBM-GEM" };
    char help_flag[]{ "--help" };
    char * argv[]{ program, help_flag };
    const int argc{ static_cast<int>(std::size(argv)) };

    EXPECT_EQ(rg::RunCommandCLI(argc, argv), 0);
}

TEST(CommandCatalogTest, RunCommandCLIReturnsFailureForMissingCommandInput)
{
    char program[]{ "RHBM-GEM" };
    char command[]{ "map_simulation" };
    char * argv[]{ program, command };
    const int argc{ static_cast<int>(std::size(argv)) };

    EXPECT_NE(rg::RunCommandCLI(argc, argv), 0);
}
