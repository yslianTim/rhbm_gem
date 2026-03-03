#include <gtest/gtest.h>

#include <string>

#include <CLI/CLI.hpp>

#include "Application.hpp"

namespace rg = rhbm_gem;

TEST(ApplicationTest, FileOnlyCommandsHideDatabaseOptionFromHelp)
{
    CLI::App app{"RHBM-GEM"};
    rg::Application controller(app);

    auto * map_simulation{ app.get_subcommand("map_simulation") };
    ASSERT_NE(map_simulation, nullptr);
    const std::string map_simulation_help{
        map_simulation->help(map_simulation->get_name(), CLI::AppFormatMode::Sub)
    };
    EXPECT_EQ(map_simulation_help.find("--database"), std::string::npos);

    auto * map_visualization{ app.get_subcommand("map_visualization") };
    ASSERT_NE(map_visualization, nullptr);
    const std::string map_visualization_help{
        map_visualization->help(map_visualization->get_name(), CLI::AppFormatMode::Sub)
    };
    EXPECT_EQ(map_visualization_help.find("--database"), std::string::npos);

    auto * potential_analysis{ app.get_subcommand("potential_analysis") };
    ASSERT_NE(potential_analysis, nullptr);
    const std::string potential_analysis_help{
        potential_analysis->help(potential_analysis->get_name(), CLI::AppFormatMode::Sub)
    };
    EXPECT_NE(potential_analysis_help.find("--database"), std::string::npos);
}
