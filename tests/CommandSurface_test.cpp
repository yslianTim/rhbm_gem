#include <gtest/gtest.h>

#include <string_view>

#include "CommandRegistry.hpp"
#include "HRLModelTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "RegisterBuiltInCommands.hpp"
#include "ResultDumpCommand.hpp"

namespace rg = rhbm_gem;

namespace {

const rg::CommandRegistry::CommandInfo * FindCommand(std::string_view name)
{
    rg::RegisterBuiltInCommands();
    const auto & commands{ rg::CommandRegistry::Instance().GetCommands() };
    for (const auto & info : commands)
    {
        if (info.name == name)
        {
            return &info;
        }
    }
    return nullptr;
}

template <typename CommandType>
void ExpectCommandInfoMatches()
{
    const auto * info{ FindCommand(CommandType::CommandName()) };
    ASSERT_NE(info, nullptr);

    const auto expected_surface{ CommandType::StaticCommandSurface() };
    EXPECT_EQ(info->description, CommandType::CommandDescription());
    EXPECT_EQ(info->surface.common_options, expected_surface.common_options);
    EXPECT_EQ(
        info->surface.deprecated_hidden_options,
        expected_surface.deprecated_hidden_options);
    EXPECT_EQ(
        info->surface.requires_database_manager,
        expected_surface.requires_database_manager);
    EXPECT_EQ(info->surface.uses_output_folder, expected_surface.uses_output_folder);
    EXPECT_EQ(info->surface.exposed_to_python, expected_surface.exposed_to_python);
}

} // namespace

TEST(CommandSurfaceTest, RegistrySurfacesMatchConcreteCommandMetadata)
{
    ExpectCommandInfoMatches<rg::PotentialAnalysisCommand>();
    ExpectCommandInfoMatches<rg::PotentialDisplayCommand>();
    ExpectCommandInfoMatches<rg::ResultDumpCommand>();
    ExpectCommandInfoMatches<rg::MapSimulationCommand>();
    ExpectCommandInfoMatches<rg::MapVisualizationCommand>();
    ExpectCommandInfoMatches<rg::PositionEstimationCommand>();
    ExpectCommandInfoMatches<rg::HRLModelTestCommand>();
}
