#include <gtest/gtest.h>

#include "BuiltInCommandCatalog.hpp"
#include "CommandRegistry.hpp"
#include "RegisterBuiltInCommands.hpp"

namespace rg = rhbm_gem;

namespace {

const rg::CommandRegistry::CommandInfo * FindCommand(rg::CommandId id)
{
    rg::RegisterBuiltInCommands();
    const auto & commands{ rg::CommandRegistry::Instance().GetCommands() };
    for (const auto & info : commands)
    {
        if (info.id == id)
        {
            return &info;
        }
    }
    return nullptr;
}

void ExpectCommandInfoMatches(const rg::CommandDescriptor & descriptor)
{
    const auto * info{ FindCommand(descriptor.id) };
    ASSERT_NE(info, nullptr);
    EXPECT_EQ(info->name, descriptor.name);
    EXPECT_EQ(info->description, descriptor.description);
    EXPECT_EQ(info->surface.common_options, descriptor.surface.common_options);
    EXPECT_EQ(
        info->surface.deprecated_hidden_options,
        descriptor.surface.deprecated_hidden_options);
    EXPECT_EQ(info->database_usage, descriptor.database_usage);
    EXPECT_EQ(info->binding_exposure, descriptor.binding_exposure);
    EXPECT_EQ(info->python_binding_name, descriptor.python_binding_name);
}

} // namespace

TEST(CommandSurfaceTest, RegistrySurfacesMatchBuiltInCommandCatalog)
{
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        ExpectCommandInfoMatches(descriptor);
    }
}

TEST(CommandSurfaceTest, BuiltInCommandMetadataMatchesExpectedSurfacePolicies)
{
    const auto & potential_analysis{
        rg::FindCommandDescriptor(rg::CommandId::PotentialAnalysis)
    };
    EXPECT_TRUE(rg::HasCommonOption(
        potential_analysis.surface.common_options, rg::CommonOption::Database));
    EXPECT_EQ(potential_analysis.database_usage, rg::DatabaseUsage::Required);
    EXPECT_EQ(potential_analysis.binding_exposure, rg::BindingExposure::PythonPublic);

    const auto & potential_display{
        rg::FindCommandDescriptor(rg::CommandId::PotentialDisplay)
    };
    EXPECT_TRUE(rg::HasCommonOption(
        potential_display.surface.common_options, rg::CommonOption::Database));
    EXPECT_EQ(potential_display.database_usage, rg::DatabaseUsage::Required);
    EXPECT_EQ(potential_display.binding_exposure, rg::BindingExposure::PythonPublic);

    const auto & result_dump{
        rg::FindCommandDescriptor(rg::CommandId::ResultDump)
    };
    EXPECT_TRUE(rg::HasCommonOption(
        result_dump.surface.common_options, rg::CommonOption::Database));
    EXPECT_EQ(result_dump.database_usage, rg::DatabaseUsage::Required);
    EXPECT_EQ(result_dump.binding_exposure, rg::BindingExposure::PythonPublic);

    const auto & map_simulation{
        rg::FindCommandDescriptor(rg::CommandId::MapSimulation)
    };
    EXPECT_FALSE(rg::HasCommonOption(
        map_simulation.surface.common_options, rg::CommonOption::Database));
    EXPECT_TRUE(rg::HasCommonOption(
        map_simulation.surface.deprecated_hidden_options, rg::CommonOption::Database));
    EXPECT_EQ(map_simulation.database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(map_simulation.binding_exposure, rg::BindingExposure::PythonPublic);

    const auto & map_visualization{
        rg::FindCommandDescriptor(rg::CommandId::MapVisualization)
    };
    EXPECT_EQ(map_visualization.database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(map_visualization.binding_exposure, rg::BindingExposure::CliOnly);

    const auto & position_estimation{
        rg::FindCommandDescriptor(rg::CommandId::PositionEstimation)
    };
    EXPECT_EQ(position_estimation.database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(position_estimation.binding_exposure, rg::BindingExposure::CliOnly);

    const auto & model_test{
        rg::FindCommandDescriptor(rg::CommandId::ModelTest)
    };
    EXPECT_EQ(model_test.database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(model_test.binding_exposure, rg::BindingExposure::CliOnly);
}
