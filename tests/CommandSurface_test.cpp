#include <gtest/gtest.h>

#include "BuiltInCommandCatalog.hpp"

namespace rg = rhbm_gem;

namespace {

const rg::CommandDescriptor * FindCommand(rg::CommandId id)
{
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        if (descriptor.id == id)
        {
            return &descriptor;
        }
    }
    return nullptr;
}

} // namespace

TEST(CommandSurfaceTest, BuiltInCommandSurfacePoliciesMatchExpectedMetadata)
{
    const auto * potential_analysis{ FindCommand(rg::CommandId::PotentialAnalysis) };
    ASSERT_NE(potential_analysis, nullptr);
    EXPECT_TRUE(rg::HasCommonOption(
        potential_analysis->surface.common_options,
        rg::CommonOption::Database));
    EXPECT_EQ(potential_analysis->database_usage, rg::DatabaseUsage::Required);
    EXPECT_EQ(potential_analysis->binding_exposure, rg::BindingExposure::PythonPublic);

    const auto * potential_display{ FindCommand(rg::CommandId::PotentialDisplay) };
    ASSERT_NE(potential_display, nullptr);
    EXPECT_TRUE(rg::HasCommonOption(
        potential_display->surface.common_options,
        rg::CommonOption::Database));
    EXPECT_EQ(potential_display->database_usage, rg::DatabaseUsage::Required);
    EXPECT_EQ(potential_display->binding_exposure, rg::BindingExposure::PythonPublic);

    const auto * result_dump{ FindCommand(rg::CommandId::ResultDump) };
    ASSERT_NE(result_dump, nullptr);
    EXPECT_TRUE(rg::HasCommonOption(
        result_dump->surface.common_options,
        rg::CommonOption::Database));
    EXPECT_EQ(result_dump->database_usage, rg::DatabaseUsage::Required);
    EXPECT_EQ(result_dump->binding_exposure, rg::BindingExposure::PythonPublic);

    const auto * map_simulation{ FindCommand(rg::CommandId::MapSimulation) };
    ASSERT_NE(map_simulation, nullptr);
    EXPECT_FALSE(rg::HasCommonOption(
        map_simulation->surface.common_options,
        rg::CommonOption::Database));
    EXPECT_EQ(map_simulation->database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(map_simulation->binding_exposure, rg::BindingExposure::PythonPublic);

    const auto * map_visualization{ FindCommand(rg::CommandId::MapVisualization) };
    ASSERT_NE(map_visualization, nullptr);
    EXPECT_FALSE(rg::HasCommonOption(
        map_visualization->surface.common_options,
        rg::CommonOption::Database));
    EXPECT_EQ(map_visualization->database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(map_visualization->binding_exposure, rg::BindingExposure::CliOnly);

    const auto * position_estimation{ FindCommand(rg::CommandId::PositionEstimation) };
    ASSERT_NE(position_estimation, nullptr);
    EXPECT_FALSE(rg::HasCommonOption(
        position_estimation->surface.common_options,
        rg::CommonOption::Database));
    EXPECT_EQ(position_estimation->database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(position_estimation->binding_exposure, rg::BindingExposure::CliOnly);

    const auto * model_test{ FindCommand(rg::CommandId::ModelTest) };
    ASSERT_NE(model_test, nullptr);
    EXPECT_FALSE(rg::HasCommonOption(
        model_test->surface.common_options,
        rg::CommonOption::Database));
    EXPECT_EQ(model_test->database_usage, rg::DatabaseUsage::NotUsed);
    EXPECT_EQ(model_test->binding_exposure, rg::BindingExposure::CliOnly);
}
