#include <gtest/gtest.h>

#include "BuiltInCommandCatalogInternal.hpp"

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
    EXPECT_TRUE(rg::UsesDatabaseAtRuntime(potential_analysis->surface));
    EXPECT_TRUE(rg::UsesOutputFolder(potential_analysis->surface));
    EXPECT_EQ(potential_analysis->python_binding_name, "PotentialAnalysisCommand");

    const auto * potential_display{ FindCommand(rg::CommandId::PotentialDisplay) };
    ASSERT_NE(potential_display, nullptr);
    EXPECT_TRUE(rg::UsesDatabaseAtRuntime(potential_display->surface));
    EXPECT_TRUE(rg::UsesOutputFolder(potential_display->surface));
    EXPECT_EQ(potential_display->python_binding_name, "PotentialDisplayCommand");

    const auto * result_dump{ FindCommand(rg::CommandId::ResultDump) };
    ASSERT_NE(result_dump, nullptr);
    EXPECT_TRUE(rg::UsesDatabaseAtRuntime(result_dump->surface));
    EXPECT_TRUE(rg::UsesOutputFolder(result_dump->surface));
    EXPECT_EQ(result_dump->python_binding_name, "ResultDumpCommand");

    const auto * map_simulation{ FindCommand(rg::CommandId::MapSimulation) };
    ASSERT_NE(map_simulation, nullptr);
    EXPECT_FALSE(rg::UsesDatabaseAtRuntime(map_simulation->surface));
    EXPECT_TRUE(rg::UsesOutputFolder(map_simulation->surface));
    EXPECT_EQ(map_simulation->python_binding_name, "MapSimulationCommand");

    const auto * map_visualization{ FindCommand(rg::CommandId::MapVisualization) };
    ASSERT_NE(map_visualization, nullptr);
    EXPECT_FALSE(rg::UsesDatabaseAtRuntime(map_visualization->surface));
    EXPECT_TRUE(rg::UsesOutputFolder(map_visualization->surface));
    EXPECT_EQ(map_visualization->python_binding_name, "MapVisualizationCommand");

    const auto * position_estimation{ FindCommand(rg::CommandId::PositionEstimation) };
    ASSERT_NE(position_estimation, nullptr);
    EXPECT_FALSE(rg::UsesDatabaseAtRuntime(position_estimation->surface));
    EXPECT_TRUE(rg::UsesOutputFolder(position_estimation->surface));
    EXPECT_EQ(position_estimation->python_binding_name, "PositionEstimationCommand");

    const auto * model_test{ FindCommand(rg::CommandId::ModelTest) };
    ASSERT_NE(model_test, nullptr);
    EXPECT_FALSE(rg::UsesDatabaseAtRuntime(model_test->surface));
    EXPECT_TRUE(rg::UsesOutputFolder(model_test->surface));
    EXPECT_EQ(model_test->python_binding_name, "HRLModelTestCommand");
}
