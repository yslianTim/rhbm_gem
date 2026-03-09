#include <gtest/gtest.h>

#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>

namespace rg = rhbm_gem;

TEST(PotentialAnalysisExperimentalBondGateTest, ExecuteRemainsReachableWhenFeatureIsEnabled)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_bond_gate" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto maps_dir{ temp_dir.path() / "maps" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };

    std::filesystem::create_directories(maps_dir);
    const auto map_path{
        command_test::GenerateMapFile(maps_dir, model_path, "bond_gate_map", "1.0")
    };

    rg::PotentialAnalysisCommand command;
    command.SetDatabasePath(database_path);
    command.SetModelFilePath(model_path);
    command.SetMapFilePath(map_path);
    command.SetSavedKeyTag("bond_gate_model");
    command.SetSamplingSize(25);

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(command.Execute());
}
