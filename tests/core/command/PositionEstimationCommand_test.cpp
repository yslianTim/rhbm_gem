#include <gtest/gtest.h>

#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/PositionEstimationCommand.hpp>

namespace rg = rhbm_gem;

TEST(PositionEstimationCommandTest, ExecuteDoesNotRequireDatabaseConfiguration)
{
    command_test::ScopedTempDir temp_dir{"position_estimation"};
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::PositionEstimationCommand command{};
    command.SetFolderPath(output_dir);
    command.SetMapFilePath(map_path);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".cmm"), 1);
}
