#include <gtest/gtest.h>

#include "CommandTestHelpers.hpp"
#include "MapVisualizationCommand.hpp"

namespace rg = rhbm_gem;

TEST(MapVisualizationCommandTest, InvalidAtomIdFailsWithoutWritingOutput)
{
    command_test::ScopedTempDir temp_dir{"map_visualization_invalid"};
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")
    };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::MapVisualizationCommand command;
    command.SetFolderPath(output_dir);
    command.SetModelFilePath(model_path);
    command.SetMapFilePath(map_path);
    command.SetAtomSerialID(999);

    EXPECT_FALSE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 0);
}

#ifdef HAVE_ROOT
TEST(MapVisualizationCommandTest, ExecuteWritesPdfToConfiguredFolder)
{
    command_test::ScopedTempDir temp_dir{"map_visualization_valid"};
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")
    };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::MapVisualizationCommand command;
    command.SetFolderPath(output_dir);
    command.SetModelFilePath(model_path);
    command.SetMapFilePath(map_path);
    command.SetAtomSerialID(1);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 1);
}
#endif
