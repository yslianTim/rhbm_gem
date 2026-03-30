#include <gtest/gtest.h>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/command/MapVisualizationCommand.hpp"

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

    rg::MapVisualizationCommand command{};
    rg::MapVisualizationRequest request{};
    request.common.folder_path = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 999;
    command.ApplyRequest(request);

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

    rg::MapVisualizationCommand command{};
    rg::MapVisualizationRequest request{};
    request.common.folder_path = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 1;
    command.ApplyRequest(request);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 1);
}
#endif
