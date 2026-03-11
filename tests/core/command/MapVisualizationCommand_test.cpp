#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/MapVisualizationCommand.hpp>

namespace rg = rhbm_gem;

TEST(MapVisualizationCommandTest, NonPositiveAtomIdBecomesParseValidationError)
{
    rg::MapVisualizationCommand command{ command_test::BuildDataIoServices() };
    command.SetAtomSerialID(0);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--atom-id",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(MapVisualizationCommandTest, NonPositiveWindowSizeBecomesParseValidationError)
{
    rg::MapVisualizationCommand command{ command_test::BuildDataIoServices() };
    command.SetWindowSize(0.0);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--window-size",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(MapVisualizationCommandTest, SamplingSizeNormalizationReportsAutoCorrectedWarning)
{
    rg::MapVisualizationCommand command{ command_test::BuildDataIoServices() };
    command.SetSamplingSize(0);

    const auto * issue{
        command_test::FindValidationIssue(
            command,
            "--sampling",
            rg::ValidationPhase::Parse,
            LogLevel::Warning)
    };
    ASSERT_NE(issue, nullptr);
    EXPECT_TRUE(issue->auto_corrected);
}

TEST(MapVisualizationCommandTest, InvalidAtomIdFailsWithoutWritingOutput)
{
    command_test::ScopedTempDir temp_dir{"map_visualization_invalid"};
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")
    };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::MapVisualizationCommand command{ command_test::BuildDataIoServices() };
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

    rg::MapVisualizationCommand command{ command_test::BuildDataIoServices() };
    command.SetFolderPath(output_dir);
    command.SetModelFilePath(model_path);
    command.SetMapFilePath(map_path);
    command.SetAtomSerialID(1);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 1);
}
#endif
