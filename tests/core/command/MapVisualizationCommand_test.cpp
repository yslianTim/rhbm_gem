#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/MapVisualizationCommand.hpp>

namespace rg = rhbm_gem;

TEST(MapVisualizationCommandTest, NonPositiveAtomIdBecomesParseValidationError)
{
    rg::MapVisualizationCommand command{};
    rg::MapVisualizationRequest request{};
    request.atom_serial_id = 0;
    command.ApplyRequest(request);

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
    rg::MapVisualizationCommand command{};
    rg::MapVisualizationRequest request{};
    request.window_size = 0.0;
    command.ApplyRequest(request);

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
    rg::MapVisualizationCommand command{};
    rg::MapVisualizationRequest request{};
    request.sampling_size = 0;
    command.ApplyRequest(request);

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
