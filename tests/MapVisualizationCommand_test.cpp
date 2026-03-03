#include <gtest/gtest.h>

#include <algorithm>

#include "CommandTestHelpers.hpp"
#include "MapVisualizationCommand.hpp"

namespace rg = rhbm_gem;

TEST(MapVisualizationCommandTest, NonPositiveAtomIdBecomesValidationErrorAtPrepare)
{
    rg::MapVisualizationCommand command;
    command.SetAtomSerialID(0);

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--atom-id";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(MapVisualizationCommandTest, NonPositiveWindowSizeBecomesValidationErrorAtPrepare)
{
    rg::MapVisualizationCommand command;
    command.SetWindowSize(0.0);

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--window-size";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(MapVisualizationCommandTest, SamplingSizeNormalizationUsesOptionsDefault)
{
    rg::MapVisualizationCommand command;
    command.SetSamplingSize(0);

    const auto & options{
        static_cast<const rg::MapVisualizationCommand::Options &>(command.GetOptions())
    };
    EXPECT_EQ(options.sampling_size, 100);

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--sampling";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Warning);
    EXPECT_TRUE(issue_iter->auto_corrected);
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
