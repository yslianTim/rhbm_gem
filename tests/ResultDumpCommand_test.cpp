#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>

#include "CommandTestHelpers.hpp"
#include "ResultDumpCommand.hpp"

namespace rg = rhbm_gem;

TEST(ResultDumpCommandTest, MapPrinterRequiresMapFile)
{
    rg::ResultDumpCommand command;
    command.SetPrinterChoice(rg::PrinterType::MAP_VALUE);
    command.SetModelKeyTagList("model");

    EXPECT_FALSE(command.PrepareForExecution());
}

TEST(ResultDumpCommandTest, InvalidPrinterChoiceBecomesValidationError)
{
    rg::ResultDumpCommand command;
    command.SetPrinterChoice(static_cast<rg::PrinterType>(999));
    command.SetModelKeyTagList("model");

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--printer";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(ResultDumpCommandTest, ExecuteTwiceDoesNotReuseStaleLoadedModels)
{
    command_test::ScopedTempDir temp_dir{"result_dump"};
    const auto database_path{ temp_dir.path() / "database.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_dir_a{ temp_dir.path() / "out_a" };
    const auto output_dir_b{ temp_dir.path() / "out_b" };

    command_test::SeedSavedModel(database_path, model_path, "key_a", "MODEL_A");
    command_test::SeedSavedModel(database_path, model_path, "key_b", "MODEL_B");

    rg::ResultDumpCommand command;
    command.SetDatabasePath(database_path);
    command.SetPrinterChoice(rg::PrinterType::ATOM_POSITION);

    command.SetFolderPath(output_dir_a);
    command.SetModelKeyTagList("key_a");
    ASSERT_TRUE(command.Execute());
    EXPECT_TRUE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_A.csv"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir_a, ".csv"), 1);

    command.SetFolderPath(output_dir_b);
    command.SetModelKeyTagList("key_b");
    ASSERT_TRUE(command.Execute());
    EXPECT_TRUE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_B.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_A.csv"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir_b, ".csv"), 1);
}
