#include <gtest/gtest.h>

#include <filesystem>

#include "CommandValidationAssertions.hpp"
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/ResultDumpCommand.hpp"

namespace rg = rhbm_gem;

TEST(ResultDumpCommandTest, MapPrinterWithoutMapFileReportsMapValidationError)
{
    rg::ResultDumpCommand command{rg::CommonOptionProfile::DatabaseWorkflow};
    rg::ResultDumpRequest request{};
    request.printer_choice = rg::PrinterType::MAP_VALUE;
    request.model_key_tag_list = "model";
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--map",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
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

    rg::ResultDumpCommand command{rg::CommonOptionProfile::DatabaseWorkflow};
    rg::ResultDumpRequest request{};
    request.common.database_path = database_path;
    request.printer_choice = rg::PrinterType::ATOM_POSITION;

    request.common.folder_path = output_dir_a;
    request.model_key_tag_list = "key_a";
    command.ApplyRequest(request);
    ASSERT_TRUE(command.Execute());
    EXPECT_TRUE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_A.csv"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir_a, ".csv"), 1);

    request.common.folder_path = output_dir_b;
    request.model_key_tag_list = "key_b";
    command.ApplyRequest(request);
    ASSERT_TRUE(command.Execute());
    EXPECT_TRUE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_B.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_A.csv"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir_b, ".csv"), 1);
}

TEST(ResultDumpCommandTest, ExecuteUsesCurrentDatabasePathAndOutputFolderOptions)
{
    command_test::ScopedTempDir temp_dir{"result_dump_context_paths"};
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto database_path_a{ temp_dir.path() / "db_a" / "database.sqlite" };
    const auto database_path_b{ temp_dir.path() / "db_b" / "database.sqlite" };
    const auto output_dir_a{ temp_dir.path() / "output_a" };
    const auto output_dir_b{ temp_dir.path() / "output_b" };

    std::filesystem::create_directories(database_path_a.parent_path());
    std::filesystem::create_directories(database_path_b.parent_path());
    command_test::SeedSavedModel(database_path_a, model_path, "shared_key", "MODEL_FROM_A");
    command_test::SeedSavedModel(database_path_b, model_path, "shared_key", "MODEL_FROM_B");

    rg::ResultDumpCommand command{rg::CommonOptionProfile::DatabaseWorkflow};
    rg::ResultDumpRequest request{};
    request.printer_choice = rg::PrinterType::ATOM_POSITION;
    request.model_key_tag_list = "shared_key";

    request.common.database_path = database_path_a;
    request.common.folder_path = output_dir_a;
    command.ApplyRequest(request);
    ASSERT_TRUE(command.Execute());
    EXPECT_TRUE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_A.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_B.csv"));

    request.common.database_path = database_path_b;
    request.common.folder_path = output_dir_b;
    command.ApplyRequest(request);
    ASSERT_TRUE(command.Execute());
    EXPECT_TRUE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_B.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_A.csv"));
}
