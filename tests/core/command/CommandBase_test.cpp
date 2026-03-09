#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include <CLI/CLI.hpp>

#include <rhbm_gem/core/CommandBase.hpp>
#include "CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct TestCommandOptions : public rg::CommandOptions
{
    bool force_invalid{ false };
};

class TestCommand final
    : public rg::CommandWithOptions<
          TestCommandOptions,
          rg::CommandId::ModelTest,
          rg::CommonOption::Threading
              | rg::CommonOption::Verbose
              | rg::CommonOption::Database
              | rg::CommonOption::OutputFolder>
{
public:
    using Options = TestCommandOptions;

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}
    int validate_count{ 0 };
    int reset_count{ 0 };
    int execute_count{ 0 };

    void SetForceInvalid(bool value)
    {
        MutateOptions([&]() { m_options.force_invalid = value; });
    }

    void ValidateOptions() override
    {
        ++validate_count;
        ResetPrepareIssues("--test");
        if (m_options.force_invalid)
        {
            AddValidationError("--test", "forced invalid config");
        }
    }

    void ResetRuntimeState() override
    {
        ++reset_count;
    }

private:
    bool ExecuteImpl() override
    {
        ++execute_count;
        return true;
    }
};

} // namespace

TEST(CommandBaseTest, SettersDoNotCreateDirectoriesUntilPrepareForExecution)
{
    command_test::ScopedTempDir temp_dir{"command_base_setters"};
    const auto database_path{ temp_dir.path() / "db" / "database.sqlite" };
    const auto folder_path{ temp_dir.path() / "out" };

    TestCommand command;
    command.SetDatabasePath(database_path);
    command.SetFolderPath(folder_path);

    EXPECT_FALSE(std::filesystem::exists(database_path.parent_path()));
    EXPECT_FALSE(std::filesystem::exists(folder_path));

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(std::filesystem::exists(database_path.parent_path()));
    EXPECT_TRUE(std::filesystem::exists(folder_path));
}

TEST(CommandBaseTest, PrepareForExecutionReportsValidationIssues)
{
    TestCommand command;
    command.SetForceInvalid(true);

    testing::internal::CaptureStderr();
    EXPECT_FALSE(command.PrepareForExecution());
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    ASSERT_FALSE(command.GetValidationIssues().empty());
    EXPECT_NE(error_output.find("Option --test: forced invalid config"), std::string::npos);
}

TEST(CommandBaseTest, BaseSettersInvalidatePreparedState)
{
    TestCommand command;

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);

    command.SetThreadSize(2);
    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_count, 1);

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(command.validate_count, 3);
    EXPECT_EQ(command.reset_count, 3);

    command.SetFolderPath("output");
    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 4);
    EXPECT_EQ(command.reset_count, 4);
    EXPECT_EQ(command.execute_count, 2);
}

TEST(CommandBaseTest, ValidationFailureSkipsFilesystemPreflight)
{
    command_test::ScopedTempDir temp_dir{"command_base_prepare_validation_failure"};
    const auto database_path{ temp_dir.path() / "db" / "database.sqlite" };
    const auto folder_path{ temp_dir.path() / "out" };

    TestCommand command;
    command.SetDatabasePath(database_path);
    command.SetFolderPath(folder_path);
    command.SetForceInvalid(true);

    ASSERT_FALSE(command.PrepareForExecution());
    EXPECT_FALSE(std::filesystem::exists(database_path.parent_path()));
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
