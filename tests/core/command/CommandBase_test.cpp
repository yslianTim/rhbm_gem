#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include <rhbm_gem/core/command/CommandBase.hpp>
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
    explicit TestCommand(const rg::DataIoServices & data_io_services) :
        CommandWithOptions{ data_io_services }
    {
    }

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}

    void SetForceInvalid(bool value)
    {
        MutateOptions([&]() { m_options.force_invalid = value; });
    }

    void ValidateOptions() override
    {
        ResetPrepareIssues("--test");
        if (m_options.force_invalid)
        {
            AddValidationError("--test", "forced invalid config");
        }
    }

    void ResetRuntimeState() override {}

private:
    bool ExecuteImpl() override
    {
        return true;
    }
};

} // namespace

TEST(CommandBaseTest, SettersDoNotCreateDirectoriesUntilPrepareForExecution)
{
    command_test::ScopedTempDir temp_dir{"command_base_setters"};
    const auto database_path{ temp_dir.path() / "db" / "database.sqlite" };
    const auto folder_path{ temp_dir.path() / "out" };

    const auto data_io_services{ command_test::BuildDataIoServices() };
    TestCommand command{ data_io_services };
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
    const auto data_io_services{ command_test::BuildDataIoServices() };
    TestCommand command{ data_io_services };
    command.SetForceInvalid(true);

    testing::internal::CaptureStderr();
    EXPECT_FALSE(command.PrepareForExecution());
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    ASSERT_FALSE(command.GetValidationIssues().empty());
    EXPECT_NE(error_output.find("Option --test: forced invalid config"), std::string::npos);
}

TEST(CommandBaseTest, ValidationFailureSkipsFilesystemPreflight)
{
    command_test::ScopedTempDir temp_dir{"command_base_prepare_validation_failure"};
    const auto database_path{ temp_dir.path() / "db" / "database.sqlite" };
    const auto folder_path{ temp_dir.path() / "out" };

    const auto data_io_services{ command_test::BuildDataIoServices() };
    TestCommand command{ data_io_services };
    command.SetDatabasePath(database_path);
    command.SetFolderPath(folder_path);
    command.SetForceInvalid(true);

    ASSERT_FALSE(command.PrepareForExecution());
    EXPECT_FALSE(std::filesystem::exists(database_path.parent_path()));
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
