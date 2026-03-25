#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "command/internal/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct TestCommandOptions : public rg::CommandOptions
{
    bool force_invalid{ false };
};

class TestCommand final
    : public rg::CommandWithOptions<TestCommandOptions>
{
public:
    using Options = TestCommandOptions;
    explicit TestCommand() :
        CommandWithOptions{
            rg::CommonOption::Threading
            | rg::CommonOption::Verbose
            | rg::CommonOption::Database
            | rg::CommonOption::OutputFolder}
    {
    }

    void SetForceInvalid(bool value)
    {
        AssignOption(m_options.force_invalid, value);
    }

    void ConfigureFilesystemOptions(
        const std::filesystem::path & database_path,
        const std::filesystem::path & folder_path)
    {
        SetDatabasePath(database_path);
        SetFolderPath(folder_path);
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

TEST(CommandBaseTest, PrepareCreatesOutputFolderButLeavesDatabaseParentToDatabaseLayer)
{
    command_test::ScopedTempDir temp_dir{"command_base_setters"};
    const auto database_path{ temp_dir.path() / "db" / "database.sqlite" };
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommand command{};
    command.ConfigureFilesystemOptions(database_path, folder_path);

    EXPECT_FALSE(std::filesystem::exists(database_path.parent_path()));
    EXPECT_FALSE(std::filesystem::exists(folder_path));

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_FALSE(std::filesystem::exists(database_path.parent_path()));
    EXPECT_TRUE(std::filesystem::exists(folder_path));
}

TEST(CommandBaseTest, PrepareForExecutionReportsValidationIssues)
{
    TestCommand command{};
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
    TestCommand command{};
    command.ConfigureFilesystemOptions(database_path, folder_path);
    command.SetForceInvalid(true);

    ASSERT_FALSE(command.PrepareForExecution());
    EXPECT_FALSE(std::filesystem::exists(database_path.parent_path()));
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
