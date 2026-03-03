#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

class TestCommand final : public rg::CommandBase
{
public:
    struct Options : public rg::CommandOptions
    {
        bool force_invalid{ false };
    };

    bool Execute() override { return true; }
    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}
    const rg::CommandOptions & GetOptions() const override { return m_options; }
    rg::CommandOptions & GetOptions() override { return m_options; }

    void SetForceInvalid(bool value)
    {
        m_options.force_invalid = value;
    }

    void ValidateOptions() override
    {
        ClearValidationIssuesForOption("--test");
        if (m_options.force_invalid)
        {
            AddValidationError("--test", "forced invalid config");
        }
    }

private:
    Options m_options{};
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
