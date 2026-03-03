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

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}
    const rg::CommandOptions & GetOptions() const override { return m_options; }
    rg::CommandOptions & GetOptions() override { return m_options; }
    rg::CommandId GetCommandId() const override { return rg::CommandId::PotentialAnalysis; }
    bool Prepared() const { return IsPreparedForExecution(); }

    void SetForceInvalid(bool value)
    {
        InvalidatePreparedState();
        m_options.force_invalid = value;
    }

    void ValidateOptions() override
    {
        ClearValidationIssues("--test", rg::ValidationPhase::Prepare);
        if (m_options.force_invalid)
        {
            AddValidationError("--test", "forced invalid config");
        }
    }

private:
    bool ExecuteImpl() override { return true; }
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

TEST(CommandBaseTest, BaseSettersInvalidatePreparedState)
{
    TestCommand command;

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(command.Prepared());

    command.SetThreadSize(2);
    EXPECT_FALSE(command.Prepared());

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(command.Prepared());

    command.SetFolderPath("output");
    EXPECT_FALSE(command.Prepared());
}
