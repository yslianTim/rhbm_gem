#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

class CommandOptionHelperCommand final : public rg::CommandBase
{
public:
    struct Options : public rg::CommandOptions
    {
        std::filesystem::path required_path;
        std::filesystem::path optional_path;
        int count{ 1 };
        int mode{ 0 };
    };

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}
    const rg::CommandOptions & GetOptions() const override { return m_options; }
    rg::CommandOptions & GetOptions() override { return m_options; }
    rg::CommandId GetCommandId() const override { return rg::CommandId::ModelTest; }

    void SetMode(int value)
    {
        MutateOptions([&]() { m_options.mode = value; });
    }

    void SetRequiredPath(const std::filesystem::path & value)
    {
        SetRequiredExistingPathOption(m_options.required_path, value, "--input", "Input file");
    }

    void SetOptionalPath(const std::filesystem::path & value)
    {
        SetOptionalExistingPathOption(m_options.optional_path, value, "--optional", "Optional file");
    }

    void SetPositiveCount(int value)
    {
        SetNormalizedScalarOption(
            m_options.count,
            value,
            "--count",
            [](int candidate) { return candidate > 0; },
            4,
            "Count must be positive, reset to default 4");
    }

    bool Prepared() const { return IsPreparedForExecution(); }
    int Count() const { return m_options.count; }

private:
    bool ExecuteImpl() override { return true; }

    Options m_options{};
};

} // namespace

TEST(CommandOptionHelperTest, MutateOptionsInvalidatesPreparedState)
{
    CommandOptionHelperCommand command;

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(command.Prepared());

    command.SetMode(3);
    EXPECT_FALSE(command.Prepared());
}

TEST(CommandOptionHelperTest, PathHelpersValidateRequiredAndOptionalInputs)
{
    command_test::ScopedTempDir temp_dir{ "command_option_helper" };
    const auto existing_file{ temp_dir.path() / "fixture.txt" };
    {
        std::ofstream output(existing_file);
        output << "fixture";
    }
    const auto missing_file{ temp_dir.path() / "missing.txt" };

    CommandOptionHelperCommand command;
    command.SetRequiredPath(existing_file);
    EXPECT_TRUE(command.GetValidationIssues().empty());

    command.SetRequiredPath(missing_file);
    const auto required_issue{
        std::find_if(
            command.GetValidationIssues().begin(),
            command.GetValidationIssues().end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--input"
                    && issue.phase == rg::ValidationPhase::Parse
                    && issue.level == LogLevel::Error;
            })
    };
    EXPECT_NE(required_issue, command.GetValidationIssues().end());

    command.SetOptionalPath({});
    const auto optional_empty_issue{
        std::find_if(
            command.GetValidationIssues().begin(),
            command.GetValidationIssues().end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--optional";
            })
    };
    EXPECT_EQ(optional_empty_issue, command.GetValidationIssues().end());

    command.SetOptionalPath(missing_file);
    const auto optional_issue{
        std::find_if(
            command.GetValidationIssues().begin(),
            command.GetValidationIssues().end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--optional"
                    && issue.phase == rg::ValidationPhase::Parse
                    && issue.level == LogLevel::Error;
            })
    };
    EXPECT_NE(optional_issue, command.GetValidationIssues().end());
}

TEST(CommandOptionHelperTest, NormalizedScalarHelperReportsAutoCorrectedWarning)
{
    CommandOptionHelperCommand command;

    command.SetPositiveCount(0);

    EXPECT_EQ(command.Count(), 4);
    ASSERT_EQ(command.GetValidationIssues().size(), 1u);
    EXPECT_EQ(command.GetValidationIssues().front().option_name, "--count");
    EXPECT_EQ(command.GetValidationIssues().front().phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(command.GetValidationIssues().front().level, LogLevel::Warning);
    EXPECT_TRUE(command.GetValidationIssues().front().auto_corrected);
}
