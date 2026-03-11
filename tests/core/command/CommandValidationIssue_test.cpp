#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include <CLI/CLI.hpp>

#include <rhbm_gem/core/command/CommandBase.hpp>
#include "CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct ValidationIssueCommandOptions : public rg::CommandOptions
{
    bool add_prepare_error{ false };
};

class ValidationIssueCommand final
    : public rg::CommandWithOptions<
          ValidationIssueCommandOptions,
          rg::CommandId::ModelTest,
          rg::CommonOption::Threading
              | rg::CommonOption::Verbose
              | rg::CommonOption::OutputFolder>
{
public:
    using Options = ValidationIssueCommandOptions;

    explicit ValidationIssueCommand(const rg::DataIoServices & data_io_services) :
        rg::CommandWithOptions<
            ValidationIssueCommandOptions,
            rg::CommandId::ModelTest,
            rg::CommonOption::Threading
                | rg::CommonOption::Verbose
                | rg::CommonOption::OutputFolder>{ data_io_services }
    {
    }

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}

    void SetProblematicValue(int value)
    {
        MutateOptions([&]()
        {
            ResetParseIssues("--problem");
            if (value <= 0)
            {
                AddNormalizationWarning("--problem", "normalized to 1");
                AddNormalizationWarning("--problem", "clamped to safe range");
            }
        });
    }

    void SetPrepareError(bool value)
    {
        MutateOptions([&]() { m_options.add_prepare_error = value; });
    }

    void ValidateOptions() override
    {
        ResetPrepareIssues("--problem");
        if (m_options.add_prepare_error)
        {
            AddValidationError("--problem", "semantic validation failed");
        }
    }

private:
    bool ExecuteImpl() override { return true; }
};

} // namespace

TEST(CommandValidationIssueTest, KeepsParseAndPrepareIssuesForSameOption)
{
    const auto data_io_services{ command_test::BuildDataIoServices() };
    ValidationIssueCommand command{ data_io_services };
    command.SetProblematicValue(0);
    command.SetPrepareError(true);

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    ASSERT_EQ(issues.size(), 3u);

    const auto parse_warning_count{
        std::count_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--problem"
                    && issue.phase == rg::ValidationPhase::Parse
                    && issue.level == LogLevel::Warning
                    && issue.auto_corrected;
            })
    };
    const auto prepare_error_count{
        std::count_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--problem"
                    && issue.phase == rg::ValidationPhase::Prepare
                    && issue.level == LogLevel::Error;
            })
    };

    EXPECT_EQ(parse_warning_count, 2);
    EXPECT_EQ(prepare_error_count, 1);
}

TEST(CommandValidationIssueTest, BaseNormalizationWarningsAreProgrammaticallyVisible)
{
    const auto data_io_services{ command_test::BuildDataIoServices() };
    ValidationIssueCommand command{ data_io_services };
    command.SetThreadSize(0);
    command.SetVerboseLevel(99);

    const auto & issues{ command.GetValidationIssues() };
    ASSERT_EQ(issues.size(), 2u);

    const auto jobs_issue{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--jobs";
            })
    };
    ASSERT_NE(jobs_issue, issues.end());
    EXPECT_EQ(jobs_issue->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(jobs_issue->level, LogLevel::Warning);
    EXPECT_TRUE(jobs_issue->auto_corrected);

    const auto verbose_issue{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--verbose";
            })
    };
    ASSERT_NE(verbose_issue, issues.end());
    EXPECT_EQ(verbose_issue->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(verbose_issue->level, LogLevel::Warning);
    EXPECT_TRUE(verbose_issue->auto_corrected);
}

TEST(CommandValidationIssueTest, ValidationIssuesUseOnlyParseAndPreparePhases)
{
    const auto data_io_services{ command_test::BuildDataIoServices() };
    ValidationIssueCommand command{ data_io_services };
    command.SetProblematicValue(0);
    command.SetPrepareError(true);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_FALSE(command.GetValidationIssues().empty());

    for (const auto & issue : command.GetValidationIssues())
    {
        EXPECT_TRUE(
            issue.phase == rg::ValidationPhase::Parse
            || issue.phase == rg::ValidationPhase::Prepare);
    }
}
