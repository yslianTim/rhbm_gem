#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>

#include "internal/command/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandEnumClass.hpp>

namespace rg = rhbm_gem;

namespace {

struct CommandOptionHelperCommandOptions
{
    std::filesystem::path required_path;
    std::filesystem::path optional_path;
    int count{ 1 };
    rg::PrinterType printer{ rg::PrinterType::GAUS_ESTIMATES };
};

class CommandOptionHelperCommand final : public rg::CommandBase
{
public:
    CommandOptionHelperCommand()
    {
        BindCommonRequest(m_common_request);
    }

    void SetRequiredPath(const std::filesystem::path & value)
    {
        m_options.required_path = value;
        ValidateRequiredPath(m_options.required_path, "--input", "Input file");
    }

    void SetOptionalPath(const std::filesystem::path & value)
    {
        m_options.optional_path = value;
        ValidateOptionalPath(m_options.optional_path, "--optional", "Optional file");
    }

    void SetPositiveCount(int value)
    {
        m_options.count = value;
        CoerceScalar(
            m_options.count,
            "--count",
            [](int candidate) { return candidate > 0; },
            4,
            LogLevel::Warning,
            "Count must be positive, reset to default 4");
    }

    int Count() const { return m_options.count; }

    void SetPrinter(rg::PrinterType value)
    {
        m_options.printer = value;
        CoerceEnum(
            m_options.printer,
            "--printer",
            rg::PrinterType::GAUS_ESTIMATES,
            "Printer choice");
    }

    rg::PrinterType Printer() const { return m_options.printer; }

    void ValidateOptions() override {}

    void ResetRuntimeState() override {}

private:
    rg::CommonCommandRequest m_common_request{};
    CommandOptionHelperCommandOptions m_options{};

    bool ExecuteImpl() override
    {
        return true;
    }
};

} // namespace

TEST(CommandOptionHelperTest, PathHelpersValidateRequiredAndOptionalInputs)
{
    command_test::ScopedTempDir temp_dir{ "command_option_helper" };
    const auto existing_file{ temp_dir.path() / "fixture.txt" };
    {
        std::ofstream output(existing_file);
        output << "fixture";
    }
    const auto missing_file{ temp_dir.path() / "missing.txt" };
    CommandOptionHelperCommand command{};
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
    CommandOptionHelperCommand command{};

    command.SetPositiveCount(0);

    EXPECT_EQ(command.Count(), 4);
    ASSERT_EQ(command.GetValidationIssues().size(), 1u);
    EXPECT_EQ(command.GetValidationIssues().front().option_name, "--count");
    EXPECT_EQ(command.GetValidationIssues().front().phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(command.GetValidationIssues().front().level, LogLevel::Warning);
    EXPECT_TRUE(command.GetValidationIssues().front().auto_corrected);
}

TEST(CommandOptionHelperTest, ValidatedEnumHelperFallsBackAndReportsParseError)
{
    CommandOptionHelperCommand command{};

    command.SetPrinter(static_cast<rg::PrinterType>(999));

    EXPECT_EQ(command.Printer(), rg::PrinterType::GAUS_ESTIMATES);
    ASSERT_EQ(command.GetValidationIssues().size(), 1u);
    EXPECT_EQ(command.GetValidationIssues().front().option_name, "--printer");
    EXPECT_EQ(command.GetValidationIssues().front().phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(command.GetValidationIssues().front().level, LogLevel::Error);
    EXPECT_FALSE(command.GetValidationIssues().front().auto_corrected);
}
