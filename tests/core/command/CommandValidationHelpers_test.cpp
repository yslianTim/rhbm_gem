#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

#include "internal/command/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandEnums.hpp>

namespace rg = rhbm_gem;

namespace {

struct ValidationHelperCommandOptions
{
    std::filesystem::path required_path;
    std::filesystem::path optional_path;
    int normalized_count{ 1 };
    rg::PrinterType printer{ rg::PrinterType::GAUS_ESTIMATES };
    double finite_positive_value{ 2.0 };
    double finite_non_negative_value{ 0.0 };
    int positive_count{ 1 };
    double command_local_value{ 1.25 };
    bool add_prepare_error{ false };
};

class ValidationHelperCommand final : public rg::CommandBase
{
public:
    ValidationHelperCommand()
    {
        BindBaseRequest(m_base_request);
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
        m_options.normalized_count = value;
        CoerceScalar(
            m_options.normalized_count,
            "--count",
            [](int candidate) { return candidate > 0; },
            4,
            LogLevel::Warning,
            "Count must be positive, reset to default 4");
    }

    int Count() const { return m_options.normalized_count; }

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

    void SetFinitePositiveValue(double value)
    {
        m_options.finite_positive_value = value;
        CoerceScalar(
            m_options.finite_positive_value,
            "--finite-positive",
            [](double candidate)
            {
                return std::isfinite(candidate) && candidate > 0.0;
            },
            2.0,
            LogLevel::Error,
            "Finite positive value must be a finite positive number.");
    }

    void SetFiniteNonNegativeValue(double value)
    {
        m_options.finite_non_negative_value = value;
        CoerceScalar(
            m_options.finite_non_negative_value,
            "--finite-non-negative",
            [](double candidate)
            {
                return std::isfinite(candidate) && candidate >= 0.0;
            },
            0.0,
            LogLevel::Error,
            "Finite non-negative value must be a finite non-negative number.");
    }

    void SetPositiveCountValue(int value)
    {
        m_options.positive_count = value;
        CoerceScalar(
            m_options.positive_count,
            "--positive-count",
            [](int candidate) { return candidate > 0; },
            1,
            LogLevel::Error,
            "Positive count must be positive.");
    }

    void SetCommandLocalValidatedValue(double value)
    {
        InvalidatePreparedState();
        ClearParseIssues("--validated");
        if (value == 3.5)
        {
            m_options.command_local_value = value;
            return;
        }

        m_options.command_local_value = 1.25;
        AddValidationError(
            "--validated",
            "Validated value must equal 3.5.",
            rg::ValidationPhase::Parse);
    }

    double CommandLocalValidatedValue() const { return m_options.command_local_value; }
    double FinitePositiveValue() const { return m_options.finite_positive_value; }
    double FiniteNonNegativeValue() const { return m_options.finite_non_negative_value; }
    int PositiveCountValue() const { return m_options.positive_count; }

    void SetProblematicValue(int value)
    {
        InvalidatePreparedState();
        ClearParseIssues("--problem");
        if (value <= 0)
        {
            AddNormalizationWarning("--problem", "normalized to 1");
            AddNormalizationWarning("--problem", "clamped to safe range");
        }
    }

    void SetPrepareError(bool value)
    {
        m_options.add_prepare_error = value;
        InvalidatePreparedState();
    }

    void SetCommonOptionsForTest(int thread_size, int verbose_level)
    {
        m_base_request.job_count = thread_size;
        m_base_request.verbosity = verbose_level;
        CoerceBaseRequest(m_base_request);
    }

    void ValidateOptions() override
    {
        ClearPrepareIssues("--problem");
        if (m_options.add_prepare_error)
        {
            AddValidationError("--problem", "semantic validation failed");
        }
    }

    void ResetRuntimeState() override {}

private:
    rg::CommandRequestBase m_base_request{};
    ValidationHelperCommandOptions m_options{};

    bool ExecuteImpl() override { return true; }
};

const rg::ValidationIssue * FindIssue(
    const ValidationHelperCommand & command,
    const char * option_name)
{
    const auto & issues{ command.GetValidationIssues() };
    const auto iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [option_name](const rg::ValidationIssue & issue)
            {
                return issue.option_name == option_name;
            })
    };
    return iter == issues.end() ? nullptr : &(*iter);
}

} // namespace

TEST(CommandValidationHelpersTest, PathHelpersValidateRequiredAndOptionalInputs)
{
    command_test::ScopedTempDir temp_dir{ "command_option_helper" };
    const auto existing_file{ temp_dir.path() / "fixture.txt" };
    {
        std::ofstream output(existing_file);
        output << "fixture";
    }
    const auto missing_file{ temp_dir.path() / "missing.txt" };
    ValidationHelperCommand command{};
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

TEST(CommandValidationHelpersTest, NormalizedScalarHelperReportsAutoCorrectedWarning)
{
    ValidationHelperCommand command{};

    command.SetPositiveCount(0);

    EXPECT_EQ(command.Count(), 4);
    ASSERT_EQ(command.GetValidationIssues().size(), 1u);
    EXPECT_EQ(command.GetValidationIssues().front().option_name, "--count");
    EXPECT_EQ(command.GetValidationIssues().front().phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(command.GetValidationIssues().front().level, LogLevel::Warning);
    EXPECT_TRUE(command.GetValidationIssues().front().auto_corrected);
}

TEST(CommandValidationHelpersTest, ValidatedEnumHelperFallsBackAndReportsParseError)
{
    ValidationHelperCommand command{};

    command.SetPrinter(static_cast<rg::PrinterType>(999));

    EXPECT_EQ(command.Printer(), rg::PrinterType::GAUS_ESTIMATES);
    ASSERT_EQ(command.GetValidationIssues().size(), 1u);
    EXPECT_EQ(command.GetValidationIssues().front().option_name, "--printer");
    EXPECT_EQ(command.GetValidationIssues().front().phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(command.GetValidationIssues().front().level, LogLevel::Error);
    EXPECT_FALSE(command.GetValidationIssues().front().auto_corrected);
}

TEST(CommandValidationHelpersTest, CommandLocalValidationPatternRejectsInvalidInputWithParseError)
{
    ValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);

    EXPECT_DOUBLE_EQ(command.CommandLocalValidatedValue(), 1.25);
    const auto * issue{ FindIssue(command, "--validated") };
    ASSERT_NE(issue, nullptr);
    EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue->level, LogLevel::Error);
    EXPECT_FALSE(issue->auto_corrected);
}

TEST(CommandValidationHelpersTest, CommandLocalValidationPatternAcceptsValidInputAndClearsPriorIssue)
{
    ValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);
    ASSERT_NE(FindIssue(command, "--validated"), nullptr);

    command.SetCommandLocalValidatedValue(3.5);

    EXPECT_DOUBLE_EQ(command.CommandLocalValidatedValue(), 3.5);
    EXPECT_EQ(FindIssue(command, "--validated"), nullptr);
}

TEST(CommandValidationHelpersTest, FinitePositiveScalarOptionRejectsZeroNegativeNanAndInfinity)
{
    const double invalid_values[]{
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };

    for (double value : invalid_values)
    {
        ValidationHelperCommand command{};
        command.SetFinitePositiveValue(value);

        EXPECT_DOUBLE_EQ(command.FinitePositiveValue(), 2.0);
        const auto * issue{ FindIssue(command, "--finite-positive") };
        ASSERT_NE(issue, nullptr);
        EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
        EXPECT_EQ(issue->level, LogLevel::Error);
        EXPECT_FALSE(issue->auto_corrected);
    }
}

TEST(CommandValidationHelpersTest, FiniteNonNegativeScalarOptionRejectsNegativeNanAndInfinity)
{
    const double invalid_values[]{
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };

    for (double value : invalid_values)
    {
        ValidationHelperCommand command{};
        command.SetFiniteNonNegativeValue(value);

        EXPECT_DOUBLE_EQ(command.FiniteNonNegativeValue(), 0.0);
        const auto * issue{ FindIssue(command, "--finite-non-negative") };
        ASSERT_NE(issue, nullptr);
        EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
        EXPECT_EQ(issue->level, LogLevel::Error);
        EXPECT_FALSE(issue->auto_corrected);
    }
}

TEST(CommandValidationHelpersTest, PositiveScalarOptionRejectsNonPositiveIntegers)
{
    ValidationHelperCommand command{};

    command.SetPositiveCountValue(0);

    EXPECT_EQ(command.PositiveCountValue(), 1);
    const auto * issue{ FindIssue(command, "--positive-count") };
    ASSERT_NE(issue, nullptr);
    EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue->level, LogLevel::Error);
    EXPECT_FALSE(issue->auto_corrected);
}

TEST(CommandValidationHelpersTest, KeepsParseAndPrepareIssuesForSameOption)
{
    ValidationHelperCommand command{};
    command.SetProblematicValue(0);
    command.SetPrepareError(true);

    EXPECT_FALSE(command.Run());
    EXPECT_FALSE(command.WasPrepared());

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

TEST(CommandValidationHelpersTest, BaseNormalizationWarningsAreProgrammaticallyVisible)
{
    ValidationHelperCommand command{};
    command.SetCommonOptionsForTest(0, 99);

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
