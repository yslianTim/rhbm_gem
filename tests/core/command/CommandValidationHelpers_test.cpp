#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string_view>
#include <vector>

#include "command/detail/CommandBase.hpp"
#include "support/CommandValidationAssertions.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandTypes.hpp>

namespace rg = rhbm_gem;

namespace {

struct ValidationHelperCommandOptions
{
    std::filesystem::path required_path;
    bool validate_required_path{ false };
    std::filesystem::path optional_path;
    bool validate_optional_path{ false };
    int normalized_count{ 1 };
    bool validate_count{ false };
    rg::PrinterType printer{ rg::PrinterType::GAUS_ESTIMATES };
    bool validate_printer{ false };
    double finite_positive_value{ 2.0 };
    bool validate_finite_positive{ false };
    double finite_non_negative_value{ 0.0 };
    bool validate_finite_non_negative{ false };
    int positive_count{ 1 };
    bool validate_positive_count{ false };
    double command_local_value{ 3.5 };
    bool validate_command_local_value{ false };
    bool add_problematic_parse_warnings{ false };
    bool add_prepare_error{ false };
};

class ValidationHelperCommand final : public rg::CommandBase<rg::CommandRequestBase>
{
public:
    void SetRequiredPath(const std::filesystem::path & value)
    {
        m_options.required_path = value;
        m_options.validate_required_path = true;
    }

    void SetOptionalPath(const std::filesystem::path & value)
    {
        m_options.optional_path = value;
        m_options.validate_optional_path = true;
    }

    void SetPositiveCount(int value)
    {
        m_options.normalized_count = value;
        m_options.validate_count = true;
    }

    int Count() const { return m_options.normalized_count; }

    void SetPrinter(rg::PrinterType value)
    {
        m_options.printer = value;
        m_options.validate_printer = true;
    }

    rg::PrinterType Printer() const { return m_options.printer; }

    void SetFinitePositiveValue(double value)
    {
        m_options.finite_positive_value = value;
        m_options.validate_finite_positive = true;
    }

    void SetFiniteNonNegativeValue(double value)
    {
        m_options.finite_non_negative_value = value;
        m_options.validate_finite_non_negative = true;
    }

    void SetPositiveCountValue(int value)
    {
        m_options.positive_count = value;
        m_options.validate_positive_count = true;
    }

    void SetCommandLocalValidatedValue(double value)
    {
        m_options.command_local_value = value;
        m_options.validate_command_local_value = true;
    }

    double CommandLocalValidatedValue() const { return m_options.command_local_value; }
    double FinitePositiveValue() const { return m_options.finite_positive_value; }
    double FiniteNonNegativeValue() const { return m_options.finite_non_negative_value; }
    int PositiveCountValue() const { return m_options.positive_count; }

    void SetProblematicValue(int value)
    {
        m_options.add_problematic_parse_warnings = value <= 0;
    }

    void SetPrepareError(bool value)
    {
        m_options.add_prepare_error = value;
    }

    void NormalizeAndValidateRequest(rg::CommandRequestBase &) override
    {
        if (m_options.validate_required_path)
        {
            ValidateRequiredPath(m_options.required_path, "--input", "Input file");
        }
        if (m_options.validate_optional_path)
        {
            ValidateOptionalPath(m_options.optional_path, "--optional", "Optional file");
        }
        if (m_options.validate_count)
        {
            CoercePositiveScalar(
                m_options.normalized_count,
                "--count",
                4,
                LogLevel::Warning,
                "Count");
        }
        if (m_options.validate_printer)
        {
            CoerceEnum(
                m_options.printer,
                "--printer",
                rg::PrinterType::GAUS_ESTIMATES,
                "Printer choice");
        }
        if (m_options.validate_finite_positive)
        {
            CoerceFinitePositiveScalar(
                m_options.finite_positive_value,
                "--finite-positive",
                2.0,
                LogLevel::Error,
                "Finite positive value");
        }
        if (m_options.validate_finite_non_negative)
        {
            CoerceFiniteNonNegativeScalar(
                m_options.finite_non_negative_value,
                "--finite-non-negative",
                0.0,
                LogLevel::Error,
                "Finite non-negative value");
        }
        if (m_options.validate_positive_count)
        {
            CoercePositiveScalar(
                m_options.positive_count,
                "--positive-count",
                1,
                LogLevel::Error,
                "Positive count");
        }
        if (m_options.validate_command_local_value && m_options.command_local_value != 3.5)
        {
            m_options.command_local_value = 1.25;
            AddParseError(
                "--validated",
                "Validated value must equal 3.5.");
        }
        if (!m_options.add_problematic_parse_warnings) return;
        AddParseNormalizationWarning("--problem", "normalized to 1");
        AddParseNormalizationWarning("--problem", "clamped to safe range");
    }

    void ValidatePreparedRequest(const rg::CommandRequestBase &) override
    {
        RequirePrepareCondition(
            !m_options.add_prepare_error,
            "--problem",
            "semantic validation failed");
    }

private:
    ValidationHelperCommandOptions m_options{};

    bool ExecuteImpl(const rg::CommandRequestBase &) override { return true; }
};

std::size_t CountDiagnostics(
    const std::vector<rg::CommandDiagnostic> & issues,
    std::string_view option_name)
{
    return static_cast<std::size_t>(
        std::count_if(
            issues.begin(),
            issues.end(),
            [option_name](const rg::CommandDiagnostic & issue)
            {
                return issue.option_name == option_name;
            }));
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
    ValidationHelperCommand required_command{};
    required_command.SetRequiredPath(existing_file);
    EXPECT_TRUE(required_command.ExecuteRequest(rg::CommandRequestBase{}).issues.empty());

    required_command.SetRequiredPath(missing_file);
    const auto required_result{ required_command.ExecuteRequest(rg::CommandRequestBase{}) };
    EXPECT_FALSE(required_result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(required_result.issues, "--input"),
        nullptr);

    ValidationHelperCommand optional_command{};
    optional_command.SetOptionalPath({});
    const auto optional_empty_result{ optional_command.ExecuteRequest(rg::CommandRequestBase{}) };
    EXPECT_TRUE(optional_empty_result.succeeded);
    EXPECT_EQ(
        command_test::FindValidationIssue(optional_empty_result.issues, "--optional"),
        nullptr);

    optional_command.SetOptionalPath(missing_file);
    const auto optional_result{ optional_command.ExecuteRequest(rg::CommandRequestBase{}) };
    EXPECT_FALSE(optional_result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(optional_result.issues, "--optional"),
        nullptr);
}

TEST(CommandValidationHelpersTest, NormalizedScalarHelperReportsAutoCorrectedWarning)
{
    ValidationHelperCommand command{};

    command.SetPositiveCount(0);

    testing::internal::CaptureStderr();
    const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    EXPECT_TRUE(result.succeeded);
    EXPECT_EQ(command.Count(), 4);
    ASSERT_EQ(result.issues.size(), 1u);
    EXPECT_EQ(result.issues.front().option_name, "--count");
    EXPECT_NE(error_output.find("[validation; auto-corrected] Option --count"), std::string::npos);
}

TEST(CommandValidationHelpersTest, ValidatedEnumHelperFallsBackAndReportsParseError)
{
    ValidationHelperCommand command{};

    command.SetPrinter(static_cast<rg::PrinterType>(999));

    const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.Printer(), rg::PrinterType::GAUS_ESTIMATES);
    ASSERT_EQ(result.issues.size(), 1u);
    EXPECT_EQ(result.issues.front().option_name, "--printer");
}

TEST(CommandValidationHelpersTest, CommandLocalValidationPatternRejectsInvalidInputWithParseError)
{
    ValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);

    const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_DOUBLE_EQ(command.CommandLocalValidatedValue(), 1.25);
    EXPECT_NE(
        command_test::FindValidationIssue(result.issues, "--validated"),
        nullptr);
}

TEST(CommandValidationHelpersTest, ExecuteRequestClearsPriorValidationIssues)
{
    ValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);
    ASSERT_NE(
        command_test::FindValidationIssue(
            command.ExecuteRequest(rg::CommandRequestBase{}).issues,
            "--validated"),
        nullptr);

    command.SetCommandLocalValidatedValue(3.5);
    const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };

    EXPECT_TRUE(result.succeeded);
    EXPECT_EQ(
        command_test::FindValidationIssue(result.issues, "--validated"),
        nullptr);
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

        const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };

        EXPECT_FALSE(result.succeeded);
        EXPECT_DOUBLE_EQ(command.FinitePositiveValue(), 2.0);
        EXPECT_NE(
            command_test::FindValidationIssue(result.issues, "--finite-positive"),
            nullptr);
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

        const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };

        EXPECT_FALSE(result.succeeded);
        EXPECT_DOUBLE_EQ(command.FiniteNonNegativeValue(), 0.0);
        EXPECT_NE(
            command_test::FindValidationIssue(result.issues, "--finite-non-negative"),
            nullptr);
    }
}

TEST(CommandValidationHelpersTest, PositiveScalarOptionRejectsNonPositiveIntegers)
{
    ValidationHelperCommand command{};

    command.SetPositiveCountValue(0);

    const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.PositiveCountValue(), 1);
    EXPECT_NE(
        command_test::FindValidationIssue(result.issues, "--positive-count"),
        nullptr);
}

TEST(CommandValidationHelpersTest, KeepsWarningsAndErrorsForSameOption)
{
    ValidationHelperCommand command{};
    command.SetProblematicValue(0);
    command.SetPrepareError(true);

    const auto result{ command.ExecuteRequest(rg::CommandRequestBase{}) };

    EXPECT_FALSE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 3u);
    EXPECT_EQ(CountDiagnostics(result.issues, "--problem"), 3u);
}

TEST(CommandValidationHelpersTest, BaseNormalizationWarningsAreProgrammaticallyVisible)
{
    ValidationHelperCommand command{};
    rg::CommandRequestBase request{};
    request.job_count = 0;
    request.verbosity = 99;

    testing::internal::CaptureStderr();
    const auto result{ command.ExecuteRequest(request) };
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    EXPECT_TRUE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 2u);
    EXPECT_NE(command_test::FindValidationIssue(result.issues, "--jobs"), nullptr);
    EXPECT_NE(command_test::FindValidationIssue(result.issues, "--verbose"), nullptr);
    EXPECT_NE(error_output.find("[validation; auto-corrected] Option --jobs"), std::string::npos);
    EXPECT_NE(error_output.find("[validation; auto-corrected] Option --verbose"), std::string::npos);
}
