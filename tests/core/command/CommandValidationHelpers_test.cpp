#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "command/detail/CommandRunner.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/CommandTypes.hpp>

using namespace rhbm_gem;
using namespace rhbm_gem::core;

namespace rhbm_gem::core {

struct ValidationHelperRequest : CommandRequestBase
{
    std::filesystem::path required_path;
    std::filesystem::path optional_path;
    int normalized_count{ 1 };
    PrinterType printer{ PrinterType::GAUS_ESTIMATES };
    double finite_positive_value{ 2.0 };
    double finite_non_negative_value{ 0.0 };
    int positive_count{ 1 };
    double command_local_value{ 3.5 };
    int problematic_value{ 1 };
    std::vector<std::string> required_list{ "entry" };
    double alpha_r{ 0.1 };
    int uncataloged_value{ 1 };
};

} // namespace rhbm_gem::core

namespace rhbm_gem::core::command_internal {

template <>
struct RequestFieldCatalog<ValidationHelperRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = ValidationHelperRequest;
        VisitFieldList(visitor,
            RequestField{ "required_path", "-i,--input", "Input file", &Self::required_path },
            RequestField{ "optional_path", "--optional", "Optional file", &Self::optional_path },
            RequestField{ "normalized_count", "-c,--count", "Count", &Self::normalized_count },
            RequestField{ "printer", "--printer", "Printer choice", &Self::printer },
            RequestField{ "finite_positive_value", "--finite-positive",
                "Finite positive value",
                &Self::finite_positive_value },
            RequestField{ "finite_non_negative_value", "--finite-non-negative",
                "Finite non-negative value",
                &Self::finite_non_negative_value },
            RequestField{ "positive_count", "--positive-count",
                "Positive count",
                &Self::positive_count },
            RequestField{ "command_local_value", "--validated",
                "Validated value",
                &Self::command_local_value },
            RequestField{ "problematic_value", "--problem",
                "Problematic value",
                &Self::problematic_value },
            RequestField{ "required_list", "--required-list",
                "Required list",
                &Self::required_list },
            RequestField{ "alpha_r", "--alpha-r", "Alpha R", &Self::alpha_r });
    }
};

} // namespace rhbm_gem::core::command_internal

namespace {

struct ValidationHelperCommandOptions
{
    bool validate_required_path{ false };
    bool validate_optional_path{ false };
    bool validate_count{ false };
    bool validate_printer{ false };
    bool validate_finite_positive{ false };
    bool validate_finite_non_negative{ false };
    bool validate_positive_count{ false };
    bool validate_command_local_value{ false };
    bool validate_required_list{ false };
    bool validate_alpha_r{ false };
    bool validate_uncataloged_value{ false };
    bool add_problematic_parse_warnings{ false };
    bool add_prepare_error{ false };
};

class ValidationHelperCommand
{
public:
    int validate_prepared_count{ 0 };
    int execute_impl_count{ 0 };

    void SetRequiredPath(const std::filesystem::path & value)
    {
        m_configured_request.required_path = value;
        m_options.validate_required_path = true;
    }

    void SetOptionalPath(const std::filesystem::path & value)
    {
        m_configured_request.optional_path = value;
        m_options.validate_optional_path = true;
    }

    void SetPositiveCount(int value)
    {
        m_configured_request.normalized_count = value;
        m_options.validate_count = true;
    }

    int Count() const { return m_last_request.normalized_count; }

    void SetPrinter(PrinterType value)
    {
        m_configured_request.printer = value;
        m_options.validate_printer = true;
    }

    PrinterType Printer() const { return m_last_request.printer; }

    void SetFinitePositiveValue(double value)
    {
        m_configured_request.finite_positive_value = value;
        m_options.validate_finite_positive = true;
    }

    void SetFiniteNonNegativeValue(double value)
    {
        m_configured_request.finite_non_negative_value = value;
        m_options.validate_finite_non_negative = true;
    }

    void SetPositiveCountValue(int value)
    {
        m_configured_request.positive_count = value;
        m_options.validate_positive_count = true;
    }

    void SetCommandLocalValidatedValue(double value)
    {
        m_configured_request.command_local_value = value;
        m_options.validate_command_local_value = true;
    }

    void SetRequiredList(std::vector<std::string> value)
    {
        m_configured_request.required_list = std::move(value);
        m_options.validate_required_list = true;
    }

    void SetAlphaR(double value)
    {
        m_configured_request.alpha_r = value;
        m_options.validate_alpha_r = true;
    }

    void SetUncatalogedValue(int value)
    {
        m_configured_request.uncataloged_value = value;
        m_options.validate_uncataloged_value = true;
    }

    double CommandLocalValidatedValue() const { return m_last_request.command_local_value; }
    double FinitePositiveValue() const { return m_last_request.finite_positive_value; }
    double FiniteNonNegativeValue() const { return m_last_request.finite_non_negative_value; }
    int PositiveCountValue() const { return m_last_request.positive_count; }

    void SetProblematicValue(int value)
    {
        m_configured_request.problematic_value = value;
        m_options.add_problematic_parse_warnings = value <= 0;
    }

    void SetPrepareError(bool value)
    {
        m_options.add_prepare_error = value;
    }

    CommandResult ExecuteConfiguredRequest()
    {
        return ExecuteRequest(m_configured_request);
    }

    CommandResult ExecuteRequest(const ValidationHelperRequest & request)
    {
        return m_runner.Run(
            request,
            [this](auto & runner, auto & prepared_request)
            {
                NormalizeAndValidateRequest(runner, prepared_request);
            },
            [this](auto & runner, const auto & prepared_request)
            {
                ValidatePreparedRequest(runner, prepared_request);
            },
            [this](const auto & prepared_request)
            {
                return ExecutePreparedRequest(prepared_request);
            });
    }

private:
    void NormalizeAndValidateRequest(
        CommandRunner<ValidationHelperRequest> & runner,
        ValidationHelperRequest & request)
    {
        if (m_options.validate_required_path)
        {
            runner.RequireExistingPath(request, &ValidationHelperRequest::required_path);
        }
        if (m_options.validate_optional_path)
        {
            runner.RequireOptionalExistingPath(request, &ValidationHelperRequest::optional_path);
        }
        if (m_options.validate_count)
        {
            runner.NormalizePositiveScalar(
                request,
                &ValidationHelperRequest::normalized_count,
                4);
        }
        if (m_options.validate_printer)
        {
            runner.RequireEnum(request, &ValidationHelperRequest::printer);
        }
        if (m_options.validate_finite_positive)
        {
            runner.RequireFinitePositiveScalar(
                request,
                &ValidationHelperRequest::finite_positive_value);
        }
        if (m_options.validate_finite_non_negative)
        {
            runner.RequireFiniteNonNegativeScalar(
                request,
                &ValidationHelperRequest::finite_non_negative_value);
        }
        if (m_options.validate_positive_count)
        {
            runner.RequirePositiveScalar(
                request,
                &ValidationHelperRequest::positive_count);
        }
        if (m_options.validate_required_list)
        {
            runner.RequireNonEmptyList(request, &ValidationHelperRequest::required_list);
        }
        if (m_options.validate_alpha_r)
        {
            runner.RequireFinitePositiveScalar(
                request,
                &ValidationHelperRequest::alpha_r);
        }
        if (m_options.validate_uncataloged_value)
        {
            runner.RequirePositiveScalar(
                request,
                &ValidationHelperRequest::uncataloged_value);
        }
        if (m_options.validate_command_local_value && request.command_local_value != 3.5)
        {
            request.command_local_value = 1.25;
            runner.AddFieldValidationError(
                &ValidationHelperRequest::command_local_value,
                "Validated value must equal 3.5.");
        }
        if (m_options.add_problematic_parse_warnings)
        {
            runner.AddFieldNormalizationWarning(
                &ValidationHelperRequest::problematic_value,
                "normalized to 1");
            runner.AddFieldNormalizationWarning(
                &ValidationHelperRequest::problematic_value,
                "clamped to safe range");
        }
        m_last_request = request;
    }

    void ValidatePreparedRequest(
        CommandRunner<ValidationHelperRequest> & runner,
        const ValidationHelperRequest &)
    {
        ++validate_prepared_count;
        runner.RequirePrepareCondition(
            !m_options.add_prepare_error,
            "semantic validation failed");
    }

    CommandRunner<ValidationHelperRequest> m_runner{};
    ValidationHelperCommandOptions m_options{};
    ValidationHelperRequest m_configured_request{};
    ValidationHelperRequest m_last_request{};

    bool ExecutePreparedRequest(const ValidationHelperRequest &)
    {
        ++execute_impl_count;
        return true;
    }
};

std::size_t CountDiagnostics(
    const std::vector<CommandDiagnostic> & issues,
    std::string_view option_name)
{
    return static_cast<std::size_t>(
        std::count_if(
            issues.begin(),
            issues.end(),
            [option_name](const CommandDiagnostic & issue)
            {
                return issue.option_name == option_name;
            }));
}

bool HasDiagnosticForOption(
    const std::vector<CommandDiagnostic> & issues,
    std::string_view option_name)
{
    return CountDiagnostics(issues, option_name) > 0;
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
    EXPECT_TRUE(required_command.ExecuteConfiguredRequest().issues.empty());

    required_command.SetRequiredPath(missing_file);
    const auto required_result{ required_command.ExecuteConfiguredRequest() };
    EXPECT_FALSE(required_result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(required_result.issues, "-i,--input"));
    ASSERT_FALSE(required_result.issues.empty());
    EXPECT_NE(required_result.issues.front().message.find("required_path"), std::string::npos);

    ValidationHelperCommand optional_command{};
    optional_command.SetOptionalPath({});
    const auto optional_empty_result{ optional_command.ExecuteConfiguredRequest() };
    EXPECT_TRUE(optional_empty_result.succeeded);
    EXPECT_FALSE(HasDiagnosticForOption(optional_empty_result.issues, "--optional"));

    optional_command.SetOptionalPath(missing_file);
    const auto optional_result{ optional_command.ExecuteConfiguredRequest() };
    EXPECT_FALSE(optional_result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(optional_result.issues, "--optional"));
}

TEST(CommandValidationHelpersTest, NormalizedScalarHelperReportsAutoCorrectedWarning)
{
    ValidationHelperCommand command{};

    command.SetPositiveCount(0);

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_TRUE(result.succeeded);
    EXPECT_EQ(command.Count(), 4);
    ASSERT_EQ(result.issues.size(), 1u);
    EXPECT_EQ(result.issues.front().option_name, "-c,--count");
    EXPECT_NE(result.issues.front().message.find("normalized_count"), std::string::npos);
}

TEST(CommandValidationHelpersTest, RequiredEnumHelperRejectsInvalidValueWithoutFallback)
{
    ValidationHelperCommand command{};

    command.SetPrinter(static_cast<PrinterType>(999));

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.Printer(), static_cast<PrinterType>(999));
    ASSERT_EQ(result.issues.size(), 1u);
    EXPECT_EQ(result.issues.front().option_name, "--printer");
}

TEST(CommandValidationHelpersTest, CommandLocalValidationPatternRejectsInvalidInputWithParseError)
{
    ValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    EXPECT_DOUBLE_EQ(command.CommandLocalValidatedValue(), 1.25);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "--validated"));
}

TEST(CommandValidationHelpersTest, ExecuteRequestClearsPriorValidationIssues)
{
    ValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);
    ASSERT_TRUE(
        HasDiagnosticForOption(
            command.ExecuteConfiguredRequest().issues,
            "--validated"));

    command.SetCommandLocalValidatedValue(3.5);
    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_TRUE(result.succeeded);
    EXPECT_FALSE(HasDiagnosticForOption(result.issues, "--validated"));
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

        const auto result{ command.ExecuteConfiguredRequest() };

        EXPECT_FALSE(result.succeeded);
        if (std::isnan(value))
        {
            EXPECT_TRUE(std::isnan(command.FinitePositiveValue()));
        }
        else
        {
            EXPECT_DOUBLE_EQ(command.FinitePositiveValue(), value);
        }
        EXPECT_TRUE(HasDiagnosticForOption(result.issues, "--finite-positive"));
        ASSERT_FALSE(result.issues.empty());
        EXPECT_EQ(result.issues.front().message.find("Using"), std::string::npos);
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

        const auto result{ command.ExecuteConfiguredRequest() };

        EXPECT_FALSE(result.succeeded);
        if (std::isnan(value))
        {
            EXPECT_TRUE(std::isnan(command.FiniteNonNegativeValue()));
        }
        else
        {
            EXPECT_DOUBLE_EQ(command.FiniteNonNegativeValue(), value);
        }
        EXPECT_TRUE(HasDiagnosticForOption(result.issues, "--finite-non-negative"));
        ASSERT_FALSE(result.issues.empty());
        EXPECT_EQ(result.issues.front().message.find("Using"), std::string::npos);
    }
}

TEST(CommandValidationHelpersTest, PositiveScalarOptionRejectsNonPositiveIntegers)
{
    ValidationHelperCommand command{};

    command.SetPositiveCountValue(0);

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.PositiveCountValue(), 0);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "--positive-count"));
}

TEST(CommandValidationHelpersTest, KeepsWarningsAndPrepareErrors)
{
    ValidationHelperCommand command{};
    command.SetProblematicValue(0);
    command.SetPrepareError(true);

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 3u);
    EXPECT_EQ(CountDiagnostics(result.issues, "--problem"), 2u);
    EXPECT_EQ(CountDiagnostics(result.issues, "request"), 1u);
}

TEST(CommandValidationHelpersTest, NonEmptyListUsesCatalogMetadata)
{
    ValidationHelperCommand command{};
    command.SetRequiredList({});

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 1u);
    EXPECT_EQ(result.issues.front().option_name, "--required-list");
    EXPECT_NE(result.issues.front().message.find("required_list"), std::string::npos);
}

TEST(CommandValidationHelpersTest, FieldMetadataUsesCliFlagsAsOptionName)
{
    ValidationHelperCommand command{};
    command.SetAlphaR(0.0);

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 1u);
    EXPECT_EQ(result.issues.front().option_name, "--alpha-r");
    EXPECT_EQ(result.issues.front().message.find("Using"), std::string::npos);
}

TEST(CommandValidationHelpersTest, FieldLevelErrorSkipsPreparedValidation)
{
    ValidationHelperCommand command{};
    command.SetPositiveCountValue(0);
    command.SetPrepareError(true);

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.validate_prepared_count, 0);
    EXPECT_EQ(command.execute_impl_count, 0);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "--positive-count"));
    EXPECT_FALSE(HasDiagnosticForOption(result.issues, "request"));
}

TEST(CommandValidationHelpersTest, MissingRequestFieldMetadataReportsLogicError)
{
    ValidationHelperCommand command{};
    command.SetUncatalogedValue(0);

    EXPECT_THROW(command.ExecuteConfiguredRequest(), std::logic_error);
}

TEST(CommandValidationHelpersTest, BaseNormalizationWarningsAreProgrammaticallyVisible)
{
    ValidationHelperCommand command{};
    ValidationHelperRequest request{};
    request.job_count = 0;
    request.verbosity = 99;

    const auto result{ command.ExecuteRequest(request) };

    EXPECT_TRUE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 2u);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "-j,--jobs"));
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "-v,--verbose"));
}
