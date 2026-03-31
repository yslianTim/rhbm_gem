#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>


#include "internal/command/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct CommandScalarValidationHelperCommandOptions
{
    double finite_positive_value{ 2.0 };
    double finite_non_negative_value{ 0.0 };
    int positive_count{ 1 };
    double command_local_value{ 1.25 };
};

class CommandScalarValidationHelperCommand final : public rg::CommandBase
{
public:
    CommandScalarValidationHelperCommand() = default;

    void SetFinitePositiveValue(double value)
    {
        m_options.finite_positive_value = value;
        NormalizeScalar(
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
        NormalizeScalar(
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
        NormalizeScalar(
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

private:
    rg::CommonCommandRequest m_common_request{};
    CommandScalarValidationHelperCommandOptions m_options{};

    rg::CommonCommandRequest & MutableCommonRequest() override { return m_common_request; }
    const rg::CommonCommandRequest & CommonRequest() const override { return m_common_request; }

    bool ExecuteImpl() override { return true; }
};

const rg::ValidationIssue * FindIssue(
    const CommandScalarValidationHelperCommand & command,
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

TEST(CommandScalarValidationHelperTest, CommandLocalValidationPatternRejectsInvalidInputWithParseError)
{
    CommandScalarValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);

    EXPECT_DOUBLE_EQ(command.CommandLocalValidatedValue(), 1.25);
    const auto * issue{ FindIssue(command, "--validated") };
    ASSERT_NE(issue, nullptr);
    EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue->level, LogLevel::Error);
    EXPECT_FALSE(issue->auto_corrected);
}

TEST(CommandScalarValidationHelperTest, CommandLocalValidationPatternAcceptsValidInputAndClearsPriorIssue)
{
    CommandScalarValidationHelperCommand command{};

    command.SetCommandLocalValidatedValue(2.0);
    ASSERT_NE(FindIssue(command, "--validated"), nullptr);

    command.SetCommandLocalValidatedValue(3.5);

    EXPECT_DOUBLE_EQ(command.CommandLocalValidatedValue(), 3.5);
    EXPECT_EQ(FindIssue(command, "--validated"), nullptr);
}

TEST(CommandScalarValidationHelperTest, FinitePositiveScalarOptionRejectsZeroNegativeNanAndInfinity)
{
    const double invalid_values[]{
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };

    for (double value : invalid_values)
    {
        CommandScalarValidationHelperCommand command{};
        command.SetFinitePositiveValue(value);

        EXPECT_DOUBLE_EQ(command.FinitePositiveValue(), 2.0);
        const auto * issue{ FindIssue(command, "--finite-positive") };
        ASSERT_NE(issue, nullptr);
        EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
        EXPECT_EQ(issue->level, LogLevel::Error);
        EXPECT_FALSE(issue->auto_corrected);
    }
}

TEST(CommandScalarValidationHelperTest, FiniteNonNegativeScalarOptionRejectsNegativeNanAndInfinity)
{
    const double invalid_values[]{
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };

    for (double value : invalid_values)
    {
        CommandScalarValidationHelperCommand command{};
        command.SetFiniteNonNegativeValue(value);

        EXPECT_DOUBLE_EQ(command.FiniteNonNegativeValue(), 0.0);
        const auto * issue{ FindIssue(command, "--finite-non-negative") };
        ASSERT_NE(issue, nullptr);
        EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
        EXPECT_EQ(issue->level, LogLevel::Error);
        EXPECT_FALSE(issue->auto_corrected);
    }
}

TEST(CommandScalarValidationHelperTest, PositiveScalarOptionRejectsNonPositiveIntegers)
{
    CommandScalarValidationHelperCommand command{};

    command.SetPositiveCountValue(0);

    EXPECT_EQ(command.PositiveCountValue(), 1);
    const auto * issue{ FindIssue(command, "--positive-count") };
    ASSERT_NE(issue, nullptr);
    EXPECT_EQ(issue->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue->level, LogLevel::Error);
    EXPECT_FALSE(issue->auto_corrected);
}
