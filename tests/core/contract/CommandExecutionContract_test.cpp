#include <gtest/gtest.h>

#include <algorithm>
#include <string_view>
#include <vector>

#include "command/detail/CommandBase.hpp"
#include <rhbm_gem/core/CommandSystem.hpp>

using namespace rhbm_gem;
using namespace rhbm_gem::core;

namespace {

struct LifecycleCommandOptions
{
    bool fail_prepare{ false };
};

struct UnsupportedRequest : CommandRequestBase
{
};

class LifecycleCommand final : public CommandBase<CommandRequestBase>
{
public:
    int validate_count{ 0 };
    int execute_impl_count{ 0 };

    void SetFailPrepare(bool value)
    {
        m_options.fail_prepare = value;
    }

    void ValidatePreparedRequest(const CommandRequestBase &) override
    {
        ++validate_count;
        RequirePrepareCondition(!m_options.fail_prepare, "prepare failed");
    }

private:
    LifecycleCommandOptions m_options{};

    bool ExecuteImpl(const CommandRequestBase &) override
    {
        ++execute_impl_count;
        return true;
    }
};

bool HasDiagnosticForOption(
    const std::vector<CommandDiagnostic> & issues,
    std::string_view option_name)
{
    return std::any_of(
        issues.begin(),
        issues.end(),
        [option_name](const CommandDiagnostic & issue)
        {
            return issue.option_name == option_name;
        });
}

} // namespace

TEST(CommandExecutionContractTest, RunValidatesBeforeExecuteImpl)
{
    LifecycleCommand command{};
    command.SetFailPrepare(true);

    const auto result{ command.ExecuteRequest(CommandRequestBase{}) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.execute_impl_count, 0);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "request"));
}

TEST(CommandExecutionContractTest, RunExecutesValidationAndExecuteOnce)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.ExecuteRequest(CommandRequestBase{}).succeeded);
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.execute_impl_count, 1);
}

TEST(CommandExecutionContractTest, RepeatedRunRecomputesPrepareIssues)
{
    LifecycleCommand command{};

    command.SetFailPrepare(true);
    const auto failed_result{ command.ExecuteRequest(CommandRequestBase{}) };
    ASSERT_FALSE(failed_result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(failed_result.issues, "request"));

    command.SetFailPrepare(false);
    const auto succeeded_result{ command.ExecuteRequest(CommandRequestBase{}) };
    ASSERT_TRUE(succeeded_result.succeeded);

    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.execute_impl_count, 1);
    EXPECT_FALSE(HasDiagnosticForOption(succeeded_result.issues, "request"));
}

TEST(CommandExecutionContractTest, RepeatedRunExecutesEachTime)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.ExecuteRequest(CommandRequestBase{}).succeeded);
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.execute_impl_count, 1);

    ASSERT_TRUE(command.ExecuteRequest(CommandRequestBase{}).succeeded);
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.execute_impl_count, 2);
}

TEST(CommandExecutionContractTest, PublicRunEntryPointReportsPreparationFailureAndValidationIssues)
{
    const auto result{ RunCommand(MapSimulationRequest{}) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "-a,--model"));
}

TEST(CommandExecutionContractTest, PublicRunEntryPointReportsUnsupportedRequestType)
{
    const auto result{ RunCommand(UnsupportedRequest{}) };

    EXPECT_FALSE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 1);
    EXPECT_EQ(result.issues.front().option_name, "request_type");
    EXPECT_EQ(result.issues.front().message, "Unsupported command request type.");
}
