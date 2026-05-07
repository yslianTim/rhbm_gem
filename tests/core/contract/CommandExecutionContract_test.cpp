#include <gtest/gtest.h>

#include "command/detail/CommandBase.hpp"
#include "support/CommandValidationAssertions.hpp"
#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandSystem.hpp>

namespace rg = rhbm_gem;

namespace {

struct LifecycleCommandOptions
{
    bool fail_prepare{ false };
};

struct UnsupportedRequest : rg::CommandRequestBase
{
};

class LifecycleCommand final : public rg::CommandBase<rg::CommandRequestBase>
{
public:
    int validate_count{ 0 };
    int execute_impl_count{ 0 };

    void SetFailPrepare(bool value)
    {
        m_options.fail_prepare = value;
    }

    void ValidatePreparedRequest() override
    {
        ++validate_count;
        RequirePrepareCondition(!m_options.fail_prepare, "--contract", "prepare failed");
    }

private:
    LifecycleCommandOptions m_options{};

    bool ExecuteImpl() override
    {
        ++execute_impl_count;
        return true;
    }
};

} // namespace

TEST(CommandExecutionContractTest, RunValidatesBeforeExecuteImpl)
{
    LifecycleCommand command{};
    command.SetFailPrepare(true);

    EXPECT_FALSE(command.Run());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.execute_impl_count, 0);
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--contract",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandExecutionContractTest, RunExecutesValidationAndExecuteOnce)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.execute_impl_count, 1);
}

TEST(CommandExecutionContractTest, RepeatedRunRecomputesPrepareIssues)
{
    LifecycleCommand command{};

    command.SetFailPrepare(true);
    ASSERT_FALSE(command.Run());
    EXPECT_NE(
        command_test::FindValidationIssue(
            command,
            "--contract",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);

    command.SetFailPrepare(false);
    ASSERT_TRUE(command.Run());

    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.execute_impl_count, 1);
    EXPECT_EQ(
        command_test::FindValidationIssue(
            command,
            "--contract",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(CommandExecutionContractTest, RepeatedRunExecutesEachTime)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.execute_impl_count, 1);

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.execute_impl_count, 2);
}

TEST(CommandExecutionContractTest, PublicRunEntryPointReportsPreparationFailureAndValidationIssues)
{
    const auto result{ rg::RunCommand(rg::MapSimulationRequest{}) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.issues,
            "--model"),
        nullptr);
}

TEST(CommandExecutionContractTest, PublicRunEntryPointReportsUnsupportedRequestType)
{
    const auto result{ rg::RunCommand(UnsupportedRequest{}) };

    EXPECT_FALSE(result.succeeded);
    ASSERT_EQ(result.issues.size(), 1);
    EXPECT_EQ(result.issues.front().option_name, "request_type");
    EXPECT_EQ(result.issues.front().message, "Unsupported command request type.");
}
