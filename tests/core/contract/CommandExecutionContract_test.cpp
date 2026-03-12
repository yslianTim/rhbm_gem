#include <gtest/gtest.h>

#include <vector>

#include <CLI/CLI.hpp>

#include <rhbm_gem/core/command/CommandBase.hpp>
#include "CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct LifecycleCommandOptions : public rg::CommandOptions
{
    bool fail_prepare{ false };
};

class LifecycleCommand final
    : public rg::CommandWithOptions<
          LifecycleCommandOptions,
          rg::CommandId::ModelTest,
          rg::CommonOption::Threading
              | rg::CommonOption::Verbose
              | rg::CommonOption::OutputFolder>
{
public:
    using Options = LifecycleCommandOptions;

    explicit LifecycleCommand() :
        rg::CommandWithOptions<
            LifecycleCommandOptions,
            rg::CommandId::ModelTest,
            rg::CommonOption::Threading
                | rg::CommonOption::Verbose
                | rg::CommonOption::OutputFolder>{}
    {
    }

    int validate_count{ 0 };
    int reset_count{ 0 };
    int execute_impl_count{ 0 };
    std::vector<int> runtime_state{};

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}

    void SetFailPrepare(bool value)
    {
        MutateOptions([&]() { m_options.fail_prepare = value; });
    }

    void ValidateOptions() override
    {
        ++validate_count;
        ResetPrepareIssues("--contract");
        if (m_options.fail_prepare)
        {
            AddValidationError("--contract", "prepare failed");
        }
    }

    void ResetRuntimeState() override
    {
        ++reset_count;
        runtime_state.clear();
    }

private:
    bool ExecuteImpl() override
    {
        ++execute_impl_count;
        if (!runtime_state.empty())
        {
            return false;
        }
        runtime_state.push_back(1);
        return true;
    }
};

} // namespace

TEST(CommandExecutionContractTest, ExecuteRunsPrepareBeforeExecuteImpl)
{
    LifecycleCommand command{};
    command.SetFailPrepare(true);

    EXPECT_FALSE(command.Execute());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);
    EXPECT_EQ(command.execute_impl_count, 0);
}

TEST(CommandExecutionContractTest, ExplicitPrepareSkipsDuplicatePreflightInsideExecute)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);
    EXPECT_EQ(command.execute_impl_count, 1);
}

TEST(CommandExecutionContractTest, RepeatedExecuteResetsRuntimeStateBetweenRuns)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.Execute());
    ASSERT_TRUE(command.Execute());

    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_impl_count, 2);
}
