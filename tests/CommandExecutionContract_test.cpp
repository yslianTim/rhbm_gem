#include <gtest/gtest.h>

#include <vector>

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

namespace rg = rhbm_gem;

namespace {

class LifecycleCommand final : public rg::CommandBase
{
public:
    struct Options : public rg::CommandOptions
    {
        bool fail_prepare{ false };
        bool fail_execute{ false };
    };

    int validate_count{ 0 };
    int reset_count{ 0 };
    int execute_impl_count{ 0 };
    std::vector<int> runtime_state{};

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}
    const rg::CommandOptions & GetOptions() const override { return m_options; }
    rg::CommandOptions & GetOptions() override { return m_options; }

    void SetFailPrepare(bool value) { m_options.fail_prepare = value; }
    void SetFailExecute(bool value) { m_options.fail_execute = value; }

    void ValidateOptions() override
    {
        ++validate_count;
        ClearValidationIssues("--contract", rg::ValidationPhase::Prepare);
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
        return !m_options.fail_execute;
    }

    Options m_options{};
};

} // namespace

TEST(CommandExecutionContractTest, ExecuteRunsPrepareBeforeExecuteImpl)
{
    LifecycleCommand command;
    command.SetFailPrepare(true);

    EXPECT_FALSE(command.Execute());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);
    EXPECT_EQ(command.execute_impl_count, 0);
}

TEST(CommandExecutionContractTest, ExplicitPrepareSkipsDuplicatePreflightInsideExecute)
{
    LifecycleCommand command;

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
    LifecycleCommand command;

    ASSERT_TRUE(command.Execute());
    ASSERT_TRUE(command.Execute());

    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_impl_count, 2);
}
