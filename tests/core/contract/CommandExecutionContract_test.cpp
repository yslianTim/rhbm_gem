#include <gtest/gtest.h>

#include <vector>

#include "internal/command/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct LifecycleCommandOptions
{
    bool fail_prepare{ false };
    bool execution_toggle{ false };
};

class LifecycleCommand final : public rg::CommandBase
{
public:
    LifecycleCommand()
    {
        BindCommonRequest(m_common_request);
    }

    int validate_count{ 0 };
    int reset_count{ 0 };
    int execute_impl_count{ 0 };
    std::vector<int> runtime_state{};

    void SetFailPrepare(bool value)
    {
        m_options.fail_prepare = value;
        InvalidatePreparedState();
    }

    void SetExecutionToggle(bool value)
    {
        m_options.execution_toggle = value;
        InvalidatePreparedState();
    }

    void ValidateOptions() override
    {
        ++validate_count;
        RequireCondition(!m_options.fail_prepare, "--contract", "prepare failed");
    }

    void ResetRuntimeState() override
    {
        ++reset_count;
        runtime_state.clear();
    }

private:
    rg::CommonCommandRequest m_common_request{};
    LifecycleCommandOptions m_options{};

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

TEST(CommandExecutionContractTest, RunValidatesAndResetsBeforeExecuteImpl)
{
    LifecycleCommand command{};
    command.SetFailPrepare(true);

    EXPECT_FALSE(command.Run());
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);
    EXPECT_EQ(command.execute_impl_count, 0);
}

TEST(CommandExecutionContractTest, RunExecutesValidationResetAndExecuteOnce)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.Run());
    EXPECT_TRUE(command.WasPrepared());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);
    EXPECT_EQ(command.execute_impl_count, 1);
}

TEST(CommandExecutionContractTest, RepeatedRunResetsRuntimeStateBetweenRuns)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.Run());
    ASSERT_TRUE(command.Run());

    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_impl_count, 2);
}

TEST(CommandExecutionContractTest, MutatingAssignedOptionBeforeSecondRunStillUsesFreshPreparation)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);
    EXPECT_EQ(command.execute_impl_count, 1);

    command.SetExecutionToggle(true);

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_impl_count, 2);
}
