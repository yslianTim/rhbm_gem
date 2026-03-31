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
    LifecycleCommand() = default;

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
        ClearPrepareIssues("--contract");
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
    rg::CommonCommandRequest m_common_request{};
    LifecycleCommandOptions m_options{};

    rg::CommonCommandRequest & MutableCommonRequest() override { return m_common_request; }
    const rg::CommonCommandRequest & CommonRequest() const override { return m_common_request; }

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

TEST(CommandExecutionContractTest, MutatingAssignedOptionAfterPrepareForcesFreshPrepareInExecute)
{
    LifecycleCommand command{};

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);

    command.SetExecutionToggle(true);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_impl_count, 1);
}
