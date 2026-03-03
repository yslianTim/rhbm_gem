#include <gtest/gtest.h>

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

namespace rg = rhbm_gem;

namespace {

class PreparedStateCommand final : public rg::CommandBase
{
public:
    struct Options : public rg::CommandOptions
    {
        int mode{ 0 };
    };

    int validate_count{ 0 };
    int reset_count{ 0 };
    int execute_count{ 0 };

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}
    const rg::CommandOptions & GetOptions() const override { return m_options; }
    rg::CommandOptions & GetOptions() override { return m_options; }
    rg::CommandId GetCommandId() const override { return rg::CommandId::ModelTest; }

    void SetMode(int value)
    {
        InvalidatePreparedState();
        m_options.mode = value;
    }

    bool Prepared() const { return IsPreparedForExecution(); }

    void ValidateOptions() override
    {
        ++validate_count;
    }

    void ResetRuntimeState() override
    {
        ++reset_count;
    }

private:
    bool ExecuteImpl() override
    {
        ++execute_count;
        return true;
    }

    Options m_options{};
};

} // namespace

TEST(PreparedStateInvalidationTest, BaseOptionMutationInvalidatesPreparedState)
{
    PreparedStateCommand command;

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(command.Prepared());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);

    command.SetThreadSize(3);
    EXPECT_FALSE(command.Prepared());

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_count, 1);
}

TEST(PreparedStateInvalidationTest, ConcreteOptionMutationInvalidatesPreparedState)
{
    PreparedStateCommand command;

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(command.Prepared());

    command.SetMode(1);
    EXPECT_FALSE(command.Prepared());

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_count, 1);
}
