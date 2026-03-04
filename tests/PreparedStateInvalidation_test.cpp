#include <gtest/gtest.h>

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

namespace rg = rhbm_gem;

namespace {

struct PreparedStateCommandOptions : public rg::CommandOptions
{
    int mode{ 0 };
};

class PreparedStateCommand final
    : public rg::CommandWithOptions<PreparedStateCommandOptions, rg::CommandId::ModelTest>
{
public:
    using Options = PreparedStateCommandOptions;

    int validate_count{ 0 };
    int reset_count{ 0 };
    int execute_count{ 0 };

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}

    void SetMode(int value)
    {
        MutateOptions([&]() { m_options.mode = value; });
    }

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
};

} // namespace

TEST(PreparedStateInvalidationTest, BaseOptionMutationInvalidatesPreparedState)
{
    PreparedStateCommand command;

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);

    command.SetThreadSize(3);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_count, 1);
}

TEST(PreparedStateInvalidationTest, ConcreteOptionMutationInvalidatesPreparedState)
{
    PreparedStateCommand command;

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_EQ(command.validate_count, 1);
    EXPECT_EQ(command.reset_count, 1);

    command.SetMode(1);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command.validate_count, 2);
    EXPECT_EQ(command.reset_count, 2);
    EXPECT_EQ(command.execute_count, 1);
}
