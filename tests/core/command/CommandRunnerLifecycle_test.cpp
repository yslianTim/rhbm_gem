#include <gtest/gtest.h>

#include <filesystem>

#include "command/detail/CommandRunner.hpp"
#include "support/CommandTestHelpers.hpp"

using namespace rhbm_gem;
using namespace rhbm_gem::core;

namespace {

struct TestCommandOptions
{
    bool force_invalid{ false };
};

class TestCommandHarness
{
public:
    int execute_impl_count{ 0 };

    void SetForceInvalid(bool value)
    {
        m_options.force_invalid = value;
    }

    void ConfigureFilesystemOptions(const std::filesystem::path & folder_path)
    {
        m_configured_request.output_dir = folder_path;
    }

    CommandResult ExecuteConfiguredRequest()
    {
        return m_runner.Run(
            m_configured_request,
            [](auto &, auto &) {},
            [this](auto & runner, const auto &)
            {
                runner.RequirePrepareCondition(
                    !m_options.force_invalid,
                    "forced invalid config");
            },
            [this](const auto &)
            {
                ++execute_impl_count;
                return true;
            });
    }

private:
    CommandRunner<CommandRequestBase> m_runner{};
    TestCommandOptions m_options{};
    CommandRequestBase m_configured_request{};
};

} // namespace

TEST(CommandRunnerLifecycleTest, RunCreatesOutputFolder)
{
    command_test::ScopedTempDir temp_dir{ "command_base_setters" };
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommandHarness command{};
    command.ConfigureFilesystemOptions(folder_path);

    EXPECT_FALSE(std::filesystem::exists(folder_path));

    ASSERT_TRUE(command.ExecuteConfiguredRequest().succeeded);
    EXPECT_EQ(command.execute_impl_count, 1);
    EXPECT_TRUE(std::filesystem::exists(folder_path));
}

TEST(CommandRunnerLifecycleTest, RunReportsValidationIssues)
{
    TestCommandHarness command{};
    command.SetForceInvalid(true);

    const auto result{ command.ExecuteConfiguredRequest() };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.execute_impl_count, 0);
    ASSERT_FALSE(result.issues.empty());
}

TEST(CommandRunnerLifecycleTest, ValidationFailureSkipsFilesystemPreflight)
{
    command_test::ScopedTempDir temp_dir{ "command_base_prepare_validation_failure" };
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommandHarness command{};
    command.ConfigureFilesystemOptions(folder_path);
    command.SetForceInvalid(true);

    ASSERT_FALSE(command.ExecuteConfiguredRequest().succeeded);
    EXPECT_EQ(command.execute_impl_count, 0);
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
