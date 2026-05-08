#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "command/detail/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"

using namespace rhbm_gem;

namespace {

struct TestCommandOptions
{
    bool force_invalid{ false };
};

class TestCommand final : public CommandBase<CommandRequestBase>
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

    void ValidatePreparedRequest(const CommandRequestBase &) override
    {
        RequirePrepareCondition(!m_options.force_invalid, "forced invalid config");
    }

private:
    TestCommandOptions m_options{};
    CommandRequestBase m_configured_request{};

    bool ExecuteImpl(const CommandRequestBase &) override
    {
        ++execute_impl_count;
        return true;
    }

public:
    CommandResult ExecuteConfiguredRequest()
    {
        return ExecuteRequest(m_configured_request);
    }
};

} // namespace

TEST(CommandBaseLifecycleTest, RunCreatesOutputFolder)
{
    command_test::ScopedTempDir temp_dir{ "command_base_setters" };
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommand command{};
    command.ConfigureFilesystemOptions(folder_path);

    EXPECT_FALSE(std::filesystem::exists(folder_path));

    ASSERT_TRUE(command.ExecuteConfiguredRequest().succeeded);
    EXPECT_EQ(command.execute_impl_count, 1);
    EXPECT_TRUE(std::filesystem::exists(folder_path));
}

TEST(CommandBaseLifecycleTest, RunReportsValidationIssues)
{
    TestCommand command{};
    command.SetForceInvalid(true);

    testing::internal::CaptureStderr();
    const auto result{ command.ExecuteConfiguredRequest() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    EXPECT_FALSE(result.succeeded);
    EXPECT_EQ(command.execute_impl_count, 0);
    ASSERT_FALSE(result.issues.empty());
    EXPECT_NE(error_output.find("Option request: forced invalid config"), std::string::npos);
}

TEST(CommandBaseLifecycleTest, ValidationFailureSkipsFilesystemPreflight)
{
    command_test::ScopedTempDir temp_dir{ "command_base_prepare_validation_failure" };
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommand command{};
    command.ConfigureFilesystemOptions(folder_path);
    command.SetForceInvalid(true);

    ASSERT_FALSE(command.ExecuteConfiguredRequest().succeeded);
    EXPECT_EQ(command.execute_impl_count, 0);
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
