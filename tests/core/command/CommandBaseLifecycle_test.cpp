#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "internal/command/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct TestCommandOptions
{
    bool force_invalid{ false };
};

class TestCommand final : public rg::CommandBase
{
public:
    TestCommand()
    {
        BindBaseRequest(m_base_request);
    }

    void SetForceInvalid(bool value)
    {
        m_options.force_invalid = value;
        InvalidatePreparedState();
    }

    void ConfigureFilesystemOptions(const std::filesystem::path & folder_path)
    {
        m_base_request.output_dir = folder_path;
        CoerceBaseRequest(m_base_request);
    }

    void ValidateOptions() override
    {
        RequireCondition(!m_options.force_invalid, "--test", "forced invalid config");
    }

    void ResetRuntimeState() override {}

private:
    rg::CommandRequestBase m_base_request{};
    TestCommandOptions m_options{};

    bool ExecuteImpl() override
    {
        return true;
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

    ASSERT_TRUE(command.Run());
    EXPECT_TRUE(command.WasPrepared());
    EXPECT_TRUE(std::filesystem::exists(folder_path));
}

TEST(CommandBaseLifecycleTest, RunReportsValidationIssues)
{
    TestCommand command{};
    command.SetForceInvalid(true);

    testing::internal::CaptureStderr();
    EXPECT_FALSE(command.Run());
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    EXPECT_FALSE(command.WasPrepared());
    ASSERT_FALSE(command.GetValidationIssues().empty());
    EXPECT_NE(error_output.find("Option --test: forced invalid config"), std::string::npos);
}

TEST(CommandBaseLifecycleTest, ValidationFailureSkipsFilesystemPreflight)
{
    command_test::ScopedTempDir temp_dir{ "command_base_prepare_validation_failure" };
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommand command{};
    command.ConfigureFilesystemOptions(folder_path);
    command.SetForceInvalid(true);

    ASSERT_FALSE(command.Run());
    EXPECT_FALSE(command.WasPrepared());
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
