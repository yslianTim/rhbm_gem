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
    TestCommand() = default;

    void SetForceInvalid(bool value)
    {
        m_options.force_invalid = value;
        InvalidatePreparedState();
    }

    void ConfigureFilesystemOptions(const std::filesystem::path & folder_path)
    {
        m_common_request.folder_path = folder_path;
        NormalizeCommonRequest(m_common_request);
    }

    void ValidateOptions() override
    {
        ClearPrepareIssues("--test");
        if (m_options.force_invalid)
        {
            AddValidationError("--test", "forced invalid config");
        }
    }

    void ResetRuntimeState() override {}

private:
    rg::CommonCommandRequest m_common_request{};
    TestCommandOptions m_options{};

    rg::CommonCommandRequest & MutableCommonRequest() override { return m_common_request; }
    const rg::CommonCommandRequest & CommonRequest() const override { return m_common_request; }

    bool ExecuteImpl() override
    {
        return true;
    }
};

} // namespace

TEST(CommandBaseTest, PrepareCreatesOutputFolder)
{
    command_test::ScopedTempDir temp_dir{"command_base_setters"};
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommand command{};
    command.ConfigureFilesystemOptions(folder_path);

    EXPECT_FALSE(std::filesystem::exists(folder_path));

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(std::filesystem::exists(folder_path));
}

TEST(CommandBaseTest, PrepareForExecutionReportsValidationIssues)
{
    TestCommand command{};
    command.SetForceInvalid(true);

    testing::internal::CaptureStderr();
    EXPECT_FALSE(command.PrepareForExecution());
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    ASSERT_FALSE(command.GetValidationIssues().empty());
    EXPECT_NE(error_output.find("Option --test: forced invalid config"), std::string::npos);
}

TEST(CommandBaseTest, ValidationFailureSkipsFilesystemPreflight)
{
    command_test::ScopedTempDir temp_dir{"command_base_prepare_validation_failure"};
    const auto folder_path{ temp_dir.path() / "out" };
    TestCommand command{};
    command.ConfigureFilesystemOptions(folder_path);
    command.SetForceInvalid(true);

    ASSERT_FALSE(command.PrepareForExecution());
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
