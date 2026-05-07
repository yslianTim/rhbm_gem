#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "command/detail/CommandBase.hpp"
#include "support/CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

struct TestCommandOptions
{
    bool force_invalid{ false };
};

class TestCommand final : public rg::CommandBase<rg::CommandRequestBase>
{
public:
    int execute_impl_count{ 0 };

    void SetForceInvalid(bool value)
    {
        m_options.force_invalid = value;
    }

    void ConfigureFilesystemOptions(const std::filesystem::path & folder_path)
    {
        rg::CommandRequestBase request{};
        request.output_dir = folder_path;
        ApplyRequest(request);
    }

    void ValidatePreparedRequest() override
    {
        RequirePrepareCondition(!m_options.force_invalid, "--test", "forced invalid config");
    }

private:
    TestCommandOptions m_options{};

    bool ExecuteImpl() override
    {
        ++execute_impl_count;
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
    EXPECT_EQ(command.execute_impl_count, 1);
    EXPECT_TRUE(std::filesystem::exists(folder_path));
}

TEST(CommandBaseLifecycleTest, RunReportsValidationIssues)
{
    TestCommand command{};
    command.SetForceInvalid(true);

    testing::internal::CaptureStderr();
    EXPECT_FALSE(command.Run());
    const std::string error_output{ testing::internal::GetCapturedStderr() };

    EXPECT_EQ(command.execute_impl_count, 0);
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
    EXPECT_EQ(command.execute_impl_count, 0);
    EXPECT_FALSE(std::filesystem::exists(folder_path));
}
