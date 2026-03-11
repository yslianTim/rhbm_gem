#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <string>

#include "CommandTestHelpers.hpp"

namespace {

int RunScaffoldDryRun(const std::filesystem::path & script_path, const char * python_bin)
{
    const std::string command{
        std::string(python_bin)
        + " \"" + script_path.string() + "\""
        + " --name ScaffoldSyncProbe --profile FileWorkflow --dry-run --wire --strict"
    };
    return std::system(command.c_str());
}

} // namespace

TEST(CommandScaffoldScriptTest, StrictDryRunWiringPassesForCurrentLayout)
{
    const auto project_root{ command_test::ProjectRootPath() };
    const auto script_path{ project_root / "scripts" / "command_scaffold.py" };
    ASSERT_TRUE(std::filesystem::exists(script_path));

    int status{ RunScaffoldDryRun(script_path, "python3") };
    if (status != 0)
    {
        status = RunScaffoldDryRun(script_path, "python");
    }

    EXPECT_EQ(status, 0);
}
