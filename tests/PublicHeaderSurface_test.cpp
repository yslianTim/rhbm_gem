#include <gtest/gtest.h>

#include <filesystem>

TEST(PublicHeaderSurfaceTest, ExperimentalBondWorkflowHeaderIsNotPublic)
{
    const auto project_root{
        std::filesystem::path(__FILE__).parent_path().parent_path()
    };
    const auto leaked_header{
        project_root / "include" / "core" / "PotentialAnalysisBondWorkflow.hpp"
    };

    EXPECT_FALSE(std::filesystem::exists(leaked_header)) << leaked_header.string();
}
