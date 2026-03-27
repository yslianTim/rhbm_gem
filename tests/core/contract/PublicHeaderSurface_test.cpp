#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/CommandMetadata.hpp>

TEST(PublicHeaderSurfaceTest, CorePublicHeadersMatchApprovedSurface) {
    const std::vector<std::string> expected{
        "core/command/CommandApi.hpp",
        "core/command/CommandList.def",
        "core/command/CommandMetadata.hpp",
        "core/command/OptionEnumClass.hpp",
        "core/command/OptionEnumTraits.hpp",
        "core/painter/AtomPainter.hpp",
        "core/painter/ComparisonPainter.hpp",
        "core/painter/DemoPainter.hpp",
        "core/painter/GausPainter.hpp",
        "core/painter/ModelPainter.hpp",
        "core/painter/PainterBase.hpp"};

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("core"), expected);
}

TEST(PublicHeaderSurfaceTest, CommandMetadataProfilesRemainStable) {
    constexpr auto file_workflow_mask{
        rhbm_gem::CommonOptionMaskForProfile(rhbm_gem::CommonOptionProfile::FileWorkflow)};
    constexpr auto database_workflow_mask{
        rhbm_gem::CommonOptionMaskForProfile(rhbm_gem::CommonOptionProfile::DatabaseWorkflow)};

    EXPECT_TRUE(rhbm_gem::HasCommonOption(file_workflow_mask, rhbm_gem::CommonOption::Threading));
    EXPECT_TRUE(rhbm_gem::HasCommonOption(file_workflow_mask, rhbm_gem::CommonOption::Verbose));
    EXPECT_TRUE(rhbm_gem::HasCommonOption(file_workflow_mask, rhbm_gem::CommonOption::OutputFolder));
    EXPECT_FALSE(rhbm_gem::HasCommonOption(file_workflow_mask, rhbm_gem::CommonOption::Database));

    EXPECT_TRUE(rhbm_gem::HasCommonOption(database_workflow_mask, rhbm_gem::CommonOption::Threading));
    EXPECT_TRUE(rhbm_gem::HasCommonOption(database_workflow_mask, rhbm_gem::CommonOption::Verbose));
    EXPECT_TRUE(rhbm_gem::HasCommonOption(database_workflow_mask, rhbm_gem::CommonOption::Database));
    EXPECT_TRUE(rhbm_gem::HasCommonOption(database_workflow_mask, rhbm_gem::CommonOption::OutputFolder));
}

TEST(PublicHeaderSurfaceTest, CommandMetadataHeaderProvidesSharedValidationAndDefaults) {
    const auto default_data_root{ rhbm_gem::GetDefaultDataRootPath() };
    const auto default_database_path{ rhbm_gem::GetDefaultDatabasePath() };

    rhbm_gem::ValidationIssue issue{
        "--example",
        rhbm_gem::ValidationPhase::Prepare,
        LogLevel::Warning,
        "example",
        true};

    EXPECT_FALSE(default_data_root.empty());
    EXPECT_EQ(default_database_path.filename(), "database.sqlite");
    EXPECT_EQ(issue.option_name, "--example");
    EXPECT_EQ(issue.phase, rhbm_gem::ValidationPhase::Prepare);
    EXPECT_EQ(issue.level, LogLevel::Warning);
    EXPECT_TRUE(issue.auto_corrected);
}

TEST(PublicHeaderSurfaceTest, CommandApiHeaderProvidesExecutionReport) {
    rhbm_gem::ExecutionReport report{};
    report.validation_issues.push_back(rhbm_gem::ValidationIssue{
        "--example",
        rhbm_gem::ValidationPhase::Prepare,
        LogLevel::Warning,
        "example",
        true});

    EXPECT_FALSE(report.prepared);
    EXPECT_FALSE(report.executed);
    ASSERT_EQ(report.validation_issues.size(), 1u);
    EXPECT_EQ(report.validation_issues.front().option_name, "--example");
    EXPECT_EQ(report.validation_issues.front().phase, rhbm_gem::ValidationPhase::Prepare);
}
