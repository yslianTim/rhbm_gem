#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/CommandContract.hpp>
#include <rhbm_gem/core/command/CommandEnumClass.hpp>

namespace rg = rhbm_gem;

namespace {

template <typename EnumType>
std::set<int> BuildValueSetFromCLI()
{
    std::set<int> values;
    for (const auto & [token, value] : rg::BuildCommandEnumCliMap<EnumType>())
    {
        static_cast<void>(token);
        values.insert(static_cast<int>(value));
    }
    return values;
}

template <typename EnumType>
std::set<int> BuildValueSetFromBindings()
{
    std::set<int> values;
    const auto binding_entries{ rg::GetCommandEnumBindingEntries<EnumType>() };
    for (const auto & entry : binding_entries)
    {
        values.insert(static_cast<int>(entry.value));
    }
    return values;
}

} // namespace

TEST(PublicHeaderSurfaceTest, CorePublicHeadersMatchApprovedSurface) {
    const std::vector<std::string> expected{
        "core/command/CommandApi.hpp",
        "core/command/CommandContract.hpp",
        "core/command/CommandEnumClass.hpp",
        "core/command/CommandList.def",
        "core/painter/AtomPainter.hpp",
        "core/painter/ComparisonPainter.hpp",
        "core/painter/DemoPainter.hpp",
        "core/painter/GausPainter.hpp",
        "core/painter/ModelPainter.hpp",
        "core/painter/PainterBase.hpp"};

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("core"), expected);
}

TEST(PublicHeaderSurfaceTest, CommandContractProfilesRemainStable) {
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

TEST(PublicHeaderSurfaceTest, CommandContractHeaderProvidesSharedValidationAndDefaults) {
    const auto default_data_root{ rhbm_gem::GetDefaultDataRootPath() };
    const auto default_database_path{ rhbm_gem::GetDefaultDatabasePath() };
    const auto catalog{ rhbm_gem::GetCommandCatalog() };

    rhbm_gem::ValidationIssue issue{
        "--example",
        rhbm_gem::ValidationPhase::Prepare,
        LogLevel::Warning,
        "example",
        true};

    EXPECT_FALSE(default_data_root.empty());
    EXPECT_EQ(default_database_path.filename(), "database.sqlite");
    EXPECT_FALSE(catalog.empty());
    EXPECT_EQ(catalog.front().id, rhbm_gem::CommandId::PotentialAnalysis);
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

TEST(PublicHeaderSurfaceTest, PainterTypeEnumMappingsStayInSync)
{
    const auto cli_map{ rg::BuildCommandEnumCliMap<rg::PainterType>() };
    EXPECT_EQ(cli_map.at("gaus"), rg::PainterType::GAUS);
    EXPECT_EQ(cli_map.at("0"), rg::PainterType::GAUS);
    EXPECT_EQ(cli_map.at("atom"), rg::PainterType::ATOM);
    const auto binding_entries{ rg::GetCommandEnumBindingEntries<rg::PainterType>() };
    EXPECT_EQ(binding_entries[0].token, "GAUS");
    EXPECT_EQ(BuildValueSetFromCLI<rg::PainterType>(), BuildValueSetFromBindings<rg::PainterType>());
}

TEST(PublicHeaderSurfaceTest, PrinterTypeEnumMappingsStayInSync)
{
    const auto cli_map{ rg::BuildCommandEnumCliMap<rg::PrinterType>() };
    EXPECT_EQ(cli_map.at("atom_out"), rg::PrinterType::ATOM_OUTLIER);
    EXPECT_EQ(cli_map.at("3"), rg::PrinterType::ATOM_OUTLIER);
    EXPECT_EQ(
        BuildValueSetFromCLI<rg::PrinterType>(),
        BuildValueSetFromBindings<rg::PrinterType>());
}

TEST(PublicHeaderSurfaceTest, PotentialModelEnumMappingsStayInSync)
{
    const auto cli_map{ rg::BuildCommandEnumCliMap<rg::PotentialModel>() };
    EXPECT_EQ(cli_map.at("five"), rg::PotentialModel::FIVE_GAUS_CHARGE);
    EXPECT_EQ(cli_map.at("1"), rg::PotentialModel::FIVE_GAUS_CHARGE);
    EXPECT_EQ(
        BuildValueSetFromCLI<rg::PotentialModel>(),
        BuildValueSetFromBindings<rg::PotentialModel>());
}

TEST(PublicHeaderSurfaceTest, PartialChargeEnumMappingsStayInSync)
{
    const auto cli_map{ rg::BuildCommandEnumCliMap<rg::PartialCharge>() };
    EXPECT_EQ(cli_map.at("amber"), rg::PartialCharge::AMBER);
    EXPECT_EQ(cli_map.at("2"), rg::PartialCharge::AMBER);
    EXPECT_EQ(
        BuildValueSetFromCLI<rg::PartialCharge>(),
        BuildValueSetFromBindings<rg::PartialCharge>());
}

TEST(PublicHeaderSurfaceTest, TesterTypeEnumMappingsStayInSync)
{
    const auto cli_map{ rg::BuildCommandEnumCliMap<rg::TesterType>() };
    EXPECT_EQ(cli_map.at("benchmark"), rg::TesterType::BENCHMARK);
    EXPECT_EQ(cli_map.at("0"), rg::TesterType::BENCHMARK);
    EXPECT_EQ(BuildValueSetFromCLI<rg::TesterType>(), BuildValueSetFromBindings<rg::TesterType>());
}

TEST(PublicHeaderSurfaceTest, SupportedEnumValueHelperMatchesDeclaredOptions)
{
    EXPECT_TRUE(rg::IsSupportedCommandEnumValue(rg::PainterType::MODEL));
    EXPECT_TRUE(rg::IsSupportedCommandEnumValue(rg::PrinterType::GAUS_ESTIMATES));
    EXPECT_TRUE(rg::IsSupportedCommandEnumValue(rg::PotentialModel::SINGLE_GAUS_USER));
    EXPECT_TRUE(rg::IsSupportedCommandEnumValue(rg::PartialCharge::PARTIAL));
    EXPECT_TRUE(rg::IsSupportedCommandEnumValue(rg::TesterType::MODEL_ALPHA_MEMBER));
}
