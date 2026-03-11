#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandMetadata.hpp>

namespace {

std::vector<std::string> CollectPublicHeadersForDomain(const std::string& domain) {
    const auto include_root{command_test::ProjectRootPath() / "include" / "rhbm_gem"};
    const auto domain_root{include_root / domain};
    std::vector<std::string> headers;

    for (const auto& entry : std::filesystem::recursive_directory_iterator(domain_root)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        headers.push_back(std::filesystem::relative(entry.path(), include_root).generic_string());
    }
    std::sort(headers.begin(), headers.end());
    return headers;
}

} // namespace

TEST(PublicHeaderSurfaceTest, CorePublicHeadersMatchApprovedSurface) {
    const std::vector<std::string> expected{
        "core/command/Application.hpp",
        "core/command/CommandBase.hpp",
        "core/command/CommandMetadata.hpp",
        "core/command/CommandOptionBinding.hpp",
        "core/command/HRLModelTestCommand.hpp",
        "core/command/MapSampling.hpp",
        "core/command/MapSimulationCommand.hpp",
        "core/command/MapVisualizationCommand.hpp",
        "core/command/OptionEnumClass.hpp",
        "core/command/OptionEnumTraits.hpp",
        "core/command/PositionEstimationCommand.hpp",
        "core/command/PotentialAnalysisCommand.hpp",
        "core/command/PotentialDisplayCommand.hpp",
        "core/command/ResultDumpCommand.hpp",
        "core/painter/AtomPainter.hpp",
        "core/painter/ComparisonPainter.hpp",
        "core/painter/DemoPainter.hpp",
        "core/painter/GausPainter.hpp",
        "core/painter/ModelPainter.hpp",
        "core/painter/PainterBase.hpp",
        "core/painter/PotentialPlotBuilder.hpp"};

    EXPECT_EQ(CollectPublicHeadersForDomain("core"), expected);
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
