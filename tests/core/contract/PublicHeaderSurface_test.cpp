#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"
#include <rhbm_gem/core/command/CommandMetadata.hpp>

namespace {

std::string ReadText(const std::filesystem::path & path)
{
    std::ifstream input(path);
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

bool PublicSectionContainsSetMethod(const std::string & header_text)
{
    std::size_t cursor{ 0 };
    while (true)
    {
        const auto public_pos{ header_text.find("public:", cursor) };
        if (public_pos == std::string::npos)
        {
            return false;
        }

        const auto private_pos{ header_text.find("private:", public_pos + 1) };
        const auto protected_pos{ header_text.find("protected:", public_pos + 1) };
        const auto next_public_pos{ header_text.find("public:", public_pos + 1) };

        const auto next_access{ std::min({ private_pos, protected_pos, next_public_pos, header_text.size() }) };
        const auto public_block{ header_text.substr(public_pos, next_access - public_pos) };
        if (public_block.find("Set") != std::string::npos)
        {
            return true;
        }
        cursor = next_access;
    }
}

} // namespace

TEST(PublicHeaderSurfaceTest, CorePublicHeadersMatchApprovedSurface) {
    const std::vector<std::string> expected{
        "core/command/Application.hpp",
        "core/command/CommandApi.hpp",
        "core/command/CommandBase.hpp",
        "core/command/CommandMetadata.hpp",
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
        "core/command/internal/CommandIdEntries.generated.inc",
        "core/painter/AtomPainter.hpp",
        "core/painter/ComparisonPainter.hpp",
        "core/painter/DemoPainter.hpp",
        "core/painter/GausPainter.hpp",
        "core/painter/ModelPainter.hpp",
        "core/painter/PainterBase.hpp",
        "core/painter/PotentialPlotBuilder.hpp"};

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

TEST(PublicHeaderSurfaceTest, CommandHeadersExposeRequestApplyInsteadOfPublicSetters)
{
    const auto root{ command_test::ProjectRootPath() / "include" / "rhbm_gem" / "core" / "command" };
    const std::vector<std::string> command_headers{
        "PotentialAnalysisCommand.hpp",
        "PotentialDisplayCommand.hpp",
        "ResultDumpCommand.hpp",
        "MapSimulationCommand.hpp",
        "MapVisualizationCommand.hpp",
        "PositionEstimationCommand.hpp",
        "HRLModelTestCommand.hpp",
    };

    for (const auto & header : command_headers)
    {
        const auto text{ ReadText(root / header) };
        EXPECT_NE(text.find("ApplyRequest("), std::string::npos) << header;
        EXPECT_FALSE(PublicSectionContainsSetMethod(text)) << header;
    }
}
