#include <gtest/gtest.h>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/command/PotentialAnalysisCommand.hpp"

namespace rg = rhbm_gem;

TEST(PotentialAnalysisExperimentalFeatureGateTest, ExecuteRemainsReachableWhenFeatureIsEnabled)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_bond_gate" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto maps_dir{ temp_dir.path() / "maps" };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };

    std::filesystem::create_directories(maps_dir);
    const auto map_path{
        command_test::GenerateMapFile(maps_dir, model_path, "bond_gate_map", "1.0")
    };

    rg::PotentialAnalysisCommand command{rg::CommonOptionProfile::DatabaseWorkflow};
    rg::PotentialAnalysisRequest request{};
    request.common.database_path = database_path;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.saved_key_tag = "bond_gate_model";
    request.sampling_size = 25;
    command.ApplyRequest(request);

    ASSERT_TRUE(command.PrepareForExecution());
    EXPECT_TRUE(command.Execute());
}
