#include <gtest/gtest.h>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/command/MapSimulationCommand.hpp"

namespace rg = rhbm_gem;

TEST(MapSimulationCommandTest, ExecuteFiltersInvalidBlurringWidths)
{
    command_test::ScopedTempDir temp_dir{"map_simulation_valid"};

    rg::MapSimulationCommand command{rg::CommonOptionProfile::FileWorkflow};
    rg::MapSimulationRequest request{};
    request.common.folder_path = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = "1.0,-2.0,3.0";
    command.ApplyRequest(request);

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 2);
}

TEST(MapSimulationCommandTest, ExecuteFailsWhenAllBlurringWidthsAreInvalid)
{
    command_test::ScopedTempDir temp_dir{"map_simulation_invalid"};

    rg::MapSimulationCommand command{rg::CommonOptionProfile::FileWorkflow};
    rg::MapSimulationRequest request{};
    request.common.folder_path = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = "-1.0,0.0";
    command.ApplyRequest(request);

    EXPECT_FALSE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 0);
}
