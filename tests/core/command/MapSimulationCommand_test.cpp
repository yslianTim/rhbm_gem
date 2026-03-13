#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/MapSimulationCommand.hpp"

namespace rg = rhbm_gem;

TEST(MapSimulationCommandTest, ExecuteFiltersInvalidBlurringWidths)
{
    command_test::ScopedTempDir temp_dir{"map_simulation_valid"};

    rg::MapSimulationCommand command{};
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

    rg::MapSimulationCommand command{};
    rg::MapSimulationRequest request{};
    request.common.folder_path = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = "-1.0,0.0";
    command.ApplyRequest(request);

    EXPECT_FALSE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 0);
}

TEST(MapSimulationCommandTest, InvalidPotentialModelBecomesValidationError)
{
    rg::MapSimulationCommand command{};
    rg::MapSimulationRequest request{};
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = "1.0";
    request.potential_model_choice = static_cast<rg::PotentialModel>(999);
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--potential-model",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(MapSimulationCommandTest, InvalidPartialChargeChoiceBecomesValidationError)
{
    rg::MapSimulationCommand command{};
    rg::MapSimulationRequest request{};
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = "1.0";
    request.partial_charge_choice = static_cast<rg::PartialCharge>(999);
    command.ApplyRequest(request);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--charge",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}
