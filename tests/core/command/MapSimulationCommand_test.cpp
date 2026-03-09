#include <gtest/gtest.h>

#include "CommandValidationAssertions.hpp"
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/core/MapSimulationCommand.hpp>

namespace rg = rhbm_gem;

TEST(MapSimulationCommandTest, ExecuteFiltersInvalidBlurringWidths)
{
    command_test::ScopedTempDir temp_dir{"map_simulation_valid"};

    rg::MapSimulationCommand command;
    command.SetFolderPath(temp_dir.path());
    command.SetMapFileName("sim_map");
    command.SetModelFilePath(command_test::TestDataPath("test_model.cif"));
    command.SetBlurringWidthList("1.0,-2.0,3.0");

    ASSERT_TRUE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 2);
}

TEST(MapSimulationCommandTest, ExecuteFailsWhenAllBlurringWidthsAreInvalid)
{
    command_test::ScopedTempDir temp_dir{"map_simulation_invalid"};

    rg::MapSimulationCommand command;
    command.SetFolderPath(temp_dir.path());
    command.SetMapFileName("sim_map");
    command.SetModelFilePath(command_test::TestDataPath("test_model.cif"));
    command.SetBlurringWidthList("-1.0,0.0");

    EXPECT_FALSE(command.Execute());
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 0);
}

TEST(MapSimulationCommandTest, InvalidPotentialModelBecomesValidationError)
{
    rg::MapSimulationCommand command;
    command.SetModelFilePath(command_test::TestDataPath("test_model.cif"));
    command.SetBlurringWidthList("1.0");
    command.SetPotentialModelChoice(static_cast<rg::PotentialModel>(999));

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
    rg::MapSimulationCommand command;
    command.SetModelFilePath(command_test::TestDataPath("test_model.cif"));
    command.SetBlurringWidthList("1.0");
    command.SetPartialChargeChoice(static_cast<rg::PartialCharge>(999));

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--charge",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}
