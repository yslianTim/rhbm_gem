#include <gtest/gtest.h>

#include <algorithm>

#include "CommandTestHelpers.hpp"
#include "MapSimulationCommand.hpp"

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

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--potential-model";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}

TEST(MapSimulationCommandTest, InvalidPartialChargeChoiceBecomesValidationError)
{
    rg::MapSimulationCommand command;
    command.SetModelFilePath(command_test::TestDataPath("test_model.cif"));
    command.SetBlurringWidthList("1.0");
    command.SetPartialChargeChoice(static_cast<rg::PartialCharge>(999));

    EXPECT_FALSE(command.PrepareForExecution());

    const auto & issues{ command.GetValidationIssues() };
    const auto issue_iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [](const rg::ValidationIssue & issue)
            {
                return issue.option_name == "--charge";
            })
    };
    ASSERT_NE(issue_iter, issues.end());
    EXPECT_EQ(issue_iter->phase, rg::ValidationPhase::Parse);
    EXPECT_EQ(issue_iter->level, LogLevel::Error);
}
