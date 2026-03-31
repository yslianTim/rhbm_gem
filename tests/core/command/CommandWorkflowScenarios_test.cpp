#include <gtest/gtest.h>

#include <filesystem>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "command/MapSimulationCommand.hpp"
#include "command/ResultDumpCommand.hpp"

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
#include "command/MapVisualizationCommand.hpp"
#include "command/PositionEstimationCommand.hpp"
#endif

namespace rg = rhbm_gem;

TEST(CommandWorkflowScenariosTest, MapSimulationGeneratesMapForEachValidBlurringWidth)
{
    command_test::ScopedTempDir temp_dir{ "map_simulation_valid" };

    rg::MapSimulationCommand command{};
    rg::MapSimulationRequest request{};
    request.output_dir = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = { 1.0, -2.0, 3.0 };
    command.ApplyRequest(request);

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 2);
}

TEST(CommandWorkflowScenariosTest, ResultDumpRerunRefreshesRuntimeStateAndUsesCurrentPaths)
{
    command_test::ScopedTempDir temp_dir{ "result_dump_context_paths" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto database_path_a{ temp_dir.path() / "db_a" / "database.sqlite" };
    const auto database_path_b{ temp_dir.path() / "db_b" / "database.sqlite" };
    const auto output_dir_a{ temp_dir.path() / "output_a" };
    const auto output_dir_b{ temp_dir.path() / "output_b" };

    std::filesystem::create_directories(database_path_a.parent_path());
    std::filesystem::create_directories(database_path_b.parent_path());
    command_test::SeedSavedModel(database_path_a, model_path, "shared_key", "MODEL_FROM_A");
    command_test::SeedSavedModel(database_path_b, model_path, "shared_key", "MODEL_FROM_B");

    rg::ResultDumpCommand command{};
    rg::ResultDumpRequest request{};
    request.printer_choice = rg::PrinterType::ATOM_POSITION;
    request.model_key_tag_list = { "shared_key" };

    request.database_path = database_path_a;
    request.output_dir = output_dir_a;
    command.ApplyRequest(request);
    ASSERT_TRUE(command.Run());
    EXPECT_TRUE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_A.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_B.csv"));

    request.database_path = database_path_b;
    request.output_dir = output_dir_b;
    command.ApplyRequest(request);
    ASSERT_TRUE(command.Run());
    EXPECT_TRUE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_B.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_A.csv"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir_b, ".csv"), 1);
}

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE

TEST(CommandWorkflowScenariosTest, MapVisualizationInvalidAtomIdFailsWithoutWritingOutput)
{
    command_test::ScopedTempDir temp_dir{ "map_visualization_invalid" };
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")
    };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::MapVisualizationCommand command{};
    rg::MapVisualizationRequest request{};
    request.output_dir = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 999;
    command.ApplyRequest(request);

    EXPECT_FALSE(command.Run());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 0);
}

#ifdef HAVE_ROOT
TEST(CommandWorkflowScenariosTest, MapVisualizationWritesPdfToConfiguredFolder)
{
    command_test::ScopedTempDir temp_dir{ "map_visualization_valid" };
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")
    };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::MapVisualizationCommand command{};
    rg::MapVisualizationRequest request{};
    request.output_dir = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 1;
    command.ApplyRequest(request);

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 1);
}
#endif

TEST(CommandWorkflowScenariosTest, PositionEstimationDoesNotRequireDatabaseConfiguration)
{
    command_test::ScopedTempDir temp_dir{ "position_estimation" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::PositionEstimationCommand command{};
    rg::PositionEstimationRequest request{};
    request.output_dir = output_dir;
    request.map_file_path = map_path;
    command.ApplyRequest(request);

    ASSERT_TRUE(command.Run());
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".cmm"), 1);
}

#endif
