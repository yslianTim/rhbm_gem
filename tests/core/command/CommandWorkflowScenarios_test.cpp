#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <string>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/command/CommandSystem.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>

namespace rg = rhbm_gem;

namespace {

std::size_t CountSubstring(const std::string & text, const std::string & target)
{
    std::size_t count{ 0 };
    std::size_t position{ 0 };
    while ((position = text.find(target, position)) != std::string::npos)
    {
        ++count;
        position += target.size();
    }
    return count;
}

} // namespace

TEST(CommandWorkflowScenariosTest, MapSimulationGeneratesMapForEachValidBlurringWidth)
{
    command_test::ScopedTempDir temp_dir{ "map_simulation_valid" };

    rg::MapSimulationRequest request{};
    request.output_dir = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = { 1.0, -2.0, 3.0 };

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto result{ rg::RunCommand(request) };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    const std::string standard_output{ testing::internal::GetCapturedStdout() };

    ASSERT_TRUE(result.succeeded);
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 2);
    EXPECT_EQ(CountSubstring(standard_output, "Building KD-Tree from"), 1);
    EXPECT_EQ(
        error_output.find("MapObject::SetMapValueArray - MapObject already has a map value array"),
        std::string::npos);
}

TEST(CommandWorkflowScenariosTest, MapSimulationEmptyModelUsesZeroOrigin)
{
    command_test::ScopedTempDir temp_dir{ "map_simulation_empty_model" };

    rg::MapSimulationRequest request{};
    request.output_dir = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model_no_atoms.cif");
    request.blurring_width_list = { 1.0 };

    ASSERT_TRUE(rg::RunCommand(request).succeeded);
    ASSERT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 1);

    std::filesystem::path map_path;
    for (const auto & entry : std::filesystem::directory_iterator(temp_dir.path()))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".map")
        {
            map_path = entry.path();
            break;
        }
    }
    ASSERT_FALSE(map_path.empty());

    auto loaded_map{ rg::ReadMap(map_path) };
    ASSERT_NE(loaded_map, nullptr);
    EXPECT_EQ(loaded_map->GetGridSize(), (std::array<int, 3>{ 1, 1, 1 }));
    EXPECT_EQ(loaded_map->GetOrigin(), (std::array<float, 3>{ 0.0f, 0.0f, 0.0f }));
    EXPECT_NE(loaded_map->GetMapValueArray(), nullptr);
}

TEST(CommandWorkflowScenariosTest, ResultDumpUsesCurrentRequestPaths)
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

    rg::ResultDumpRequest request{};
    request.printer_choice = rg::PrinterType::ATOM_POSITION;
    request.model_key_tag_list = { "shared_key" };

    request.database_path = database_path_a;
    request.output_dir = output_dir_a;
    ASSERT_TRUE(rg::RunCommand(request).succeeded);
    EXPECT_TRUE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_A.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_B.csv"));

    request.database_path = database_path_b;
    request.output_dir = output_dir_b;
    ASSERT_TRUE(rg::RunCommand(request).succeeded);
    EXPECT_TRUE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_B.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_A.csv"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir_b, ".csv"), 1);
}

TEST(CommandWorkflowScenariosTest, ResultDumpDatabaseLoadFailureRetainsCommandContext)
{
    command_test::ScopedTempDir temp_dir{ "result_dump_missing_model" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto database_path{ temp_dir.path() / "db" / "database.sqlite" };

    std::filesystem::create_directories(database_path.parent_path());
    command_test::SeedSavedModel(database_path, model_path, "present_key", "MODEL_PRESENT");

    rg::ResultDumpRequest request{};
    request.database_path = database_path;
    request.model_key_tag_list = { "missing_key" };
    request.printer_choice = rg::PrinterType::ATOM_POSITION;

    testing::internal::CaptureStderr();
    EXPECT_FALSE(rg::RunCommand(request).succeeded);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    EXPECT_NE(error_output.find("ResultDumpCommand::BuildDataObjectList"), std::string::npos);
    EXPECT_NE(error_output.find("Failed to load dump inputs"), std::string::npos);
    EXPECT_NE(error_output.find("missing_key"), std::string::npos);
}

TEST(CommandWorkflowScenariosTest, MapSimulationFileLoadFailureRetainsCommandContext)
{
    command_test::ScopedTempDir temp_dir{ "map_simulation_wrong_input_type" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto wrong_input_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    rg::MapSimulationRequest request{};
    request.output_dir = temp_dir.path() / "out";
    request.model_file_path = wrong_input_path;
    request.map_file_name = "sim_map";
    request.blurring_width_list = { 1.0 };

    testing::internal::CaptureStderr();
    EXPECT_FALSE(rg::RunCommand(request).succeeded);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    EXPECT_NE(error_output.find("MapSimulationCommand::BuildDataObject"), std::string::npos);
    EXPECT_NE(error_output.find("Failed to load model file"), std::string::npos);
    EXPECT_NE(error_output.find("ModelObject"), std::string::npos);
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

    rg::MapVisualizationRequest request{};
    request.output_dir = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 999;

    EXPECT_FALSE(rg::RunCommand(request).succeeded);
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

    rg::MapVisualizationRequest request{};
    request.output_dir = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 1;

    ASSERT_TRUE(rg::RunCommand(request).succeeded);
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

    rg::PositionEstimationRequest request{};
    request.output_dir = output_dir;
    request.map_file_path = map_path;

    ASSERT_TRUE(rg::RunCommand(request).succeeded);
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".cmm"), 1);
}

#endif
