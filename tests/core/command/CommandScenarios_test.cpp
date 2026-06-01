#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>

using namespace rhbm_gem;
using namespace rhbm_gem::core;

namespace {

bool HasDiagnosticForOption(
    const std::vector<CommandDiagnostic> & issues,
    std::string_view option_name)
{
    return std::any_of(
        issues.begin(),
        issues.end(),
        [option_name](const CommandDiagnostic & issue)
        {
            return issue.option_name == option_name;
        });
}

MapObject MakeConstantMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        values[i] = 1.0f;
    }
    return MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

std::filesystem::path WriteConstantMapFixture(const std::filesystem::path & output_dir)
{
    std::filesystem::create_directories(output_dir);
    const auto map_path{ output_dir / "constant.map" };
    auto map_object{ MakeConstantMapObject() };
    WriteMap(map_path, map_object);
    return map_path;
}

std::filesystem::path WriteCarbonHydrogenModelFixture(const std::filesystem::path & output_dir)
{
    std::filesystem::create_directories(output_dir);
    const auto model_path{ output_dir / "carbon_hydrogen.cif" };
    std::ofstream outfile{ model_path };
    outfile
        << "data_carbon_hydrogen\n"
        << "#\n"
        << "loop_\n"
        << "_atom_type.symbol\n"
        << "C\n"
        << "H\n"
        << "#\n"
        << "loop_\n"
        << "_atom_site.group_PDB\n"
        << "_atom_site.id\n"
        << "_atom_site.type_symbol\n"
        << "_atom_site.label_atom_id\n"
        << "_atom_site.label_alt_id\n"
        << "_atom_site.label_comp_id\n"
        << "_atom_site.label_asym_id\n"
        << "_atom_site.label_seq_id\n"
        << "_atom_site.Cartn_x\n"
        << "_atom_site.Cartn_y\n"
        << "_atom_site.Cartn_z\n"
        << "_atom_site.occupancy\n"
        << "_atom_site.B_iso_or_equiv\n"
        << "_atom_site.pdbx_PDB_model_num\n"
        << "ATOM 1 C CA . ALA A 1 0.0 0.0 0.0 1.0 0.0 1\n"
        << "ATOM 2 H H . ALA A 1 1.0 0.0 0.0 1.0 0.0 1\n"
        << "#\n";
    return model_path;
}

PotentialAnalysisRequest MakeNormalizationScenarioRequest(
    const std::filesystem::path & temp_dir,
    const std::filesystem::path & map_path)
{
    PotentialAnalysisRequest request{};
    request.database_path = temp_dir / "analysis.sqlite";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.map_file_path = map_path;
    return request;
}

std::string CapturePotentialAnalysisStderr(const PotentialAnalysisRequest & request)
{
    testing::internal::CaptureStderr();
    static_cast<void>(RunCommand(request));
    return testing::internal::GetCapturedStderr();
}

std::string CapturePotentialAnalysisCliStderr(std::vector<std::string> args)
{
    std::vector<char *> argv;
    argv.reserve(args.size());
    for (auto & arg : args)
    {
        argv.push_back(arg.data());
    }

    testing::internal::CaptureStderr();
    static_cast<void>(RunCommandCLI(static_cast<int>(argv.size()), argv.data()));
    return testing::internal::GetCapturedStderr();
}

constexpr const char * kMapNormalizationWarning{ "MapObject::MapValueArrayNormalization" };

} // namespace

TEST(CommandScenariosTest, MapSimulationRejectsAllInvalidBlurringWidthsAtPrepare)
{
    MapSimulationRequest request{};
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = { -1.0, 0.0 };

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "request"));
}

TEST(CommandScenariosTest, PotentialAnalysisRequiresPositiveResolutionWhenSimulationEnabled)
{
    PotentialAnalysisRequest request{};
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.map_file_path = command_test::TestDataPath("test_model.cif");
    request.simulation_flag = true;
    request.simulated_map_resolution = 0.0;

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "request"));
}

TEST(CommandScenariosTest, PotentialAnalysisDefaultsMapNormalizationOn)
{
    PotentialAnalysisRequest request{};

    EXPECT_TRUE(request.map_normalization_flag);
}

TEST(CommandScenariosTest, PotentialAnalysisDefaultsToIncludingHydrogen)
{
    PotentialAnalysisRequest request{};

    EXPECT_FALSE(request.exclude_hydrogen);
}

TEST(CommandScenariosTest, PotentialAnalysisDefaultsToFibonacciSamplingMethod)
{
    PotentialAnalysisRequest request{};

    EXPECT_EQ(
        SphereSamplingMethod::FibonacciDeterministic,
        request.sampling_method);
}

TEST(CommandScenariosTest, PotentialAnalysisRejectsInvalidSamplingMethodAtParse)
{
    PotentialAnalysisRequest request{};
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.map_file_path = command_test::TestDataPath("test_model.cif");
    request.sampling_method = static_cast<SphereSamplingMethod>(99);

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "--sampling-method"));
}

TEST(CommandScenariosTest, PotentialAnalysisRejectsEmptySavedKeyAtParse)
{
    PotentialAnalysisRequest request{};
    request.saved_key_tag = "";

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "-k,--save-key"));
}

TEST(CommandScenariosTest, PotentialAnalysisNormalizesMapByDefault)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_normalization_default" };
    const auto map_path{ WriteConstantMapFixture(temp_dir.path() / "maps") };
    auto request{ MakeNormalizationScenarioRequest(temp_dir.path(), map_path) };

    const auto error_output{ CapturePotentialAnalysisStderr(request) };

    EXPECT_NE(error_output.find(kMapNormalizationWarning), std::string::npos);
}

TEST(CommandScenariosTest, PotentialAnalysisSkipsMapNormalizationForSimulationDefault)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_normalization_simulation" };
    const auto map_path{ WriteConstantMapFixture(temp_dir.path() / "maps") };
    auto request{ MakeNormalizationScenarioRequest(temp_dir.path(), map_path) };
    request.simulation_flag = true;
    request.simulated_map_resolution = 1.5;

    const auto error_output{ CapturePotentialAnalysisStderr(request) };

    EXPECT_EQ(error_output.find(kMapNormalizationWarning), std::string::npos);
}

TEST(CommandScenariosTest, PotentialAnalysisCliMapNormalizationTrueStillSkipsSimulationNormalization)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_normalization_cli_true" };
    const auto map_path{ WriteConstantMapFixture(temp_dir.path() / "maps") };

    const auto error_output{
        CapturePotentialAnalysisCliStderr({
            "RHBM-GEM",
            "potential_analysis",
            "--database",
            (temp_dir.path() / "analysis.sqlite").string(),
            "--model",
            command_test::TestDataPath("test_model.cif").string(),
            "--map",
            map_path.string(),
            "--simulation",
            "true",
            "--sim-resolution",
            "1.5",
            "--map-normalization",
            "true",
        })
    };

    EXPECT_EQ(error_output.find(kMapNormalizationWarning), std::string::npos);
}

TEST(CommandScenariosTest, PotentialAnalysisCliMapNormalizationFalseSkipsNormalization)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_normalization_cli_false" };
    const auto map_path{ WriteConstantMapFixture(temp_dir.path() / "maps") };

    const auto error_output{
        CapturePotentialAnalysisCliStderr({
            "RHBM-GEM",
            "potential_analysis",
            "--database",
            (temp_dir.path() / "analysis.sqlite").string(),
            "--model",
            command_test::TestDataPath("test_model.cif").string(),
            "--map",
            map_path.string(),
            "--map-normalization",
            "false",
        })
    };

    EXPECT_EQ(error_output.find(kMapNormalizationWarning), std::string::npos);
}

TEST(CommandScenariosTest, PotentialAnalysisCliAcceptsExcludeHydrogen)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_exclude_hydrogen_cli" };
    const auto model_path{ WriteCarbonHydrogenModelFixture(temp_dir.path() / "models") };
    const auto map_path{
        command_test::GenerateMapFile(temp_dir.path() / "maps", model_path)
    };

    std::vector<std::string> args{
        "RHBM-GEM",
        "potential_analysis",
        "--database",
        (temp_dir.path() / "analysis.sqlite").string(),
        "--model",
        model_path.string(),
        "--map",
        map_path.string(),
        "--save-key",
        "exclude_hydrogen_cli",
        "--exclude-hydrogen",
        "true",
    };
    std::vector<char *> argv;
    argv.reserve(args.size());
    for (auto & arg : args)
    {
        argv.push_back(arg.data());
    }

    const auto result{ RunCommandCLI(static_cast<int>(argv.size()), argv.data()) };

    EXPECT_EQ(result, 0);
}

TEST(CommandScenariosTest, PotentialAnalysisExcludeHydrogenFiltersOnlyHydrogenAtoms)
{
    command_test::ScopedTempDir temp_dir{ "potential_analysis_exclude_hydrogen_behavior" };
    const auto model_path{ WriteCarbonHydrogenModelFixture(temp_dir.path() / "models") };
    const auto map_path{
        command_test::GenerateMapFile(temp_dir.path() / "maps", model_path)
    };
    const auto database_path{ temp_dir.path() / "analysis.sqlite" };

    PotentialAnalysisRequest include_request{};
    include_request.database_path = database_path;
    include_request.model_file_path = model_path;
    include_request.map_file_path = map_path;
    include_request.saved_key_tag = "include_hydrogen";
    include_request.exclude_hydrogen = false;
    ASSERT_TRUE(RunCommand(include_request).succeeded);

    PotentialAnalysisRequest exclude_request{ include_request };
    exclude_request.saved_key_tag = "exclude_hydrogen";
    exclude_request.exclude_hydrogen = true;
    ASSERT_TRUE(RunCommand(exclude_request).succeeded);

    DataRepository repository{ database_path };
    auto include_model{ repository.LoadModel("include_hydrogen") };
    auto exclude_model{ repository.LoadModel("exclude_hydrogen") };
    ASSERT_NE(include_model, nullptr);
    ASSERT_NE(exclude_model, nullptr);

    EXPECT_EQ(include_model->GetSelectedAtomCount(), 2u);
    ASSERT_EQ(exclude_model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(exclude_model->GetSelectedAtoms().front()->GetElement(), Element::CARBON);
}

TEST(CommandScenariosTest, RHBMTestRejectsInvertedFitRangeAtPrepare)
{
    RHBMTestRequest request{};
    request.fit_range_min = 2.0;
    request.fit_range_max = 1.0;

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "request"));
}

TEST(CommandScenariosTest, PotentialDisplayRejectsMalformedReferenceGroups)
{
    PotentialDisplayRequest request{};
    request.painter_choice = PainterType::GAUS;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups[""] = { "invalid" };

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "--ref-group"));
}

TEST(CommandScenariosTest, PotentialDisplayAllowsWellFormedReferenceGroupsPastPrepare)
{
    PotentialDisplayRequest request{};
    request.painter_choice = PainterType::GAUS;
    request.model_key_tag_list = { "model_key" };
    request.reference_model_groups = {
        { "with_charge", { "ref_a", "ref_b" } },
        { "no_charge", { "ref_c" } },
    };

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_FALSE(HasDiagnosticForOption(result.issues, "--ref-group"));
}

TEST(CommandScenariosTest, ResultDumpRequiresMapFileForMapPrinterAtPrepare)
{
    ResultDumpRequest request{};
    request.printer_choice = PrinterType::MAP_VALUE;
    request.model_key_tag_list = { "model" };

    const auto result{ RunCommand(request) };

    EXPECT_FALSE(result.succeeded);
    EXPECT_TRUE(HasDiagnosticForOption(result.issues, "request"));
}

TEST(CommandScenariosTest, MapSimulationGeneratesMapForEachValidBlurringWidth)
{
    command_test::ScopedTempDir temp_dir{ "map_simulation_valid" };

    MapSimulationRequest request{};
    request.output_dir = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = { 1.0, -2.0, 3.0 };

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto result{ RunCommand(request) };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    static_cast<void>(testing::internal::GetCapturedStdout());

    ASSERT_TRUE(result.succeeded);
    EXPECT_EQ(command_test::CountFilesWithExtension(temp_dir.path(), ".map"), 2);
    EXPECT_EQ(
        error_output.find("MapObject::SetMapValueArray - MapObject already has a map value array"),
        std::string::npos);

    std::filesystem::path generated_map_path;
    for (const auto & entry : std::filesystem::directory_iterator(temp_dir.path()))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".map")
        {
            generated_map_path = entry.path();
            break;
        }
    }
    ASSERT_FALSE(generated_map_path.empty());

    auto loaded_map{ ReadMap(generated_map_path) };
    ASSERT_NE(loaded_map, nullptr);
    ASSERT_NE(loaded_map->GetMapValueArray(), nullptr);
    bool has_non_zero_value{ false };
    for (size_t i = 0; i < loaded_map->GetMapValueArraySize(); i++)
    {
        if (loaded_map->GetMapValue(i) != 0.0f)
        {
            has_non_zero_value = true;
            break;
        }
    }
    EXPECT_TRUE(has_non_zero_value);
}

TEST(CommandScenariosTest, MapSimulationEmptyModelUsesZeroOrigin)
{
    command_test::ScopedTempDir temp_dir{ "map_simulation_empty_model" };

    MapSimulationRequest request{};
    request.output_dir = temp_dir.path();
    request.map_file_name = "sim_map";
    request.model_file_path = command_test::TestDataPath("test_model_no_atoms.cif");
    request.blurring_width_list = { 1.0 };

    ASSERT_TRUE(RunCommand(request).succeeded);
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

    auto loaded_map{ ReadMap(map_path) };
    ASSERT_NE(loaded_map, nullptr);
    EXPECT_EQ(loaded_map->GetGridSize(), (std::array<int, 3>{ 1, 1, 1 }));
    EXPECT_EQ(loaded_map->GetOrigin(), (std::array<float, 3>{ 0.0f, 0.0f, 0.0f }));
    EXPECT_NE(loaded_map->GetMapValueArray(), nullptr);
}

TEST(CommandScenariosTest, ResultDumpUsesCurrentRequestPaths)
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

    ResultDumpRequest request{};
    request.printer_choice = PrinterType::ATOM_POSITION;
    request.model_key_tag_list = { "shared_key" };

    request.database_path = database_path_a;
    request.output_dir = output_dir_a;
    ASSERT_TRUE(RunCommand(request).succeeded);
    EXPECT_TRUE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_A.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_a / "atom_position_list_MODEL_FROM_B.csv"));

    request.database_path = database_path_b;
    request.output_dir = output_dir_b;
    ASSERT_TRUE(RunCommand(request).succeeded);
    EXPECT_TRUE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_B.csv"));
    EXPECT_FALSE(std::filesystem::exists(output_dir_b / "atom_position_list_MODEL_FROM_A.csv"));
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir_b, ".csv"), 1);
}

TEST(CommandScenariosTest, ResultDumpDatabaseLoadFailureRetainsCommandContext)
{
    command_test::ScopedTempDir temp_dir{ "result_dump_missing_model" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto database_path{ temp_dir.path() / "db" / "database.sqlite" };

    std::filesystem::create_directories(database_path.parent_path());
    command_test::SeedSavedModel(database_path, model_path, "present_key", "MODEL_PRESENT");

    ResultDumpRequest request{};
    request.database_path = database_path;
    request.model_key_tag_list = { "missing_key" };
    request.printer_choice = PrinterType::ATOM_POSITION;

    testing::internal::CaptureStderr();
    EXPECT_FALSE(RunCommand(request).succeeded);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    EXPECT_NE(error_output.find("ResultDumpCommand::BuildDataObjectList"), std::string::npos);
    EXPECT_NE(error_output.find("Failed to load dump inputs"), std::string::npos);
    EXPECT_NE(error_output.find("missing_key"), std::string::npos);
}

TEST(CommandScenariosTest, MapSimulationFileLoadFailureRetainsCommandContext)
{
    command_test::ScopedTempDir temp_dir{ "map_simulation_wrong_input_type" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto wrong_input_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    MapSimulationRequest request{};
    request.output_dir = temp_dir.path() / "out";
    request.model_file_path = wrong_input_path;
    request.map_file_name = "sim_map";
    request.blurring_width_list = { 1.0 };

    testing::internal::CaptureStderr();
    EXPECT_FALSE(RunCommand(request).succeeded);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    EXPECT_NE(error_output.find("MapSimulationCommand::BuildDataObject"), std::string::npos);
    EXPECT_NE(error_output.find("Failed to load model file"), std::string::npos);
    EXPECT_NE(error_output.find("ModelObject"), std::string::npos);
}

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE

TEST(CommandScenariosTest, MapVisualizationInvalidAtomIdFailsWithoutWritingOutput)
{
    command_test::ScopedTempDir temp_dir{ "map_visualization_invalid" };
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")
    };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    MapVisualizationRequest request{};
    request.output_dir = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 999;

    EXPECT_FALSE(RunCommand(request).succeeded);
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 0);
}

#ifdef HAVE_ROOT
TEST(CommandScenariosTest, MapVisualizationWritesPdfToConfiguredFolder)
{
    command_test::ScopedTempDir temp_dir{ "map_visualization_valid" };
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")
    };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    MapVisualizationRequest request{};
    request.output_dir = output_dir;
    request.model_file_path = model_path;
    request.map_file_path = map_path;
    request.atom_serial_id = 1;

    ASSERT_TRUE(RunCommand(request).succeeded);
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".pdf"), 1);
}
#endif

TEST(CommandScenariosTest, PositionEstimationDoesNotRequireDatabaseConfiguration)
{
    command_test::ScopedTempDir temp_dir{ "position_estimation" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_dir{ temp_dir.path() / "map" };
    const auto output_dir{ temp_dir.path() / "out" };
    const auto map_path{ command_test::GenerateMapFile(map_dir, model_path, "fixture_map") };

    PositionEstimationRequest request{};
    request.output_dir = output_dir;
    request.map_file_path = map_path;

    ASSERT_TRUE(RunCommand(request).succeeded);
    EXPECT_EQ(command_test::CountFilesWithExtension(output_dir, ".cmm"), 1);
}

#endif
