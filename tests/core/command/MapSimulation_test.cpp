#include <gtest/gtest.h>
#include <boost/json.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <numbers>

#include "support/CommandTestHelpers.hpp"
#include "command/detail/MapSimulation.hpp"
#include "command/detail/SimulationManifest.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>

using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace simulation = rhbm_gem::core::simulation;
namespace json = boost::json;

namespace {
std::string ReadBytes(const std::filesystem::path & path)
{
    std::ifstream input(path, std::ios::binary);
    return { std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>() };
}

json::object ReadManifest(const std::filesystem::path & map_path)
{
    return json::parse(ReadBytes(map_path.string() + ".simulation.json")).as_object();
}

std::filesystem::path FindMap(const std::filesystem::path & directory)
{
    for (const auto & entry : std::filesystem::directory_iterator(directory))
        if (entry.path().extension() == ".map") return entry.path();
    throw std::runtime_error("No simulation map found.");
}

// Independent scalar reference, including the existing single-Gaussian charge cutoff.
double SingleGaussianReference(Element element, double distance, double width, double charge)
{
    if (element == Element::OXYGEN) width *= 0.8;
    else if (element == Element::NITROGEN) width *= 0.9;
    const double variance{ width * width };
    double charge_term{};
    if (distance < 1.0e-5) charge_term = charge * std::sqrt(2.0 / std::numbers::pi) / width;
    else if (distance <= 2.5) charge_term = charge / distance * std::erf(distance / width / std::sqrt(2.0));
    return static_cast<double>(element) * std::pow(2.0 * std::numbers::pi * variance, -1.5)
        * std::exp(-(distance * distance) / (2.0 * variance)) + charge_term;
}

double ReferenceVoxel(const std::array<double, 3> & position,
    const std::vector<simulation::SimulationAtom> & atoms, double width, double cutoff,
    PotentialModel model = PotentialModel::SINGLE_GAUS)
{
    ElectricPotential potential;
    potential.SetModelChoice(static_cast<int>(model));
    potential.SetBlurringWidth(width);
    double sum{};
    for (const auto & atom : atoms)
    {
        double distance_square{};
        for (size_t axis = 0; axis < 3; ++axis)
        {
            const double delta{ position[axis] - atom.position[axis] };
            distance_square += delta * delta;
        }
        if (distance_square > cutoff * cutoff) continue;
        const double distance{ std::sqrt(distance_square) };
        sum += model == PotentialModel::SINGLE_GAUS
            ? SingleGaussianReference(atom.element, distance, width, atom.charge_used)
            : potential.GetPotentialValue(atom.element, distance, atom.charge_used);
    }
    return sum;
}
} // namespace

class MapSimulationTest : public ::testing::Test
{
protected:
    command_test::ScopedTempDir model_directory{ "simulation_models" };

    MapSimulationRequest MakeRequest(const std::filesystem::path & output_dir)
    {
        MapSimulationRequest request;
        request.model_file_path = model_directory.path() / "overlapping.cif";
        std::ofstream model_file(request.model_file_path);
        model_file << R"cif(data_simulation
    #
    loop_
    _database_2.database_id
    _database_2.database_code
    PDB 1SIM
    #
    loop_
    _entity.id
    _entity.type
    _entity.pdbx_number_of_molecules
    1 polymer 1
    #
    loop_
    _struct_asym.id
    _struct_asym.entity_id
    A 1
    B 1
    #
    loop_
    _struct_conf.id
    _struct_conf.conf_type_id
    _struct_conf.beg_label_asym_id
    _struct_conf.beg_label_seq_id
    _struct_conf.end_label_asym_id
    _struct_conf.end_label_seq_id
    HELX1 HELX_P A 2 A 2
    #
    loop_
    _struct_sheet_range.sheet_id
    _struct_sheet_range.id
    _struct_sheet_range.beg_label_asym_id
    _struct_sheet_range.beg_label_seq_id
    _struct_sheet_range.end_label_asym_id
    _struct_sheet_range.end_label_seq_id
    S1 1 A 3 A 3
    #
    loop_
    _atom_type.symbol
    C
    O
    H
    #
    loop_
    _atom_site.group_PDB
    _atom_site.id
    _atom_site.type_symbol
    _atom_site.label_atom_id
    _atom_site.label_alt_id
    _atom_site.label_comp_id
    _atom_site.label_asym_id
    _atom_site.label_seq_id
    _atom_site.Cartn_x
    _atom_site.Cartn_y
    _atom_site.Cartn_z
    _atom_site.occupancy
    _atom_site.B_iso_or_equiv
    _atom_site.pdbx_PDB_model_num
    ATOM 1 O O . ALA A 1 0.0 0.0 0.0 1.0 0.0 1
    ATOM 2 O O . ALA A 2 0.5 0.0 0.0 1.0 0.0 1
    ATOM 3 C CB . VAL A 3 0.0 0.5 0.0 1.0 0.0 1
    ATOM 4 H H . ALA A 1 0.0 0.0 0.5 1.0 0.0 1
    ATOM 5 C Q . UNK B 4 -0.5 0.0 0.0 1.0 0.0 1
    #
    )cif";
        model_file.close();
        request.output_dir = output_dir;
        request.potential_model_choice = PotentialModel::SINGLE_GAUS;
        request.blurring_width_list = { 1.0 };
        request.grid_spacing = 0.5;
        request.cutoff_distance = 3.5;
        request.verbosity = 0;
        return request;
    }
};

TEST_F(MapSimulationTest, PreparedAtomsDistinguishTableChargesFallbacksAndNeutralMode)
{
    const auto request{ MakeRequest({}) };
    auto model{ ReadModel(request.model_file_path) };
    model->SelectAllAtoms();
    const auto atoms{ simulation::PrepareSimulationAtomList(*model, request) };
    ASSERT_EQ(atoms.atom_list.size(), 5u);
    EXPECT_DOUBLE_EQ(atoms.atom_list[0].charge_used, -0.673);
    EXPECT_DOUBLE_EQ(atoms.atom_list[1].charge_used, -0.701);
    EXPECT_EQ(atoms.atom_list[0].charge_lookup->table, ChargeTable::Buried);
    EXPECT_EQ(atoms.atom_list[1].charge_lookup->table, ChargeTable::Helix);
    EXPECT_DOUBLE_EQ(atoms.atom_list[2].charge_used, 0.0);
    EXPECT_EQ(simulation::ChargeStatusText(atoms.atom_list[2]), "found");
    EXPECT_EQ(simulation::ChargeStatusText(atoms.atom_list[3]), "unsupported_spot");
    EXPECT_EQ(simulation::ChargeStatusText(atoms.atom_list[4]), "unsupported_residue");
    for (size_t i = 3; i < 5; ++i)
    {
        EXPECT_DOUBLE_EQ(atoms.atom_list[i].charge_used, 0.0);
        EXPECT_FALSE(atoms.atom_list[i].charge_lookup->charge);
    }
    auto neutral{ request };
    neutral.partial_charge_choice = PartialCharge::NEUTRAL;
    for (const auto & atom : simulation::PrepareSimulationAtomList(*model, neutral).atom_list)
    {
        EXPECT_DOUBLE_EQ(atom.charge_used, 0.0);
        EXPECT_EQ(simulation::ChargeStatusText(atom), "neutral_mode");
        EXPECT_FALSE(atom.charge_lookup);
    }
    // The preparation vector does not collapse atoms with the same serial identifier.
    model->GetAtomList()[1]->SetSerialID(model->GetAtomList()[0]->GetSerialID());
    const auto duplicate_serial{ simulation::PrepareSimulationAtomList(*model, request) };
    EXPECT_EQ(duplicate_serial.atom_list[0].serial_id, duplicate_serial.atom_list[1].serial_id);
    EXPECT_NE(duplicate_serial.atom_list[0].charge_used, duplicate_serial.atom_list[1].charge_used);
}

TEST_F(MapSimulationTest, OverlappingVoxelsMatchSerialReferenceAndRepeatBitForBitAcrossTeams)
{
    auto request{ MakeRequest({}) };
    auto model{ ReadModel(request.model_file_path) };
    model->SelectAllAtoms();
    const auto atoms{ simulation::PrepareSimulationAtomList(*model, request) };
    for (const auto potential_model : { PotentialModel::SINGLE_GAUS, PotentialModel::FIVE_GAUS_CHARGE })
    {
        request.potential_model_choice = potential_model;
        std::vector<std::uint64_t> reference_bits;
        for (const int jobs : { 1, 2, 4, 4, 2, 1 })
        {
            request.job_count = jobs;
            auto map{ simulation::CreateMapObject(request, atoms) };
            const auto actual{ simulation::PopulateMapValueArray(*map, atoms, request, 1.0) };
#ifdef USE_OPENMP
            EXPECT_EQ(actual, jobs);
#else
            EXPECT_EQ(actual, 1);
#endif
            for (size_t i = 0; i < map->GetMapValueArraySize(); ++i)
            {
                const double value{ map->GetMapValue(i) };
                const auto bits{ std::bit_cast<std::uint64_t>(value) };
                if (reference_bits.size() <= i) reference_bits.push_back(bits);
                else EXPECT_EQ(bits, reference_bits[i]);
                EXPECT_NEAR(value, ReferenceVoxel(map->GetGridPosition(i), atoms.atom_list,
                    1.0, request.cutoff_distance, potential_model), 1.0e-12);
            }
        }
    }
}

TEST_F(MapSimulationTest, AtomCenterAndBothCutoffBoundariesArePreserved)
{
    auto request{ MakeRequest({}) };
    request.job_count = 4;
    simulation::SimulationAtomPreparationResult atoms;
    atoms.atom_list.push_back(simulation::SimulationAtom{
        .element = Element::CARBON, .position = { 0.0, 0.0, 0.0 }, .charge_used = -0.673 });
    // Includes the center, both signs of the outer cutoff, and voxels beyond it.
    MapObject map({ 17, 1, 1 }, { 0.5, 0.5, 0.5 }, { -4.0, 0.0, 0.0 });
    simulation::PopulateMapValueArray(map, atoms, request, 1.0);
    for (size_t i = 0; i < map.GetMapValueArraySize(); ++i)
        EXPECT_DOUBLE_EQ(map.GetMapValue(i), ReferenceVoxel(map.GetGridPosition(i), atoms.atom_list, 1.0, 3.5));
    EXPECT_DOUBLE_EQ(map.GetMapValue(0), 0.0);
    EXPECT_NE(map.GetMapValue(1), 0.0);
    EXPECT_NE(map.GetMapValue(15), 0.0);
    EXPECT_DOUBLE_EQ(map.GetMapValue(16), 0.0);
}

TEST_F(MapSimulationTest, ManifestReconstructsSavedVoxelsAndRecordsExactSettings)
{
    command_test::ScopedTempDir temp("simulation_manifest");
    auto request{ MakeRequest(temp.path()) };
    request.job_count = 4;
    request.blurring_width_list = { 1.1234567890123457, -2.0, 2.0 };
    ASSERT_TRUE(RunCommand(request).succeeded);
    EXPECT_EQ(command_test::CountFilesWithExtension(temp.path(), ".map"), 2u);
    EXPECT_EQ(command_test::CountFilesWithExtension(temp.path(), ".json"), 2u);
    for (const auto & entry : std::filesystem::directory_iterator(temp.path()))
    {
        if (entry.path().extension() != ".map") continue;
        const auto manifest{ ReadManifest(entry.path()) };
        EXPECT_EQ(manifest.at("schema_version").as_int64(), 2);
        EXPECT_EQ(manifest.at("source").at("model_sha256").as_string(), simulation::FileSha256(request.model_file_path));
        EXPECT_EQ(manifest.at("output").at("map_sha256").as_string(), simulation::FileSha256(entry.path()));
        EXPECT_EQ(manifest.at("output").at("map_file").as_string(), entry.path().filename().string());
        EXPECT_EQ(manifest.at("atom_count").as_int64(), 5);
        EXPECT_EQ(manifest.at("fallback_charge_count").as_int64(), 2);
        for (const auto * key : { "source_sha256", "configuration_sha256", "build_sha256" })
            EXPECT_EQ(manifest.at("generator").at(key).as_string().size(), 64u);
        const auto & settings{ manifest.at("settings") };
        const auto widths{ json::value_to<std::vector<double>>(settings.at("blurring_width_list")) };
        ASSERT_EQ(widths.size(), 2u);
        EXPECT_EQ(std::bit_cast<std::uint64_t>(widths[0]), std::bit_cast<std::uint64_t>(request.blurring_width_list[0]));
        const double width{ json::value_to<double>(settings.at("blurring_width")) };
        EXPECT_TRUE(width == widths[0] || width == widths[1]);
        EXPECT_DOUBLE_EQ(json::value_to<double>(manifest.at("kernel").at("charge_term_cutoff")), 2.5);
        EXPECT_FALSE(manifest.at("kernel").as_object().contains("effective_charge_width"));
        EXPECT_EQ(manifest.at("support").at("version"), "sphere-fma-v1");
        EXPECT_EQ(settings.at("missing_charge_policy").as_string(), "zero_with_status");
        EXPECT_EQ(manifest.at("execution").at("requested_job_count").as_int64(), 4);
#ifdef USE_OPENMP
        EXPECT_TRUE(manifest.at("execution").at("openmp_enabled").as_bool());
        EXPECT_EQ(manifest.at("execution").at("actual_job_count").as_int64(), 4);
#else
        EXPECT_FALSE(manifest.at("execution").at("openmp_enabled").as_bool());
        EXPECT_EQ(manifest.at("execution").at("actual_job_count").as_int64(), 1);
#endif
        std::vector<simulation::SimulationAtom> reconstructed;
        for (const auto & record : manifest.at("atoms").as_array())
        {
            EXPECT_EQ(json::value_to<size_t>(record.at("preparation_index")), reconstructed.size());
            const int element{ json::value_to<int>(record.at("element")) };
            const double effective{ width * (element == 8 ? 0.8 : element == 7 ? 0.9 : 1.0) };
            EXPECT_DOUBLE_EQ(json::value_to<double>(record.at("effective_gaussian_width")), effective);
            EXPECT_DOUBLE_EQ(json::value_to<double>(record.at("effective_charge_width")), effective);
            reconstructed.push_back(simulation::SimulationAtom{
                .element = static_cast<Element>(record.at("element").as_int64()),
                .position = json::value_to<std::array<double, 3>>(record.at("position")),
                .charge_used = json::value_to<double>(record.at("charge_used")) });
        }
        EXPECT_EQ(manifest.at("atoms").as_array()[0].at("table").as_string(), "buried");
        EXPECT_EQ(manifest.at("atoms").as_array()[1].at("table").as_string(), "helix");
        EXPECT_EQ(manifest.at("atoms").as_array()[2].at("lookup_status").as_string(), "found");
        EXPECT_TRUE(manifest.at("atoms").as_array()[3].at("lookup_charge").is_null());
        const auto map{ ReadMap(entry.path()) };
        const auto origin{ json::value_to<std::array<double, 3>>(settings.at("origin")) };
        const auto spacing{ json::value_to<std::array<double, 3>>(settings.at("grid_spacing")) };
        EXPECT_EQ(map->GetGridSize(), (json::value_to<std::array<int, 3>>(settings.at("grid_size"))));
        for (size_t i = 0; i < map->GetMapValueArraySize(); ++i)
        {
            const auto grid_index{ map->GetGridIndex(i) };
            std::array<double, 3> position{};
            for (size_t axis = 0; axis < 3; ++axis)
                position[axis] = origin[axis] + static_cast<double>(grid_index[axis]) * spacing[axis];
            const double expected{ ReferenceVoxel(position, reconstructed, width,
                json::value_to<double>(settings.at("cutoff_distance"))) };
            EXPECT_EQ(std::bit_cast<std::uint32_t>(static_cast<float>(map->GetMapValue(i))),
                std::bit_cast<std::uint32_t>(static_cast<float>(expected)));
        }
    }
}

TEST_F(MapSimulationTest, SavedMapsAreIdenticalAcrossRunsAndChargeModeHistory)
{
    command_test::ScopedTempDir temp("simulation_repeat");
    auto request{ MakeRequest(temp.path()) };
    std::string partial_map_bytes;
    for (const int jobs : { 1, 2, 4, 4 })
    {
        request.job_count = jobs;
        request.partial_charge_choice = PartialCharge::AMBER;
        ASSERT_TRUE(RunCommand(request).succeeded);
        request.partial_charge_choice = PartialCharge::PARTIAL;
        ASSERT_TRUE(RunCommand(request).succeeded);
        const auto path{ FindMap(temp.path()) };
        const auto bytes{ ReadBytes(path) };
        if (partial_map_bytes.empty()) partial_map_bytes = bytes;
        else EXPECT_EQ(bytes, partial_map_bytes);
        const auto manifest{ ReadManifest(path) };
        EXPECT_EQ(manifest.at("settings").at("charge_mode").as_string(), "partial");
        EXPECT_EQ(manifest.at("output").at("map_sha256").as_string(), simulation::FileSha256(path));
    }
    EXPECT_EQ(std::distance(std::filesystem::directory_iterator(temp.path()), std::filesystem::directory_iterator()), 2);
}

TEST_F(MapSimulationTest, ManifestEscapesAtomIdentityWithoutLosingCharacters)
{
    command_test::ScopedTempDir temp("simulation_strings");
    const auto request{ MakeRequest(temp.path()) };
    auto model{ ReadModel(request.model_file_path) };
    model->SelectAllAtoms();
    auto atoms{ simulation::PrepareSimulationAtomList(*model, request) };
    const std::string identity{ "quoted\"\\line\n\t\r" };
    atoms.atom_list[0].atom_id = identity;
    auto map{ simulation::CreateMapObject(request, atoms) };
    const auto actual{ simulation::PopulateMapValueArray(*map, atoms, request, 1.0) };
    const auto output{ temp.path() / "escaped.map" };
    simulation::WriteSimulationArtifacts(output, *map, atoms, request, 1.0, actual,
        { model->GetPdbID(), simulation::FileSha256(request.model_file_path) });
    EXPECT_EQ(ReadManifest(output).at("atoms").as_array()[0].at("atom_id").as_string(), identity);
}

TEST_F(MapSimulationTest, SelectionAndEmptyModelsHaveMatchingManifests)
{
    command_test::ScopedTempDir temp("simulation_selection");
    auto request{ MakeRequest(temp.path()) };
    request.partial_charge_choice = PartialCharge::NEUTRAL;
    request.exclude_hydrogen = true;
    request.only_backbone = true;
    ASSERT_TRUE(RunCommand(request).succeeded);
    const auto selected{ ReadManifest(FindMap(temp.path())) };
    EXPECT_EQ(selected.at("atom_count").as_int64(), 2);
    EXPECT_EQ(selected.at("fallback_charge_count").as_int64(), 0);
    for (const auto & atom : selected.at("atoms").as_array())
    {
        EXPECT_EQ(atom.at("lookup_status").as_string(), "neutral_mode");
        EXPECT_TRUE(atom.at("table").is_null());
        EXPECT_DOUBLE_EQ(json::value_to<double>(atom.at("charge_used")), 0.0);
    }
    request.output_dir = temp.path() / "empty";
    request.model_file_path = command_test::TestDataPath("test_model_no_atoms.cif");
    ASSERT_TRUE(RunCommand(request).succeeded);
    const auto empty_map{ FindMap(request.output_dir) };
    const auto empty{ ReadManifest(empty_map) };
    EXPECT_EQ(empty.at("atom_count").as_int64(), 0);
    EXPECT_TRUE(empty.at("atoms").as_array().empty());
    EXPECT_DOUBLE_EQ(ReadMap(empty_map)->GetMapValue(0), 0.0);
}

TEST_F(MapSimulationTest, FilenameCollisionsAndNonfiniteOutputsFailWithoutArtifacts)
{
    command_test::ScopedTempDir temp("simulation_invalid");
    auto request{ MakeRequest(temp.path()) };
    for (const auto & widths : { std::vector<double>{ 1.001, 1.004 }, std::vector<double>{ 1.0, 1.0 } })
    {
        request.blurring_width_list = widths;
        EXPECT_FALSE(RunCommand(request).succeeded);
        EXPECT_TRUE(std::filesystem::is_empty(temp.path()));
    }
    request.job_count = 4;
    for (const double width : { 1.0e-20, 1.0e-200 })
    {
        request.blurring_width_list = { width };
        const auto result{ RunCommand(request) };
        EXPECT_FALSE(result.succeeded);
        ASSERT_FALSE(result.issues.empty());
        EXPECT_TRUE(std::filesystem::is_empty(temp.path()));
    }
    request.blurring_width_list = { 1.0 };
    request.potential_model_choice = PotentialModel::SINGLE_GAUS_USER;
    EXPECT_FALSE(RunCommand(request).succeeded);
    EXPECT_TRUE(std::filesystem::is_empty(temp.path()));
}

TEST_F(MapSimulationTest, WorkerExceptionsLeaveThePreviousMapUnchanged)
{
    auto request{ MakeRequest({}) };
    request.job_count = 4;
    request.potential_model_choice = PotentialModel::FIVE_GAUS_CHARGE;
    simulation::SimulationAtomPreparationResult atoms;
    atoms.atom_list.push_back(simulation::SimulationAtom{ .element = Element::HELIUM });
    MapObject map({ 3, 3, 3 }, { 0.5, 0.5, 0.5 }, { 0.0, 0.0, 0.0 });
    auto previous{ std::make_unique<double[]>(27) };
    std::fill_n(previous.get(), 27, 7.0);
    map.SetMapValueArray(std::move(previous));
    EXPECT_THROW(simulation::PopulateMapValueArray(map, atoms, request, 1.0), std::out_of_range);
    for (size_t i = 0; i < 27; ++i) EXPECT_DOUBLE_EQ(map.GetMapValue(i), 7.0);
}

TEST_F(MapSimulationTest, OutputFailurePreservesExistingMapAndCleansStaging)
{
    command_test::ScopedTempDir temp("simulation_write_failure");
    auto request{ MakeRequest(temp.path()) };
    const auto output{ temp.path() / "sim_map_1SIM_bw1.00.map" };
    { std::ofstream file(output); file << "previous map"; }
    const auto blocker{ output.string() + ".simulation.json" };
    std::filesystem::create_directory(blocker);
    EXPECT_FALSE(RunCommand(request).succeeded);
    EXPECT_EQ(ReadBytes(output), "previous map");
    std::filesystem::remove(blocker);
    MapObject invalid_map;
    EXPECT_THROW(simulation::WriteSimulationArtifacts(output, invalid_map, {}, request, 1.0, 1, {}), std::runtime_error);
    EXPECT_EQ(ReadBytes(output), "previous map");
    EXPECT_EQ(std::distance(std::filesystem::directory_iterator(temp.path()), std::filesystem::directory_iterator()), 1);
}

TEST_F(MapSimulationTest, FileHashesMatchStandardSha256Vectors)
{
    command_test::ScopedTempDir temp("simulation_sha256");
    const auto path{ temp.path() / "hash_input" };
    { std::ofstream file(path); }
    EXPECT_EQ(simulation::FileSha256(path), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    { std::ofstream file(path); file << "abc"; }
    EXPECT_EQ(simulation::FileSha256(path), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
}
