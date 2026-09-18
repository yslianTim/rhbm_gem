#include "SimulationManifest.hpp"
#include "SimulationBuildInfo.hpp"

// Compile Boost.JSON once inside the library; the existing Boost dependency supplies its headers.
#include <boost/json/src.hpp>
#include <boost/hash2/sha2.hpp>

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>

#include <atomic>
#include <chrono>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <system_error>

namespace rhbm_gem::core::simulation {
namespace {
namespace json = boost::json;

json::value OptionalNumber(const std::optional<double> & value)
{
    return value ? json::value(*value) : json::value(nullptr);
}

json::value TableName(const std::optional<ChargeTable> & table)
{
    if (!table) return nullptr;
    switch (*table)
    {
    case ChargeTable::Buried: return "buried";
    case ChargeTable::Helix: return "helix";
    case ChargeTable::Sheet: return "sheet";
    case ChargeTable::Amber95: return "amber95";
    }
    throw std::logic_error("Invalid charge table.");
}

const char * PotentialModelName(PotentialModel model)
{
    switch (model)
    {
    case PotentialModel::SINGLE_GAUS: return "single_gaus";
    case PotentialModel::FIVE_GAUS_CHARGE: return "five_gaus_charge";
    case PotentialModel::SINGLE_GAUS_USER: return "single_gaus_user";
    }
    throw std::logic_error("Invalid potential model.");
}

const char * ChargeModeName(PartialCharge mode)
{
    switch (mode)
    {
    case PartialCharge::NEUTRAL: return "neutral";
    case PartialCharge::PARTIAL: return "partial";
    case PartialCharge::AMBER: return "amber";
    }
    throw std::logic_error("Invalid charge mode.");
}

json::object MakeManifest(const std::filesystem::path & output, const std::string & map_sha256,
    const MapObject & map, const SimulationAtomPreparationResult & atoms,
    const MapSimulationRequest & request, double blurring_width, int actual_job_count,
    const SimulationSource & source)
{
    ElectricPotential potential;
    potential.SetModelChoice(static_cast<int>(request.potential_model_choice));
    potential.SetBlurringWidth(blurring_width);
    const auto kernel{ potential.GetKernelSettings() };
    json::array atom_records;
    atom_records.reserve(atoms.atom_list.size());
    size_t fallback_count{ 0 };
    for (const auto & atom : atoms.atom_list)
    {
        const auto widths{ potential.GetEffectiveWidths(atom.element) };
        if (atom.charge_lookup && !atom.charge_lookup->charge) ++fallback_count;
        atom_records.emplace_back(json::object{
            { "preparation_index", atom_records.size() },
            { "serial_id", atom.serial_id }, { "sequence_id", atom.sequence_id },
            { "chain_id", atom.chain_id }, { "component_id", atom.component_id },
            { "atom_id", atom.atom_id }, { "alternate_indicator", atom.alternate_indicator },
            { "element", static_cast<int>(atom.element) },
            { "residue", static_cast<int>(atom.residue) },
            { "spot", static_cast<int>(atom.spot) },
            { "structure", static_cast<int>(atom.structure) },
            { "position", json::value_from(atom.position) },
            { "table", atom.charge_lookup ? TableName(atom.charge_lookup->table) : json::value(nullptr) },
            { "lookup_status", std::string(ChargeStatusText(atom)) },
            { "lookup_charge", atom.charge_lookup ? OptionalNumber(atom.charge_lookup->charge) : json::value(nullptr) },
            { "charge_used", atom.charge_used },
            { "effective_gaussian_width", OptionalNumber(widths.gaussian) },
            { "effective_charge_width", OptionalNumber(widths.charge) }
        });
    }
#ifdef USE_OPENMP
    constexpr bool openmp_enabled{ true };
#else
    constexpr bool openmp_enabled{ false };
#endif
    return {
        { "schema_version", 2 },
        { "generator", json::object{
            { "version", RHBM_GEM_SIMULATION_VERSION },
            { "source_sha256", RHBM_GEM_SIMULATION_SOURCE_SHA256 },
            { "configuration_sha256", RHBM_GEM_SIMULATION_CONFIG_SHA256 },
            { "build_sha256", RHBM_GEM_SIMULATION_BUILD_SHA256 }
        } },
        { "source", json::object{
            { "model_path", std::filesystem::absolute(request.model_file_path).string() },
            { "pdb_id", source.pdb_id }, { "model_sha256", source.model_sha256 }
        } },
        { "output", json::object{
            { "map_file", output.filename().string() }, { "map_sha256", map_sha256 },
            { "format", "ccp4" }, { "calculation_precision", "float64" }, { "storage_precision", "float32" }
        } },
        { "settings", json::object{
            { "potential_model", PotentialModelName(request.potential_model_choice) },
            { "potential_model_code", static_cast<int>(request.potential_model_choice) },
            { "charge_mode", ChargeModeName(request.partial_charge_choice) },
            { "charge_mode_code", static_cast<int>(request.partial_charge_choice) },
            { "blurring_width", blurring_width },
            { "blurring_width_list", json::value_from(request.blurring_width_list) },
            { "grid_spacing", json::value_from(map.GetGridSpacing()) },
            { "grid_size", json::value_from(map.GetGridSize()) },
            { "origin", json::value_from(map.GetOrigin()) },
            { "cutoff_distance", request.cutoff_distance },
            { "coordinate_unit", "angstrom" }, { "exclude_hydrogen", request.exclude_hydrogen },
            { "only_backbone", request.only_backbone }, { "missing_charge_policy", "zero_with_status" },
            { "occupancy_applied", false }, { "temperature_factor_applied", false },
            { "normalization_applied", false }
        } },
        { "kernel", json::object{
            { "version", std::string(PotentialModelName(request.potential_model_choice)) + "-v1" },
            { "width_policy", request.potential_model_choice == PotentialModel::SINGLE_GAUS
                ? json::value(json::object{ { "id", "element-scaled-v1" },
                    { "oxygen", 0.8 }, { "nitrogen", 0.9 }, { "other", 1.0 } })
                : json::value(json::object{ { "id", "model-specific-v1" } }) },
            { "charge_term_cutoff", OptionalNumber(kernel.charge_term_cutoff) },
            { "near_zero_distance", OptionalNumber(kernel.near_zero_distance) },
            { "minimum_charge_width", OptionalNumber(kernel.minimum_charge_width) }
        } },
        { "support", json::object{
            { "version", "sphere-fma-v1" }, { "outer_cutoff", request.cutoff_distance },
            { "comparison", "squared_distance<=cutoff*cutoff" },
            { "coordinates", "fma(index,spacing,origin)" },
            { "squared_distance", "fma(dz,dz,fma(dy,dy,dx*dx))" }
        } },
        { "execution", json::object{
            { "openmp_enabled", openmp_enabled }, { "requested_job_count", request.job_count },
            { "actual_job_count", actual_job_count }, { "accumulation", "z_planes_in_preparation_order" }
        } },
        { "charge_semantics", "argument_passed_to_electric_potential" },
        { "atom_count", atoms.atom_list.size() }, { "fallback_charge_count", fallback_count },
        { "atoms", std::move(atom_records) }
    };
}

class StagingDirectory
{
public:
    std::filesystem::path path;
    bool preserve{ false };

    explicit StagingDirectory(const std::filesystem::path & parent)
    {
        static std::atomic<unsigned long long> sequence{ 0 };
        const auto stamp{ std::chrono::steady_clock::now().time_since_epoch().count() };
        path = parent / (".simulation-" + std::to_string(stamp) + "-" + std::to_string(sequence++));
        if (!std::filesystem::create_directory(path))
            throw std::runtime_error("Cannot create simulation staging directory: " + path.string());
    }
    ~StagingDirectory()
    {
        if (preserve) return;
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
};

} // namespace

std::string FileSha256(const std::filesystem::path & path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read file for SHA-256: " + path.string());
    boost::hash2::sha2_256 hash;
    std::array<char, 65536> buffer{};
    while (input.read(buffer.data(), static_cast<std::streamsize>(buffer.size())) || input.gcount() > 0)
        hash.update(buffer.data(), static_cast<size_t>(input.gcount()));
    if (!input.eof() || input.bad())
        throw std::runtime_error("Failed while hashing file: " + path.string());
    return boost::hash2::to_string(hash.result());
}

void WriteSimulationArtifacts(const std::filesystem::path & output, const MapObject & map,
    const SimulationAtomPreparationResult & atoms, const MapSimulationRequest & request,
    double blurring_width, int actual_job_count, const SimulationSource & source)
{
    const std::array<std::filesystem::path, 2> destinations{ output, output.string() + ".simulation.json" };
    for (const auto & destination : destinations)
    {
        if (std::filesystem::exists(destination) && !std::filesystem::is_regular_file(destination))
            throw std::runtime_error("Simulation output is not a regular file: " + destination.string());
    }
    StagingDirectory staging(output.parent_path());
    const std::array<std::filesystem::path, 2> temporary{ staging.path / "new.map", staging.path / "new.json" };
    const std::array<std::filesystem::path, 2> backup{ staging.path / "previous.map", staging.path / "previous.json" };
    WriteMap(temporary[0], map);
    const auto manifest{ MakeManifest(output, FileSha256(temporary[0]), map, atoms,
        request, blurring_width, actual_job_count, source) };
    std::ofstream json_file;
    json_file.exceptions(std::ios::badbit | std::ios::failbit);
    json_file.open(temporary[1], std::ios::binary);
    json_file << json::serialize(manifest) << '\n';
    json_file.close();

    std::array<bool, 2> backed_up{}, published{};
    try
    {
        for (size_t i = 0; i < destinations.size(); ++i)
        {
            if (std::filesystem::exists(destinations[i]))
            {
                std::filesystem::rename(destinations[i], backup[i]);
                backed_up[i] = true;
            }
        }
        // The manifest is published last; its map hash also detects interruption between renames.
        for (size_t i = 0; i < destinations.size(); ++i)
        {
            std::filesystem::rename(temporary[i], destinations[i]);
            published[i] = true;
        }
    }
    catch (...)
    {
        const auto original_error{ std::current_exception() };
        for (size_t i = 0; i < destinations.size(); ++i)
        {
            std::error_code error;
            if (published[i]) std::filesystem::remove(destinations[i], error);
            if (error) staging.preserve = true;
            if (backed_up[i]) std::filesystem::rename(backup[i], destinations[i], error);
            if (error) staging.preserve = true;
        }
        if (staging.preserve)
            throw std::runtime_error("Simulation publication failed; recovery files retained in " + staging.path.string());
        std::rethrow_exception(original_error);
    }
}

} // namespace rhbm_gem::core::simulation
