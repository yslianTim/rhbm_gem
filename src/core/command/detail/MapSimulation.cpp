#include "MapSimulation.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <map>
#include <stdexcept>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core::simulation {
namespace {

void LogMapSummary(const MapObject & map_object)
{
    std::ostringstream oss;
    oss << "MapObject Summary:\n";
    oss << " o=====================================================o\n";
    oss << " |  Map Object  |   X-axis   |   Y-axis   |   Z-axis   |\n";
    oss << " o=====================================================o\n";
    oss << " | Grid size    | ";
    oss << std::setw(10) << map_object.GetGridSize().at(0) << " | "
        << std::setw(10) << map_object.GetGridSize().at(1) << " | "
        << std::setw(10) << map_object.GetGridSize().at(2) << " |\n";
    oss << " | Grid Spacing | ";
    oss << std::setw(10) << map_object.GetGridSpacing().at(0) << " | "
        << std::setw(10) << map_object.GetGridSpacing().at(1) << " | "
        << std::setw(10) << map_object.GetGridSpacing().at(2) << " |\n";
    oss << " | Origin (A)   | ";
    oss << std::setw(10) << map_object.GetOrigin().at(0) << " | "
        << std::setw(10) << map_object.GetOrigin().at(1) << " | "
        << std::setw(10) << map_object.GetOrigin().at(2) << " |\n";
    oss << " | Map Length(A)| ";
    oss << std::setw(10)
        << static_cast<double>(map_object.GetGridSize().at(0)) * map_object.GetGridSpacing().at(0)
        << " | "
        << std::setw(10)
        << static_cast<double>(map_object.GetGridSize().at(1)) * map_object.GetGridSpacing().at(1)
        << " | "
        << std::setw(10)
        << static_cast<double>(map_object.GetGridSize().at(2)) * map_object.GetGridSpacing().at(2)
        << " |\n";
    oss << " |-----------------------------------------------------|\n";
    oss << " | Map value min  | " << std::setw(34) << map_object.GetMapValueMin() << " |\n";
    oss << " | Map value max  | " << std::setw(34) << map_object.GetMapValueMax() << " |\n";
    oss << " | Map value mean | " << std::setw(34) << map_object.GetMapValueMean() << " |\n";
    oss << " | Map value s.d. | " << std::setw(34) << map_object.GetMapValueSD() << " |\n";
    oss << " o=====================================================o\n";
    Logger::Log(LogLevel::Info, oss.str());
}

} // namespace

std::unique_ptr<MapObject> CreateMapObject(
    const MapSimulationRequest & request,
    const SimulationAtomPreparationResult & result)
{
    ScopeTimer timer("MapSimulationCommand::CreateMapObject");
    std::array<double, 3> grid_spacing{
        request.grid_spacing,
        request.grid_spacing,
        request.grid_spacing
    };
    std::array<double, 3> origin{ 0.0, 0.0, 0.0 };
    std::array<int, 3> grid_size{ 1, 1, 1 };
    if (!result.atom_list.empty())
    {
        for (size_t axis = 0; axis < 3; ++axis)
        {
            const auto cells{ std::ceil((result.range_max[axis] - result.range_min[axis]) / grid_spacing[axis]) };
            if (!std::isfinite(result.range_min[axis]) || !std::isfinite(result.range_max[axis])
                || !std::isfinite(cells) || cells < 1.0
                || cells > static_cast<double>(std::numeric_limits<int>::max()))
                throw std::runtime_error("Simulation grid dimensions are not representable.");
        }
        origin = {
            std::floor(result.range_min[0] / grid_spacing[0]) * grid_spacing[0],
            std::floor(result.range_min[1] / grid_spacing[1]) * grid_spacing[1],
            std::floor(result.range_min[2] / grid_spacing[2]) * grid_spacing[2]
        };
        grid_size = {
            static_cast<int>(std::ceil((result.range_max[0] - result.range_min[0]) / grid_spacing[0])),
            static_cast<int>(std::ceil((result.range_max[1] - result.range_min[1]) / grid_spacing[1])),
            static_cast<int>(std::ceil((result.range_max[2] - result.range_min[2]) / grid_spacing[2]))
        };
    }
    auto map_object{ std::make_unique<MapObject>(grid_size, grid_spacing, origin) };
    map_object->ClearMapValueArray();

    return map_object;
}

std::string_view ChargeStatusText(const SimulationAtom & atom)
{
    if (!atom.charge_lookup) return "neutral_mode";
    switch (atom.charge_lookup->status)
    {
    case ChargeLookupStatus::Found: return "found";
    case ChargeLookupStatus::UnsupportedResidue: return "unsupported_residue";
    case ChargeLookupStatus::UnsupportedStructure: return "unsupported_structure";
    case ChargeLookupStatus::UnsupportedSpot: return "unsupported_spot";
    case ChargeLookupStatus::TableDataMismatch: return "table_data_mismatch";
    }
    throw std::logic_error("Invalid charge lookup status.");
}

SimulationAtomPreparationResult PrepareSimulationAtomList(
    const ModelObject & model, const MapSimulationRequest & request)
{
    SimulationAtomPreparationResult result;
    result.atom_list.reserve(model.GetSelectedAtomCount());
    std::map<std::string_view, size_t> failures;
    for (const auto * atom : model.GetSelectedAtoms())
    {
        SimulationAtom prepared{
            .serial_id = atom->GetSerialID(),
            .sequence_id = atom->GetSequenceID(),
            .chain_id = atom->GetChainID(),
            .component_id = atom->GetComponentID(),
            .atom_id = atom->GetAtomID(),
            .alternate_indicator = atom->GetIndicator(),
            .element = atom->GetElement(),
            .residue = atom->GetResidue(),
            .spot = atom->GetSpot(),
            .structure = atom->GetStructure(),
            .position = atom->GetPosition()
        };
        switch (request.partial_charge_choice)
        {
        case PartialCharge::NEUTRAL: break;
        case PartialCharge::PARTIAL:
        case PartialCharge::AMBER:
            prepared.charge_lookup = ComponentHelper::LookupPartialCharge(
                prepared.residue, prepared.spot, prepared.structure,
                request.partial_charge_choice == PartialCharge::AMBER);
            prepared.charge_used = prepared.charge_lookup->charge.value_or(0.0);
            if (!prepared.charge_lookup->charge) ++failures[ChargeStatusText(prepared)];
            break;
        default: throw std::invalid_argument("Invalid simulation charge mode.");
        }
        for (size_t axis = 0; axis < 3; ++axis)
        {
            if (!std::isfinite(prepared.position[axis]))
                throw std::runtime_error("Non-finite simulation atom position.");
            result.range_min[axis] = std::min(result.range_min[axis], prepared.position[axis]);
            result.range_max[axis] = std::max(result.range_max[axis], prepared.position[axis]);
        }
        result.atom_list.emplace_back(std::move(prepared));
    }
    if (!result.atom_list.empty())
    {
        for (size_t axis = 0; axis < 3; ++axis)
        {
            result.range_min[axis] -= request.cutoff_distance;
            result.range_max[axis] += request.cutoff_distance;
        }
    }
    for (const auto & [reason, count] : failures)
    {
        Logger::Log(LogLevel::Warning, "Simulation charge lookup: " + std::string(reason)
            + " for " + std::to_string(count) + " atoms; using zero charge (recorded in simulation JSON).");
    }
    Logger::Log(LogLevel::Info, "Number of selected atoms to be simulated = "
        + std::to_string(result.atom_list.size()) + " / "
        + std::to_string(model.GetNumberOfAtom()) + " atoms.");
    return result;
}

int PopulateMapValueArray(MapObject & map, const SimulationAtomPreparationResult & atoms,
    const MapSimulationRequest & request, double blurring_width)
{
    ScopeTimer timer("MapSimulationCommand::PopulateMapValueArray");
    ElectricPotential potential;
    potential.SetBlurringWidth(blurring_width);
    potential.SetModelChoice(static_cast<int>(request.potential_model_choice));
    const auto grid_size{ map.GetGridSize() };
    const auto spacing{ map.GetGridSpacing() };
    const auto origin{ map.GetOrigin() };
    const auto radius_square{ request.cutoff_distance * request.cutoff_distance };
    const auto voxel_count{ map.GetMapValueArraySize() };
    auto values{ std::make_unique<double[]>(voxel_count) };

    struct Bounds { std::array<int, 3> lower{}, upper{}; };
    std::vector<Bounds> bounds(atoms.atom_list.size());
    std::vector<std::vector<size_t>> plane_atoms(static_cast<size_t>(grid_size[2]));
    for (size_t i = 0; i < atoms.atom_list.size(); ++i)
    {
        const auto & center{ atoms.atom_list[i].position };
        for (size_t axis = 0; axis < 3; ++axis)
        {
            const auto lower{ std::floor((center[axis] - request.cutoff_distance - origin[axis]) / spacing[axis]) };
            const auto upper{ std::floor((center[axis] + request.cutoff_distance - origin[axis]) / spacing[axis]) };
            bounds[i].lower[axis] = static_cast<int>(std::clamp(lower, 0.0, static_cast<double>(grid_size[axis])));
            bounds[i].upper[axis] = static_cast<int>(std::clamp(upper, -1.0, static_cast<double>(grid_size[axis] - 1)));
        }
        for (int z = bounds[i].lower[2]; z <= bounds[i].upper[2]; ++z)
            plane_atoms[static_cast<size_t>(z)].push_back(i);
    }

    int actual_job_count{ 1 };
    size_t completed_planes{ 0 };
    std::vector<std::exception_ptr> errors(plane_atoms.size());
#ifdef USE_OPENMP
    #pragma omp parallel num_threads(request.job_count)
#endif
    {
#ifdef USE_OPENMP
        #pragma omp single
        actual_job_count = omp_get_num_threads();
        #pragma omp for schedule(static)
#endif
        for (int z = 0; z < grid_size[2]; ++z)
        {
            try
            {
                for (const auto i : plane_atoms[static_cast<size_t>(z)])
                {
                    const auto & atom{ atoms.atom_list[i] };
                    for (int y = bounds[i].lower[1]; y <= bounds[i].upper[1]; ++y)
                    {
                        for (int x = bounds[i].lower[0]; x <= bounds[i].upper[0]; ++x)
                        {
                            const std::array<double, 3> position{
                                origin[0] + static_cast<double>(x) * spacing[0],
                                origin[1] + static_cast<double>(y) * spacing[1],
                                origin[2] + static_cast<double>(z) * spacing[2]
                            };
                            const auto dx{ position[0] - atom.position[0] };
                            const auto dy{ position[1] - atom.position[1] };
                            const auto dz{ position[2] - atom.position[2] };
                            if (dx * dx + dy * dy + dz * dz > radius_square) continue;
                            const auto index{ static_cast<size_t>(x) + static_cast<size_t>(grid_size[0])
                                * (static_cast<size_t>(y) + static_cast<size_t>(grid_size[1]) * static_cast<size_t>(z)) };
                            const auto distance{ array_helper::ComputeNorm(atom.position, map.GetGridPosition(index)) };
                            const auto contribution{ potential.GetPotentialValue(atom.element, distance, atom.charge_used) };
                            values[index] += contribution;
                            if (!std::isfinite(contribution) || !std::isfinite(values[index]))
                                throw std::runtime_error("Non-finite simulation potential at voxel " + std::to_string(index)
                                    + ", atom serial " + std::to_string(atom.serial_id));
                        }
                    }
                }
#ifdef USE_OPENMP
                #pragma omp critical(rhbm_simulation_progress)
#endif
                {
                    ++completed_planes;
                    Logger::ProgressPercent(completed_planes, plane_atoms.size());
                }
            }
            catch (...)
            {
                // Each plane owns its error slot; exceptions must not escape an OpenMP region.
                errors[static_cast<size_t>(z)] = std::current_exception();
            }
        }
    }
    for (const auto & error : errors) if (error) std::rethrow_exception(error);
    for (size_t i = 0; i < voxel_count; ++i)
    {
        if (std::abs(values[i]) > static_cast<double>(std::numeric_limits<float>::max()))
            throw std::runtime_error("Simulation voxel cannot be stored as finite float32: " + std::to_string(i));
    }
    map.ClearMapValueArray();
    map.SetMapValueArray(std::move(values));
    LogMapSummary(map);
    return actual_job_count;
}

} // namespace rhbm_gem::core::simulation
