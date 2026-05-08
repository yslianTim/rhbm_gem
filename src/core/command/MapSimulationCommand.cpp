#include "detail/CommandBase.hpp"
#include "detail/DataObjectSummaryLog.hpp"
#include "data/detail/MapSpatialIndex.hpp"

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

class MapSimulationCommand final : public CommandBase<MapSimulationRequest>
{
public:
    MapSimulationCommand();

private:
    void NormalizeAndValidateRequest(MapSimulationRequest & request) override;
    void ValidatePreparedRequest(const MapSimulationRequest & request) override;
    bool ExecuteImpl(const MapSimulationRequest & request) override;
};

namespace {

struct SimulationAtomPreparationResult
{
    std::vector<AtomObject *> atom_list;
    std::unordered_map<int, double> atom_charge_map;
    std::array<float, 3> range_min{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() };
    std::array<float, 3> range_max{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() };
    bool has_atom{ false };
};

double CalculateAtomChargeForSimulation(
    const AtomObject & atom,
    PartialCharge partial_charge_choice)
{
    switch (partial_charge_choice)
    {
    case PartialCharge::NEUTRAL:
        return 0.0;
    case PartialCharge::PARTIAL:
        return ComponentHelper::GetPartialCharge(
            atom.GetResidue(),
            atom.GetSpot(),
            atom.GetStructure());
    case PartialCharge::AMBER:
        return ComponentHelper::GetPartialCharge(
            atom.GetResidue(),
            atom.GetSpot(),
            atom.GetStructure(),
            true);
    default:
        Logger::Log(LogLevel::Error,
            "PrepareSimulationAtomList reached invalid partial-charge choice: "
            + std::to_string(static_cast<int>(partial_charge_choice)));
        return 0.0;
    }
}

SimulationAtomPreparationResult PrepareSimulationAtomList(
    ModelObject & model_object,
    const MapSimulationRequest & request)
{
    SimulationAtomPreparationResult result;
    result.atom_list.reserve(model_object.GetNumberOfAtom());
    for (auto & atom : model_object.GetAtomList())
    {
        result.atom_list.emplace_back(atom.get());
        result.atom_charge_map.emplace(
            atom->GetSerialID(),
            CalculateAtomChargeForSimulation(*atom, request.partial_charge_choice));

        const auto & atom_position{ atom->GetPositionRef() };
        result.range_min[0] = std::min(result.range_min[0], atom_position[0]);
        result.range_min[1] = std::min(result.range_min[1], atom_position[1]);
        result.range_min[2] = std::min(result.range_min[2], atom_position[2]);
        result.range_max[0] = std::max(result.range_max[0], atom_position[0]);
        result.range_max[1] = std::max(result.range_max[1], atom_position[1]);
        result.range_max[2] = std::max(result.range_max[2], atom_position[2]);
    }
    result.has_atom = !result.atom_list.empty();
    if (result.has_atom)
    {
        for (size_t i = 0; i < result.range_min.size(); ++i)
        {
            result.range_min[i] -= static_cast<float>(request.cutoff_distance);
            result.range_max[i] += static_cast<float>(request.cutoff_distance);
        }
    }

    Logger::Log(LogLevel::Info,
        "Number of selected atoms to be simulated = "
        + std::to_string(result.atom_list.size()) +" / "
        + std::to_string(model_object.GetNumberOfAtom()) + " atoms.");
    return result;
}

std::unique_ptr<MapObject> CreateMapObject(
    const MapSimulationRequest & request,
    const SimulationAtomPreparationResult & result)
{
    ScopeTimer timer("MapSimulationCommand::CreateMapObject");
    std::array<float, 3> grid_spacing{
        static_cast<float>(request.grid_spacing),
        static_cast<float>(request.grid_spacing),
        static_cast<float>(request.grid_spacing)
    };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    std::array<int, 3> grid_size{ 1, 1, 1 };
    if (result.has_atom)
    {
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

void PopulateMapValueArray(
    MapObject * map_object,
    const MapSpatialIndex & spatial_index,
    const SimulationAtomPreparationResult & atom_list,
    const MapSimulationRequest & request,
    double blurring_width,
    int thread_size)
{
    ScopeTimer timer("MapSimulationCommand::PopulateMapValueArray");
    Logger::Log(LogLevel::Info,
        " /- Start map value array production with blurring width = "+
        string_helper::ToStringWithPrecision<double>(blurring_width, 2)
    );

    auto electric_potential{ std::make_unique<ElectricPotential>() };
    electric_potential->SetBlurringWidth(blurring_width);
    electric_potential->SetModelChoice(static_cast<int>(request.potential_model_choice));

    auto voxel_size{ map_object->GetMapValueArraySize() };
    auto map_value_array{ std::make_unique<float[]>(voxel_size) };
    std::fill_n(map_value_array.get(), voxel_size, 0.0f);

    auto atom_size{ atom_list.atom_list.size() };
    size_t atom_count{ 0 };
    std::vector<size_t> in_range_grid_index_list;

    Logger::Log(LogLevel::Info,
        " /- Total number of atoms to be processed: "+ std::to_string(atom_size) + " atoms.");
    
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(thread_size) private(in_range_grid_index_list)
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ atom_list.atom_list[i] };
        auto charge{ atom_list.atom_charge_map.at(atom->GetSerialID()) };
        auto element{ atom->GetElement() };
        auto atom_position{ atom->GetPosition() };
        spatial_index.CollectGridIndicesInRange(
            atom_position, static_cast<float>(request.cutoff_distance), in_range_grid_index_list);

        for (const auto grid_index : in_range_grid_index_list)
        {
            auto distance{
                array_helper::ComputeNorm(atom_position, map_object->GetGridPosition(grid_index))
            };
            map_value_array[grid_index] += static_cast<float>(
                electric_potential->GetPotentialValue(element, distance, charge)
            );
        }

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_size);
        }
    }

    map_object->ClearMapValueArray();
    map_object->SetMapValueArray(std::move(map_value_array));
    command_detail::LogMapSummary(*map_object);
}
} // namespace

MapSimulationCommand::MapSimulationCommand() : CommandBase<MapSimulationRequest>{}
{
}

void MapSimulationCommand::NormalizeAndValidateRequest(MapSimulationRequest & request)
{
    RequireExistingPath(request, &MapSimulationRequest::model_file_path);
    RequireEnum(request, &MapSimulationRequest::potential_model_choice);
    RequireEnum(request, &MapSimulationRequest::partial_charge_choice);
    NormalizeFinitePositiveScalar(request, &MapSimulationRequest::cutoff_distance, 5.0);
    NormalizeFinitePositiveScalar(request, &MapSimulationRequest::grid_spacing, 0.5);

    std::vector<double> filtered_widths;
    filtered_widths.reserve(request.blurring_width_list.size());
    for (const auto width : request.blurring_width_list)
    {
        if (!numeric_validation::IsFinitePositive(width))
        {
            AddFieldNormalizationWarning(&MapSimulationRequest::blurring_width_list,
                "Blurring width must be a finite positive value, dropping current setting: "
                    + std::to_string(width));
            continue;
        }
        filtered_widths.push_back(width);
    }
    request.blurring_width_list = std::move(filtered_widths);
}

bool MapSimulationCommand::ExecuteImpl(const MapSimulationRequest & request)
{
    std::unique_ptr<ModelObject> model_object;
    try
    {
        model_object = ReadModel(request.model_file_path);
        model_object->SetKeyTag("model");
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "MapSimulationCommand::BuildDataObject : Failed to load model file '"
                + request.model_file_path.string() + "' as ModelObject: " + std::string(e.what()));
        return false;
    }

    auto atom_list{ PrepareSimulationAtomList(*model_object, request) };
    Logger::Log(LogLevel::Info,
        "Total number of blurring width sets to be simulated: "
        + std::to_string(request.blurring_width_list.size()));

    auto map_object{ CreateMapObject(request, atom_list) };
    MapSpatialIndex spatial_index(*map_object, request.job_count);
    spatial_index.Build();
    for (auto & blurring_width : request.blurring_width_list)
    {
        auto map_key_tag{
            model_object->GetPdbID() + "_bw" +
            string_helper::ToStringWithPrecision<double>(blurring_width, 2)
        };
        PopulateMapValueArray(
            map_object.get(), spatial_index, atom_list, request, blurring_width, request.job_count);
        auto output{ request.output_dir / (request.map_file_name + "_" + map_key_tag + ".map") };
        WriteMap(output, *map_object);
    }
    return true;
}

void MapSimulationCommand::ValidatePreparedRequest(const MapSimulationRequest & request)
{
    RequirePrepareCondition(
        !request.blurring_width_list.empty(),
        "At least one positive blurring width is required.");
}

namespace command_internal {

CommandResult ExecuteMapSimulationCommand(const MapSimulationRequest & request)
{
    MapSimulationCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem
