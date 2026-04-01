#include "MapSimulationCommand.hpp"
#include "detail/DataObjectSummaryLog.hpp"
#include "data/object/MapSpatialIndex.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace {
constexpr std::string_view kModelKey{ "model" };
constexpr std::string_view kModelOption{ "--model" };
constexpr std::string_view kPotentialModelOption{ "--potential-model" };
constexpr std::string_view kChargeOption{ "--charge" };
constexpr std::string_view kCutoffOption{ "--cut-off" };
constexpr std::string_view kGridSpacingOption{ "--grid-spacing" };
constexpr std::string_view kBlurringWidthOption{ "--blurring-width" };

struct SimulationAtomPreparationResult
{
    std::vector<rhbm_gem::AtomObject *> atom_list;
    std::unordered_map<int, double> atom_charge_map;
    std::array<float, 3> range_minimum{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() };
    std::array<float, 3> range_maximum{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() };
    bool has_atom{ false };
};

double CalculateAtomChargeForSimulation(
    const rhbm_gem::AtomObject & atom,
    rhbm_gem::PartialCharge partial_charge_choice)
{
    switch (partial_charge_choice)
    {
    case rhbm_gem::PartialCharge::NEUTRAL:
        return 0.0;
    case rhbm_gem::PartialCharge::PARTIAL:
        return ComponentHelper::GetPartialCharge(
            atom.GetResidue(),
            atom.GetSpot(),
            atom.GetStructure());
    case rhbm_gem::PartialCharge::AMBER:
        return ComponentHelper::GetPartialCharge(
            atom.GetResidue(),
            atom.GetSpot(),
            atom.GetStructure(),
            true);
    default:
        Logger::Log(
            LogLevel::Error,
            "PrepareSimulationAtoms reached invalid partial-charge choice: "
            + std::to_string(static_cast<int>(partial_charge_choice)));
        return 0.0;
    }
}

SimulationAtomPreparationResult PrepareSimulationAtoms(
    rhbm_gem::ModelObject & model_object,
    rhbm_gem::PartialCharge partial_charge_choice,
    bool include_unknown_atoms)
{
    SimulationAtomPreparationResult result;
    result.atom_list.reserve(model_object.GetNumberOfAtom());
    for (auto & atom : model_object.GetAtomList())
    {
        if (!include_unknown_atoms && atom->IsUnknownAtom())
        {
            continue;
        }
        result.atom_list.emplace_back(atom.get());
        result.atom_charge_map.emplace(
            atom->GetSerialID(),
            CalculateAtomChargeForSimulation(*atom, partial_charge_choice));

        const auto & atom_position{ atom->GetPositionRef() };
        result.range_minimum[0] = std::min(result.range_minimum[0], atom_position[0]);
        result.range_minimum[1] = std::min(result.range_minimum[1], atom_position[1]);
        result.range_minimum[2] = std::min(result.range_minimum[2], atom_position[2]);
        result.range_maximum[0] = std::max(result.range_maximum[0], atom_position[0]);
        result.range_maximum[1] = std::max(result.range_maximum[1], atom_position[1]);
        result.range_maximum[2] = std::max(result.range_maximum[2], atom_position[2]);
    }
    result.has_atom = !result.atom_list.empty();
    return result;
}
}

namespace rhbm_gem {

MapSimulationCommand::MapSimulationCommand() :
    CommandWithRequest<MapSimulationRequest>{},
    m_selected_atom_list{}, m_atom_charge_map{}, m_model_object{ nullptr },
    m_atom_range_minimum{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() },
    m_atom_range_maximum{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() }
{
}

void MapSimulationCommand::NormalizeRequest()
{
    auto & request{ MutableRequest() };
    ValidateRequiredPath(request.model_file_path, kModelOption, "Model file");
    CoerceEnum(
        request.potential_model_choice,
        kPotentialModelOption,
        PotentialModel::FIVE_GAUS_CHARGE,
        "Potential model");
    CoerceEnum(
        request.partial_charge_choice,
        kChargeOption,
        PartialCharge::PARTIAL,
        "Partial charge choice");
    CoerceScalar(
        request.cutoff_distance,
        kCutoffOption,
        [](double candidate) { return candidate > 0.0; },
        5.0,
        LogLevel::Warning,
        "Cutoff distance must be positive, reset to default 5.0");
    CoerceScalar(
        request.grid_spacing,
        kGridSpacingOption,
        [](double candidate) { return candidate > 0.0; },
        0.5,
        LogLevel::Warning,
        "Grid spacing must be positive, reset to default 0.5");

    InvalidatePreparedState();
    ClearParseIssues(kBlurringWidthOption);
    std::vector<double> filtered_widths;
    filtered_widths.reserve(request.blurring_width_list.size());
    for (const auto width : request.blurring_width_list)
    {
        if (width <= 0.0)
        {
            AddNormalizationWarning(
                kBlurringWidthOption,
                "Blurring width must be positive, dropping current setting: "
                    + std::to_string(width));
            continue;
        }
        filtered_widths.push_back(width);
    }
    request.blurring_width_list = std::move(filtered_widths);
}

bool MapSimulationCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    RunMapSimulation();
    return true;
}

void MapSimulationCommand::ValidateOptions()
{
    const auto & request{ RequestOptions() };
    RequireCondition(
        !request.blurring_width_list.empty(),
        kBlurringWidthOption,
        "At least one positive blurring width is required.");
}

void MapSimulationCommand::ResetRuntimeState()
{
    m_selected_atom_list.clear();
    m_atom_charge_map.clear();
    m_model_object.reset();
    m_atom_range_minimum = {
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max()
    };
    m_atom_range_maximum = {
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest()
    };
}

bool MapSimulationCommand::BuildDataObject()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer("MapSimulationCommand::BuildDataObject");
    try
    {
        m_model_object = LoadInputFile<ModelObject>(
            request.model_file_path,
            std::string(kModelKey));
        BuildAtomList(m_model_object.get());
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "MapSimulationCommand::BuildDataObject : Failed to load model file '"
                + request.model_file_path.string() + "' as ModelObject: " + std::string(e.what()));
        return false;
    }
    return true;
}

void MapSimulationCommand::RunMapSimulation()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer("MapSimulationCommand::RunMapSimulation");
    
    Logger::Log(LogLevel::Info,
        "Total number of blurring width sets to be simulated: "
        + std::to_string(request.blurring_width_list.size()));
    
    auto map_object{ CreateMapObject() };
    for (auto & blurring_width : request.blurring_width_list)
    {
        auto map_key_tag{
            m_model_object->GetPdbID() + "_bw" +
            StringHelper::ToStringWithPrecision<double>(blurring_width, 2)
        };
        PopulateMapValueArray(map_object.get(), blurring_width);
        const auto output_file_name{
            BuildOutputPath(request.map_file_name + "_" + map_key_tag, ".map")
        };
        WriteMap(output_file_name, *map_object);
    }
}

void MapSimulationCommand::BuildAtomList(ModelObject * model_object)
{
    const auto & request{ RequestOptions() };
    if (model_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "MapSimulationCommand::BuildAtomList(): model object is null.");
        return;
    }

    auto result{ PrepareSimulationAtoms(*model_object, request.partial_charge_choice, true) };

    m_selected_atom_list = result.atom_list;
    m_atom_charge_map = result.atom_charge_map;
    if (result.has_atom)
    {
        m_atom_range_minimum = result.range_minimum;
        m_atom_range_maximum = result.range_maximum;
        m_atom_range_minimum[0] -= static_cast<float>(request.cutoff_distance);
        m_atom_range_minimum[1] -= static_cast<float>(request.cutoff_distance);
        m_atom_range_minimum[2] -= static_cast<float>(request.cutoff_distance);
        m_atom_range_maximum[0] += static_cast<float>(request.cutoff_distance);
        m_atom_range_maximum[1] += static_cast<float>(request.cutoff_distance);
        m_atom_range_maximum[2] += static_cast<float>(request.cutoff_distance);
    }

    Logger::Log(LogLevel::Info,
        "Number of selected atoms to be simulated = "
        + std::to_string(m_selected_atom_list.size()) +" / "
        + std::to_string(model_object->GetNumberOfAtom()) + " atoms.");
}

std::unique_ptr<MapObject> MapSimulationCommand::CreateMapObject()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer("MapSimulationCommand::CreateMapObject");
    std::array<float, 3> grid_spacing{
        static_cast<float>(request.grid_spacing),
        static_cast<float>(request.grid_spacing),
        static_cast<float>(request.grid_spacing)
    };

    auto origin{ CalculateOrigin(grid_spacing) };
    auto grid_size{ CalculateGridSize(grid_spacing) };
    auto map_object{ std::make_unique<MapObject>(grid_size, grid_spacing, origin) };
    map_object->SetThreadSize(ThreadSize());

    auto voxel_size{ map_object->GetMapValueArraySize() };
    auto map_value_array{ std::make_unique<float[]>(voxel_size) };
    std::fill_n(map_value_array.get(), voxel_size, 0.0f);
    map_object->SetMapValueArray(std::move(map_value_array));
    map_object->GetSpatialIndex().Build(ThreadSize());

    return map_object;
}

void MapSimulationCommand::PopulateMapValueArray(MapObject * map_object, double blurring_width)
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer("MapSimulationCommand::PopulateMapValueArray");
    Logger::Log(LogLevel::Info,
        " /- Start map value array production with blurring width = "+
        StringHelper::ToStringWithPrecision<double>(blurring_width, 2)
    );

    auto electric_potential{ std::make_unique<ElectricPotential>() };
    electric_potential->SetBlurringWidth(blurring_width);
    electric_potential->SetModelChoice(static_cast<int>(request.potential_model_choice));

    auto voxel_size{ map_object->GetMapValueArraySize() };
    auto map_value_array{ std::make_unique<float[]>(voxel_size) };
    std::fill_n(map_value_array.get(), voxel_size, 0.0f);

    auto atom_size{ m_selected_atom_list.size() };
    map_object->GetSpatialIndex().Build(ThreadSize());
    auto kd_tree_root{ map_object->GetSpatialIndex().GetRoot() };
    size_t atom_count{ 0 };
    std::vector<GridNode*> in_range_grid_node_list;

    Logger::Log(LogLevel::Info,
        " /- Total number of atoms to be processed: "+ std::to_string(atom_size) + " atoms."
    );
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(ThreadSize()) private(in_range_grid_node_list)
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ m_selected_atom_list[i] };
        auto charge{ m_atom_charge_map.at(atom->GetSerialID()) };
        auto element{ atom->GetElement() };
        auto atom_position{ atom->GetPosition() };
        MapObject query_map_object({1, 1, 1}, {1.0, 1.0, 1.0}, atom_position);
        GridNode query_node(0, &query_map_object);

        KDTreeAlgorithm<GridNode>::RangeSearch(
            kd_tree_root, &query_node, request.cutoff_distance, in_range_grid_node_list);
        
        for (auto & grid_node : in_range_grid_node_list)
        {
            auto distance{
                ArrayStats<float>::ComputeNorm(atom_position, grid_node->GetPosition())
            };
            map_value_array[grid_node->GetGridIndex()] += static_cast<float>(
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

    map_object->SetMapValueArray(std::move(map_value_array));
    command_detail::LogMapSummary(*map_object);
}

std::array<int, 3> MapSimulationCommand::CalculateGridSize(
    const std::array<float, 3> & grid_spacing) const
{
    auto selected_atom_size{ m_selected_atom_list.size() };
    if (selected_atom_size == 0)
    {
        Logger::Log(LogLevel::Warning, "No atoms selected. Grid size is set to [1,1,1].");
        return { 1, 1, 1 };
    }

    return {
        static_cast<int>(
            std::ceil((m_atom_range_maximum[0] - m_atom_range_minimum[0]) / grid_spacing[0])),
        static_cast<int>(
            std::ceil((m_atom_range_maximum[1] - m_atom_range_minimum[1]) / grid_spacing[1])),
        static_cast<int>(
            std::ceil((m_atom_range_maximum[2] - m_atom_range_minimum[2]) / grid_spacing[2]))
    };
}

std::array<float, 3> MapSimulationCommand::CalculateOrigin(
    const std::array<float, 3> & grid_spacing) const
{
    return {
        std::floor(m_atom_range_minimum[0] / grid_spacing[0]) * grid_spacing[0],
        std::floor(m_atom_range_minimum[1] / grid_spacing[1]) * grid_spacing[1],
        std::floor(m_atom_range_minimum[2] / grid_spacing[2]) * grid_spacing[2]
    };
}

} // namespace rhbm_gem
